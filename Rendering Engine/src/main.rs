extern crate sdl2;
use sdl2::ttf;
use sdl2::pixels::Color;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::rect::Rect;
use sdl2::surface::Surface;
use sdl2::pixels::PixelFormatEnum;
use sdl2::mouse::{MouseWheelDirection,MouseState};
use sdl2::keyboard::KeyboardState;

use rand::{random};

use noise::{Perlin, Seedable};

use std::f32::INFINITY;
use std::collections::HashMap;
use std::time::Duration;
use std::time::Instant;
use std::str::FromStr;
use std::io;
use std::thread;
use std::f32::consts::PI;

mod modules;
use crate::modules::rendu::{DataRendu,Mesh,Vec3D,Camera,convert_rad,Triangle,cree_cube,regle_angle_deg,rasterisation_monocoeur};
use crate::modules::geo::*;
use crate::modules::chunks::{Chunks,GenParamTerre,BlockType,GenParamBase};

const LARGEUR_IMAGE:f32 = 1280.0;//1280.0
const HAUTEUR_IMAGE:f32 = 720.0;//720.0
const LARGEUR_IMAGE_USIZE:usize = LARGEUR_IMAGE as usize;
const HAUTEUR_IMAGE_USIZE:usize = HAUTEUR_IMAGE as usize;
const LARGEUR_IMAGE_I32:i32 = LARGEUR_IMAGE as i32;
const HAUTEUR_IMAGE_I32:i32 = HAUTEUR_IMAGE as i32;
const FROTTEMENTS:f32 = 0.8;
const LIMITE_FROTTEMENTS:f32 = 0.1;
const LUX_MIN:(u8,u8,u8) = (80,80,80);
const CHUNK_SIZE:usize = 10;
const CHUNK_SIZE_I64:i64 = CHUNK_SIZE as i64;
const CHUNK_SIZE_POW3:usize = CHUNK_SIZE.pow(3);
const RENDER_DISTANCE:i64 = 16;

struct MouseData {
    left:bool,
    right:bool,
    middle:bool,
    scroll:i32,
} 

// trait that defines what a user tool does
trait Outil {
    fn interaction_monde(&mut self, monde:&mut World, mode:ModeOutil, symmetrie:bool) -> SortieOutil;
    fn interaction_souris_cam(&mut self, souris:&MouseData, cam:&Camera,monde:&mut World, touches:&KeyboardState, symmetrie:bool) -> (SortieOutil, String);
}

enum ModeOutil {
    Pose((i64,i64,i64)),
    Casse(Point3D),
    PoseLux,
    EnleveDerniereLux
}

enum SortieOutil {
    ToucheMonde((i64,i64,i64)),
    ChangeCurseur([u16;6]),
    ToucheMondeVec(Vec<(i64,i64,i64)>),
    Rien,
}

impl SortieOutil {
    fn unwrap_pos(&self) -> (i64,i64,i64) {
        // à n'utiliser que si on sait que la sortie est ToucheMonde
        match self {
            &Self::ToucheMonde(coord) => return coord,
            _ => ()
        }
        (0,0,0)
    }
    fn unwrap_vecpos(&self) -> Vec<(i64,i64,i64)> {
        // à n'utiliser que si on sait que la sortie est ToucheMondeVec
        let mut out = Vec::new();
        match &self {
            &Self::ToucheMondeVec(coords) => for coord in coords {
                out.push(*coord)
            },
            _ => ()
        }
        out
    }
}
struct Bloc {
    textures:[u16;6],
}
fn lux_u8_a_f32((r,g,b):(u8,u8,u8)) -> (f32,f32,f32) {
    (r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
} 
struct VoxelBuilder {
    bloc_choisi:usize,
    blocs_dispo:Vec<Bloc>,
}

impl Outil for VoxelBuilder {
    fn interaction_monde(&mut self, monde:&mut World, mode:ModeOutil, symmetrie:bool) -> SortieOutil {
        match mode {
            ModeOutil::Casse(impact) => {
                let point_avant = impact.get_co_int();
                let mut lux_avant:(f32,f32,f32) = lux_u8_a_f32(monde.chunks.lumiere_a(point_avant));
                let mut max_lux = lux_u8_a_f32(LUX_MIN);
                let decalages:[(i64,i64,i64) ; 6] = [(1,0,0),(-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)];
                let mut plus_eclaire = point_avant;
                for decalage in decalages.iter() {
                    let lux_a = lux_u8_a_f32(monde.chunks.lumiere_a((point_avant.0 + decalage.0, point_avant.1 + decalage.1, point_avant.2 + decalage.2)));
                    if lux_avant < lux_a {
                        lux_avant = lux_a;
                        plus_eclaire = (point_avant.0 + decalage.0, point_avant.1 + decalage.1, point_avant.2 + decalage.2);
                    }
                }
                for decalage in decalages.iter() {
                    let lux_a = lux_u8_a_f32(monde.chunks.lumiere_a((plus_eclaire.0 + decalage.0, plus_eclaire.1 + decalage.1, plus_eclaire.2 + decalage.2)));
                    if max_lux < lux_a {
                        max_lux = lux_a;
                    }
                }
                let mut case = monde.chunks.get_float_mut(&impact);
                case.rempli = false;
                //case.interaction = None;
                if max_lux < lux_avant {
                    case.lux.0 = (lux_avant.0 * (max_lux.0/lux_avant.0) * 255.0) as u8;
                    case.lux.1 = (lux_avant.1 * (max_lux.1/lux_avant.1) * 255.0) as u8;
                    case.lux.2 = (lux_avant.2 * (max_lux.2/lux_avant.2) * 255.0) as u8;
                }
                else {
                    case.lux.0 = (lux_avant.0 * (lux_avant.0/max_lux.0) * 255.0) as u8;
                    case.lux.1 = (lux_avant.1 * (lux_avant.1/max_lux.1) * 255.0) as u8;
                    case.lux.2 = (lux_avant.2 * (lux_avant.2/max_lux.2) * 255.0) as u8;
                }
                monde.chunks.get_chunk_at(point_avant).num_remp -= 1;
                return SortieOutil::ToucheMonde(point_avant)
            },
            ModeOutil::Pose(pos) => {
                monde.chunks.get_mut_at(pos).rempli = true;
                monde.chunks.get_mut_at(pos).block_type = self.bloc_choisi as u16;
                //monde.chunks.get_mut_at(pos).interaction = self.blocs_dispo[self.bloc_choisi].interaction.clone();
                monde.chunks.get_chunk_at(pos).num_remp += 1;
                return SortieOutil::ToucheMonde(pos)
            },
            _ => SortieOutil::Rien,
        }
    }
    fn interaction_souris_cam(&mut self, souris:&MouseData,cam:&Camera, monde:&mut World, touches:&KeyboardState,symmetrie:bool) -> (SortieOutil, String) {
        let mut sortie = (SortieOutil::Rien, String::from("Pas d'interaction"));
        if souris.left || souris.right {
            let rayon = cam.get_raydir();
            match monde.chunks.ray_traverse(&rayon, 100.0) {
                Some(coef) => {
                    let point = rayon.at(coef);
                    let point_avant = rayon.at(coef - 0.1);
                    if point.dist(&monde.ents[0].pos) < 15.0 {
                        
                        if souris.left {
                            if symmetrie {
                                sortie = (SortieOutil::ToucheMondeVec(vec![
                                self.interaction_monde(monde, ModeOutil::Casse(point),symmetrie).unwrap_pos(),
                                self.interaction_monde(monde, ModeOutil::Casse(Point3D::new(-point.x, point.y, -point.z)),symmetrie).unwrap_pos()
                                ])
                                ,String::from("Pas d'interaction"));
                            }
                            else {
                                sortie = (self.interaction_monde(monde, ModeOutil::Casse(point),symmetrie),String::from("Pas d'interaction"));
                            }
                            thread::sleep(Duration::from_millis(50));
                        }
                        else if souris.right {
                            if symmetrie {
                                sortie = (SortieOutil::ToucheMondeVec(vec![
                                self.interaction_monde(monde, ModeOutil::Pose(point_avant.get_co_int()),symmetrie).unwrap_pos(),
                                self.interaction_monde(monde, ModeOutil::Pose(Point3D::new(-point_avant.x, point_avant.y, -point_avant.z).get_co_int()),symmetrie).unwrap_pos()
                                ])
                                ,String::from("Pas d'interaction"));
                            }
                            else {
                                sortie = (self.interaction_monde(monde, ModeOutil::Pose(point_avant.get_co_int()),symmetrie),String::from("Pas d'interaction"));
                            }
                            thread::sleep(Duration::from_millis(50));
                        }
                    }
                },
                None => ()
            }
        }
        if souris.scroll != 0 {
            if souris.scroll > 0 {
                self.bloc_choisi += (souris.scroll as usize).clamp(0, self.blocs_dispo.len() - self.bloc_choisi - 1);
            }
            else {
                self.bloc_choisi -= (souris.scroll.abs() as usize).clamp(0, self.bloc_choisi);
            }
            sortie = (SortieOutil::ChangeCurseur(self.blocs_dispo[self.bloc_choisi].textures),String::from("A"))
        }
        sortie
    }
}

struct Lumiere {
    pos:(i64,i64,i64),
    puissance:usize,
    col:(u8,u8,u8),
}
struct LuxBuilder {
    niveau:usize,
    lumieres_placees:Vec<Lumiere>,
    col:(u8,u8,u8),
}

impl Outil for LuxBuilder {
    fn interaction_monde(&mut self, monde:&mut World, mode:ModeOutil,symmetie:bool) -> SortieOutil {
        match mode {
            ModeOutil::Pose(pos) => {
                let mut peut = true;
                for lux in &self.lumieres_placees {
                    if lux.pos == pos {
                        peut = false;
                        break
                    }
                }
                if peut {
                    self.lumieres_placees.push(Lumiere{puissance:self.niveau, pos,col:self.col})
                }
                SortieOutil::Rien
            },
            ModeOutil::EnleveDerniereLux => {
                if self.lumieres_placees.len() > 0 {
                    let lux = &self.lumieres_placees[self.lumieres_placees.len() - 1];
                    let mut lux_update = Vec::<(i64,i64,i64)>::new();
                    let rayon_check = (lux.puissance + 10).div_euclid(CHUNK_SIZE) as i64;
                    let chunk_lux = (lux.pos.0.div_euclid(CHUNK_SIZE_I64), lux.pos.1.div_euclid(CHUNK_SIZE_I64), lux.pos.2.div_euclid(CHUNK_SIZE_I64));
                    for xc in chunk_lux.0-rayon_check..chunk_lux.0+rayon_check {
                        for yc in chunk_lux.1-rayon_check..chunk_lux.1+rayon_check {
                            for zc in chunk_lux.2-rayon_check..chunk_lux.2+rayon_check {
                                lux_update.push((xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64))
                            }
                        }
                    }
                    monde.chunks.enleve_lumiere(lux.pos, lux.puissance);
                    self.lumieres_placees.remove(self.lumieres_placees.len() - 1);
                    return SortieOutil::ToucheMondeVec(lux_update)
                }
                SortieOutil::Rien
            },
            ModeOutil::PoseLux => {
                let mut lux_update = Vec::<(i64,i64,i64)>::new();
                for lux in &self.lumieres_placees {
                    let rayon_check = (lux.puissance + 10).div_euclid(CHUNK_SIZE) as i64;
                    let chunk_lux = (lux.pos.0.div_euclid(CHUNK_SIZE_I64), lux.pos.1.div_euclid(CHUNK_SIZE_I64), lux.pos.2.div_euclid(CHUNK_SIZE_I64));
                    for xc in chunk_lux.0-rayon_check..chunk_lux.0+rayon_check {
                        for yc in chunk_lux.1-rayon_check..chunk_lux.1+rayon_check {
                            for zc in chunk_lux.2-rayon_check..chunk_lux.2+rayon_check {
                                let mut est_libre = true;
                                for chunk_lux in &lux_update {
                                    if (xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64) == *chunk_lux {
                                        est_libre = false;
                                        break
                                    }
                                }
                                if est_libre {
                                    lux_update.push((xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64))
                                }
                            }
                        }
                    }
                    monde.chunks.rajoute_lumiere(lux.pos, lux.puissance, lux.col);
                }
                SortieOutil::ToucheMondeVec(lux_update)
            }
            _ => SortieOutil::Rien,
        }
    }
    fn interaction_souris_cam(&mut self, souris:&MouseData, cam:&Camera, monde:&mut World, touches:&KeyboardState,symmetrie:bool) -> (SortieOutil, String) {
        if souris.left || souris.right || souris.middle {
            let rayon = cam.get_raydir();
            match monde.chunks.ray_traverse(&rayon, 100.0) {
                Some(coef) => {
                    if souris.left {
                        if symmetrie {
                            let mut data_out = self.interaction_monde(monde, ModeOutil::EnleveDerniereLux,symmetrie).unwrap_vecpos();
                            for pos in self.interaction_monde(monde, ModeOutil::EnleveDerniereLux,symmetrie).unwrap_vecpos() {
                                data_out.push(pos);
                            }
                            return (SortieOutil::ToucheMondeVec(data_out), format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2))
                        }
                        else {
                            return (self.interaction_monde(monde, ModeOutil::EnleveDerniereLux,symmetrie), format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2));
                        }
                    }
                    else if souris.right {
                        std::thread::sleep(Duration::from_secs_f32(0.1));
                        return (self.interaction_monde(monde, ModeOutil::Pose(rayon.at(coef - 0.1).get_co_int()),symmetrie), format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2));
                    }
                    else if souris.middle {
                        if symmetrie {
                            let mut data_out = self.interaction_monde(monde, ModeOutil::PoseLux,symmetrie).unwrap_vecpos();
                            for pos in self.interaction_monde(monde, ModeOutil::PoseLux,symmetrie).unwrap_vecpos() {
                                data_out.push(pos);
                            }
                            return (SortieOutil::ToucheMondeVec(data_out), format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2))
                        }
                        else {
                            return (self.interaction_monde(monde, ModeOutil::PoseLux,symmetrie), format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2));
                        }
                    }
                },
                None => ()
            }
        }
        if souris.scroll != 0 {
            if souris.scroll > 0 {
                self.niveau += (souris.scroll as usize).clamp(0, 10);
            }
            else {
                self.niveau -= (souris.scroll.abs() as usize).clamp(0, self.niveau);
            }
        }
        if touches.is_scancode_pressed(sdl2::keyboard::Scancode::X) {
            self.col = (demande_input_recurs::<u8>("Red ?"),demande_input_recurs::<u8>("Green ?"),demande_input_recurs::<u8>("Blue ?"));
        }
        (SortieOutil::Rien, format!("lights:{}, power:{}, col:({},{},{})", self.lumieres_placees.len(), self.niveau, self.col.0, self.col.1, self.col.2))
    }
}


struct Vaisseau {
    ent:usize,
    carburant:i32,
    carb_regen:i32,
    max_carburant:i32,
    health:i32,
    max_health:i32,
    armor:i32,
    equipe:u8,
    munitions:i32,
    munitions_max:i32,
}

impl Vaisseau {
    fn new(ent:usize,max_health:i32, max_carb:i32,carb_regen:i32,armor:i32,equipe:u8,munitions_max:i32) -> Vaisseau {
        Vaisseau { ent, carburant: max_carb, carb_regen, max_carburant: max_carb, health: max_health, max_health, armor, equipe, munitions: munitions_max, munitions_max}
    }
}

#[derive(Clone, Copy)]
struct EntityController {
    pos:Point3D,
    model:usize,
    taille:f32,
    vitesse:Point3D,
    orientation:(f32,f32),
}

impl EntityController {
    fn new(pos:Point3D, model:usize,taille:f32, orientation:(f32,f32)) -> EntityController {
        EntityController{pos,model,taille,vitesse:Point3D::new(0.0, 0.0, 0.0), orientation}
    }
    fn frottements(&mut self) {
        self.vitesse.x *= FROTTEMENTS;
        self.vitesse.y *= FROTTEMENTS;
        self.vitesse.z *= FROTTEMENTS;
        if self.vitesse.z.abs() < LIMITE_FROTTEMENTS && self.vitesse.y.abs() < LIMITE_FROTTEMENTS && self.vitesse.x.abs() < LIMITE_FROTTEMENTS {
            self.vitesse.z = 0.0;
            self.vitesse.x = 0.0;
            self.vitesse.y = 0.0;
        }
    }
    fn collisions_volume(&mut self,volume:&mut Chunks,vitesse:&Point3D) -> bool {
        
        if volume.plein_a(((self.pos.x + vitesse.x) as i64, (self.pos.y + vitesse.y) as i64, (self.pos.z + vitesse.z) as i64)) {
            true
        } 
        else {
            false
        }
    }
    fn collisions_comp(&mut self,volume:&mut Chunks) -> bool {
        // if speed is smaller than 0.5, only test at the next position, else test at multiple positions until then
        let totvitesse = Point3D::new(0.0,0.0,0.0).dist(&self.vitesse);
        if totvitesse <= 0.5 {
            let dupli_vitesse = self.vitesse;
            if self.collisions_volume(volume, &dupli_vitesse) {
                self.vitesse = Point3D::new(0.0,0.0,0.0);
                return true
            }
        }
        else {
            let mut coef = 0.0;
            while coef < totvitesse {
                let vitesse_test = Point3D::new(self.vitesse.x * coef, self.vitesse.y * coef, self.vitesse.z * coef);
                if self.collisions_volume(volume, &vitesse_test) {
                    self.vitesse = Point3D::new(0.0,0.0,0.0);
                    return true
                }
                coef += 0.1;
            }
        }
        false
    }
    fn mouvement(&mut self,volume:&mut Chunks) -> bool {
        if self.collisions_comp(volume) {
            true
        }
        else {
            self.pos.x += self.vitesse.x;
            self.pos.y += self.vitesse.y;
            self.pos.z += self.vitesse.z;
            self.frottements();
            false
        }
    }
}
struct PlayerController {
    cam:Camera,
    ent:usize,
    chunks_charges:HashMap<(i64,i64,i64), Vec<Triangle>>,
    outils:Vec<Box<dyn Outil>>,
    outil_choisi:usize,
    vaisseau:usize,
    pseudo:String,
    symmetrie:bool,
    render_distance:i64,
}

impl PlayerController {
    fn new(pseudo:String,ent:usize,vaisseau:usize) -> PlayerController {
        PlayerController{ent, cam:Camera{co: Vec3D {x:5.2235, y:5.224, z:5.2342, col:(0.0, 0.0)}, orient_h:0.1, orient_v:0.1,roll:0.1},chunks_charges:HashMap::new(), outil_choisi:0, vaisseau,
        outils:vec![
            Box::new(VoxelBuilder{bloc_choisi:0, 
                blocs_dispo:vec![
                Bloc{textures:[0 ; 6]},
                Bloc{textures:[1 ; 6]},
                Bloc{textures:[2 ; 6]},
                Bloc{textures:[3 ; 6]},
                Bloc{textures:[4 ; 6]},
                Bloc{textures:[5 ; 6]},
                Bloc{textures:[6 ; 6]},
                Bloc{textures:[7 ; 6]},
                ] 
            
            }),
            Box::new(LuxBuilder{niveau:0, lumieres_placees:Vec::new(),col:(255,255,255)})
            ],
        pseudo,
        symmetrie:false,render_distance:12}
    }
    fn gere_outils(&mut self, souris:&MouseData, monde:&mut World, meshes: &mut Vec<Vec<Mesh>>, touches:&KeyboardState) -> String {
        // handles using the current tool at the disposal of the user, and putting the 3D cursor where they are looking
        let propre_outil = &mut self.outils[self.outil_choisi];
        let ray = self.cam.get_raydir();
        match monde.chunks.ray_traverse(&ray, 100.0) {
            Some(coef) => {
                let impact = ray.at(coef);
                let mut mesh = &mut meshes[1][0];
                mesh.pos.x = impact.x * 10.0;
                mesh.pos.y = impact.y * 10.0;
                mesh.pos.z = impact.z * 10.0;
            },
            None => (),
        }
        let resultat = propre_outil.interaction_souris_cam(souris, &self.cam, monde,touches,self.symmetrie);
        match resultat.0 {
            SortieOutil::Rien => (),
            SortieOutil::ChangeCurseur(texture) => {
                meshes[1][0].triangles.clear();
                cree_cube(2.0, [(255,255,255);6], texture, &mut meshes[1][0].triangles)
            },
            SortieOutil::ToucheMonde(pos) => self.update_static(&mut monde.chunks, pos),
            SortieOutil::ToucheMondeVec(positions) => for pos in positions {self.update_static(&mut monde.chunks, pos)},
        }
        resultat.1
    }

    fn update_coord(&mut self, monde:&World) {
        self.cam.co = Vec3D::new(monde.ents[self.ent].pos.x * 10.0, monde.ents[self.ent].pos.y * 10.0,monde.ents[self.ent].pos.z * 10.0,(0.0,0.0));
    }
    fn genere_static(&mut self,chunks:&mut Chunks) -> Vec<Triangle> {
        let mut triangles = Vec::new();
        let (xc, yc, zc) = ((self.cam.co.x as i64).div_euclid(CHUNK_SIZE_I64*10), (self.cam.co.y as i64).div_euclid(CHUNK_SIZE_I64*10), (self.cam.co.z as i64).div_euclid(CHUNK_SIZE_I64*10));
        for x in xc - RENDER_DISTANCE..xc + RENDER_DISTANCE {
            for y in yc - RENDER_DISTANCE..yc + RENDER_DISTANCE {
                for z in zc - RENDER_DISTANCE..zc + RENDER_DISTANCE {
                    match self.chunks_charges.get(&(x,y,z)) {
                        Some(triangles_chunk) => {
                            for triangle in triangles_chunk {
                                triangles.push(*triangle);
                            }
                        },
                        None => {
                            let mut triangles_chunk = Vec::new();
                            chunks.get_triangles_chunk(&mut triangles_chunk, (x,y,z));
                            for triangle in &triangles_chunk {
                                triangles.push(*triangle);
                            }
                            self.chunks_charges.insert((x,y,z), triangles_chunk);
                        }
                    }
                }
            }
        }
        triangles
    }

    fn genere_static_cone(&mut self,chunks:&mut Chunks,render_distance:i64) -> Vec<Triangle> {
        // this is supposed to generate the meshes of all chunks in a frustrum in front of the camera, but it only works properly when looking horizontally
        let mut triangles = Vec::with_capacity(render_distance as usize * CHUNK_SIZE_POW3 as usize * 10);
        let (xc, yc, zc) = ((self.cam.co.x as i64).div_euclid(CHUNK_SIZE_I64*10), (self.cam.co.y as i64).div_euclid(CHUNK_SIZE_I64*10), (self.cam.co.z as i64).div_euclid(CHUNK_SIZE_I64*10));
        let pointcam = Point3D::new( self.cam.co.x / (CHUNK_SIZE as f32 * 10.0), self.cam.co.y / (CHUNK_SIZE as f32 * 10.0), self.cam.co.z / (CHUNK_SIZE as f32 * 10.0));
        let rayhautgauche = Rayon::new_orient(pointcam, convert_rad(self.cam.orient_h - 45.0), convert_rad(self.cam.orient_v + 45.0));
        let raybasgauche = Rayon::new_orient(pointcam, convert_rad(self.cam.orient_h - 45.0), convert_rad(self.cam.orient_v));
        let rayhautdroite = Rayon::new_orient(pointcam, convert_rad(self.cam.orient_h + 45.0), convert_rad(self.cam.orient_v + 45.0));
        let mut affiches = HashMap::<(i64,i64,i64), bool>::new();
        match self.chunks_charges.get(&(xc,yc,zc)) {
            Some(triangles_chunk) => {
                for triangle in triangles_chunk {
                    triangles.push(*triangle);
                }
            },
            None => {
                let mut triangles_chunk = Vec::with_capacity(2000);
                chunks.get_triangles_chunk(&mut triangles_chunk, (xc,yc,zc));
                for triangle in &triangles_chunk {
                    triangles.push(*triangle);
                }
                self.chunks_charges.insert((xc,yc,zc), triangles_chunk);
            }
        }
        affiches.insert((xc,yc,zc), true);
        let mut avancement = 0.0;
        while avancement <= render_distance as f32 {
            let hautgauche = rayhautgauche.at(avancement);
            let basgauche = raybasgauche.at(avancement);
            let hautdroite = rayhautdroite.at(avancement);
            let dist_vert = hautgauche.dist(&basgauche);
            let dist_horiz = hautgauche.dist(&hautdroite);
            let ray_vert = Rayon::new_orient_vers(hautgauche, basgauche);
            let mut av_vert = 0.0;
            while av_vert <= dist_vert {
                let mut ray_horiz = Rayon::new_orient_vers(hautgauche, hautdroite);
                ray_horiz.orig = ray_vert.at(av_vert);
                let mut av_horiz = 0.0;
                while av_horiz <= dist_horiz {
                    let point_chunk = ray_horiz.at(av_horiz).co_chunk();
                    //if point_chunk.1 <= 300 && point_chunk.1 >= -300 {
                    match affiches.get(&point_chunk) {
                        Some(_) => (),
                        None => {
                            match self.chunks_charges.get(&point_chunk) {
                                Some(triangles_chunk) => {
                                    for triangle in triangles_chunk {
                                        triangles.push(*triangle);
                                    }
                                },
                                None => {
                                    let mut triangles_chunk = Vec::with_capacity(2000);
                                    chunks.get_triangles_chunk(&mut triangles_chunk, point_chunk);
                                    for triangle in &triangles_chunk {
                                        triangles.push(*triangle);
                                    }
                                    self.chunks_charges.insert(point_chunk, triangles_chunk);
                                }
                            }
                            affiches.insert(point_chunk, true);
                        }
                    }
                    //}
                    av_horiz += 0.5
                }
                av_vert += 0.5
            }
            avancement += 0.5
        }
        triangles
    }

    fn update_static(&mut self, chunks:&mut Chunks, coord:(i64,i64,i64)) {
        // when modifying the world, this function handles updating the mesh of a given chunk and all the ones around it
        let (xc,yc,zc) = (coord.0.div_euclid(CHUNK_SIZE_I64), coord.1.div_euclid(CHUNK_SIZE_I64), coord.2.div_euclid(CHUNK_SIZE_I64));
        let (xv,yv,zv) = (coord.0 - xc * CHUNK_SIZE_I64, coord.0 - yc * CHUNK_SIZE_I64, coord.0 - zc * CHUNK_SIZE_I64);
        let mut liste_a_tester = vec![(xc, yc, zc)];
        liste_a_tester.push((xc-1,yc,zc));
        liste_a_tester.push((xc+1,yc,zc));
        liste_a_tester.push((xc,yc-1,zc));
        liste_a_tester.push((xc,yc+1,zc));
        liste_a_tester.push((xc,yc,zc-1));
        liste_a_tester.push((xc,yc,zc+1));
        for test in liste_a_tester {
            match self.chunks_charges.get_mut(&test) {
                Some(triangles_chunk) => {
                    triangles_chunk.clear();
                    chunks.get_triangles_chunk(triangles_chunk, test);
                },
                None => ()
            }
        }
    }
}

struct World {
    ents:Vec<EntityController>,
    chunks:Chunks,
}

impl World {
    fn new(chunks:Chunks) -> World {
        World{ents:Vec::new(), chunks}
    }

    fn update_spd(&mut self) {
        for ent in &mut self.ents {
            ent.mouvement(&mut self.chunks);
        }
    }
    fn update_tout(&mut self) {
        self.update_spd();
    }
}





fn demande_input(texte_input: &str) -> String {
    let mut input = String::new();
    println!("{}", texte_input);
    match io::stdin().read_line(&mut input) {
        Ok(_) => {
            ()
        },
        Err(e) => println!("There was an error : {}", e)
    }
    input
}

fn demande_input_recurs<T:FromStr>(text_input:&str) -> T {
    // asks for a given type until given a valid answer
    match T::from_str(&demande_input(text_input).trim()) {
        Ok(output) => return output,
        Err(_) => return demande_input_recurs(text_input),
    }
}


fn clear_buffer(buffer:&mut Vec<u8>, couleur:(u8,u8,u8)) {
    // clears a framebuffer
    let mut i = 0;
    while i < LARGEUR_IMAGE_USIZE*HAUTEUR_IMAGE_USIZE*3 {
        buffer[i] = couleur.0;
        buffer[i+1] = couleur.1;
        buffer[i+2] = couleur.2;
        i += 3;
    }
}

fn clear_zbuffer(buffer:&mut Vec<f32>) {
    // clears a zbuffer
    let mut i = 0;
    while i < LARGEUR_IMAGE_USIZE*HAUTEUR_IMAGE_USIZE*3 {
        buffer[i] = INFINITY;
        i += 3;
    }
}

fn get_carte_valide() -> Chunks {
    // asks for a new map name until a valid one is given
    match Chunks::charge_fichier(demande_input("Name of the map ?").trim()) {
        Ok(carte) => carte,
        Err(_) => get_carte_valide(),
    }
}

fn orient_souris_simple(souris:(i32,i32),sensi:f32, (orient_h, orient_v):(f32,f32)) -> (f32,f32) {
    (
        regle_angle_deg(orient_v - (((souris.1 - HAUTEUR_IMAGE_I32/2) / 2) as f32 * sensi / (HAUTEUR_IMAGE)) * 45.0),
        regle_angle_deg(orient_h - (((-souris.0 + LARGEUR_IMAGE_I32/2) / 2) as f32 * sensi / (LARGEUR_IMAGE)) * 29.5)
    )
}

fn gere_event(player:&mut PlayerController, monde:&mut World,event:Event,nb_quit:&mut usize, scroll:&mut i32, a_focus:&mut bool,data_rendu:&DataRendu) -> bool { // stop ?
    // this fuction handles part of user input
    let camera = &mut player.cam;
    let mut ent_cam = &mut monde.ents[player.ent];
    match event { 
        Event::Quit {..} |
        Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
            *nb_quit += 1;
            if *nb_quit > 50 {
                return true;
            }
        },
        Event::KeyDown { keycode: Some(Keycode::Space), .. } => {
            ent_cam.vitesse.y += 0.5;
        },
        Event::KeyDown { keycode: Some(Keycode::LCtrl), .. } => {
            ent_cam.vitesse.y -= 0.5;
        },
        Event::KeyDown { keycode: Some(Keycode::Z), .. } => {
            let raydir = camera.get_raydir();
            ent_cam.vitesse.x += raydir.dir.x * 0.5;
            ent_cam.vitesse.y += raydir.dir.y * 0.5;
            ent_cam.vitesse.z += raydir.dir.z * 0.5;
        },
        Event::KeyDown { keycode: Some(Keycode::S), .. } => {
            let raydir = camera.get_raydir();
            ent_cam.vitesse.x -= raydir.dir.x * 0.5;
            ent_cam.vitesse.y -= raydir.dir.y * 0.5;
            ent_cam.vitesse.z -= raydir.dir.z * 0.5;
        },
        Event::KeyDown { keycode: Some(Keycode::T), .. } => {
            if player.outil_choisi < player.outils.len() - 1 {
                player.outil_choisi += 1;
            }
        },
        Event::KeyDown { keycode: Some(Keycode::G), .. } => {
            if player.outil_choisi != 0 {
                player.outil_choisi -= 1;
            }
        },
        Event::KeyDown { keycode: Some(Keycode::U), .. } => {
            if player.render_distance < 64 {
                player.render_distance += 1;
            }
        },
        Event::KeyDown { keycode: Some(Keycode::J), .. } => {
            if player.render_distance > 1 {
                player.render_distance -= 1;
            }
        },
        Event::KeyDown { keycode: Some(Keycode::Y), .. } => {
            if player.symmetrie {
                player.symmetrie = false;
            }
            else {
                player.symmetrie = true;
            }
        },
        Event::KeyDown {keycode:Some(Keycode::B), ..} => {
            let mut perlin = Perlin::new();
            let seed = random();
            perlin = perlin.set_seed(seed);
            let mut perlin_temp = Perlin::new();
            perlin_temp = perlin_temp.set_seed(seed);
            monde.chunks = Chunks::new((0,0,0), Some(Box::new(GenParamTerre{seed, freq_3d:0.05, max_height:100.0,perlin,perlin_temp,freq_temp:0.005,num_text:monde.chunks.num_text, frequency:0.005,niveau_mer:27.0,profondeur:4})),monde.chunks.num_text,monde.chunks.block_types.clone());
            player.chunks_charges.clear();
        }

        Event::KeyDown {keycode:Some(Keycode::N), ..} => {
            let mut perlin = Perlin::new();
            let seed = random();
            perlin = perlin.set_seed(seed);
            monde.chunks = Chunks::new((0,0,0), Some(Box::new(GenParamBase{seed, perlin,num_text:monde.chunks.num_text,min_plein1:0.0,min_plein2:0.0})),monde.chunks.num_text,monde.chunks.block_types.clone());
            player.chunks_charges.clear();
        }
        Event::KeyDown {keycode:Some(Keycode::V), ..} => {
            monde.chunks = Chunks::new((0,0,0), None,monde.chunks.num_text,monde.chunks.block_types.clone());
            monde.chunks.get_mut_at((0,0,0)).block_type = 1;
            monde.chunks.get_mut_at((0,1,0)).rempli = false;
            ent_cam.pos = Point3D::new(0.0,1.0,0.0);
            player.chunks_charges.clear();
        }
        Event::KeyDown {keycode:Some(Keycode::C), ..} => {
            monde.chunks = get_carte_valide();
            ent_cam.pos = Point3D::new(0.0,1.0,0.0);
            player.chunks_charges.clear();
        }

        Event::KeyDown { keycode: Some(Keycode::L), .. } => {
            monde.chunks.sauvegarde_fichier(demande_input("name of the file ?").trim(), "CHUNKS DE TEST 1", data_rendu.get_textures());
        },
        Event::MouseWheel {y:posy, direction: MouseWheelDirection::Normal, ..} => {
            *scroll = posy;
        },
        Event::Window {win_event:event,..} => {
            match event {
                sdl2::event::WindowEvent::FocusGained => *a_focus = true,
                sdl2::event::WindowEvent::FocusLost => *a_focus = false,
                _ => (),
            }
        }
        _ => (),
    }
    false
}
fn introduction() {
    println!("Hello, this is my software rendering engine");
    println!("It is currently single-threaded");
    println!("The keys currently used for interaction are : ");
    println!(" - Z -> forwards / S -> backwards");
    println!(" - Space -> Up / left control -> down");
    println!(" - L -> Saves the current map / C -> loads a map");
    println!(" - B -> generates a new map with 'continental' generation / N -> generates a new map with 3D perlin noise generation / V -> generates a completely full map");
    println!(" - T -> scrolls up one tool / G -> scrolls down one tool (there are currently only 2 tools)");
    println!(" - Escape -> stops the program when held down for long enough");
    println!(" - U -> increases render distance (up to 64) / J -> reduces render distance (down to 1)");
    println!("There are currently 2 tools : ");
    println!(" - a 'Voxel Builder' which allows the user to place and destroy blocks from a given selection left click destroys, right click places, scrollwheel changes the held block");
    println!(" - a 'Light Builder' which allows the user to place and remove lights with a reach given by the scroll wheel. Left click places invisible lights, right click removes them, and scroll wheel click computes the spread of the lights.");
    println!("   X allows you to pick the R,G,B color of the light.");
    println!("Focusing on another window will keep this program from warping the mouse to the center of its window");
    println!("Have fun !");
}

fn mode_constru() -> Result<(),String> {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    
    let window = video_subsystem.window("This is a window", LARGEUR_IMAGE as u32, HAUTEUR_IMAGE as u32)
        .position_centered()
        .build()
        .unwrap();
    let mut canvas = window.into_canvas().build().unwrap();
    let mut event_pump = sdl_context.event_pump().unwrap();
    // initialisation du côtté graphique
    
    canvas.set_draw_color(Color::RGB(0, 255, 255));
    canvas.clear();
    canvas.present();
    let mut player = PlayerController::new("Pseudo ?".to_string(),0,0);
    player.cam.co.y = 100.0;

    let mut data_rendu = DataRendu::new();
    data_rendu.charge_textures(vec!["eau_prof.png"], 1);
    data_rendu.charge_textures(vec!["eau.png"],1);
    data_rendu.charge_textures(vec!["sable.png"],1);
    data_rendu.charge_textures(vec!["terre_herbe.png"],1);
    data_rendu.charge_textures(vec!["terre_cail.png"],1);
    data_rendu.charge_textures(vec!["terre.png"],1);
    data_rendu.charge_textures(vec!["roche.png"],1);
    data_rendu.charge_textures(vec!["neige.png"],1);

    let mut perlin = Perlin::new();
    let seed = 14123;
    perlin = perlin.set_seed(seed);
    let mut monde = World::new(Chunks::new((0,0,0), Some(Box::new(GenParamBase{seed,perlin,num_text:8,min_plein1:0.0,min_plein2:0.0})), 8,
    vec![
        BlockType{texture:0},
        BlockType{texture:1},
        BlockType{texture:2},
        BlockType{texture:3},
        BlockType{texture:4},
        BlockType{texture:5},
        BlockType{texture:6},
        BlockType{texture:7},
    ]));
    monde.ents.push(EntityController::new(Point3D::new(player.cam.co.x,player.cam.co.y,player.cam.co.z), 1, 1.0,(0.0,0.0)));

    let ttf_context = ttf::init().map_err(|e| e.to_string())?;

    let mut font = ttf_context.load_font("bitmap fonts/VCR_OSD_MONO_1.001.ttf", 20)?;
    let texture_creator = canvas.texture_creator();

    let mut souris = event_pump.mouse_state();
    let mut scroll = 0;
    let mut string_outil = String::from("AAAAA");

    let mut a_focus = false;
    let mut framebuffer_fini = vec![0;LARGEUR_IMAGE_USIZE*HAUTEUR_IMAGE_USIZE*3];
    let mut zbuffer = vec![INFINITY; LARGEUR_IMAGE_USIZE*HAUTEUR_IMAGE_USIZE*3];
    let mut nb_quit = 0;

    let mut heure = 0.0;

    'running: loop  {
        
        canvas.clear();
        event_pump.pump_events();
        let clavier = KeyboardState::new(&event_pump);
        souris = MouseState::new(&event_pump);
        if a_focus {
            player.cam.orient_souris((souris.x(), souris.y()), 2.5);
            sdl_context.mouse().warp_mouse_in_window(canvas.window(), LARGEUR_IMAGE_I32/2, HAUTEUR_IMAGE_I32/2);
        }
        let mousedata = MouseData {left:souris.left(), right:souris.right(), middle:souris.middle(), scroll};
        if scroll != 0 {
            scroll = 0;
        }
        string_outil = player.gere_outils(&mousedata, &mut monde, &mut data_rendu.meshes, &clavier);
        for event in event_pump.poll_iter() {
            if gere_event(&mut player, &mut monde, event, &mut nb_quit, &mut scroll, &mut a_focus, &data_rendu) {
                break 'running;
            }
        }
        let debut = Instant::now();
        data_rendu.meshes[0][0].triangles = player.genere_static_cone(&mut monde.chunks,player.render_distance);
        
        rasterisation_monocoeur(data_rendu.get_meshes(), &player.cam, data_rendu.get_textures(), &mut framebuffer_fini, &mut zbuffer,heure);
        
        
        let surface = font
        .render(format!("FPS : {} ; TRIS : {} ; RENDER DISTANCE : {}", 1.0/debut.elapsed().as_secs_f32(), data_rendu.meshes[0][0].triangles.len(), player.render_distance).as_str())
        .blended(Color::RGBA(255, 0, 0, 255))
        .map_err(|e| e.to_string())?;
        let texture_text = texture_creator
        .create_texture_from_surface(&surface)
        .map_err(|e| e.to_string())?;
        

        let surface_outil = font
        .render(string_outil.as_str())
        .blended(Color::RGBA(255, 0, 0, 255))
        .map_err(|e| e.to_string())?;
        

        let texture_text_outil = texture_creator
        .create_texture_from_surface(&surface_outil)
        .map_err(|e| e.to_string())?;
        

        let long_framebuffer = framebuffer_fini.len();
        let surf_image = Surface::from_data(&mut framebuffer_fini[0..long_framebuffer] ,LARGEUR_IMAGE as u32, HAUTEUR_IMAGE as u32, (LARGEUR_IMAGE * 3.0)as u32 ,PixelFormatEnum::RGB24).unwrap();
        match texture_creator.create_texture_from_surface(surf_image) {
            Ok(texture) => {
                canvas.copy(&texture, None, None);
                canvas.copy(&texture_text, None, Rect::new(0,0,surface.width(),surface.height()));
                canvas.copy(&texture_text_outil, None, Rect::new(0,surface.height() as i32,surface_outil.width(),surface_outil.height()));
            },
            Err(err) => println!("grosse erreur {}", err)
        };

        canvas.present();

        
        data_rendu.update_texture_st();
        monde.update_tout();
        player.update_coord(&monde);

        let heure_soleil = heure + PI/2.0;
        
        clear_buffer(&mut framebuffer_fini, (0,(80.0 + 70.0 * heure_soleil.cos()) as u8,(128.0 + 100.0 * heure_soleil.cos()) as u8));
        clear_zbuffer(&mut zbuffer);
        heure += 0.01;
    }
    Ok(())
}

fn demande_mode() {
    introduction();
    match demande_input("Do you want to start (Y/N) ?").trim() {
        "Y" => {mode_constru();},
        _ => println!("Ok goodbye !")
    }
}

fn main() {

    demande_mode();
}
