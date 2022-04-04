extern crate image;

use image::io::Reader as ImageReader;
use image::ImageBuffer;
use image::Rgb;
use sdl2::rect::Point;

use std::f32::INFINITY;
use std::f32::consts::PI;

use super::geo::*;

const NEAR_CLIPPING_PLANE:f32 = 1.0;
const LARGEUR_IMAGE:f32 = 1280.0;//1280.0
const HAUTEUR_IMAGE:f32 = 720.0;//720.0
const MOIT_LARGEUR_IMAGE:f32 = LARGEUR_IMAGE/2.0;
const MOIT_HAUTEUR_IMAGE:f32 = HAUTEUR_IMAGE/2.0;
const LARGEUR_IMAGE_USIZE:usize = LARGEUR_IMAGE as usize;
const HAUTEUR_IMAGE_USIZE:usize = HAUTEUR_IMAGE as usize;
const LARGEUR_IMAGE_I32:i32 = LARGEUR_IMAGE as i32;
const HAUTEUR_IMAGE_I32:i32 = HAUTEUR_IMAGE as i32;
const ASPECT_RATIO:f32 = LARGEUR_IMAGE/HAUTEUR_IMAGE;
const PI_SUR_180:f32 = PI/180.0;
const DRAW_DISTANCE:f32 = 1.0/8000.0;

#[derive(Copy, Clone)]
pub struct Vec3D {
    pub x:f32,
    pub y:f32,
    pub z:f32,
    pub col:(f32, f32),
}

impl Vec3D {
    pub fn new(x:f32, y:f32, z:f32, col:(f32, f32)) -> Vec3D {
        Vec3D{x,y,z,col}
    }
}

// textures have their own type that allow for dynamic animation
// detatick is the number of frames between each step of a texture animation
// data is a Vec<> containing tuples of a Vec<u8> representing a texture, and its average color
pub struct RasterTexture {
    largeur:f32,
    hauteur:f32,
    data:Vec<(Vec<(u8, u8, u8)>, (u8,u8,u8))>,
    index:usize,
    tick:usize,
    pub deltatick:usize,
    pub noms:Vec<String>,
}
impl RasterTexture {
    fn avance_index(&mut self) {
        self.tick += 1;
        if self.tick >= self.deltatick {
            self.index += 1;
            self.tick = 0;
        }
        if self.index >= self.data.len() {
            self.index = 0;
        }
    }
    fn get_data(&self) -> &Vec<(u8,u8,u8)> {
        &self.data[self.index].0
    }
    fn get_moy(&self) -> (u8,u8,u8) {
        self.data[self.index].1
    }
}
pub fn regle_angle_deg(ang: f32) -> f32 {
    let mut nouv_ang:f32 = ang;
    if ang < 0.0 {
        nouv_ang = 360.0 + ang;
    }
    else if ang > 360.0 {
        nouv_ang = ang - 360.0;
    }
    nouv_ang
}

#[derive(Clone,Copy)]
pub struct Triangle {
    V1:Vec3D,
    V2:Vec3D,
    V3:Vec3D,
    texture_id:u16,
    lumiere:(u8,u8,u8),
}

impl Triangle {
    pub fn new(V1:Vec3D, V2:Vec3D, V3:Vec3D, texture_id:u16, lumiere:(u8,u8,u8)) -> Triangle {
        Triangle{V1, V2, V3, texture_id, lumiere}
    }
    fn calcule_normale(&self) -> Point3D {
        let u = Point3D::new(self.V2.x - self.V1.x,self.V2.y - self.V1.y,self.V2.z - self.V1.z);
        let v = Point3D::new(self.V3.x - self.V1.x,self.V3.y - self.V1.y,self.V3.z - self.V1.z);
        Point3D::new(
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x
        )

    }
    fn raster_points(&self, camera:&Vec3D, sinang:f32, cosang:f32, cosangv:f32, sinangv:f32) -> Triangle {
        Triangle::new(self.V1.to_raster(camera, sinang, cosang, cosangv, sinangv), self.V2.to_raster(camera, sinang, cosang, cosangv, sinangv), self.V3.to_raster(camera, sinang, cosang, cosangv, sinangv), self.texture_id, self.lumiere)
    }
    fn range_points(&self) -> (&Vec3D, &Vec3D, &Vec3D) {
        // sorts points of a triangle from top to bottom, meant to be used on triangles already prohjected onto the view frustrum
        let top:&Vec3D;
        let bot:&Vec3D;
        let autre:&Vec3D;
        if self.V1.y < self.V2.y && self.V1.y < self.V3.y {
            top = &self.V1;
            if self.V2.y > self.V3.y {
                bot = &self.V2;
                autre = &self.V3;
            }
            else {
                bot = &self.V3;
                autre = &self.V2;
            }
        }
        else if self.V2.y < self.V1.y && self.V2.y < self.V3.y {
            top = &self.V2;
            if self.V1.y > self.V3.y {
                bot = &self.V1;
                autre = &self.V3;
            }
            else {
                bot = &self.V3;
                autre = &self.V1;
            }
        }
        else {
            top = &self.V3;
            if self.V1.y > self.V2.y {
                bot = &self.V1;
                autre = &self.V2;
            }
            else {
                bot = &self.V2;
                autre = &self.V1;
            }
        }
        
        (top,bot,autre)
    }
    

    #[inline(always)]
    fn dans_triangle_aire_garanti(&self, point:&Vec2D, aire:&f32, larg:&f32, long:&f32, unsurv1z:&f32,unsurv2z:&f32,unsurv3z:&f32) -> (f32,f32,f32) {
        let w0 = edge_function(&self.V2, &self.V3, point) * aire;
        let w1 = edge_function(&self.V3, &self.V1, point) * aire;
        let w2 = edge_function(&self.V1, &self.V2, point) * aire;
        let z = 1.0/(w0  * self.V1.z + w1  * self.V2.z + w2  * self.V3.z );


        ( 
            (self.V1.col.0 * w0 + self.V2.col.0 * w1 + self.V3.col.0 * w2) * z * larg,
            (self.V1.col.1 * w0 + self.V2.col.1 * w1 + self.V3.col.1 * w2) * z * long,
            w0 * unsurv1z + 
            w1 * unsurv2z + 
            w2 * unsurv3z
            
        )
    }
    #[inline(always)]
    fn dans_triangle_aire_garanti_que_d(&self, point:&Vec2D, aire:&f32,unsurv1z:&f32,unsurv2z:&f32,unsurv3z:&f32) -> f32 {
        let w0 = edge_function(&self.V2, &self.V3, point) * aire;
        let w1 = edge_function(&self.V3, &self.V1, point) * aire;
        let w2 = edge_function(&self.V1, &self.V2, point) * aire;


        w0 * unsurv1z + w1 * unsurv2z + w2 * unsurv3z
        
    }
}

impl Vec3D {
    fn to_raster(&self,camera:&Vec3D,sinang:f32, cosang:f32, cosangv:f32, sinangv:f32) -> Vec3D {
        let mut zcam = self.z - camera.z;
        let mut xcam = self.x - camera.x;
        let xcam1 = xcam;
        xcam = xcam * cosang - zcam * sinang;
        zcam = zcam * cosang + xcam1 * sinang;

        let mut ycam = self.y - camera.y;
        
        let ycam1 = ycam;
        
        ycam = ycam * cosangv - zcam * sinangv;
        zcam = zcam * cosangv + ycam1 * sinangv;

        
        let z = 1.0/zcam;
        Vec3D {
            x:(1.0 + (NEAR_CLIPPING_PLANE * xcam * z)) * MOIT_LARGEUR_IMAGE,
            y:(1.0 - (NEAR_CLIPPING_PLANE * ycam * z)) * MOIT_HAUTEUR_IMAGE * ASPECT_RATIO,
            z,
            col: (self.col.0 * z, self.col.1 * z),
        }
    }
}

struct Vec2D {
    x:f32,
    y:f32,
}

fn charge_image(nom: &str) -> Result<ImageBuffer<Rgb<u8>, Vec<u8>>, ()> {
    // tries to load an image in the "textures" folder
    let chemin = format!("textures/{}", nom);
    match ImageReader::open(chemin.as_str()) {
        Ok(texture) => {
            let image_buffer = texture.decode().unwrap().to_rgb8();
            return Ok(image_buffer);
        },
        Err(err) => { println!("texture : {} introuvable dans le dossier cartes... erreur {}", nom, err); return Err(())}
    }
}

fn charge_image_sur(nom:&str) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
    // if it can't load the given image, load the default "arbre.png" image
    match charge_image(nom) {
        Ok(image_buffer) => return image_buffer,
        Err(()) => return charge_image("arbre.png").unwrap(),
    }
}


#[inline(always)]
fn edge_function(p1:&Vec3D, p2:&Vec3D, p3:&Vec2D) -> f32 { // p3 est la cible
    (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x)
}


fn edge_function_2D(p1:&Vec3D, p2:&Vec3D, p3:&Vec3D) -> f32 { // p3 est la cible
    (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x)
}

fn aire_raster(points:&Triangle) -> f32 {
    edge_function_2D(&points.V1, &points.V2, &points.V3)
}

pub fn cree_cube(taille:f32, lumiere:[(u8,u8,u8);6], textures:[u16;6], triangles: &mut Vec<Triangle>) {
    //creates a cube around the origin of a mesh

    //z - 1
    triangles.push(Triangle::new( 
        Vec3D::new(-1.0 * taille,-1.0 * taille, -1.0 * taille, (0.0,1.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, -1.0 * taille, (1.0,1.0)),
        Vec3D::new( -1.0 * taille, 1.0 * taille, -1.0 * taille, (0.0,0.0)),
        textures[5],
        lumiere[5]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(-1.0 * taille,1.0 * taille, -1.0 * taille, (0.0,0.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, -1.0 * taille, (1.0,1.0)),
        Vec3D::new( 1.0 * taille, 1.0 * taille, -1.0 * taille, (1.0,0.0)),
        textures[5],
        lumiere[5]
    ));

    //z + 1
    triangles.push(Triangle::new(
        Vec3D::new(-1.0 * taille,1.0 * taille, 1.0 * taille, (0.0,0.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( -1.0 * taille, -1.0 * taille, 1.0 * taille, (0.0,1.0)),
        textures[4],
        lumiere[4]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,1.0 * taille, 1.0 * taille, (1.0,0.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( -1.0 * taille, 1.0 * taille, 1.0 * taille, (0.0,0.0)),
        textures[4],
        lumiere[4]
    ));

    // x - 1
    triangles.push(Triangle::new(
        Vec3D::new(-1.0 * taille,1.0 * taille, -1.0 * taille, (0.0,0.0)),
        Vec3D::new(-1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( -1.0 * taille, -1.0 * taille, -1.0 * taille, (0.0,1.0)),
        textures[3],
        lumiere[3]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(-1.0 * taille,1.0 * taille, 1.0 * taille, (1.0,0.0)),
        Vec3D::new(-1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( -1.0 * taille, 1.0 * taille, -1.0 * taille, (0.0,0.0)),
        textures[3],
        lumiere[3]
    ));

    // x + 1
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,-1.0 * taille, -1.0 * taille, (0.0,1.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( 1.0 * taille, 1.0 * taille, -1.0 * taille, (0.0,0.0)),
        textures[2],
        lumiere[2]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,1.0 * taille, -1.0 * taille, (0.0,0.0)),
        Vec3D::new(1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new( 1.0 * taille, 1.0 * taille, 1.0 * taille, (1.0,0.0)),
        textures[2],
        lumiere[2]
    ));

    // y + 1
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,1.0 * taille, -1.0 * taille, (1.0,0.0)),
        Vec3D::new(-1.0 * taille, 1.0 * taille, 1.0 * taille, (0.0,1.0)),
        Vec3D::new( -1.0 * taille, 1.0 * taille, -1.0 * taille, (0.0,0.0)),
        textures[1],
        lumiere[1]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,1.0 * taille, 1.0 * taille, (1.0,1.0)),
        Vec3D::new(-1.0 * taille, 1.0 * taille, 1.0 * taille, (0.0,1.0)),
        Vec3D::new( 1.0 * taille, 1.0 * taille, -1.0 * taille, (1.0,0.0)),
        textures[1],
        lumiere[1]
    ));

    // y - 1
    triangles.push(Triangle::new(
        Vec3D::new(-1.0 * taille,-1.0 * taille, -1.0 * taille, (0.0,0.0)),
        Vec3D::new(-1.0 * taille, -1.0 * taille, 1.0 * taille, (0.0,1.0)),
        Vec3D::new( 1.0 * taille, -1.0 * taille, -1.0 * taille, (1.0,0.0)),
        textures[0],
        lumiere[0]
    ));
    triangles.push(Triangle::new(
        Vec3D::new(1.0 * taille,-1.0 * taille, -1.0 * taille, (1.0,0.0)),
        Vec3D::new(-1.0 * taille, -1.0 * taille, 1.0 * taille, (0.0,1.0)),
        Vec3D::new( 1.0 * taille, -1.0 * taille, 1.0 * taille, (1.0,1.0)),
        textures[0],
        lumiere[0]
    ));
}

// this struct houses data directly related to rendering
pub struct DataRendu {
    pub meshes:Vec<Vec<Mesh>>,
    pub textures:Vec<RasterTexture>
}

impl DataRendu {
    pub fn new() -> DataRendu {
        let mut triangles_curseur = Vec::<Triangle>::new();
        cree_cube(2.0, [(255,255,255);6],[1,1,1,1,1,1], &mut triangles_curseur);
        DataRendu {
            meshes:vec![
            vec![Mesh::new(Point3D::new(0.0, 0.0, 0.0),0.0,0.0,Vec::new())], // monde
            vec![Mesh::new(Point3D::new(0.0,0.0,0.0), 0.0, 0.0,triangles_curseur)], // personnage
            Vec::new(), // entités
            Vec::new(), // projectiles
            ],
            textures:Vec::new()
        }
    }
    pub fn update_texture_st(&mut self) {
        for texture in &mut self.textures {
            texture.avance_index();
        }
    }
    pub fn charge_textures(&mut self, noms:Vec<&str>, deltatick:usize) {
        // loads multiples textures and adds them to self's Vec<RasterTexture> as a singular struct
        let image_buffer = charge_image_sur(noms[0]);
        let (w,h) = image_buffer.dimensions();
        let mut i = 0;
        let mut texture = RasterTexture{largeur:w as f32, hauteur: h as f32, data:Vec::new(),index:0, tick:0, deltatick,noms:Vec::new()};
        for nom in noms {
            let image_buffer = charge_image_sur(nom);
            texture.data.push((Vec::new(),(0,0,0)));
            for y in 0..h {
                for x in 0..w {
                    let pixel = image_buffer.get_pixel(x, y);
                    texture.data[i].0.push((pixel.0[0], pixel.0[1], pixel.0[2]));
                }
            }
            // this part calculates the average color of the texture
            let mut nb_moy:(u64,u64,u64) = (0,0,0);
            for col in &texture.data[i].0 {
                nb_moy.0 += col.0 as u64;
                nb_moy.1 += col.1 as u64;
                nb_moy.2 += col.2 as u64;
            }
            nb_moy.0 /= (w*h) as u64;
            nb_moy.1 /= (w*h) as u64;
            nb_moy.2 /= (w*h) as u64;
            texture.data[i].1 = (nb_moy.0 as u8, nb_moy.1 as u8, nb_moy.2 as u8);


            i += 1;
            texture.noms.push(nom.to_string());
        }
        self.textures.push(texture);
    }
    pub fn get_meshes(&self) -> &Vec<Vec<Mesh>> {
        &self.meshes
    }
    pub fn get_textures(&self) -> &Vec<RasterTexture> {
        &self.textures
    }
}

pub struct Mesh {
    pub triangles:Vec<Triangle>,
    pub pos:Point3D,
    orientation_h:f32,
    orientation_v:f32,
    cosangh:f32,
    sinangh:f32,
    cosangv:f32,
    sinangv:f32,
    visible:bool,
}

impl Mesh {
    fn new(pos:Point3D,orienth:f32,orientv:f32,triangles: Vec<Triangle>) -> Mesh {
        Mesh { triangles, pos, orientation_h:orienth, cosangh: orienth.cos(), sinangh: orienth.sin(), orientation_v:orientv, sinangv:orientv.sin(), cosangv:orientv.cos(),visible:true}
    }
    fn new_avec_vis(pos:Point3D,orienth:f32,orientv:f32,triangles: Vec<Triangle>, visible:bool) -> Mesh {
        Mesh { triangles, pos, orientation_h:orienth, cosangh: orienth.cos(), sinangh: orienth.sin(), orientation_v:orientv, sinangv:orientv.sin(), cosangv:orientv.cos(), visible}
    }
    fn update_orient(&mut self, orienth:f32, orientv:f32) {
        // updates orientation of the mesh and precalculates sine and cosine
        self.cosangh = orienth.cos();
        self.sinangh = orienth.sin();
        self.orientation_h = orienth;
        self.cosangv = orientv.cos();
        self.sinangv = orientv.sin();
        self.orientation_v = orientv;
    }
    fn get_wt_triangle_raster_at(&self,index:usize) -> Triangle {
        // rotates the model around its own origin point before adding the mesh world coords to get the triangle in the world
        // but the rotation doesn't really work well

        let triangle = self.triangles[index];
        
        let world_x_1 = (triangle.V1.x * self.cosangh - triangle.V1.z * self.sinangh) + self.pos.x;
        let mut world_z_1 = triangle.V1.z * self.cosangh + triangle.V1.x * self.sinangh;

        let world_x_2 = (triangle.V2.x * self.cosangh - triangle.V2.z * self.sinangh) + self.pos.x;
        let mut world_z_2 = triangle.V2.z * self.cosangh + triangle.V2.x * self.sinangh;

        let world_x_3 = (triangle.V3.x * self.cosangh - triangle.V3.z * self.sinangh) + self.pos.x;
        let mut world_z_3 = triangle.V3.z * self.cosangh + triangle.V3.x * self.sinangh;

        let world_y_1 = (triangle.V1.y * self.cosangv - world_z_1 * self.sinangv) + self.pos.y;
        world_z_1 = (world_z_1 * self.cosangv + triangle.V1.y * self.sinangv) + self.pos.z;

        let world_y_2 = (triangle.V2.y * self.cosangv - world_z_2 * self.sinangv) + self.pos.y;
        world_z_2 = (world_z_2 * self.cosangv + triangle.V2.y * self.sinangv) + self.pos.z;

        let world_y_3 = (triangle.V3.y * self.cosangv - world_z_3 * self.sinangv) + self.pos.y;
        world_z_3 = (world_z_3 * self.cosangv + triangle.V3.y * self.sinangv) + self.pos.z;

        Triangle::new(Vec3D::new(world_x_1,world_y_1,world_z_1,triangle.V1.col),
        Vec3D::new(world_x_2,world_y_2,world_z_2,triangle.V2.col),
        Vec3D::new(world_x_3,world_y_3,world_z_3,triangle.V3.col),
        triangle.texture_id,
        triangle.lumiere
        )
    }
    fn duplique(&self) -> Mesh {
        let mut triangles_out = Vec::new();
        for triangle in &self.triangles {
            triangles_out.push(*triangle)
        }
        Mesh { triangles:triangles_out, pos: self.pos, orientation_h: self.orientation_h, cosangh: self.cosangh, sinangh: self.sinangh, visible: true, orientation_v:self.orientation_v, cosangv:self.cosangv, sinangv:self.sinangv}
    }
    fn duplique_avec_orient_pos(&self,pos:Point3D ,orientation_h:f32) -> Mesh {
        let mut triangles_out = Vec::new();
        for triangle in &self.triangles {
            triangles_out.push(*triangle)
        }
        Mesh { triangles:triangles_out, pos, orientation_h, cosangh: orientation_h.cos(), sinangh: orientation_h.sin(), visible: true, orientation_v:self.orientation_v, cosangv:self.cosangv, sinangv:self.sinangv }
    }
}

pub struct Camera {
    pub co:Vec3D,
    pub orient_h:f32,
    pub orient_v:f32,
    pub roll:f32, // unused
}
impl Camera {
    pub fn get_orient_h_centre(&self) -> f32 {
        self.orient_h * PI_SUR_180
    }
    pub fn get_orient_v_centre(&self) -> f32 {
        (self.orient_v + 23.5) * PI_SUR_180
    }
    pub fn get_orient_rad(&self) -> f32 {
        self.orient_h * PI_SUR_180
    }
    pub fn get_orient_v_rad(&self) -> f32 {
        self.orient_v * PI_SUR_180
    }
    pub fn get_roll_rad(&self) -> f32 {
        self.orient_v * PI_SUR_180
    }
    pub fn orient_souris(&mut self, souris:(i32,i32),sensi:f32) {
        self.orient_v = regle_angle_deg(self.orient_v - (((souris.1 - HAUTEUR_IMAGE_I32/2) / 2) as f32 * sensi / (HAUTEUR_IMAGE)) * 45.0);
        self.orient_h = regle_angle_deg(self.orient_h - (((-souris.0 + LARGEUR_IMAGE_I32/2) / 2) as f32 * sensi / (LARGEUR_IMAGE)) * 29.5);
    }
    pub fn get_raydir(&self) -> Rayon {
        // a function to get a ray from the center of the FOV
        let angv = self.get_orient_v_centre();
        Rayon {orig:Point3D::new(self.co.x / 10.0, self.co.y / 10.0, self.co.z / 10.0), dir:Point3D{x:self.get_orient_h_centre().sin() * angv.cos(), y:angv.sin(), z:self.get_orient_h_centre().cos() * angv.cos()}}
    }
}

pub fn convert_rad(ang:f32) -> f32 {
    ang * PI_SUR_180
}

fn convert_heure_orient(heure:f32) -> (f32,f32) {
    (
        convert_rad(30.0),
        heure,
    )
}

pub fn rasterisation_monocoeur(vecmeshes:&Vec<Vec<Mesh>>, camera:&Camera,textures:&Vec<RasterTexture>, framebuffer:&mut Vec<u8>, zbuffer:&mut Vec<f32>,heure:f32) {
    // this function handles rasterisation of all meshes

    let orient_h_rad = camera.get_orient_rad();
    let cosang = orient_h_rad.cos();
    let sinang = orient_h_rad.sin();

    
    let orient_v_rad = camera.get_orient_v_rad();
    let sinangv = orient_v_rad.sin();
    let cosangv = orient_v_rad.cos();
    //let normale_cam = Point3D::get_vec_orient(camera.get_orient_h_centre(), camera.get_orient_v_centre()).normalise();

    let orient_heure = convert_heure_orient(heure);
    let normale_soleil = Point3D::get_vec_orient(orient_heure.0, orient_heure.1);

    let mut point:Vec2D;

    let mut texture_data:&Vec<(u8,u8,u8)>;
    let mut largeur_texture:f32;
    let mut largeur_texture_usize:usize;
    let mut hauteur_texture:f32;


    let mut posbuf:usize;
    let mut index_texture:usize;

    let mut col:(u8,u8,u8);

    let mut y1:i32;
    let mut y2:i32;
    let mut y3:i32;
    let mut x1:f32;
    let mut x2:f32;
    let mut dxim:f32;
    let mut dyim:f32;
    let mut deltad:f32;
    let mut unsuraymty:f32;
    let mut unsurbymty:f32;
    let mut unsurbymay:f32;
    let mut x2mx1:f32;

    let mut unsurv1z:f32;
    let mut unsurv2z:f32;
    let mut unsurv3z:f32;

    let mut ymy1:f32;
    let mut ymy2:f32;

    let mut texture_len:usize;
    let mut zbuffer_a:&mut f32;

    let mut col_src:[u8;3];

    let mut aire:f32;
    let mut prod_scal:f32;
    let mut triangle_lumiere:(f32,f32,f32);

    // this is using the scanline method, as i could not manage to properly optimise the edge function method

    for meshes in vecmeshes.iter() {
        for mesh in meshes {
            if mesh.visible { // meshes can be made invisible
                let mut index = 0;
                while index < mesh.triangles.len() {
                    let triangle = mesh.get_wt_triangle_raster_at(index);
                    let triangle_raster = triangle.raster_points(&camera.co, sinang, cosang, cosangv, sinangv);
                    
                    
                    if triangle_raster.V1.z > 0.0 && triangle_raster.V2.z > 0.0 && triangle_raster.V3.z > 0.0 {
                        aire = aire_raster(&triangle_raster);
                        if aire > 0.0 {
                            
                            triangle_lumiere = lux_u8_a_f32(triangle.lumiere);
                            let normale_triangle = triangle.calcule_normale().normalise();
                            prod_scal = normale_triangle.produit_scalaire(&normale_soleil).clamp(0.2, 1.0);
                            triangle_lumiere.0 *= prod_scal;
                            triangle_lumiere.1 *= prod_scal;
                            triangle_lumiere.2 *= prod_scal;
                            // 1/z is precalculated, which is why the comparison is > for closer triangles
                            // this comparison lets us have 2 levels of LOD for textures, close and far
                            if triangle_raster.V1.z > DRAW_DISTANCE && triangle_raster.V2.z > DRAW_DISTANCE && triangle_raster.V3.z > DRAW_DISTANCE {
                                aire = 1.0/aire;
                                texture_data = textures[triangle.texture_id as usize].get_data();
                                largeur_texture = textures[triangle.texture_id as usize].largeur;
                                largeur_texture_usize = largeur_texture as usize;
                                hauteur_texture = textures[triangle.texture_id as usize].hauteur;
                                let (top,bot,autre) = triangle_raster.range_points();
                                // "range_points()" gives us the highest point of the triangle, lowest, and the middle one, with their UV coords
                                y1 = (top.y as i32).clamp(0,HAUTEUR_IMAGE_I32);
                                y2 = autre.y as i32;
                                y3 = (bot.y as i32).clamp(0,HAUTEUR_IMAGE_I32);

                                point = Vec2D{x:0.0, y:y1 as f32 + 0.5};


                                // all possible division is precalculated
                                unsurv1z = 1.0/triangle_raster.V1.z;
                                unsurv2z = 1.0/triangle_raster.V2.z;
                                unsurv3z = 1.0/triangle_raster.V3.z;

                                unsuraymty = (autre.x - top.x)/(autre.y - top.y);
                                unsurbymty = (bot.x - top.x)/(bot.y - top.y);
                                unsurbymay = (bot.x - autre.x)/(bot.y - autre.y);
                            
                                ymy2 = if autre.y < 0.0 {-autre.y} else {0.0};
                                ymy1 = if top.y < 0.0 {-top.y} else {0.0};

                                texture_len = texture_data.len() - 1;
                                'loopy:for y in y1..y3 {
                                    if y < y2 {
                                        x1 = ymy1 * unsuraymty + top.x;
                                    }
                                    else {
                                        x1 = ymy2 * unsurbymay  + autre.x;
                                        ymy2 += 1.0;
                                    }
                                    x2 = ymy1 * unsurbymty + top.x;
                                    ymy1 += 1.0;
                                    if x1 > x2 {
                                        std::mem::swap(&mut x1, &mut x2);
                                    }
                                    // here texture coords and distance to camera are calculated
                                    point.x = x1;
                                    let (mut xi,mut yi,mut d1) = triangle_raster.dans_triangle_aire_garanti(&point, &aire, &largeur_texture, &hauteur_texture,&unsurv1z,&unsurv2z,&unsurv3z);
                                    point.x = x2;
                                    let (xf, yf, d2) = triangle_raster.dans_triangle_aire_garanti(&point, &aire, &largeur_texture, &hauteur_texture,&unsurv1z,&unsurv2z,&unsurv3z);

                                    x2mx1 = 1.0/(x2 - x1);

                                    // calculating what to increment each texture coord and distance with as we draw the line
                                    dxim = (xf - xi) * x2mx1;
                                    dyim = (yf - yi) * x2mx1;
                                    deltad = (d2 - d1) * x2mx1;
                                    // if the line begins outside of the screen, increment accordingly before clamping to the start of the screen
                                    if x1 < 0.0 {
                                        xi += dxim*(-x1);
                                        yi += dyim*(-x1);
                                        x1 = 0.0;
                                    }

                                    // none of what is unsafe in this block is proven in any way, everything just worked with its safe counterpart before making it unsafe
                                    unsafe {
                                        posbuf = (y as usize * LARGEUR_IMAGE_USIZE + x1.to_int_unchecked::<usize>()) * 3;

                                        if x2 > LARGEUR_IMAGE {
                                            x2 = LARGEUR_IMAGE
                                        }
                                        for x in x1.to_int_unchecked::<i32>()..x2.to_int_unchecked::<i32>() {
                                            zbuffer_a = zbuffer.get_unchecked_mut(posbuf);
                                            if zbuffer_a > &mut d1 {
                                                *zbuffer_a = d1;
                                                index_texture = std::cmp::min(yi.to_int_unchecked::<usize>() * largeur_texture_usize + xi.to_int_unchecked::<usize>(),texture_len);
                                                col = *texture_data.get_unchecked(index_texture);
                                                framebuffer[posbuf..posbuf+3].copy_from_slice(&[(col.0 as f32 * triangle_lumiere.0).to_int_unchecked(), (col.1 as f32 * triangle_lumiere.1).to_int_unchecked(), (col.2 as f32 * triangle_lumiere.2).to_int_unchecked()]);
                                            }
                                            posbuf += 3;
                                            xi += dxim;
                                            yi += dyim;
                                            d1 += deltad;
                                        }
                                    }
                                    point.y += 1.0;
                                }
                            
                            }
                            // this is the level of LOD for far away triangles, colored using only the average color of their texture (which has been precalcultated when loading the texture)
                            // it doesn't give that much of a performance boost but it can be noticeable
                            else {
                                // none of what is unsafe in this block is proven in any way, everything just worked with its safe counterpart before making it unsafe
                                unsafe {
                                aire = 1.0/aire;
                                col = textures[triangle.texture_id as usize].get_moy();
                                col = (
                                    (col.0 as f32 * triangle_lumiere.0).to_int_unchecked(),
                                    (col.1 as f32 * triangle_lumiere.1).to_int_unchecked(),
                                    (col.2 as f32 * triangle_lumiere.2).to_int_unchecked()
                                );
                                let (top,bot,autre) = triangle_raster.range_points();
                                y1 = top.y.to_int_unchecked::<i32>().clamp(0, HAUTEUR_IMAGE_I32);
                                y2 = autre.y.to_int_unchecked();
                                y3 = bot.y.to_int_unchecked::<i32>().clamp(0,HAUTEUR_IMAGE_I32);

                                point = Vec2D{x:0.0, y:y1 as f32 + 0.5};

                                unsurv1z = 1.0/triangle_raster.V1.z;
                                unsurv2z = 1.0/triangle_raster.V2.z;
                                unsurv3z = 1.0/triangle_raster.V3.z;

                                unsuraymty = (autre.x - top.x)/(autre.y - top.y);
                                unsurbymty = (bot.x - top.x)/(bot.y - top.y);
                                unsurbymay = (bot.x - autre.x)/(bot.y - autre.y);
                            
                                ymy2 = if autre.y < 0.0 {-autre.y} else {0.0};
                                ymy1 = if top.y < 0.0 {-top.y} else {0.0};

                                col_src = [col.0,col.1,col.2];

                                'loopy:for y in y1..y3 {
                                    //ligne = (y * LARGEUR_IMAGE_I32) as usize;
                                    //posbuf = (ligne + xmin)*3;
                                    if y < y2 {
                                        x1 = ymy1 * unsuraymty + top.x;
                                    }
                                    else {
                                        x1 = ymy2 * unsurbymay  + autre.x;
                                        ymy2 += 1.0;
                                    }
                                    x2 = ymy1 * unsurbymty + top.x;
                                    ymy1 += 1.0;
                                    if x1 > x2 {
                                        std::mem::swap(&mut x1, &mut x2);
                                    }
                                    point.x = x1;
                                    let mut d1 = triangle_raster.dans_triangle_aire_garanti_que_d(&point, &aire, &unsurv1z,&unsurv2z,&unsurv3z);
                                    point.x = x2;
                                    let d2 = triangle_raster.dans_triangle_aire_garanti_que_d(&point, &aire, &unsurv1z,&unsurv2z,&unsurv3z);
                                    deltad = (d2 - d1) / (x2 - x1);
                                    if x1 < 0.0 {
                                        x1 = 0.0;
                                    }
                                    
                                        posbuf = (y as usize * LARGEUR_IMAGE_USIZE + x1.to_int_unchecked::<usize>()) * 3;

                                        if x2 > LARGEUR_IMAGE {
                                            x2 = LARGEUR_IMAGE
                                        }
                                        for x in x1.to_int_unchecked::<i32>()..x2.to_int_unchecked::<i32>() {
                                            zbuffer_a = zbuffer.get_unchecked_mut(posbuf);
                                            if zbuffer_a > &mut d1 {
                                                *zbuffer_a = d1;
                                                framebuffer[posbuf..posbuf+3].copy_from_slice(&col_src);
                                            }
                                            posbuf += 3;
                                            d1 += deltad;
                                        }
                                    }
                                    point.y += 1.0;
                                }
                            }   
                        }
                    }
                    index += 1;
                    
                }
            }
        }
    }
}

fn lux_u8_a_f32((r,g,b):(u8,u8,u8)) -> (f32,f32,f32) {
    (r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
} 