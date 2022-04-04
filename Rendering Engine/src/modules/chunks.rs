// this file contains everything regarding the voxel chunks, and how to interact with them


use std::collections::HashMap;
use std::str::FromStr;
use std::fs::{File, create_dir};
use std::io::prelude::*;
use std::io::BufReader;


use noise::{NoiseFn, Perlin, Seedable};

use super::rendu::{Triangle, Vec3D, RasterTexture,DataRendu};
use super::geo::*;

const CHUNK_SIZE:usize = 10;
pub const CHUNK_SIZE_I64:i64 = CHUNK_SIZE as i64;
const CHUNK_SIZE_POW3:usize = CHUNK_SIZE.pow(3);
const LUX_MIN:(u8,u8,u8) = (50,50,50);

pub struct Chunk {
    volume:Volume,
    // 2 unused fields supposed to store ids for entities and ids meshes in a given chunk
    ents:Vec<usize>, 
    modeles:Vec<usize>,
    pub num_remp:usize,
    num_lux:usize,
}

impl Chunk {
    fn new_gen(decalage:(i64,i64,i64),cote:usize,gen_aleat:&Box<dyn GenParam>,num_text:usize,ents:Vec<usize>,modeles:Vec<usize>) -> Chunk {
        let result_gen = Volume::new_avec_genparam(decalage,cote,gen_aleat,num_text);
        Chunk {volume:result_gen.0, ents, modeles, num_remp:result_gen.1, num_lux:0}
    }
    fn new_vide(cote:usize,num_text:usize,ents:Vec<usize>,modeles:Vec<usize>) -> Chunk {
        Chunk {volume:Volume::new_vide(cote,cote,cote,num_text), ents, modeles, num_remp:0, num_lux:0}
    }
    fn new_remp_avec_text(cote:usize,num_text:usize,ents:Vec<usize>,modeles:Vec<usize>,texture:u16) -> Chunk {
        Chunk {volume:Volume::new_remp_avec_text(cote,cote,cote,num_text,texture), ents, modeles, num_remp:CHUNK_SIZE_POW3, num_lux:0}
    }
}

// trait for types that can serve as generation parameters for a Chunks struct
pub trait GenParam {
    fn get_case_a(&self,coord:(i64,i64,i64)) -> Case;
    fn seed(&self) -> u32;
}

#[derive(Clone, Copy)]
pub struct GenParamTerre {
    pub seed:u32,
    pub perlin:Perlin,
    pub perlin_temp:Perlin,
    pub freq_temp:f64,

    pub max_height:f64,
    pub niveau_mer:f64,
    pub num_text:usize,
    pub frequency:f64,
    pub freq_3d:f64,
    pub profondeur:u8,
}

impl GenParam for GenParamTerre {
    fn seed(&self) -> u32 {
        self.seed
    }
    fn get_case_a(&self,(x,y,z):(i64,i64,i64)) -> Case {
        let mut case = Case::new(false, 1);
        let mut alt = 0.0;
        for i in 0..self.profondeur {
            let pow2 = 2_i32.pow(i as u32) as f64;
            alt += self.perlin.get([x as f64 * (self.frequency * pow2), z as f64 * (self.frequency * pow2)])/pow2;
        }
        alt = (((alt * self.max_height).clamp(-self.max_height, self.max_height) + self.max_height) / 2.0) * ((self.perlin_temp.get([x as f64 * self.freq_temp, z as f64 * self.freq_temp]) + 1.0)/2.0);

        if y < alt as i64 || y <= self.niveau_mer as i64 {
            case.rempli = true;
            //else if (self.perlin.get([x as f64 * self.freq_3d,y as f64 * self.freq_3d, z as f64 * self.freq_3d]) + 1.0) / 2.0 > 0.35 {
            //    case.rempli = true;
            //}
        }
        if case.rempli {
            let rapport_alt = alt/self.max_height;
                case.block_type = (rapport_alt * (self.num_text - 1) as f64) as u16;
        }
        else {
            case.rempli = false;
            case.lux = (255,255,255)
        }
        case
    }

}

#[derive(Clone, Copy)]
pub struct GenParamBase {
    pub seed:u32,
    pub perlin:Perlin,
    pub min_plein1:f32,
    pub min_plein2:f32,
    pub num_text:usize,
}

impl GenParam for GenParamBase {
    fn seed(&self) -> u32 {
        self.seed
    }
    fn get_case_a(&self,(x,y,z):(i64,i64,i64)) -> Case {
        let mut case = Case::new(false, 1);
        case.lux = (255,255,255);
        if ((self.perlin.get([x as f64/50.0, y as f64/50.0, z as f64/50.0]) + 1.0) / 2.0) * ((self.perlin.get([x as f64/25.0, y as f64/25.0, z as f64/25.0]) + 1.0) / 2.0)  > 0.26 {
            case.rempli = true;
            case.block_type = (((self.perlin.get([x as f64/5.0, (y as f64 - 1.0)/5.0, z as f64/5.0]) + 1.0) / 2.0) * (self.num_text - 1) as f64) as u16;
        }
        else {
            case.rempli = false;
        }
        case
    }

}
#[derive(Clone, Copy)]
pub struct BlockType {
    pub texture:u16,

}

pub struct Chunks {
    pub data:HashMap<(i64, i64, i64), Chunk>,
    pub gen_aleat:Option<Box<dyn GenParam>>,
    pub num_text:usize,
    pub block_types:Vec<BlockType>
}

impl Chunks {
    pub fn new(start:(i64, i64, i64), gen_aleat:Option<Box<dyn GenParam>>, num_text:usize,block_types:Vec<BlockType>) -> Chunks {
        let mut map = HashMap::new();
        let mut chunks_out = Chunks {data:map, gen_aleat, num_text,block_types};
        chunks_out.get_at(start);
        chunks_out
        
    }

    pub fn get_float_mut(&mut self, coord:&Point3D) -> &mut Case {
        
        self.get_mut_at((coord.x.round() as i64, coord.y.round() as i64, coord.z.round() as i64))
    }
    pub fn get_mut_at(&mut self, coord:(i64, i64, i64)) -> &mut Case {
        // all functions for getting a block in a chunk guarantee finding the target, it is slower but saved me a lot of time on note having to handle edge cases
        let mut bon = true;
        let (xc, yc, zc) = (coord.0.div_euclid(CHUNK_SIZE_I64), coord.1.div_euclid(CHUNK_SIZE_I64), coord.2.div_euclid(CHUNK_SIZE_I64));
        match self.data.get(&(xc, yc, zc)) {
            Some(chunk) => (),
            None => bon = false,
        }
        if bon {
            self.data.get_mut(&(xc, yc, zc))
            .unwrap()
            .volume.get_mut((coord.0 - xc * CHUNK_SIZE_I64) as usize, (coord.1 - yc * CHUNK_SIZE_I64) as usize, (coord.2 - zc * CHUNK_SIZE_I64) as usize)
        }
        else {
            match &self.gen_aleat {
                Some(param) => {
                    self.data.insert(
                (xc,yc,zc),
                Chunk::new_gen( 
                    (xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64),
                    CHUNK_SIZE,&param,self.num_text,Vec::new(),Vec::new()
                    )
                    );
                },
                None => {
                    self.data.insert(
                    (xc,yc,zc),
                    Chunk::new_remp_avec_text(
                    CHUNK_SIZE,self.num_text,Vec::new(),Vec::new(), 0
                    ));
                }
            }
            self.get_mut_at(coord)
        }
    }
    pub fn get_chunk_at(&mut self, coord:(i64,i64,i64)) -> &mut Chunk {
        let mut bon = true;
        let (xc, yc, zc) = (coord.0.div_euclid(CHUNK_SIZE_I64), coord.1.div_euclid(CHUNK_SIZE_I64), coord.2.div_euclid(CHUNK_SIZE_I64));
        match self.data.get(&(xc, yc, zc)) {
            Some(chunk) => (),
            None => bon = false,
        }
        if bon {
            self.data.get_mut(&(xc, yc, zc))
            .unwrap()
        }
        else {
            match &self.gen_aleat {
                Some(param) => {
                    self.data.insert(
                (xc,yc,zc),
                Chunk::new_gen( 
                    (xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64),
                    CHUNK_SIZE,&param,self.num_text,Vec::new(),Vec::new()
                    )
                    );
                },
                None => {
                    self.data.insert(
                    (xc,yc,zc),
                    Chunk::new_remp_avec_text(
                    CHUNK_SIZE,self.num_text,Vec::new(),Vec::new(), 0
                    ));
                }
            }
            self.get_chunk_at(coord)
        }
    }

    pub fn get_at(&mut self, coord:(i64, i64, i64)) -> &Case {
        let mut bon = true;
        let (xc, yc, zc) = (coord.0.div_euclid(CHUNK_SIZE_I64), coord.1.div_euclid(CHUNK_SIZE_I64), coord.2.div_euclid(CHUNK_SIZE_I64));
        match self.data.get(&(xc, yc, zc)) {
            Some(chunk) => (),
            None => bon = false,
        }
        if bon {
            self.data.get(&(xc, yc, zc))
            .unwrap()
            .volume.get((coord.0 - xc * CHUNK_SIZE_I64) as usize, (coord.1 - yc * CHUNK_SIZE_I64) as usize, (coord.2 - zc * CHUNK_SIZE_I64) as usize)
        }
        else {
            match &self.gen_aleat {
                Some(param) => {
                    self.data.insert(
                (xc,yc,zc),
                Chunk::new_gen( 
                    (xc * CHUNK_SIZE_I64,yc * CHUNK_SIZE_I64,zc * CHUNK_SIZE_I64),
                    CHUNK_SIZE,&param,self.num_text,Vec::new(),Vec::new()
                    )
                    );
                },
                None => {
                    self.data.insert(
                    (xc,yc,zc),
                    Chunk::new_remp_avec_text(
                    CHUNK_SIZE,self.num_text,Vec::new(),Vec::new(), 0
                    ));
                }
            }
            self.get_at(coord)
        }
    }
    pub fn get_triangles_chunk(&mut self, triangles: &mut Vec<Triangle>,coord:(i64, i64, i64)) {
        self.get_at((coord.0*10, coord.1*10, coord.2*10));
        let chunk = self.data.get(&coord).unwrap();
        if chunk.num_remp != 0 && chunk.num_remp != CHUNK_SIZE_POW3 {
            for x in coord.0*10..coord.0*10 + CHUNK_SIZE_I64 {
                for y in coord.1*10..coord.1*10 + CHUNK_SIZE_I64 {
                    for z in coord.2*10..coord.2*10 + CHUNK_SIZE_I64 {
                        if self.plein_a((x,y,z)) {
                            self.get_triangles_at((x,y,z), triangles)
                        }
                    }
                }
            }
        }
        else if chunk.num_remp != 0 {
            let (checkx1, checkx2, checky1, checky2, checkz1, checkz2) = (coord.0*10, coord.0*10 + CHUNK_SIZE_I64 - 1,coord.1*10, coord.1*10 + CHUNK_SIZE_I64 - 1,coord.2*10, coord.2*10 + CHUNK_SIZE_I64 - 1);
            for x in coord.0*10..coord.0*10 + CHUNK_SIZE_I64 {
                for y in coord.1*10..coord.1*10 + CHUNK_SIZE_I64 {
                    for z in coord.2*10..coord.2*10 + CHUNK_SIZE_I64 {
                        if x == checkx1 || x == checkx2 || y == checky1 || y == checky2 || z == checkz1 || z == checkz2 {
                            if self.plein_a((x,y,z)) {
                                self.get_triangles_at((x,y,z), triangles)
                            }
                        }
                    }
                }
            }
        }
    }
    pub fn get_triangle_at_nocheck(&mut self,(x,y,z):(i64,i64,i64), triangles:&mut Vec<Triangle>) {
        let (xf, yf, zf) = (x as f32, y as f32, z as f32);
        let textures = self.texture_a((x,y,z));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y-1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                textures,
                self.lumiere_a((x, y-1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y+1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y+1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x+1, y,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                textures,
                self.lumiere_a((x+1, y,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x-1, y,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x-1, y,z))
            ));
        
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y,z+1))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y,z+1))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y,z-1))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (1.0,0.0)),
                textures,
                self.lumiere_a((x, y,z-1))
            ));
    }
    pub fn get_triangles_at(&mut self,(x,y,z):(i64,i64,i64), triangles:&mut Vec<Triangle>) {
        let (xf, yf, zf) = (x as f32, y as f32, z as f32);
        let textures = self.texture_a((x,y,z));
        if !self.plein_a((x, y-1,z)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y-1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                textures,
                self.lumiere_a((x, y-1,z))
            ));
        }
        if !self.plein_a((x, y+1,z)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y+1,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y+1,z))
            ));
        }
        if !self.plein_a((x+1, y,z)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x+1, y,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                textures,
                self.lumiere_a((x+1, y,z))
            ));
        }
        if !self.plein_a((x-1, y,z)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x-1, y,z))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x-1, y,z))
            ));
        }
        if !self.plein_a((x,y,z+1)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (0.0,1.0)),
                textures,
                self.lumiere_a((x, y,z+1))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (1.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 + 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 + 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y,z+1))
            ));
        }
        if !self.plein_a((x,y,z-1)) {
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (0.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                textures,
                self.lumiere_a((x, y,z-1))
            ));
            triangles.push(Triangle::new(
                Vec3D::new(xf * 10.0 - 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (0.0,0.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 - 5.0, zf * 10.0 - 5.0, (1.0,1.0)),
                Vec3D::new(xf * 10.0 + 5.0, yf * 10.0 + 5.0, zf * 10.0 - 5.0, (1.0,0.0)),
                textures,
                self.lumiere_a((x, y,z-1))
            ));
        }
    }
    pub fn lumiere_a(&mut self, coord:(i64, i64, i64)) -> (u8,u8,u8) {
        self.get_at(coord).lux
    }
    pub fn plein_a(&mut self, coord:(i64, i64, i64)) -> bool {
        self.get_at(coord).rempli
    }
    pub fn texture_a(&mut self, coord:(i64, i64, i64)) -> u16 {
        let bloc_type = self.get_at(coord).block_type as usize;
        self.block_types[bloc_type].texture
    }
    pub fn enleve_lumiere(&mut self, co_lux:(i64,i64,i64),puissance:usize) {
        let mut etapes = Vec::<EtapeLumiere>::new();
        let decalages:[(i64,i64,i64) ; 6] = [(1,0,0),(-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)];
        etapes.push(EtapeLumiere{etape:0, coord:co_lux});
        //self.get_mut_at((co_lux.0, co_lux.1, co_lux.2)).source_lux = true;
        let mut etapes_a_explorer:usize = 1;
        for i in 0..puissance {
            let mut a_rajouter = Vec::<EtapeLumiere>::new();
            for j in (etapes.len() - etapes_a_explorer)..etapes.len()  {
                let etape = &etapes[j];
                let (x,y,z) = etape.coord;
                for decalage in decalages.iter() {
                    if !self.plein_a(((x + decalage.0), (y + decalage.1), (z + decalage.2))) {
                    let mut bon = true;
                    
                    for k in (etapes.len() - etapes_a_explorer)..etapes.len() {
                        let etape_2 = &etapes[k];
                        if  etape_2.meme_coord((x + decalage.0, y + decalage.1, z + decalage.2)) {
                            bon = false;
                            break;
                        }
                    }
                    if bon {
                        a_rajouter.push(EtapeLumiere { etape: i, coord:(x + decalage.0, y + decalage.1, z + decalage.2) });
                    }
                    }
                }
            }
            etapes_a_explorer = 0;
            for rajout in &a_rajouter {
                let mut bon = true;
                for etape in &etapes {
                    if etape.meme_coord(rajout.coord) {
                        bon = false;
                        break;
                    }
                }
                if bon {
                    etapes_a_explorer += 1;
                    etapes.push(*rajout);
                }
            }
        }
        for etape in &etapes {
            let mut change = false;
            {
                let mut case = self.get_mut_at((etape.coord.0, etape.coord.1, etape.coord.2));
                if case.lux != LUX_MIN {
                    change = true;
                }
                case.lux = LUX_MIN;
            }
            if change {
                self.get_chunk_at(etape.coord).num_lux -= 1;
            }
        }
    }
    pub fn rajoute_lumiere(&mut self, co_lux:(i64,i64,i64),puissance:usize,col:(u8,u8,u8)) {
        //this algorithm to add lights start from a center point, then spreads out with increasingly less strength, recording every valid block along the way.
        // once all valid blocks are available, their respective light levels are updated
        // this becomes incredibly slow as puissance gets bigger, but i haven't spent the time coming up with or researching a better solution
        let mut etapes = Vec::<EtapeLumiere>::new();
        let decalages:[(i64,i64,i64) ; 6] = [(1,0,0),(-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)];
        etapes.push(EtapeLumiere{etape:0, coord:co_lux});
        let mut etapes_a_explorer:usize = 1;
        for i in 0..puissance {
            let mut a_rajouter = Vec::<EtapeLumiere>::new();
            for j in (etapes.len() - etapes_a_explorer)..etapes.len()  {
                let etape = &etapes[j];
                let (x,y,z) = etape.coord;
                for decalage in decalages.iter() {
                    if !self.plein_a(((x + decalage.0), (y + decalage.1), (z + decalage.2))) {
                    let mut bon = true;
                    
                    for k in (etapes.len() - etapes_a_explorer)..etapes.len() {
                        let etape_2 = &etapes[k];
                        if  etape_2.meme_coord((x + decalage.0, y + decalage.1, z + decalage.2)) {
                            bon = false;
                            break;
                        }
                    }
                    if bon {
                        a_rajouter.push(EtapeLumiere { etape: i, coord:(x + decalage.0, y + decalage.1, z + decalage.2) });
                    }
                    }
                }
            }
            etapes_a_explorer = 0;
            for rajout in &a_rajouter {
                let mut bon = true;
                for etape in &etapes {
                    if etape.meme_coord(rajout.coord) {
                        bon = false;
                        break;
                    }
                }
                if bon {
                    etapes_a_explorer += 1;
                    etapes.push(*rajout);
                }
            }
        }
        for etape in &etapes {
            let mut change = false;
            {
                let mut case = self.get_mut_at((etape.coord.0, etape.coord.1, etape.coord.2));
                if case.lux.0 < (((puissance - etape.etape) as f32 / puissance as f32) * col.0 as f32) as u8 {
                    change = true;
                    case.lux.0 = (((puissance - etape.etape) as f32 / puissance as f32) * col.0 as f32) as u8;
                }
                if case.lux.1 < (((puissance - etape.etape) as f32 / puissance as f32) * col.1 as f32) as u8 {
                    change = true;
                    case.lux.1 = (((puissance - etape.etape) as f32 / puissance as f32) * col.1 as f32) as u8;
                }
                if case.lux.2 < (((puissance - etape.etape) as f32 / puissance as f32) * col.2 as f32) as u8{
                    change = true;
                    case.lux.2 = (((puissance - etape.etape) as f32 / puissance as f32) * col.2 as f32) as u8;
                }
            }
            if change {
                self.get_chunk_at(etape.coord).num_lux += 1;
            }
        }
    }
    pub fn ray_traverse(&mut self, ray:&Rayon,lim:f32) -> Option<f32> {
        let mut coef = 0.0;
        let mut pointtest:Point3D = ray.at(coef);
        while coef < lim {
            if self.plein_a((pointtest.x.round() as i64, pointtest.y.round() as i64,pointtest.z.round() as i64)) {
                return Some(coef)
            }
            else {
                coef += 0.1;
                pointtest = ray.at(coef);
            }
        }
        None
    }
    pub fn sauvegarde_fichier(&self, nom:&str,nom_comp:&str, textures:&Vec<RasterTexture>) {
        match create_dir(format!("chunks/{}",nom)) {
            Ok(()) => {
                match File::create(format!("chunks/{}/{}_main.txt", nom, nom)) {
                    Ok(fic) => {
                        let mut scripteur = fic;
                        let mut data_fichier = String::new();
                        data_fichier.push_str(format!("{}\n", nom_comp).trim_matches('\t'));
                        match &self.gen_aleat {
                            Some(opt) => {
                                data_fichier.push_str(format!("{} {} {} {} {} {} {}\n", CHUNK_SIZE, LUX_MIN.0,LUX_MIN.1,LUX_MIN.2, self.num_text, self.data.len(), opt.seed()).trim_matches('\t'));
                            }
                            None => {
                                data_fichier.push_str(format!("{} {} {} {} {} {}\n", CHUNK_SIZE,  LUX_MIN.0,LUX_MIN.1,LUX_MIN.2, self.num_text, self.data.len()).trim_matches('\t'));
                            }
                        }

                        for i in 0..self.num_text {
                            let mut ligne_text = String::new();
                            let texture = &textures[i];
                            ligne_text.push_str(format!("{}",texture.deltatick).trim());
                            for nom in &texture.noms {
                                ligne_text.push_str(format!(";{}", nom.trim()).trim())
                            }
                            data_fichier.push_str(format!("{}\n",ligne_text.trim()).trim_matches('\t'));
                        }
                        let mut total_lumiere = 0;
                        for (coord, chunk) in &self.data {
                            if (chunk.num_remp != 0 || chunk.num_lux != 0) && chunk.num_remp != CHUNK_SIZE_POW3 {
                                data_fichier.push_str(chunk.volume.genere_lignes_chunk(coord).trim_matches('\t'));
                            }
                            else if chunk.num_remp == CHUNK_SIZE_POW3 {
                                let mut plein_1 = true;
                                let texture_test = self.block_types[chunk.volume.data[0].block_type as usize].texture;
                                for case in &chunk.volume.data {
                                    if  self.block_types[case.block_type as usize].texture != texture_test {
                                        plein_1 = false;
                                        break;
                                    }
                                }
                                if plein_1 {
                                    data_fichier.push_str(format!("CHNK {} {} {} P {} {} {} {} {} {}\n", coord.0, coord.1, coord.2, texture_test, texture_test, texture_test, texture_test, texture_test, texture_test).trim_matches('\t'));
                                }
                                else {
                                    data_fichier.push_str(chunk.volume.genere_lignes_chunk(coord).trim_matches('\t'));
                                }
                            }
                            else {
                                data_fichier.push_str(format!("CHNK {} {} {} V\n", coord.0, coord.1, coord.2).trim_matches('\t'));
                            }
                            total_lumiere += chunk.num_lux;
                        }
                        write!(&mut scripteur, "{}", data_fichier.trim_matches('\t')).unwrap();
                        println!("fichier de chunks sauvegardé ! nombre de lumière {}", total_lumiere)
                    },
                    Err(_err) => println!("FQFQSF"),
                
                }
            },
            Err(_err) => println!("nom invalide")
        }
    }
    pub fn charge_fichier(nom:&str) -> Result<Chunks,()> {
        // i still need to write down the spec for the chunk file format
            match File::open(format!("chunks/{0}/{0}_main.txt",nom)) {
                Ok(fic) => {
                    let mut data_rendu = DataRendu::new();
                    let mut chunks_sortie = Chunks::new((0,0,0), None, 0,Vec::new());
                    let mut liste_lignes = Vec::<String>::new();
                    let reader = BufReader::new(fic);
                    let lignes = reader.lines();
                    for ligne in lignes {
                        let ligne_str = ligne.unwrap();
                        liste_lignes.push(ligne_str);
                    }
                    let nom = liste_lignes[0].clone();
                    if liste_lignes[1].len() >= 6 {
                        let ligne_info: Vec<&str> = liste_lignes[1].split_whitespace().collect();
                        let taille_chunk = usize::from_str(ligne_info[0]).unwrap();
                        let min_lux_t = (f32::from_str(ligne_info[1]).unwrap() as u8, f32::from_str(ligne_info[2]).unwrap() as u8,f32::from_str(ligne_info[3]).unwrap() as u8);
                        let num_text = usize::from_str(ligne_info[4]).unwrap();
                        for i in 2..2+num_text {
                            let mut ligne:Vec<&str> = liste_lignes[i].split(";").collect();
                            let deltatick = usize::from_str(ligne[0]).unwrap();
                            ligne.remove(0);
                            data_rendu.charge_textures(ligne,  deltatick);
                        }
                        let mut dans_chunk:usize = 0;
                        let mut num_remp:usize = 0;
                        let mut coord:(i64,i64,i64) = (0,0,0);
                        let mut data = Vec::new();
                        let mut num_lux_tot = 0;
                        let mut num_chunk_charg = 0;
                        for index in 2+num_text..liste_lignes.len() {
                            if dans_chunk == 0 {
                                let ligne_coord: Vec<&str> = liste_lignes[index].split_whitespace().collect();
                                if ligne_coord.len() > 0 && ligne_coord[0] == "CHNK" {
                                    coord = (i64::from_str(ligne_coord[1]).unwrap(),i64::from_str(ligne_coord[2]).unwrap(),i64::from_str(ligne_coord[3]).unwrap());
                                    data = Vec::new();
                                    num_remp = 0;
                                    if ligne_coord.len() <= 4 {
                                        dans_chunk += 1;
                                    }
                                }
                                else if ligne_coord.len() > 4 && ligne_coord[4] == "V" {
                                    coord = (i64::from_str(ligne_coord[1]).unwrap(),i64::from_str(ligne_coord[2]).unwrap(),i64::from_str(ligne_coord[3]).unwrap());
                                    chunks_sortie.data.insert(coord, Chunk::new_vide(taille_chunk, num_text, Vec::new(), Vec::new()));
                                }
                                else if ligne_coord.len() > 4 && ligne_coord[4] == "P" {
                                    coord = (i64::from_str(ligne_coord[1]).unwrap(),i64::from_str(ligne_coord[2]).unwrap(),i64::from_str(ligne_coord[3]).unwrap());
                                    chunks_sortie.data.insert(coord, Chunk::new_remp_avec_text(taille_chunk, num_text, Vec::new(), Vec::new(),
                                        u16::from_str(ligne_coord[5]).unwrap(),
                                    ));
                                }
                            }
                            else if dans_chunk <= taille_chunk.pow(3) {
                                let case = Case::lit_ligne_fic(liste_lignes[index].trim(), min_lux_t);
                                if case.rempli {
                                        num_remp += 1;
                                }
                                if case.lux != min_lux_t {
                                    num_lux_tot += 1;
                                }
                                data.push(case);
                                dans_chunk += 1;
                            }
                            if dans_chunk == taille_chunk.pow(3) + 1 {
                                let mut chunk = Chunk::new_vide(taille_chunk, num_text, Vec::new(), Vec::new());
                                chunk.num_remp = num_remp;
                                let mut data_copy = Vec::<Case>::with_capacity(taille_chunk.pow(3));
                                for case in &data { data_copy.push(case.clone())}

                                chunk.volume.data = data_copy;
                                chunks_sortie.data.insert((coord.0,coord.1,coord.2), chunk);
                                dans_chunk = 0;
                                num_chunk_charg += 1;
                            }
                            
                        }
                        let gen_aleat:Option<Box<dyn GenParam>>;
                        let mut perlin = Perlin::new();
                        perlin.set_seed(u32::from_str(ligne_info[6]).unwrap());
                        if ligne_info.len() == 7 {
                            gen_aleat = Some(Box::new(GenParamBase{seed:u32::from_str(ligne_info[6]).unwrap(), min_plein1:0.0, min_plein2:0.0,perlin,num_text}));
                        }
                        else {
                            gen_aleat = None;
                        }
                        chunks_sortie.gen_aleat = gen_aleat;
                        chunks_sortie.num_text = num_text;
                        println!("num lux charg: {}", num_lux_tot);
                        println!("num chunk charg: {}", num_chunk_charg);
                        println!("num lignes = {}", liste_lignes.len());
                        return Ok(chunks_sortie)
                    }
                },
                Err(err) => { println!("Fichier inexistant ! {}", err); return Err(()) }
            }
            Err(())
    }
}

#[derive(Clone)]
pub struct Case {
    pub rempli:bool,
    pub lux:(u8,u8,u8),
    pub block_type:u16,
    //interaction:Option<Box<dyn InteractionBloc>> 
    // this last commented field isn't used anymore, and its functionality will be replaced by block types
}

impl Case {
    fn new(rempli:bool, block_type:u16) -> Case {
        Case { rempli, lux:LUX_MIN, block_type}
    }
    fn genere_ligne_fic(&self, min_lux:(u8,u8,u8)) -> String {
        let mut ligne = String::new();
        if self.lux != min_lux {
            ligne.push_str(format!("L {} {} {} {};", self.lux.0, self.lux.1, self.lux.2, "N").trim());
        }
        if self.rempli {
            ligne.push_str(format!("T {};",self.block_type).trim());
        }
        ligne
    }
    fn lit_ligne_fic(ligne:&str, min_lux:(u8,u8,u8)) -> Case {
        let composants:Vec<&str> = ligne.split(";").collect();
        let mut block_type = 0;
        let mut lux = (min_lux, false);
        let mut rempli = false;
        for composant in composants {
            let details:Vec<&str> = composant.split_whitespace().collect();
            
            if details.len() > 0 {
                match details[0] {
                    "T" => {
                        block_type = u16::from_str(details[1]).unwrap();
                        rempli = true
                    },
                    "L" => {
                        lux = (((f32::from_str(details[1]).unwrap() * 255.0) as u8,(f32::from_str(details[2]).unwrap() * 255.0) as u8,(f32::from_str(details[3]).unwrap() * 255.0) as u8), if details[4] == "S" {true} else {false})
                    }
                    _ => (),
                }
            }
        }
        Case {block_type, rempli, lux:lux.0}//, interaction}
    }
}

struct Volume {
    longueur:usize,
    largeur:usize,
    hauteur:usize,
    data:Vec<Case>,
    min_lux:(u8,u8,u8),
    num_text:usize,
}

impl Volume {
    fn new(longueur:usize, largeur:usize, hauteur:usize, num_text:usize) -> Volume {
        Volume { longueur, largeur, hauteur, data: vec![Case::new(true, 0); longueur * largeur * hauteur], num_text, min_lux:LUX_MIN}
    }
    fn new_vide(longueur:usize, largeur:usize, hauteur:usize, num_text:usize) -> Volume {
        Volume { longueur, largeur, hauteur, data: vec![Case::new(false, 0); longueur * largeur * hauteur], num_text, min_lux:LUX_MIN}
    }
    fn new_remp_avec_text(longueur:usize, largeur:usize, hauteur:usize, num_text:usize, texture:u16) -> Volume {
        Volume { longueur, largeur, hauteur, data: vec![Case::new(true, texture); longueur * largeur * hauteur], num_text, min_lux:LUX_MIN}
    }
    fn genere_lignes_chunk(&self,coord:&(i64,i64,i64)) -> String {
        let mut data_fichier = String::new();
        data_fichier.push_str(format!("CHNK {} {} {}\n", coord.0,coord.1,coord.2).trim_matches('\t'));
        //data_fichier.push_str(format!("{} {} {} {}\n", self.hauteur, self.longueur, self.largeur, self.min_lux).trim_matches('\t'));//, self.num_text).unwrap();

        for case in &self.data {
            data_fichier.push_str(format!("{}\n",case.genere_ligne_fic(self.min_lux).trim()).trim_matches('\t'));
        }
        data_fichier
    }
    fn new_avec_genparam(decalage:(i64,i64,i64),cote:usize,genparam:&Box<dyn GenParam>,num_text:usize) -> (Volume,usize) {
        let (longueur, largeur, hauteur) = (cote, cote, cote);
        let mut volume = Volume::new(longueur, largeur, hauteur, num_text);
        let mut num_remp = 0;
        for x in decalage.0..decalage.0+longueur as i64 {
            for y in decalage.1..decalage.1+hauteur as i64 {
                for z in decalage.2..decalage.2+largeur as i64 {
                    let (xg, yg, zg) = ((x - decalage.0) as usize, (y - decalage.1) as usize, (z - decalage.2) as usize);
                    let case = genparam.get_case_a((x,y,z));
                    if case.rempli {
                        num_remp += 1;
                    }
                    *volume.get_mut(xg, yg, zg) = case;
                }
            }
        }
        (volume,num_remp)
    }
    fn get(&self, x:usize, y:usize, z:usize) -> &Case {
        &self.data[(y * self.longueur * self.largeur) + z * (self.largeur) + x]
    }
    fn get_mut(&mut self, x:usize, y:usize, z:usize) -> &mut Case {
        &mut self.data[(y * self.longueur * self.largeur) + z * (self.largeur) + x]
    }

}

#[derive(Clone, Copy)]
struct EtapeLumiere {
    etape:usize,
    coord:(i64, i64, i64)
}

impl EtapeLumiere {
    fn meme_coord(&self, (x,y,z):(i64, i64, i64)) -> bool {
        if self.coord.0 == x && self.coord.1 == y && self.coord.2 == z {
            true
        }
        else {
            false
        }
    }
}