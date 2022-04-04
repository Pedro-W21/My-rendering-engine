use super::chunks::CHUNK_SIZE_I64;

// this file is supposed to contain everything regarding geometry

#[derive(Clone, Copy)]
pub struct Point3D {
    pub x:f32,
    pub y:f32,
    pub z:f32,
}

impl Point3D {
    pub fn new(x:f32, y:f32, z:f32) -> Point3D {
        Point3D {x,y,z}
    }
    pub fn get_co_int(&self) -> (i64, i64, i64) {
        (self.x.round() as i64, self.y.round() as i64, self.z.round() as i64)
    }
    pub fn dist(&self, autre:&Point3D) -> f32 {
        ((autre.x - self.x).powi(2) + (autre.y - self.y).powi(2) + (autre.z - self.z).powi(2)).sqrt()
    }
    pub fn co_chunk(&self) -> (i64,i64,i64) {
        (((self.x * 100.0) as i64).div_euclid(CHUNK_SIZE_I64 * 10), ((self.y * 100.0) as i64).div_euclid(CHUNK_SIZE_I64 * 10), ((self.z * 100.0) as i64).div_euclid(CHUNK_SIZE_I64 * 10))
    }
    pub fn dist_signe(&self, taille1:f32, taille2:f32, pos2:&Point3D) -> f32 {
        ((pos2.x - self.x).powi(2) + (pos2.y - self.y).powi(2) + (pos2.z - self.z).powi(2)).sqrt() - (taille1 + taille2)
    }
    pub fn get_orient_vers(&self, cible:&Point3D) -> (f32,f32) {
        let dist_horiz = ((cible.x - self.x).powi(2) + (cible.z - self.z).powi(2)).sqrt();
        ((cible.x - self.x).atan2(cible.z - self.z), (cible.y - self.y).atan2(dist_horiz))
    }
    pub fn normalise(&self) -> Point3D {
        let inv_dist = 1.0/Point3D::new(0.0, 0.0, 0.0).dist(&self);
        Point3D::new(
            self.x * inv_dist,
            self.y * inv_dist, 
            self.z * inv_dist
        )
    }
    pub fn get_vec_orient(angh:f32,angv:f32) -> Point3D {
        Point3D{x:angh.sin() * angv.cos(), y:angv.sin(), z:angh.cos() * angv.cos()}
    }
    pub fn produit_scalaire(&self, autre:&Point3D) -> f32 {
        self.x * autre.x + self.y * autre.y + self.z * autre.z
    }
}
pub struct Rayon {
    pub orig:Point3D,
    pub dir:Point3D,
}
impl Rayon {
    pub fn at(&self, coef:f32) -> Point3D {
        Point3D::new(self.orig.x + self.dir.x * coef, self.orig.y + self.dir.y * coef, self.orig.z + self.dir.z * coef)
    }
    pub fn new_orient(orig:Point3D, angh:f32, angv:f32) -> Rayon {
        Rayon {orig, dir:Point3D{x:angh.sin() * angv.cos(), y:angv.sin(), z:angh.cos() * angv.cos()}}
    }
    pub fn new_orient_vers(orig:Point3D, cible:Point3D) -> Rayon {
        let dist_horiz = ((cible.x - orig.x).powi(2) + (cible.z - orig.z).powi(2)).sqrt();
        let angh = (cible.x - orig.x).atan2(cible.z - orig.z);
        let angv = (cible.y - orig.y).atan2(dist_horiz);
        Rayon {orig, dir:Point3D{x:angh.sin() * angv.cos(), y:angv.sin(), z:angh.cos() * angv.cos()}}
    }
}