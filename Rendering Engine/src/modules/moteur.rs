use super::geo::*;
use super::chunks::*;

use std::cell::{RefCell, RefMut, Ref};

// none of this file is used anywhere, it's WIP

const LIMITE_FROTTEMENTS:f32 = 0.1;
const FROTTEMENTS:f32 = 0.8;

trait Collider {
    fn get_pos(&self) -> Point3D;
    fn set_pos(&mut self,pos:Point3D);
    fn get_speed(&self) -> Point3D;
    fn set_speed(&mut self,speed:Point3D);
    fn get_orient(&self) -> (f32,f32);
    fn set_orient(&mut self,orient:(f32,f32));
    fn test_collide(&self,autre:&Box<dyn Collider>);
    fn tick(&mut self, time:f32);
    fn point_dans(&self, cible:&Point3D);

}

struct AABB {
    max:Point3D,
    min:Point3D,
}

impl AABB {
    fn new(point1:Point3D,point2:Point3D) -> AABB {
        let mut max = Point3D::new(0.0,0.0,0.0);
        let mut min = Point3D::new(0.0, 0.0, 0.0);
        if point1.x > point2.x {
            max.x = point1.x;
            min.x = point2.x;
        }
        else {
            max.x = point2.x;
            min.x = point1.x;
        }
        if point1.y > point2.y {
            max.y = point1.y;
            min.y = point2.y;
        }
        else {
            max.y = point2.y;
            min.y = point1.y;
        }
        if point1.z > point2.z {
            max.z = point1.z;
            min.z = point2.z;
        }
        else {
            max.z = point2.z;
            min.z = point1.z;
        }
        AABB { max, min }
    }
    fn collision_point(&self,cible:&Point3D) -> bool {
        return (cible.x >= self.min.x && cible.x <= self.max.x) && (cible.y >= self.min.y && cible.y <= self.max.y) && (cible.z >= self.min.z && cible.z <= self.max.z)
    }
    fn collision_aabb(&self,cible:&AABB) -> bool {
        return (cible.max.x >= self.min.x && cible.min.x <= self.max.x) && (cible.max.y >= self.min.y && cible.min.y <= self.max.y) && (cible.max.z >= self.min.z && cible.min.z <= self.max.z)
    }
    fn update_avec_spd(&mut self, speed:Point3D) {
        self.max.x += speed.x;
        self.max.y += speed.y;
        self.max.z += speed.z;
        self.min.x += speed.x;
        self.min.y += speed.y;
        self.min.z += speed.z;
    }
}

#[derive(Clone, Copy)]
struct EntityController {
    pos:Point3D,
    taille:f32,
    vitesse:Point3D,
    orientation:(f32,f32),
    equipe:u8,
}

impl EntityController {
    fn new(pos:Point3D, taille:f32, orientation:(f32,f32),equipe:u8) -> EntityController {
        EntityController{pos,taille,vitesse:Point3D::new(0.0, 0.0, 0.0), orientation,equipe}
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
    fn collisions_volume(&self,volume:&mut Chunks,vitesse:&Point3D) -> bool {
        
        if volume.plein_a(((self.pos.x + vitesse.x) as i64, (self.pos.y + vitesse.y) as i64, (self.pos.z + vitesse.z) as i64)) {
            true
        } 
        else {
            false
        }
    }
    fn collisions_comp(&self,volume:&mut Chunks) -> bool {
        let totvitesse = Point3D::new(0.0,0.0,0.0).dist(&self.vitesse);
        if totvitesse <= 0.5 {
            let dupli_vitesse = self.vitesse;
            if self.collisions_volume(volume, &dupli_vitesse) {
                return true
            }
        }
        else {
            let mut coef = 0.0;
            while coef < totvitesse {
                let vitesse_test = Point3D::new(self.vitesse.x * coef, self.vitesse.y * coef, self.vitesse.z * coef);
                if self.collisions_volume(volume, &vitesse_test) {
                    return true
                }
                coef += 0.1;
            }
        }
        false
    }
    fn mouvement(&self,id:usize,volume:&mut Chunks) -> WorldEvent {
        if self.collisions_comp(volume) {
            WorldEvent::Stop(id, self.pos)
        }
        else {
            WorldEvent::Mouv(id)
        }
    }
}


// nearly all of the ECS stuff is copied from https://ianjk.com/ecs-in-rust/

trait ComponentVec {
    fn as_any(&self) -> &dyn std::any::Any;
    fn as_any_mut(&mut self) -> &mut dyn std::any::Any;
    fn push_none(&mut self);
    fn remove_at(&mut self, id:usize);
    /* we'll add more functions here in a moment */
}

impl<T: 'static> ComponentVec for RefCell<Vec<Option<T>>> {
    // Same as before
    fn as_any(&self) -> &dyn std::any::Any {
        self as &dyn std::any::Any
    }

    // Same as before
    fn as_any_mut(&mut self) -> &mut dyn std::any::Any {
        self as &mut dyn std::any::Any
    }

    fn push_none(&mut self) {
        // `&mut self` already guarantees we have
        // exclusive access to self so can use `get_mut` here
        // which avoids any runtime checks.
        self.get_mut().push(None)
    }
    
    fn remove_at(&mut self, id:usize) {
        self.get_mut().remove(id);
    }
}

enum WorldEvent {
    Mouv(usize),//vitesse
    Stop(usize,Point3D),//endroit où l'entité s'arrête
}

struct World {
    // We'll use `entities_count` to assign each Entity a unique ID.
    entities_count: usize,
    controllers: Vec<EntityController>,
    component_vecs: Vec<Box<dyn ComponentVec>>,
    chunks:Chunks,
}

impl World {
    fn new() -> Self {
        Self {
            entities_count: 0,
            component_vecs: Vec::new(),
            controllers:Vec::new(),
            chunks:Chunks::new((0,0,0),None,0,Vec::new())
        }
    }
    
    fn new_entity(&mut self) -> usize {
        let entity_id = self.entities_count;
        for component_vec in self.component_vecs.iter_mut() {
            component_vec.push_none();
        }
        self.controllers.push(EntityController::new(Point3D::new(0.0, 0.0, 0.0), 0.0, (0.0,0.0),  0));
        self.entities_count += 1;
        entity_id
    }
    fn calcul_mouvement(&mut self, events:&mut Vec<WorldEvent>) {
        for i in 0..self.controllers.len() {
            events.push(self.controllers[i].mouvement(i,&mut self.chunks));
        }
    }

    fn applique_events(&mut self, events:&mut Vec<WorldEvent>) {
        for event in events {
            match event {
                WorldEvent::Mouv(id) => {
                    self.controllers[*id].pos.x += self.controllers[*id].vitesse.x;
                    self.controllers[*id].pos.y += self.controllers[*id].vitesse.y;
                    self.controllers[*id].pos.z += self.controllers[*id].vitesse.z;
                    self.controllers[*id].frottements();
                }
                WorldEvent::Stop(id, pos) => {
                    
                }
            }
        }
    }
    
    fn add_component_vec<ComponentType: 'static>(
        &mut self,
    ) {
        let mut new_component_vec: Vec<Option<ComponentType>> =
            Vec::with_capacity(self.entities_count);

        for _ in 0..self.entities_count {
            new_component_vec.push(None);
        }
        // Here we create a `RefCell` before inserting into `component_vecs`
        self.component_vecs
            .push(Box::new(RefCell::new(new_component_vec)));
    }
    fn add_component_to_entity<ComponentType: 'static>(
        &mut self,
        entity: usize,
        component: ComponentType,
    ) {
        for component_vec in self.component_vecs.iter_mut() {
            // The `downcast_mut` type here is changed to `RefCell<Vec<Option<ComponentType>>`
            if let Some(component_vec) = component_vec
                .as_any_mut()
                .downcast_mut::<RefCell<Vec<Option<ComponentType>>>>()
            {
                // add a `get_mut` here. Once again `get_mut` bypasses
                // `RefCell`'s runtime checks if accessing through a `&mut` reference.
                component_vec.get_mut()[entity] = Some(component);
                return;
            }
        }

        let mut new_component_vec: Vec<Option<ComponentType>> =
            Vec::with_capacity(self.entities_count);

        for _ in 0..self.entities_count {
            new_component_vec.push(None);
        }

        new_component_vec[entity] = Some(component);

        // Here we create a `RefCell` before inserting into `component_vecs`
        self.component_vecs
            .push(Box::new(RefCell::new(new_component_vec)));
    }
    fn borrow_component_vec_mut<ComponentType: 'static>(
    &self,
    ) -> Option<RefMut<Vec<Option<ComponentType>>>> {
        for component_vec in self.component_vecs.iter() {
            if let Some(component_vec) = component_vec
                .as_any()
                .downcast_ref::<RefCell<Vec<Option<ComponentType>>>>()
            {
                // Here we use `borrow_mut`. 
                // If this `RefCell` is already borrowed from this will panic.
                return Some(component_vec.borrow_mut());
            }
        }
        None
    }
    fn borrow_component_vec<ComponentType: 'static>(
    &self,
    ) -> Option<Ref<Vec<Option<ComponentType>>>> {
        for component_vec in self.component_vecs.iter() {
            if let Some(component_vec) = component_vec
                .as_any()
                .downcast_ref::<RefCell<Vec<Option<ComponentType>>>>()
            {
                // Here we use `borrow_mut`. 
                // If this `RefCell` is already borrowed from this will panic.
                return Some(component_vec.borrow());
            }
        }
        None
    }
}
