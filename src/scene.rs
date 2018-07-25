use common::*;

use tracer::{Hit, Transition};

#[derive(Clone, Debug)]
pub struct HDRITexture {
    pub pixels: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

impl HDRITexture {
    pub fn new(pixels: Vec<f32>, width: usize, height: usize) -> HDRITexture {
        HDRITexture {
            pixels: pixels,
            width: width,
            height: height,
        }
    }

    pub fn sampler<'a>(&'a self) -> HDRITextureSampler<'a> {
        HDRITextureSampler {
            texture: self,
        }
    }
}

fn repeat_f32(x: f32) -> f32 {
    if x < 0.0 {
        x - (x as i32) as f32 + 1.0
    } else {
        x - (x as i32) as f32
    }
}

pub struct HDRITextureSampler<'a> {
    texture: &'a HDRITexture,
}

impl<'a>HDRITextureSampler<'a> {
    pub fn sample(&self, tex_coord: Vec2) -> Vec3 {
        let tex_coord = Vec2::new(
            repeat_f32(tex_coord.x),
            repeat_f32(tex_coord.y),
        );

        let pix_coord_x = (tex_coord.x*self.texture.width  as f32) as usize;
        let pix_coord_y = (tex_coord.y*self.texture.height as f32) as usize;
        let index = 3*(pix_coord_y*self.texture.width + pix_coord_x);

        Vec3::new(
            self.texture.pixels[index + 0],
            self.texture.pixels[index + 1],
            self.texture.pixels[index + 2],
        )
    }

    pub fn sample_equirectangular(&self, direction: Vec3) -> Vec3 {
        let uv = Vec2::new(
            f32::atan2(direction.z, direction.x),
            f32::acos(direction.y),
        ) / Vec2::new(2.0*PI, PI);
        self.sample(uv)
    }
}

#[derive(Clone, Debug)]
pub enum Sky {
    Constant(Vec3),
    HDRI(String, Option<HDRITexture>),
}

#[derive(Clone, Debug)]
pub struct PBRParameters {
    pub reflectivity: Vec3,
    pub roughness: f32,
    pub metalness: f32,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
pub enum Material {
    None,
    Emissive(Vec3),
    Mirror,
    Translucent(f32),
    Physically(PBRParameters),
}

// @TODO: Don't use this plane in the camera
#[derive(Clone, Debug, new)]
pub struct Plane {
    pub origin: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
}

#[derive(Clone, Debug, new)]
pub struct Sphere {
    pub origin: Vec3,
    pub radius: f32,
    pub material: Material,
}

#[derive(Debug, new)]
pub struct Scene {
    pub sky: Sky,
    pub spheres: Vec<Sphere>,
    pub planes: Vec<Plane>,
}

fn intersect_sphere<'a>(sphere: &'a Sphere, ray: &Ray) -> Option<Hit<'a>> {
    let to_center = sphere.origin - ray.origin;
    let projection = ray.direction.dot(to_center);
    if projection < 0.0 {
        return None;
    }

    let on_ray_to_center = projection*ray.direction;
    let to_inner_hit = to_center - on_ray_to_center;
    let inner_hit_distance = to_inner_hit.length();
    if inner_hit_distance > sphere.radius {
        return None;
    }

    let on_ray_in_sphere = f32::sqrt(sphere.radius*sphere.radius - inner_hit_distance*inner_hit_distance);
    let t1 = projection - on_ray_in_sphere;
    let t2 = projection + on_ray_in_sphere;

    let mut parameter = -2.0;
    if t1 > parameter && t1 > 0.0 { parameter = t1; }
    if t2 < parameter && t2 > 0.0 { parameter = t2; }
    if parameter < -1.0 {
        return None;
    }

    let position = ray.origin + parameter*ray.direction;
    let material = &sphere.material;
    let mut normal = (position - sphere.origin).normalize();
    let mut transition = Transition::In;
    if ray.direction.dot(normal) > 0.0 { 
        normal = -normal;
        transition = Transition::Out;
    }
    Some(Hit::new(parameter, position, normal, material, transition))
}

fn intersect_plane<'a>(plane: &'a Plane, ray: &Ray) -> Option<Hit<'a>> {
    let n = plane.u.cross(plane.v).normalize();
    let s = plane.origin;
    
    let p = ray.origin;
    let d = ray.direction;

    let num = n.x*(s.x - p.x) + n.y*(s.y - p.y) + n.z*(s.z - p.z);
    let denum = n.x*d.x + n.y*d.y + n.z*d.z;

    if f32::abs(denum) < 0.0000001 { // @TODO: Set a reasonable epsilon
        return None;
    }

    let parameter = num / denum;
    if parameter < 0.0 {
        return None;
    }

    let position = p + parameter*d;
    let material = &plane.material;
    let (transition, normal) = if ray.direction.dot(n) < 0.0 { 
        (Transition::In, n)
    } else { 
        (Transition::Out, -n)
    };


    Some(Hit::new(parameter, position, normal, material, transition))
}

pub fn find_scene_hit<'a>(ray: &Ray, scene: &'a Scene) -> Option<Hit<'a>> {
    let mut nearest_hit: Option<Hit> = None;

    for sphere in &scene.spheres {
        if let Some(hit) = intersect_sphere(sphere, &ray) {
            nearest_hit = if let Some(nearest_hit) = nearest_hit {
                if hit.parameter < nearest_hit.parameter {
                    Some(hit)
                } else {
                    Some(nearest_hit)
                }
            } else {
                Some(hit)
            }
        }
    }

    for plane in &scene.planes {
        if let Some(hit) = intersect_plane(plane, &ray) {
            nearest_hit = if let Some(nearest_hit) = nearest_hit {
                if hit.parameter < nearest_hit.parameter {
                    Some(hit)
                } else {
                    Some(nearest_hit)
                }
            } else {
                Some(hit)
            }
        }
    }

    nearest_hit
}
