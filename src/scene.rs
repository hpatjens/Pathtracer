use common::*;

use tracer::{Hit, Transition, Camera, ToneMapping};

use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct HDRITexture {
    pub pixels: Vec<f32>,
    pub width: usize,
    pub height: usize,
    last_pixel_index: usize, // Index of the last pixel to clamp the index while accessing the array
}

impl HDRITexture {
    pub fn new(pixels: Vec<f32>, width: usize, height: usize) -> HDRITexture {
        HDRITexture {
            pixels: pixels,
            width: width,
            height: height,
            last_pixel_index: width*height*3 - 3,
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
        let index = usize::min(self.texture.last_pixel_index, usize::max(0, 3*(pix_coord_y*self.texture.width + pix_coord_x)));

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
    HDRI(String, Option<Arc<HDRITexture>>),
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

// @TODO: Distinguish finite and infinite planes
#[derive(Clone, Debug)]
pub struct Plane {
    pub origin: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
    pub normal: Vec3,
}

impl Plane {
    pub fn new(origin: Vec3, u: Vec3, v: Vec3, material: Material) -> Self {
        let normal = u.cross(v).normalize();
        Plane {
            origin: origin,
            u: u,
            v: v,
            material: material,
            normal: normal,
        }
    }

    pub fn sample(&self, position: Vec3) -> Option<Vec3> {
        if (self.origin - position).dot(self.normal) < 0.0 {
            None
        } else {
            let r1 = random32();
            let r2 = random32();

            Some(self.origin + r1*self.u + r2*self.v)
        }
    }
}

#[derive(Clone, Debug, new)]
pub struct Sphere {
    pub origin: Vec3,
    pub radius: f32,
    pub material: Material,
}

impl Sphere {
    pub fn sample(&self, _position: Vec3) -> Option<Vec3> {
        let r1 = random32();
        let r2 = random32();

        let theta = 2.0*PI*r1;
        let phi = f32::acos(2.0*r2 - 1.0);

        Some(Vec3::new(
            f32::sin(theta)*f32::cos(phi),
            f32::cos(theta),
            f32::sin(theta)*f32::sin(phi),
        ))
    }
}

#[derive(Debug, new)]
pub struct Scene {
    pub camera: Camera,
    pub sky: Sky,
    pub spheres: Vec<Sphere>,
    pub planes: Vec<Plane>,
    pub emissive_spheres: Vec<Sphere>,
    pub emissive_planes: Vec<Plane>,
}

impl Default for Scene {
    fn default() -> Self {
        let camera = Camera::new(Vec3::new(0.0, 2.0, 20.0), Vec3::zero(), Vec3::new(0.0, 1.0, 0.0), 4.0, 4.0, 10.0, ToneMapping::Exposure(1.0), 100.0);
        Scene::new(camera, Sky::Constant(Vec3::new(1.0, 1.0, 1.0)), Vec::new(), Vec::new(), Vec::new(), Vec::new())
    }
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
    let n = plane.normal;
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
    let in_plane = position - s;
    let on_u = in_plane.dot(plane.u.normalize());
    if on_u > plane.u.length() || on_u < 0.0 {
        return None;
    }
    let on_v = in_plane.dot(plane.v.normalize());
    if on_v > plane.v.length() || on_v < 0.0 {
        return None;
    }

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

    for sphere in scene.spheres.iter().chain(scene.emissive_spheres.iter()) {
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

    for plane in scene.planes.iter().chain(scene.emissive_planes.iter()) {
        if let Some(hit) = intersect_plane(plane, &ray) {
            if hit.transition == Transition::In {
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
    }

    nearest_hit
}
