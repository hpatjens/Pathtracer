use common::*;

use tracer::Hit;

#[derive(Clone, Debug)]
#[allow(dead_code)]
pub enum Material {
    None,
    Emissive(Vec3),
    Mirror,
    Glass,
    Phyiscally {
        reflectivity: Vec3,
        roughness: f32,
        metalness: f32,
    },
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
    let normal = (position - sphere.origin).normalize();
    let material = &sphere.material;
    Some(Hit::new(parameter, position, normal, material))
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
    let normal = n;
    let material = &plane.material;

    Some(Hit::new(parameter, position, normal, material))
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
