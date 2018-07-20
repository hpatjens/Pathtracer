use common::*;

use scene::{Scene, Plane, Material, find_scene_hit};

use std::sync::{Arc, RwLock};
use std::cell::UnsafeCell;

#[derive(Clone, Debug, new)]
pub struct WorkTile {
    pub tile_index: Vec2u,
    pub position: Vec2u,
    pub size: Vec2u,
}

#[derive(Debug, Clone, new)]
pub struct Hit<'a> {
    pub parameter: f32,
    pub position: Vec3,
    pub normal: Vec3,
    pub material: &'a Material,
}

#[derive(Debug, new)]
pub struct Camera {
    projection_plane: Plane,
    eye: Vec3,
}

pub fn look_at(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32) -> Camera {
    let projection_plane = {
        let z = (position - target).normalize();
        let x = up.cross(z).normalize();
        let y = z.cross(x).normalize();

        let half_width = width / 2.0;
        let half_height = height / 2.0;

        let to_plane_center = -z_near*z;

        let origin = position + to_plane_center - x*half_width - y*half_height;

        let u = x*width;
        let v = y*height;
        Plane::new(origin, u, v, Material::None)
    };
    Camera::new(projection_plane, position)
}

pub struct Backbuffer {
    pub width: u32,
    pub height: u32,
    // The worker threads are all processing distinct tiles within the Vec which enables
    // synchronization-free writes. At the same time the main thread should not have to
    // wait on mutexes to copy the data from the backbuffer to the window as this would
    // most likely result in stuttering. Therefore the pixel data is stored in an UnsafeCell
    // and the Backbuffer struct is passed around immutably.
    //
    // A different solution to this problem would be to split the Vec with chunks_mut.
    // This however would add a lot of bookkeeping to pass the right references to the 
    // workers that process the respective tile. Consider that one worker needs multiple
    // references to the individual lines of the tile!
    pub pixels: UnsafeCell<Vec<Pixel>>,
}
// UnsafeCell does not implement Sync and therefore Backbuffer could not be passed to the
// worker threads without this implementation.
unsafe impl Sync for Backbuffer {}

impl Backbuffer {
    pub fn new(width: u32, height: u32) -> Backbuffer {
        Backbuffer {
            width: width,
            height: height,
            pixels: {
                let mut pixels = Vec::new();
                let num_pixels = (width * height) as usize;
                pixels.resize(num_pixels, Pixel(0, 0, 0));
                UnsafeCell::new(pixels)
            },
        }
    }

    fn set_pixel_unsafe(&self, x: u32, y: u32, pixel: Pixel) {
        let index = (y*self.width + x) as usize;
        unsafe {
            (*self.pixels.get())[index] = pixel;
        }
    }
}

fn reflect(incoming: Vec3, n: Vec3) -> Vec3 {
    incoming - 2.0*incoming.dot(n)*n
}

fn construct_coordinate_system(normal: Vec3) -> (Vec3, Vec3, Vec3) {
    const EPS: f32 = 0.9999;
    let other = if normal.y > EPS || normal.y < -EPS {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let y = normal;
    let x = y.cross(other);
    let z = y.cross(x);
    (x, y, z)
}

fn sample_hemisphere_cos(xi: Vec2) -> Vec3 {
    let r = f32::sqrt(xi.x);
    let theta = 2.0*PI*xi.y;
 
    let x = r * f32::cos(theta);
    let z = r * f32::sin(theta);
 
    Vec3::new(x, f32::sqrt(f32::max(0.0, 1.0 - xi.x)), z)
}

fn brdf_lambert(reflectivity: Vec3, radiance: Vec3) -> Vec3 {
    reflectivity*radiance
}

fn trace_radiance(ray: &Ray, scene: &Scene, depth: u8) -> Vec3 {
    if depth <= 0 {
        return Vec3::zero();
    }

    let nearest_hit = find_scene_hit(ray, scene);

    if let Some(nearest_hit) = nearest_hit {
        let outwards_shifted_position = nearest_hit.position + 0.0001*nearest_hit.normal; // @TODO: Find a good factor

        match nearest_hit.material {
            Material::None => Vec3::one(),
            Material::Emissive(ref color) => color.clone(),
            Material::Mirror => {
                let reflection_direction = reflect(ray.direction, nearest_hit.normal);
                trace_radiance(&Ray::new(outwards_shifted_position, reflection_direction), scene, depth - 1)
            },
            Material::Phyiscally{ ref reflectivity, .. } => {
                let (ax, ay, az) = construct_coordinate_system(nearest_hit.normal);
                let xi = Vec2::new(random32(), random32());
                let h = sample_hemisphere_cos(xi);
                let direction = h.x*ax + h.y*ay + h.z*az;
                let ray = Ray::new(outwards_shifted_position, direction);

                let cos_theta_reflection = direction.dot(nearest_hit.normal);

                let reflection = trace_radiance(&ray, scene, depth - 1);

                brdf_lambert(*reflectivity, reflection)*cos_theta_reflection*PI
            },
        }
    } else {
        // @TODO: Make a nice gradient for the sky of sample an equirectangular projection.
        let theta = f32::acos(ray.direction.y);
        let t = f32::powf(theta / PI, 2.0);
        let intensity = 1.0 - 2.0*t;
        intensity*Vec3::new(0.6, 0.6, 0.8)
    }
}

#[allow(dead_code)]
fn tone_map_reinhard(radiance: Vec3) -> Vec3 {
    radiance / (Vec3::one() + radiance)
}

fn tone_map_clamp(radiance: Vec3) -> Vec3 {
    Vec3::new(
        saturatef32(radiance.x),
        saturatef32(radiance.y),
        saturatef32(radiance.z),
    )
}

pub fn render(work_tile: WorkTile, backbuffer: &Arc<Backbuffer>, camera: &Camera, scene: Arc<RwLock<Scene>>) {
    let camera_u = camera.projection_plane.u / backbuffer.width as f32;
    let camera_v = camera.projection_plane.v / backbuffer.height as f32;

    let (x0, x1) = (work_tile.position.x, work_tile.position.x + work_tile.size.x);
    let (y0, y1) = (work_tile.position.y, work_tile.position.y + work_tile.size.y);
    
    let scene = scene.read().unwrap(); // @TODO: Handle the unwrap

    for y in y0..y1 {
        for x in x0..x1 {
            let ray = {
                let origin = {
                    let du = x as f32*camera_u;
                    let dv = y as f32*camera_v;
                    camera.projection_plane.origin + du + dv
                };
                let direction = (origin - camera.eye).normalize();
                Ray::new(origin, direction)
            };

            let hdr_radiance = {
                const N: usize = 1;
                let mut hdr_radiance = Vec3::zero();
                for _ in 0..N {
                    hdr_radiance += trace_radiance(&ray, &*scene, 2);
                }
                hdr_radiance / N as f32
            };
            let ldr_radiance = tone_map_clamp(hdr_radiance);
            let color = Pixel::from_unit(ldr_radiance);
            backbuffer.set_pixel_unsafe(x, y, color);
        }
    }
}