#[macro_use] extern crate derive_new;

extern crate glium;
extern crate glutin;
extern crate hmath;
extern crate rand;
extern crate time;

mod worker;
use worker::*;

use glium::glutin::dpi::LogicalSize;

use hmath::*;

use std::cell::UnsafeCell;
use std::sync::{Arc, RwLock};

type Vec2 = Vector2<f32>;
type Vec3 = Vector3<f32>;

type Vec2u = Vector2<u32>;

#[link(name = "opengl32")]
extern "C" {
    fn glDrawPixels(width: u32, height: u32, format: i32, component_type: i32, data: *const u8);
}

const GL_RGB: i32 = 0x1907;
const GL_UNSIGNED_BYTE: i32 = 0x1401;

const PI: f32 = std::f32::consts::PI;

// Source: https://de.wikipedia.org/wiki/Xorshift
static mut X32: u32 = 314159265;
fn xorshift32() -> u32 {
    unsafe { 
        X32 ^= X32 << 13;
        X32 ^= X32 >> 17;
        X32 ^= X32 << 5;
        X32
    }
}

fn random32() -> f32 {
    let r = xorshift32();
    r as f32 / std::u32::MAX as f32
}

fn clampf32(min: f32, max: f32, x: f32) -> f32 {
    if x > max {
        max
    } else if x < min {
        min
    } else {
        x
    }
}

fn saturatef32(x: f32) -> f32 {
    clampf32(0.0, 1.0, x)
}

#[derive(Clone)]
struct Pixel(u8, u8, u8);

impl Pixel {
    fn from_unit(color: Vec3) -> Self {
        Pixel((color.x*255.0) as u8, (color.y*255.0) as u8, (color.z*255.0) as u8)
    }
}

#[derive(Clone, Debug, new)]
struct WorkTile {
    tile_index: Vec2u,
    position: Vec2u,
    size: Vec2u,
}

struct Backbuffer {
    width: u32,
    height: u32,
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
    pixels: UnsafeCell<Vec<Pixel>>,
}
// UnsafeCell does not implement Sync and therefore Backbuffer could not be passed to the
// worker threads without this implementation.
unsafe impl Sync for Backbuffer {}

impl Backbuffer {
    fn new(width: u32, height: u32) -> Backbuffer {
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

fn look_at(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32) -> Camera {
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

fn update(scene: &Arc<RwLock<Scene>>, frame_index: usize) {
    let material1 = Material::Phyiscally{
        reflectivity: Vec3::new(0.7, 0.73, 0.72),
        roughness: 1.0,
        metalness: 0.0,
    };
    let material2 = Material::Phyiscally{
        reflectivity: Vec3::new(0.5, 0.5, 0.5),
        roughness: 1.0,
        metalness: 0.0,
    };

    let x = frame_index as f32 / 40.0;
    let position1 = Vec3::new(f32::sin(x), f32::cos(x), f32::cos(x));
    let position2 = Vec3::new(f32::sin(1.12*x + 0.124), f32::cos(1.45*x + 0.7567), f32::cos(0.923*x + 0.2345));
    let light_position = Vec3::new(0.0, 3.5, -1.0);

    let mut scene = scene.write().unwrap(); // @TODO: Handle the unwrap
    *scene = Scene::new(
        vec![
            Sphere::new(light_position, 2.0, Material::Emissive(Vec3::new(1.0, 1.0, 1.0))),
            Sphere::new(position1, 1.0, material1.clone()),
            Sphere::new(position2, 1.0, Material::Emissive(Vec3::new(0.0, 1.0, 0.0))),
        ],
        vec![
            Plane::new(Vec3::new(0.0, -2.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0), material2.clone())
        ]
    )
}

fn main() {
    let width: u32 = 512;
    let height: u32 = 512;

    let logical_size = LogicalSize::new(width as f64, height as f64);

    let mut events_loop = glium::glutin::EventsLoop::new();
    let window = glium::glutin::WindowBuilder::new()
        .with_dimensions(logical_size)
        .with_title("Pathtracer");
    let context = glium::glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();

    let mut frame_index = 0;

    let camera = {
        let x = frame_index as f32 / 40.0;
        const D: f32 = 20.0;
        let position = Vec3::new(D*f32::cos(x), 2.0 + f32::cos(x), D*f32::sin(x));
        Arc::new(look_at(position, Vec3::zero(), Vec3::new(0.0, 1.0, 0.0), 4.0, 4.0, 10.0))
    };

    let backbuffer = Arc::new(Backbuffer::new(width, height));

    let scene = Arc::new(RwLock::new(Scene::new(Vec::new(), Vec::new())));

    let worker_pool = {
        const NUM_WORKER_THREADS: usize = 4;
        let backbuffer2 = backbuffer.clone();
        let camera = camera.clone();
        let scene2 = scene.clone();
        WorkerPool::new(NUM_WORKER_THREADS, Box::new(move |work_tile| {
            render(work_tile, &backbuffer2, &camera, scene2.clone());
        }))
    };

    let mut running = true;
    while running {
        let frame_time_start = time::precise_time_ns();

        update(&scene, frame_index);

        for _ in 0..1 {
            const TILE_SIZE: u32 = 32;
            let tile_size = Vec2u::new(TILE_SIZE, TILE_SIZE);
            let num_tiles_x = (backbuffer.width + TILE_SIZE - 1) / TILE_SIZE;
            let num_tiles_y = (backbuffer.height + TILE_SIZE - 1) / TILE_SIZE;
            for y in 0..num_tiles_y {
                for x in 0..num_tiles_x {
                    let tile_index = Vec2u::new(x, y);
                    let tile_position = Vec2u::new(x*TILE_SIZE, y*TILE_SIZE);
                    let work_tile = WorkTile::new(tile_index, tile_position, tile_size);
                    worker_pool.process(work_tile);
                }
            }
        }

        worker_pool.wait();
        assert!(worker_pool.queue_len() == 0);

        let target = display.draw();

        unsafe {
            let raw = &(*backbuffer.pixels.get())[0].0 as *const u8;
            glDrawPixels(backbuffer.width,
                         backbuffer.height,
                         GL_RGB,
                         GL_UNSIGNED_BYTE,
                         raw);
        };

        target.finish().unwrap();

        events_loop.poll_events(|ev| {
            match ev {
                glutin::Event::WindowEvent { event, .. } => match event {
                    glutin::WindowEvent::CloseRequested => running = false,
                    _ => (),
                },
                _ => (),
            }
        });

        frame_index += 1;

        let frame_time_end = time::precise_time_ns();
        println!("frame_time = {} ms", (frame_time_end - frame_time_start) as f64 / 1_000_000.0);
    }
}

#[derive(Clone, Debug, new)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
enum Material {
    None,
    Emissive(Vec3),
    Mirror,
    Phyiscally {
        reflectivity: Vec3,
        roughness: f32,
        metalness: f32,
    },
}

#[derive(Clone, Debug, new)]
struct Plane {
    origin: Vec3,
    u: Vec3,
    v: Vec3,
    material: Material,
}

#[derive(Clone, Debug, new)]
struct Sphere {
    origin: Vec3,
    radius: f32,
    material: Material,
}

#[derive(Debug, new)]
struct Camera {
    projection_plane: Plane,
    eye: Vec3,
}

#[derive(Debug, new)]
struct Scene {
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
}

#[derive(Debug, Clone, new)]
struct Hit<'a> {
    parameter: f32,
    position: Vec3,
    normal: Vec3,
    material: &'a Material,
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

    let parameter = if t1 < t2 { t1 } else { t2 };
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

fn find_scene_hit<'a>(ray: &Ray, scene: &'a Scene) -> Option<Hit<'a>> {
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

fn render(work_tile: WorkTile, backbuffer: &Arc<Backbuffer>, camera: &Camera, scene: Arc<RwLock<Scene>>) {
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