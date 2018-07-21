#[macro_use] extern crate derive_new;

extern crate glium;
extern crate glutin;
extern crate hmath;
extern crate rand;
extern crate time;

mod common;
use common::*;

mod scene;
mod worker;
mod tracer;

use scene::{Scene, Plane, Sphere, PBRParameters, Material};
use tracer::Camera;

use glium::glutin::dpi::LogicalSize;

use std::sync::{Arc, RwLock};

#[link(name = "opengl32")]
extern "C" {
    fn glDrawPixels(width: u32, height: u32, format: i32, component_type: i32, data: *const u8);
}

const GL_RGB: i32 = 0x1907;
const GL_UNSIGNED_BYTE: i32 = 0x1401;

fn update(scene: &Arc<RwLock<Scene>>, camera: &Arc<RwLock<tracer::Camera>>, frame_index: usize) {
    {
        let mut camera = camera.write().unwrap(); // @TODO: Handle the unwrap

        let x = frame_index as f32 / 80.0;
        const D: f32 = 25.0;
        let position = Vec3::new(D*f32::cos(x), 2.0 + f32::cos(x), D*f32::sin(x));
        camera.look_at(position, Vec3::zero(), Vec3::new(0.0, 1.0, 0.0))
    }

    {
        let material_sphere = Material::Phyiscally(PBRParameters {
            reflectivity: Vec3::new(0.45, 0.05, 0.25),
            roughness: 0.3,
            metalness: 0.0,
        });
        let material_plane = Material::Phyiscally(PBRParameters {
            reflectivity: Vec3::new(0.91, 0.92, 0.92),
            roughness: 0.6,
            metalness: 1.0,
        });

        let x = frame_index as f32 / 40.0;
        let position1 = Vec3::new(2.3*f32::sin(x), f32::cos(x), 2.1*f32::cos(x));
        let position2 = Vec3::new(1.2*f32::sin(1.12*x + 0.124), f32::cos(1.45*x + 0.7567), 1.6*f32::cos(0.923*x + 0.2345));
        let position3 = Vec3::new(f32::sin(1.43*x + 0.224), f32::cos(1.76*x + 0.2134), f32::cos(0.123*x + 0.6346));
        let light_position = Vec3::new(0.0, 3.5, -1.0);

        let mut scene = scene.write().unwrap(); // @TODO: Handle the unwrap
        *scene = Scene::new(
            vec![
                Sphere::new(light_position, 2.0, Material::Emissive(Vec3::new(1.0, 1.0, 1.0))),
                Sphere::new(position1, 1.0, material_sphere.clone()),
                Sphere::new(position2, 0.5, Material::Emissive(Vec3::new(0.0, 5.0, 0.0))),
                Sphere::new(position3, 1.2, Material::Glass),
            ],
            vec![
                Plane::new(Vec3::new(0.0, -2.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0), material_plane.clone())
            ]
        )
    }
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

    let backbuffer = Arc::new(tracer::Backbuffer::new(width, height));

    let scene = Arc::new(RwLock::new(Scene::new(Vec::new(), Vec::new())));
    let camera = Arc::new(RwLock::new(Camera::new(Vec3::one(), Vec3::zero(), Vec3::new(0.0, 1.0, 0.0), 4.0, 4.0, 10.0)));

    let worker_pool = {
        const NUM_WORKER_THREADS: usize = 8;
        let backbuffer2 = backbuffer.clone();
        let camera2 = camera.clone();
        let scene2 = scene.clone();
        worker::WorkerPool::new(NUM_WORKER_THREADS, Box::new(move |work_tile| {
            tracer::render(work_tile, &backbuffer2, camera2.clone(), scene2.clone());
        }))
    };

    let mut running = true;
    while running {
        let frame_time_start = time::precise_time_ns();

        update(&scene, &camera, frame_index);

        for _ in 0..1 {
            const TILE_SIZE: u32 = 32;
            let tile_size = Vec2u::new(TILE_SIZE, TILE_SIZE);
            let num_tiles_x = (backbuffer.width + TILE_SIZE - 1) / TILE_SIZE;
            let num_tiles_y = (backbuffer.height + TILE_SIZE - 1) / TILE_SIZE;
            for y in 0..num_tiles_y {
                for x in 0..num_tiles_x {
                    let tile_index = Vec2u::new(x, y);
                    let tile_position = Vec2u::new(x*TILE_SIZE, y*TILE_SIZE);
                    let work_tile = tracer::WorkTile::new(tile_index, tile_position, tile_size);
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
