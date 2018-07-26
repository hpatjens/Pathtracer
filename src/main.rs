#[macro_use] extern crate derive_new;

extern crate glium;
extern crate glutin;
extern crate hmath;
extern crate rand;
extern crate time;
extern crate notify;
extern crate stb_image;

mod common;
use common::*;

mod scene;
mod worker;
mod tracer;
mod parser;
mod content;

use content::Content;

use scene::{Scene, Sky, HDRITexture};
use tracer::Camera;

use glium::glutin::dpi::LogicalSize;
use notify::Watcher;

use std::time::Duration;
use std::io::prelude::*;
use std::fs::File;
use std::sync::{Arc, RwLock};
use std::sync::mpsc::channel;

#[link(name = "opengl32")]
extern "C" {
    fn glDrawPixels(width: u32, height: u32, format: i32, component_type: i32, data: *const u8);
}

const GL_RGB: i32 = 0x1907;
const GL_UNSIGNED_BYTE: i32 = 0x1401;

fn update(camera: &Arc<RwLock<tracer::Camera>>, frame_index: usize) {

    // @TODO: The camera should be specified by the scene file. Parameters like t could be used too.
    return;

    let mut camera = camera.write().unwrap(); // @TODO: Handle the unwrap

    let x = frame_index as f32 / 80.0;
    const D: f32 = 25.0;
    let position = Vec3::new(D*f32::cos(x), 2.0, D*f32::sin(x));
    camera.look_at(position, Vec3::zero(), Vec3::new(0.0, 1.0, 0.0))
}

fn load_scene(filename: &str) -> Result<Scene, parser::ParseError> {
    let content = {
        let mut content = String::new();
        let mut file = File::open(filename).expect("File not found.");
        file.read_to_string(&mut content).expect("Could not read the file.");
        content
    };
    parser::parse_scene(&*content)
}

fn main() {
    let (sender, receiver) = channel();

    // @TODO: Search for a file in the current folder and add a command line argument.
    const SCENE_FILENAME: &str = "scenes/sample/sample.scene";

    let mut watcher: notify::RecommendedWatcher = match notify::Watcher::new(sender, Duration::from_millis(100)) {
        Ok(watcher) => watcher,
        Err(_) => panic!("Could not create the watcher for the scene file."),
    };
    match watcher.watch(SCENE_FILENAME, notify::RecursiveMode::NonRecursive) {
        Ok(_) => (),
        Err(_) => panic!("Could not start the watcher for the scene file."),
    };

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

    let scene = match load_scene(SCENE_FILENAME) {
        Ok(scene) => Arc::new(RwLock::new(scene)),
        Err(err) => {
            println!("Could not load the scene. Error: {:?}", err);
            Arc::new(RwLock::new(Scene::new(Sky::Constant(Vec3::new(1.0, 1.0, 1.0)), Vec::new(), Vec::new())))
        },
    };
    let camera = Arc::new(RwLock::new(Camera::new(Vec3::new(0.0, 2.0, 20.0), Vec3::zero(), Vec3::new(0.0, 1.0, 0.0), 4.0, 4.0, 10.0)));

    let worker_pool = {
        const NUM_WORKER_THREADS: usize = 8;
        let backbuffer2 = backbuffer.clone();
        let camera2 = camera.clone();
        let scene2 = scene.clone();
        worker::WorkerPool::new(NUM_WORKER_THREADS, Box::new(move |work_tile| {
            tracer::render(work_tile, &backbuffer2, camera2.clone(), scene2.clone());
        }))
    };

    let content = Content::new();

    let mut running = true;
    while running {
        let frame_time_start = time::precise_time_ns();

        {
            let mut scene = scene.write().expect("Could not get writing access to the scene for updating the sky texture.");
            if let Sky::HDRI(ref path, ref mut option_texture) = scene.sky {
                if let None = option_texture {
                    *option_texture = match content.get_hdri_texture(path) {
                        Some(hdri_texture) => {
                            backbuffer.clear();
                            Some(hdri_texture)
                        },
                        None => None,
                    };
                }
            }
        }

        match receiver.try_recv() {
            Ok(event) => {
                if let notify::DebouncedEvent::Write(_) = event {
                    let mut scene = scene.write().unwrap(); // @TODO: Handle the unwrap
                    match load_scene(SCENE_FILENAME) {
                        Ok(loaded_scene) => {
                            *scene = loaded_scene;
                            backbuffer.clear();
                        },
                        Err(err) => {
                            println!("Could not load the scene. Error: {:?}", err);
                        },
                    };
                }
            },
            Err(_) => (), // @TODO: If disconnected panic
        }

        update(&camera, frame_index);

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
            let raw = &(*backbuffer.pixels8.get())[0].0 as *const u8;
            glDrawPixels(backbuffer.width,
                         backbuffer.height,
                         GL_RGB,
                         GL_UNSIGNED_BYTE,
                         raw);
            *backbuffer.num_samples.get() += 1;
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
