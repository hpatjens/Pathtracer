use std;
use std::time::Duration;
use std::fs;
use std::sync::{Arc, RwLock};
use std::sync::mpsc::channel;

use time;
use glutin;
use glium;
use notify;
use notify::Watcher;

use common::*;
use worker;
use content;
use scene;
use tracer;
use parser;

#[link(name = "opengl32")]
extern "C" {
    // glium only provides functionality of the 'modern' OpenGL APIs. glDrawPixels is deprecated
    // and therefore not available. However, to date it is commonly still implemented and can be
    // used without much worries.
    fn glDrawPixels(width: u32, height: u32, format: i32, component_type: i32, data: *const u8);
}

fn load_scene(filename: &str) -> Result<scene::Scene, String> {
    match fs::read_to_string(filename) {
        Ok(content) => match parser::parse_scene(&*content) {
            Ok(scene) => Ok(scene),
            Err(err) => Err(format!("Could not parse the file \"{}\". Error: \"{:?}\".", filename, err)), // @TODO: Make the ParseError printable.
        },
        Err(_) => Err(format!("Could not read the file \"{}\".", filename)),
    }
}

pub fn start() {
    //
    // READING THE COMMAND LINE ARGUMENTS
    //
    // All parameters of the image as well as the rendering process should be specified
    // within the scene file. 
    //
    let scene_file_path = match {
        const ARG: &'static str = "--path=";
        std::env::args()
            .find(|arg| arg.starts_with(ARG))
            .map(|arg| String::from(arg.split_at(ARG.len()).1))
    } {
        Some(path) => path,
        None => String::from("scenes/sample/sample.scene"), // Let's try to load this one when no other was provided via the arguments.
    };

    //
    // TRYING TO LOAD THE SCENE FOR THE FIRST TIME
    //
    let scene = Arc::new(RwLock::new(load_scene(&*scene_file_path)
        .unwrap_or_else(|err| {
            println!("{}", err);
            println!("Using the default scene.");
            scene::Scene::default()
        })));
    
    let width = scene.read().unwrap().image_settings.width as u32;
    let height = scene.read().unwrap().image_settings.height as u32;

    //
    // SETTING UP THE HOT-RELOADING OF THE SCENE FILE
    //
    // @TODO: This should not be done with panics. Not being able to load the a file
    //        from disk is something that can happen in normal use.
    //
    let (sender, receiver) = channel();
    let mut watcher: notify::RecommendedWatcher = match notify::Watcher::new(sender, Duration::from_millis(100)) {
        Ok(watcher) => watcher,
        Err(_) => panic!("Could not create the watcher for the scene file."),
    };
    match watcher.watch(&*scene_file_path, notify::RecursiveMode::NonRecursive) {
        Ok(_) => (),
        Err(_) => panic!("Could not start the watcher for the scene file."),
    };

    //
    // SETTING UP WINDOW AND EVENT LOOP
    //
    // @TODO: The backbuffer should resize when the window is resized. This is
    //        not implemented right not but also questionable as the resolution
    //        of the image should be specified by the scene file. Maybe it's
    //        best to have an internal resolution with which to export images and
    //        a resolution for displaying the image on the window.
    //
    let mut events_loop = glium::glutin::EventsLoop::new();
    let window = glium::glutin::WindowBuilder::new()
        .with_dimensions(glium::glutin::dpi::LogicalSize::new(width as f64, height as f64))
        .with_title("Pathtracer");
    let display = glium::Display::new(window, glium::glutin::ContextBuilder::new(), &events_loop).unwrap();

    //
    // SETTING UP THE BACKBUFFER
    //
    // @TODO: Implement resizing. View comment on the window stuff.
    //
    let backbuffer = Arc::new(tracer::Backbuffer::new(width, height));

    //
    // STARTING THE WORKER POOL
    //
    let worker_pool = {
        const NUM_WORKER_THREADS: usize = 8;
        let backbuffer2 = backbuffer.clone();
        let scene2 = scene.clone();
        worker::WorkerPool::new(NUM_WORKER_THREADS, Box::new(move |work_tile| {
            tracer::render(work_tile, &backbuffer2, scene2.clone());
        }))
    };

    //
    // SETTING UP CONTENT MANAGEMENT
    //
    let content = content::Content::new();

    //
    // MAIN LOOP
    //
    let mut running = true;
    while running {
        let frame_time_start = time::precise_time_ns();

        // ASSIGNING THE HDRI TEXTURE FOR THE SKY
        {
            let mut scene = scene.write().expect("Could not get writing access to the scene for updating the sky texture.");
            if let scene::Sky::HDRI(ref path, ref mut option_texture) = scene.sky {
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

        // RELOADING THE SCENE WHEN THE FILE WAS UPDATED
        match receiver.try_recv() {
            Ok(event) => {
                if let notify::DebouncedEvent::Write(_) = event {
                    let mut scene = scene.write().unwrap(); // @TODO: Handle the unwrap
                    match load_scene(&*scene_file_path) {
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

        // GENERATING WORK ITEMS FOR THE WORKER THREADS TO DO THE PATH TRACING
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

        //
        // WAITING UNTIL THE WORKER POOL HAS FINISHED THE PATH TRACING
        //
        // There are two reasons for waiting on the path tracing result:
        //   1. All pixels should be written to the backbuffer before it is
        //      displayed on the window. This however is not really important
        //      considering the amount of noise. This would only be noticeable
        //      when the camera changes every frame which could result in 
        //      complete tiles not being ready.
        //   2. There must be a mechanism that ensures that the work items for
        //      the tiles don't pile up in the queue for the worker threads.
        //      Only one work item per tile should be active at a time. Waiting
        //      for all tiles to be ready ensures this.
        //
        worker_pool.wait();
        assert!(worker_pool.queue_len() == 0);

        // WRITING THE BACKBUFFER TO THE WINDOW
        let target = display.draw();
        unsafe {
            // Constants from gl.h
            const GL_RGB: i32 = 0x1907;
            const GL_UNSIGNED_BYTE: i32 = 0x1401;

            let raw = &(*backbuffer.pixels8.get())[0].0 as *const u8;
            glDrawPixels(backbuffer.width,
                         backbuffer.height,
                         GL_RGB,
                         GL_UNSIGNED_BYTE,
                         raw);
            *backbuffer.num_samples.get() += 1;
        };
        target.finish().unwrap();

        // HANDLE THE EVENTS PROVIDED BY THE WINDOW
        events_loop.poll_events(|ev| {
            match ev {
                glutin::Event::WindowEvent { event, .. } => match event {
                    glutin::WindowEvent::CloseRequested => running = false,
                    _ => (),
                },
                _ => (),
            }
        });

        // MEASURING THE FRAME TIME
        let frame_time_end = time::precise_time_ns();
        println!("frame_time = {} ms", (frame_time_end - frame_time_start) as f64 / 1_000_000.0);
    }
}
