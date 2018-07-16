extern crate glium;
extern crate glutin;

use glium::glutin::dpi::LogicalSize;

#[link(name = "opengl32")]
extern "C" {
    fn glDrawPixels(width: u32, height: u32, format: i32, component_type: i32, data: *const u8);
}

const GL_RGB: i32 = 0x1907;
const GL_UNSIGNED_BYTE: i32 = 0x1401;

#[derive(Clone)]
struct Pixel(u8, u8, u8);

struct Backbuffer {
    width: u32,
    height: u32,
    pixels: Vec<Pixel>,
}

impl Backbuffer {
    fn new(width: u32, height: u32) -> Backbuffer {
        Backbuffer {
            width: width,
            height: height,
            pixels: {
                let mut pixels = Vec::new();
                let num_pixels = (width * height) as usize;
                pixels.resize(num_pixels, Pixel(0, 0, 0));
                pixels
            },
        }
    }
}

fn main() {
    let width: u32 = 1024;
    let height: u32 = 768;

    let logical_size = LogicalSize::new(width as f64, height as f64);

    let mut events_loop = glium::glutin::EventsLoop::new();
    let window = glium::glutin::WindowBuilder::new()
        .with_dimensions(logical_size)
        .with_title("Pathtracer");
    let context = glium::glutin::ContextBuilder::new();
    let display = glium::Display::new(window, context, &events_loop).unwrap();

    let backbuffer = Backbuffer::new(width, height);

    let mut running = true;
    while running {
        let target = display.draw();

        // render();

        unsafe {
            let raw = &backbuffer.pixels[0].0 as *const u8;
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
    }
}