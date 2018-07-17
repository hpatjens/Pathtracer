#[macro_use] extern crate derive_new;

extern crate glium;
extern crate glutin;
extern crate hmath;

use glium::glutin::dpi::LogicalSize;

use hmath::*;

type Vec3 = Vector3<f32>;

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

    fn set(&mut self, x: u32, y: u32, pixel: Pixel) {
        let index = (y*self.width + x) as usize;
        self.pixels[index] = pixel;
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

    let mut backbuffer = Backbuffer::new(width, height);

    let camera = {
        let projection_plane = {
            let origin = Vec3::new(-2.0, -2.0, -5.0);
            let u = Vec3::new(4.0 / width as f32, 0.0, 0.0);
            let v = Vec3::new(0.0, 4.0 / height as f32, 0.0);
            Plane::new(origin, u, v)
        };
        let eye = Vec3::new(0.0, 0.0, -10.0);
        Camera::new(projection_plane, eye)
    };
    let scene = Scene::new(vec![Sphere::new(Vec3::zero(), 1.0)]);

    let mut running = true;
    while running {
        let target = display.draw();

        render(&mut backbuffer, &camera, &scene);

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

#[derive(Clone, Debug, new)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

#[derive(Clone, Debug, new)]
struct Plane {
    origin: Vec3,
    u: Vec3,
    v: Vec3,
}

#[derive(Clone, Debug, new)]
struct Sphere {
    origin: Vec3,
    radius: f32,
}

#[derive(Debug, new)]
struct Camera {
    projection_plane: Plane,
    eye: Vec3,
}

#[derive(Debug, new)]
struct Scene {
    spheres: Vec<Sphere>,
}
fn intersect(sphere: &Sphere, ray: &Ray) -> bool {
    let w = sphere.origin - ray.origin;
    let e = ray.direction.dot(w)*ray.direction;
    let r = w - e;
       r.length() < sphere.radius
}

fn render(backbuffer: &mut Backbuffer, camera: &Camera, scene: &Scene) {
    for y in 0..backbuffer.height {
        for x in 0..backbuffer.width {
            let ray = {
                let origin = {
                    let du = x as f32*camera.projection_plane.u;
                    let dv = y as f32*camera.projection_plane.v;
                    camera.projection_plane.origin + du + dv
                };
                let direction = (origin - camera.eye).normalize();
                Ray::new(origin, direction)
            };

            for sphere in &scene.spheres {
                if intersect(sphere, &ray) {
                    backbuffer.set(x, y, Pixel(255, 65, 21));
                } else {
                    backbuffer.set(x, y, Pixel(34, 34, 34));
                }
            }
        }
    }
}