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

    fn from_signed_unit(color: Vec3) -> Self {
        Self::from_unit(Vec3::new(0.5, 0.5, 0.5) + 0.5*color)
    }
}

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
    let width: u32 = 256;
    let height: u32 = 256;

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
            Plane::new(origin, u, v, Material::None)
        };
        let eye = Vec3::new(0.0, 0.0, -20.0);
        Camera::new(projection_plane, eye)
    };

    let mut frame_index = 0;

    let mut running = true;
    while running {
        let target = display.draw();

        let scene = {
            let x = frame_index as f32 / 40.0;
            let position1 = Vec3::new(f32::sin(x), f32::cos(x), f32::cos(x));
            let position2 = Vec3::new(f32::sin(1.12*x + 0.124), f32::cos(1.45*x + 0.7567), f32::cos(0.923*x + 0.2345));
            Scene::new(
                vec![
                    Sphere::new(position1, 1.0, Material::Mirror),
                    Sphere::new(position2, 1.0, Material::Color(Vec3::new(0.0, 1.0, 0.0))),
                ],
                vec![
                    Plane::new(Vec3::new(0.0, -1.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0), Material::Color(Vec3::new(0.2, 0.3, 0.4)))
                ]
            )
        };
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

        frame_index += 1;
    }
}

#[derive(Clone, Debug, new)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

#[derive(Clone, Debug)]
enum Material {
    None,
    Color(Vec3),
    Mirror,
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

fn trace_radiance(ray: &Ray, scene: &Scene, depth: u8) -> Vec3 {
    if depth <= 0 {
        return Vec3::zero();
    }

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

    if let Some(nearest_hit) = nearest_hit {
        let reflection = {
            let reflection_direction = reflect(ray.direction, nearest_hit.normal);
            let outwards_shifted_position = nearest_hit.position + 0.0001*nearest_hit.normal; // @TODO: Find a good factor
            trace_radiance(&Ray::new(outwards_shifted_position, reflection_direction), scene, depth - 1)
        };

        let cos_theta = ray.direction.dot(nearest_hit.normal);

        match nearest_hit.material {
            Material::None => Vec3::one(),
            Material::Color(ref color) => color.clone(),
            Material::Mirror => reflection,
        }
    } else {
        Vec3::new(0.8, 0.8, 0.9)
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

            let hdr_radiance = trace_radiance(&ray, scene, 2);
            let ldr_radiance = tone_map_clamp(hdr_radiance);
            let color = Pixel::from_unit(ldr_radiance);
            backbuffer.set(x, y, color);
        }
    }
}