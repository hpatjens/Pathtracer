use common::*;

use scene::{Scene, Plane, PBRParameters, Material, find_scene_hit};

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

#[derive(Debug)]
pub struct Camera {
    projection_plane: Plane,
    position: Vec3,
    width: f32,
    height: f32,
    z_near: f32,
}

impl Camera {
    pub fn new(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32) -> Self {
        Camera {
            projection_plane: Self::construct_projection_plane(position, target, up, width, height, z_near),
            position: position,
            width: width,
            height: height,
            z_near: z_near,
        }
    }

    fn construct_projection_plane(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32) -> Plane {
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
    }

    pub fn look_at(&mut self, position: Vec3, target: Vec3, up: Vec3) {
        self.projection_plane = Self::construct_projection_plane(position, target, up, self.width, self.height, self.z_near);
        self.position = position;
    }

    fn sample(&self, backbuffer_width: u32, backbuffer_height: u32) -> CameraSampler {
        CameraSampler {
            camera: self,
            u: self.projection_plane.u / backbuffer_width as f32,
            v: self.projection_plane.v / backbuffer_height as f32,
        }
    }
}

struct CameraSampler<'a> {
    camera: &'a Camera,
    u: Vec3,
    v: Vec3,
}

// @TODO: Make a trait when more camera types are added. Different camera models can
//        precompute different things.
impl<'a> CameraSampler<'a> {
    fn pinhole_ray(&self, x: u32, y: u32) -> Ray {
        let origin = {
            let du = x as f32*self.u;
            let dv = y as f32*self.v;

            self.camera.projection_plane.origin + du + dv
        };
        let direction = (origin - self.camera.position).normalize();
        Ray::new(origin, direction)
    }

    fn thin_lens_ray(&self, x: u32, y: u32) -> Ray {
        let Ray{ ref origin, ref direction } = self.pinhole_ray(x, y);

        // In this thin-lens model a focus distance and an aperture is defined. To 
        // determine the ray with which to sample, the point on the focus plane is
        // computed which would be hit by the ideal ray through the pixel. Then the
        // circle of confusion is computed which would be produced by the apterture.
        // A point within the circle of confusion is determined randomly as the new
        // origin point of the ray. The direction is computed by the difference of
        // the point of the focus plane and the point within the circle of confusion
        // on the projection plane.

        const FOCUS_DISTANCE: f32 = 10.0;

        // @TODO: This results in a spherical focus "plane" which is not really ideal
        //        Use the unnormalized direction here dividing by the distance to the
        //        projection plane to get the distance right.
        let point_on_focus_plane = origin + FOCUS_DISTANCE*direction;

        // To sample a circle correctly (with a normal distribution), the angle and
        // radius have to be computed as follows:
        //
        //   alpha = 2*PI*R*sqrt(r1)
        //   r = R*sqrt(r2)
        //
        // where 
        //   R = maximum radius (circle of confusion in this case)
        //   r1, r2 = random normally distributed numbers in the interval [0, 1]
        //
        // This however is not really what we want in this case as the probability
        // of a ray hitting the sensor at a given position is not normally distributed
        // within the circle of confusion.

        const LENS_APERTURE: f32 = 4.0;

        // The factor 22 means that an apterture of 22 is perfectly sharp and that 
        // lower apertures > 0 produce larger confusion circles.
        let circle_of_confusion = 22.0/LENS_APERTURE;

        let r1 = random32();
        let r2 = random32();

        let alpha = 2.0*PI*r1;
        let r = circle_of_confusion*r2;

        let l_u = r*f32::sin(alpha);
        let l_v = r*f32::cos(alpha);

        let lens_origin = origin + l_u*self.u + l_v*self.v;
        let lens_direction = (point_on_focus_plane - lens_origin).normalize();

        Ray::new(lens_origin, lens_direction)
    }
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

fn reflect(incoming: Vec3, normal: Vec3) -> Vec3 {
    // Bear in mind that incoming is directed at the surface!
    incoming - 2.0*incoming.dot(normal)*normal
}

fn refract(incoming: Vec3, normal: Vec3, n1: f32, n2: f32) -> Vec3 {
    // Bear in mind that incoming is directed at the surface!

    // Pythagorean trigonometric identity: sin^2(a) + cos^2(a) = 1
    // Since cos(a) is easy to compute, sin(a) can be computed by
    // sin(a) = sqrt(1 - cos^2(a))
    let cos_theta1 = -incoming.dot(normal);
    // @TODO: Test whether this is really faster.
    // @TODO: Implement a faster but worse sqrt.
    // Look at this article: https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
    let sin_theta1 = f32::sqrt(1.0 - cos_theta1*cos_theta1);

    // Snell's law: sin(theta_1) / sin(theta_2) = n_2 / n_1
    // where theta_1: angle between the incoming ray and normal
    //       theta_2: angle between the outgoing ray and -normal
    //       n1: index of refraction for the medium above the surface
    //       n2: index of refraction for the medium below the surface
    let sin_theta2 = (sin_theta1*n1)/n2;
    
    // Corresponds to the projection of 'incoming' onto the ground plane.
    let p = incoming + cos_theta1*normal;
    
    -normal + p.normalize()*sin_theta2
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

#[allow(dead_code)]
fn brdf_lambert(pbr_parameters: &PBRParameters) -> Vec3 {
    let &PBRParameters{ ref reflectivity, .. } = pbr_parameters;
    reflectivity/PI
}

fn normal_distribution_ggx(normal: Vec3, half: Vec3, alpha: f32) -> f32 {
    // The normal distribution function computes how much the microfacets of the surface
    // contribute to a reflection whose half vector is 'half'.
    
    // More in-depth explanation: The microfacet theory assumes that the relevant reflection 
    // characteristics of every surface can be modelled by an in infinite amount of
    // infinitesimal small mirrors where the orientation and visibility of these mirrors are 
    // only defined statistically. The normal distribution function represents the statistical
    // orientation of the microfacets in the sense that it returns the probability of a
    // microfacet being oriented in the 'half' direction. Alternatively this can be seen as
    // determining the surface area of the microfacets being oriented in the half direction
    // relativ to the surface area of the macro surface. Therefore, when the intensity of
    // a reflection from direction d_1 to direction d_2 has to be determined, the half vector
    // between both directions can be determined and the normal distribution function returns
    // the strength of the reflection as those microfacets oriented in the 'half' direction
    // contribute to exactly that reflection.

    let n = normal; // normal of the surface
    let h = half; // half vector (normal of microfacets that contribute to the reflection)
    let a = alpha; // alpha value that is based on the roughness of the surface (e.g. roughness^2)

    let a2 = a*a;
    let n_dot_h = f32::max(n.dot(h), 0.0);
    let n_dot_h2 = n_dot_h*n_dot_h;

    let d = n_dot_h2 * (a2 - 1.0) + 1.0;
	
    a2 / (PI*d*d)
}

fn geometry_smith(normal: Vec3, view: Vec3, light: Vec3, k: f32) -> f32 {
    // This function computes the probability of light being blocked either by 
    //   1. the microfacet being shadowed by another microfacet when the light is incoming or
    //   2. the microfacet being masked by another microfacet when the light is outgoing.
    // Another way to lock at is would be: This function reduces the amount of reflected energy
    // by considering the proportional surface area of microfacets that is shadowed or masked
    // by other microfacets relativ to the surface area of the macro surface.

    let n = normal; // normal of the surface
    let v = view; // direction to the viewer
    let l = light; // direction to the light source

    let geometry_schlick_ggx = |a: Vec3| -> f32 {
        let n_dot_a = f32::max(0.0, n.dot(a));
        let nom = n_dot_a;
        let denom = n_dot_a * (1.0 - k) + k;
        nom / denom
    };
    geometry_schlick_ggx(l)*geometry_schlick_ggx(v)
}

fn fresnel_schlick(cos_theta: f32, f0: Vec3) -> Vec3 {
    f0 + (Vec3::one() - f0)*f32::powf(1.0 - cos_theta, 5.0)
}

fn brdf_cook_torrance(view: Vec3, light: Vec3, normal: Vec3, pbr_parameters: &PBRParameters) -> Vec3 {
    // Bear in mind that 'view' as well as 'light' are pointing away from the surface!

    let v = view;
    let l = light;
    let n = normal;

    let &PBRParameters{ reflectivity, roughness, metalness } = pbr_parameters;

    let alpha = roughness*roughness;
    let h = (v + l).normalize();
    let cos_theta = h.dot(v);

    // @TODO: There are multiple definitions for this. What's the best?
    // let k = (alpha + 1.0)*(alpha + 1.0) / 8.0;
    let k = alpha / 2.0; // Was/is used in the Unreal Engine 4 according to Brian Karis' blog.

    const R: f32 = 0.04;
    let f0 = mix_vec3(Vec3::new(R, R, R), reflectivity, metalness);

    let num = normal_distribution_ggx(n, h, alpha)*geometry_smith(n, v, l, k)*fresnel_schlick(cos_theta, f0);
    let denum = 4.0*n.dot(l)*n.dot(v);
    num / denum
}

fn trace_radiance(ray: &Ray, scene: &Scene, depth: u8) -> Vec3 {
    // @TODO: Find all the places where NANs can be generated and fix as many as it makes sense.

    if depth <= 0 {
        return Vec3::zero();
    }

    let nearest_hit = find_scene_hit(ray, scene);

    if let Some(nearest_hit) = nearest_hit {
        const SHIFT_AMOUNT: f32 = 0.0001; // @TODO: Find a good factor and maybe make it dependent on the slope
        let outwards_shifted_position = ||{ nearest_hit.position + SHIFT_AMOUNT*nearest_hit.normal };
        let inwards_shifted_position  = ||{ nearest_hit.position - SHIFT_AMOUNT*nearest_hit.normal };

        match nearest_hit.material {
            Material::None => Vec3::one(),
            Material::Emissive(ref color) => color.clone(),
            Material::Mirror => {
                let reflection_direction = reflect(ray.direction, nearest_hit.normal);
                trace_radiance(&Ray::new(outwards_shifted_position(), reflection_direction), scene, depth - 1)
            },
            Material::Glass => {
                // @TODO: Here, the direction has to be taken into account (inside-out or outside-in).
                let refraction_direction = refract(ray.direction, nearest_hit.normal, 1.0, 1.5);
                trace_radiance(&Ray::new(inwards_shifted_position(), refraction_direction), scene, depth - 1)
            },
            Material::Phyiscally(ref pbr_parameters) => {
                let (ax, ay, az) = construct_coordinate_system(nearest_hit.normal);
                let xi = Vec2::new(random32(), random32());
                let h = sample_hemisphere_cos(xi);
                // @TODO: Here, importance sampling is needed. Otherwise the reflection becomes really 
                //        dark with low roughness values.
                let reflection_direction = h.x*ax + h.y*ay + h.z*az;
                let reflection_ray = Ray::new(outwards_shifted_position(), reflection_direction);

                let cos_theta_reflection = reflection_direction.dot(nearest_hit.normal);

                let reflection = trace_radiance(&reflection_ray, scene, depth - 1);

                // Specular reflection:
                let view = -ray.direction;
                let light = reflection_direction;
                let normal = nearest_hit.normal;
                
                // @TODO: Does this have to be multiplied by PI or 2*PI?
                let l_spec = brdf_cook_torrance(view, light, normal, pbr_parameters)*reflection*cos_theta_reflection*PI;

                // Diffuse reflection:
                // @TODO: Does this have to be multiplied by PI or 2*PI?
                let l_diff = brdf_lambert(pbr_parameters)*reflection*cos_theta_reflection*PI;

                // @TODO: Make this energy conserving
                l_spec + l_diff
            },
        }
    } else {
        // @TODO: Make a nice gradient for the sky of sample an equirectangular projection.
        let theta = f32::acos(ray.direction.y);
        let t = f32::powf(theta / PI, 2.0);
        let intensity = 1.0 - 2.0*t;
        2.0*intensity*Vec3::new(0.6, 0.6, 0.8)
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

pub fn render(work_tile: WorkTile, backbuffer: &Arc<Backbuffer>, camera: Arc<RwLock<Camera>>, scene: Arc<RwLock<Scene>>) {
    let scene = scene.read().unwrap(); // @TODO: Handle the unwrap
    let camera = camera.read().unwrap(); // @TODO: Handle the unwrap

    let sampler = camera.sample(backbuffer.width, backbuffer.height);

    let (x0, x1) = (work_tile.position.x, work_tile.position.x + work_tile.size.x);
    let (y0, y1) = (work_tile.position.y, work_tile.position.y + work_tile.size.y);

    for y in y0..y1 {
        for x in x0..x1 {
            // Using different camera models has quite a substantial impact
            // in performance. This might go down as the scene description
            // gets more complicated. Test this from time to time and maybe
            // optimze.
            let ray = sampler.thin_lens_ray(x, y);

            let hdr_radiance = {
                const N: usize = 1;
                let mut hdr_radiance = Vec3::zero();
                for _ in 0..N {
                    hdr_radiance += trace_radiance(&ray, &*scene, 4);
                }
                hdr_radiance / N as f32
            };
            let ldr_radiance = tone_map_clamp(hdr_radiance);
            let color = Pixel::from_unit(ldr_radiance);
            backbuffer.set_pixel_unsafe(x, y, color);
        }
    }
}