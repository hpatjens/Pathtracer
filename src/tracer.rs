use common::*;

use scene::{Scene, Sky, PBRParameters, Material, find_scene_hit};

use std::sync::{Arc, RwLock};
use std::cell::UnsafeCell;

#[derive(Clone, Debug, new)]
pub struct WorkTile {
    pub tile_index: Vec2u,
    pub position: Vec2u,
    pub size: Vec2u,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Transition { In, Out }

#[derive(Debug, Clone, new)]
pub struct Hit<'a> {
    pub parameter: f32,
    pub position: Vec3,
    pub normal: Vec3,
    pub material: &'a Material,
    pub transition: Transition,
}

#[derive(Debug, Clone, new)]
pub struct Meta {
    pixel_position: Vec2u,
    sample_index: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct Basis(pub Vec3, pub Vec3, pub Vec3);

#[derive(Clone, Debug)]
pub enum ToneMapping {
    Clamp,
    Reinhard,
    Exposure(f32),
}

#[derive(Clone, Debug, new)]
struct CameraPlane {
    origin: Vec3,
    u: Vec3,
    v: Vec3,
}

#[derive(Debug)]
pub struct Camera {
    projection_plane: CameraPlane,
    position: Vec3,
    width: f32,
    height: f32,
    z_near: f32,
    tone_mapping: ToneMapping,
    iso: f32,
}

impl Camera {
    pub fn new(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32, tone_mapping: ToneMapping, iso: f32) -> Self {
        Camera {
            projection_plane: Self::construct_projection_plane(position, target, up, width, height, z_near),
            position: position,
            width: width,
            height: height,
            z_near: z_near,
            tone_mapping: tone_mapping,
            iso: iso,
        }
    }

    fn construct_projection_plane(position: Vec3, target: Vec3, up: Vec3, width: f32, height: f32, z_near: f32) -> CameraPlane {
        let z = (position - target).normalize();
        let x = up.cross(z).normalize();
        let y = z.cross(x).normalize();

        let half_width = width / 2.0;
        let half_height = height / 2.0;

        let to_plane_center = -z_near*z;

        let origin = position + to_plane_center - x*half_width - y*half_height;

        let u = x*width;
        let v = y*height;
        CameraPlane::new(origin, u, v)
    }

    #[allow(dead_code)]
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
    #[allow(dead_code)]
    fn pinhole_ray(&self, x: u32, y: u32) -> Ray {
        let origin = {
            let ru = random32();
            let rv = random32();

            let du = (x as f32 + ru)*self.u;
            let dv = (y as f32 + rv)*self.v;

            self.camera.projection_plane.origin + du + dv
        };
        let direction = (origin - self.camera.position).normalize();
        Ray::new(origin, direction)
    }

    #[allow(dead_code)]
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
    pub num_samples: UnsafeCell<usize>,
    pub pixels32: UnsafeCell<Vec<Pixel32>>,
    pub pixels8: UnsafeCell<Vec<Pixel8>>,
}
// UnsafeCell does not implement Sync and therefore Backbuffer could not be passed to the
// worker threads without this implementation.
unsafe impl Sync for Backbuffer {}

impl Backbuffer {
    pub fn new(width: u32, height: u32) -> Backbuffer {
        let num_pixels = (width * height) as usize;
        Backbuffer {
            width: width,
            height: height,
            num_samples: UnsafeCell::new(1),
            pixels32: {
                let mut pixels = Vec::new();
                pixels.resize(num_pixels, Pixel32(0, 0, 0));
                UnsafeCell::new(pixels)
            },
            pixels8: {
                let mut pixels = Vec::new();
                pixels.resize(num_pixels, Pixel8(0, 0, 0));
                UnsafeCell::new(pixels)
            },

        }
    }

    pub fn clear(&self) {
        unsafe {
            let ref mut pixels32 = *self.pixels32.get();
            for i in 0..pixels32.len() {
                pixels32[i] = Pixel32(0, 0, 0);
            }

            let ref mut pixels8 = *self.pixels8.get();
            for i in 0..pixels8.len() {
                pixels8[i] = Pixel8(0, 0, 0);
            }

            *self.num_samples.get() = 1;
        }
    }

    fn add_pixel32_unsafe(&self, x: u32, y: u32, pixel: Pixel32) {
        let index = (y*self.width + x) as usize;
        unsafe {
            let old_pixel = (*self.pixels32.get())[index].clone();
            let new_pixel = Pixel32(
                old_pixel.0 + pixel.0,
                old_pixel.1 + pixel.1,
                old_pixel.2 + pixel.2,
            );
            (*self.pixels32.get())[index] = new_pixel;
        }
    }

    fn assign_pixel8_unsafe(&self, x: u32, y: u32, iso_factor: f32) {
        let index = (y*self.width + x) as usize;
        unsafe {
            let num_samples = *self.num_samples.get() as u32;
            let ref pixel32 = (*self.pixels32.get())[index];
            (*self.pixels8.get())[index] = Pixel8(
                f32::min(255.0, iso_factor*pixel32.0 as f32 / num_samples as f32) as u8,
                f32::min(255.0, iso_factor*pixel32.1 as f32 / num_samples as f32) as u8,
                f32::min(255.0, iso_factor*pixel32.2 as f32 / num_samples as f32) as u8,
            );
        }
    }
}

fn reflect(incoming: Vec3, normal: Vec3) -> Vec3 {
    // Bear in mind that incoming is directed at the surface!
    incoming - 2.0*incoming.dot(normal)*normal
}

fn refract(incoming: Vec3, normal: Vec3, n1: f32, n2: f32) -> Vec3 {
    // Bear in mind that incoming is directed at the surface!

    let v = -incoming;
    let n = normal;

    // Pythagorean trigonometric identity: sin^2(a) + cos^2(a) = 1
    // Since cos(a) is easy to compute, sin(a) can be computed by
    // sin(a) = sqrt(1 - cos^2(a))
    let cos_theta1 = v.dot(n);

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
    let cos_theta2 = f32::sqrt(1.0 - sin_theta2*sin_theta2);
    
    // Corresponds to the projection of 'incoming' onto the ground plane.
    let p = incoming + cos_theta1*normal;
    
    -normal*cos_theta2 + p.normalize()*sin_theta2
}

fn construct_coordinate_system(normal: Vec3) -> Basis {
    const EPS: f32 = 0.9999;
    let other = if normal.y > EPS || normal.y < -EPS {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let y = normal;
    let x = y.cross(other);
    let z = y.cross(x);
    Basis(x, y, z)
}

#[allow(dead_code)]
fn brdf_lambert(pbr_parameters: &PBRParameters) -> Vec3 {
    let &PBRParameters{ ref reflectivity, .. } = pbr_parameters;
    reflectivity/PI
}

#[allow(dead_code)]
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

fn brdf_cook_torrance(view: Vec3, light: Vec3, normal: Vec3, reflectivity: Vec3, roughness: f32, f0: Vec3) -> Vec3 {
    // Bear in mind that 'view' as well as 'light' are pointing away from the surface!

    let v = view;
    let l = light;
    let n = normal;

    let alpha = roughness*roughness;
    let h = (v + l).normalize();
    let cos_theta = f32::max(0.0, h.dot(v));

    // @TODO: There are multiple definitions for this. What's the best?
    // let k = (alpha + 1.0)*(alpha + 1.0) / 8.0;
    let k = alpha / 2.0; // Was/is used in the Unreal Engine 4 according to Brian Karis' blog.

    // The normal distribution function is canceled out as it functions as the probability 
    // density function for the monte carlo integration.
    let num = fresnel_schlick(cos_theta, f0)*geometry_smith(n, v, l, k);
    let denum = 4.0*f32::max(0.0, n.dot(v));
    (num / denum)
}

// @TODO: Ensure that the reflection direction is not below the horizon.
#[allow(dead_code)]
fn importance_sample_ggx(xi: Vec2, roughness: f32) -> Vec3 {
    let a = roughness*roughness;

    let phi = 2.0*PI*xi.x;
    let sin_phi = f32::sin(phi);
    let cos_phi = f32::cos(phi);

    let cos_theta = f32::sqrt((1.0 - xi.y) / (xi.y*(a*a - 1.0) + 1.0));
    let sin_theta = f32::sqrt(1.0 - cos_theta*cos_theta);

    Vec3::new(sin_theta*cos_phi, cos_theta, sin_theta*sin_phi)
}

#[allow(dead_code)]
fn importance_sample_cos(xi: Vec2) -> Vec3 {
    let r = f32::sqrt(xi.x);
    let theta = 2.0*PI*xi.y;
 
    let x = r * f32::cos(theta);
    let z = r * f32::sin(theta);
 
    Vec3::new(x, f32::sqrt(f32::max(0.0, 1.0 - xi.x)), z)
}

pub fn to_basis(basis: Basis, v: Vec3) -> Vec3 {
    let Basis(x, y, z) = basis;
    x*v.x + y*v.y + z*v.z
}

fn trace_radiance(meta: &Meta, ray: &Ray, scene: &Scene, depth: u8) -> Vec3 {
    // @TODO: Find all the places where NANs can be generated and fix as many as it makes sense.

    if depth <= 0 {
        let num_light_sources = scene.emissive_spheres.len() + scene.emissive_planes.len();        
        let r = xorshift32() as usize % num_light_sources;

        let (sample_point, light_radiance) = if r < scene.emissive_spheres.len() {
            let ref light = scene.emissive_spheres[r];
            let radiance = match light.material {
                Material::Emissive(radiance) => radiance,
                _ => panic!("Non-emissive material in emissive material Vec."),
            };
            (light.sample(ray.origin), radiance)
        } else {
            let ref light = scene.emissive_planes[r];
            let radiance = match light.material {
                Material::Emissive(radiance) => radiance,
                _ => panic!("Non-emissive material in emissive material Vec."),
            };
            (light.sample(ray.origin), radiance)
        };

        let sample_point = if let Some(point) = sample_point {
            point
        } else {
            return Vec3::zero()
        };

        let light_ray = sample_point - ray.origin;
        let light_distance = light_ray.length();

        let ray = Ray::new(ray.origin, light_ray.normalize());
        return if let Some(nearest_hit) = find_scene_hit(&ray, scene) { // @TODO: Make a function for finding any hit.
            // @TODO: This should not be done by distance but by object id
            if nearest_hit.parameter > light_distance - 0.0001 {
                light_radiance
            } else {
                Vec3::zero()
            }
        } else {
            Vec3::zero()
        }
    }

    let nearest_hit = find_scene_hit(ray, scene);

    if let Some(nearest_hit) = nearest_hit {
        const SHIFT_AMOUNT: f32 = 0.0001; // @TODO: Find a good factor and maybe make it dependent on the slope
        let outwards_shifted_position = ||{ nearest_hit.position + SHIFT_AMOUNT*nearest_hit.normal };
        let inwards_shifted_position  = ||{ nearest_hit.position - SHIFT_AMOUNT*nearest_hit.normal };

        const R: f32 = 0.04;

        match nearest_hit.material {
            Material::None => Vec3::one(),
            Material::Emissive(ref color) => {
                color.clone()
            },
            Material::Mirror => {
                let reflection_direction = reflect(ray.direction, nearest_hit.normal);
                trace_radiance(meta, &Ray::new(outwards_shifted_position(), reflection_direction), scene, depth - 1)
            },
            Material::Translucent(ior) => {
                let cos_theta = nearest_hit.normal.dot(-ray.direction);
                let fresnel = fresnel_schlick(cos_theta, Vec3::new(R, R, R));

                const IOR_AIR: f32 = 1.0;
                let (n1, n2) = match nearest_hit.transition {
                    Transition::In  => (IOR_AIR, *ior),
                    Transition::Out => (*ior, IOR_AIR),
                };

                let k_refl = fresnel;
                let reflection_direction = reflect(ray.direction, nearest_hit.normal);
                let l_refl = trace_radiance(meta, &Ray::new(outwards_shifted_position(), reflection_direction), scene, depth - 1);

                let k_refr = Vec3::one() - fresnel;
                let refraction_direction = refract(ray.direction, nearest_hit.normal, n1, n2);                
                let l_refr = trace_radiance(meta, &Ray::new(inwards_shifted_position(), refraction_direction), scene, depth - 1);

                k_refl*l_refl + k_refr*l_refr
            },
            Material::Physically(ref pbr_parameters) => {
                let tangent_space = construct_coordinate_system(nearest_hit.normal);

                let view = -ray.direction;
                let normal = nearest_hit.normal;

                let PBRParameters{ reflectivity, roughness, metalness } = *pbr_parameters;
                let f0 = mix_vec3(Vec3::new(R, R, R), reflectivity, metalness);

                // Decide whether this reflection is specular or diffuse
                let specular_reflection = random32() < f0.as_array()[(xorshift32() % 3) as usize];

                if specular_reflection {
                    let xi = Vec2::new(random32(), random32());

                    // Half vector for the reflection
                    let h = to_basis(tangent_space, importance_sample_ggx(xi, pbr_parameters.roughness));

                    // The reflection direction is called light. Not to be confused  with a ray towards a light source.
                    let light = reflect(ray.direction, h);
                   
                    let light_ray = Ray::new(outwards_shifted_position(), light);
                    let light_cos_theta = light.dot(normal);

                    let light_radiance = trace_radiance(meta, &light_ray, scene, depth - 1);

                    // The cos(theta) was canceled out as it is in the denominator of the Cook-Torrance BRDF.
                    PI*brdf_cook_torrance(view, light, normal, reflectivity, roughness, f0)*light_radiance
                } else {
                    let xi = Vec2::new(random32(), random32());

                    let light = to_basis(tangent_space, importance_sample_cos(xi));
                    let light_ray = Ray::new(outwards_shifted_position(), light);
                    let light_cos_theta = light.dot(normal);

                    let light_radiance = trace_radiance(meta, &light_ray, scene, depth - 1);

                    PI*brdf_lambert(pbr_parameters)*light_radiance*light_cos_theta
                }
            },
        }
    } else {
        match &scene.sky {
            &Sky::Constant(radiance) => radiance,
            &Sky::HDRI(ref _path, ref option_texture) => {
                if let Some(texture) = option_texture {
                    let sampler = texture.sampler();
                    sampler.sample_equirectangular(ray.direction)
                } else {
                    Vec3::one() // Returning white until the sky is loaded
                }
            },
        }
    }
}

#[allow(dead_code)]
fn tone_map_reinhard(radiance: Vec3) -> Vec3 {
    radiance / (Vec3::one() + radiance)
}

#[allow(dead_code)]
fn tone_map_clamp(radiance: Vec3) -> Vec3 {
    Vec3::new(
        saturatef32(radiance.x),
        saturatef32(radiance.y),
        saturatef32(radiance.z),
    )
}

#[allow(dead_code)]
fn tone_map_exposure(radiance: Vec3, exposure: f32) -> Vec3 {
    Vec3::new(
        1.0 - f32::exp(-exposure*radiance.x),
        1.0 - f32::exp(-exposure*radiance.y),
        1.0 - f32::exp(-exposure*radiance.z),
    )
}

fn gamma_correction(radiance: Vec3) -> Vec3 {
    const GAMMA: f32 = 1.0/2.2;
    Vec3::new(
        f32::powf(radiance.x, GAMMA),
        f32::powf(radiance.y, GAMMA),
        f32::powf(radiance.z, GAMMA),
    )
}

pub fn render(work_tile: WorkTile, backbuffer: &Arc<Backbuffer>, scene: Arc<RwLock<Scene>>) {
    let scene = scene.read().unwrap(); // @TODO: Handle the unwrap

    let iso_factor = scene.camera.iso/100.0;

    let sampler = scene.camera.sample(backbuffer.width, backbuffer.height);

    let (x0, x1) = (work_tile.position.x, work_tile.position.x + work_tile.size.x);
    let (y0, y1) = (work_tile.position.y, work_tile.position.y + work_tile.size.y);

    for y in y0..y1 {
        for x in x0..x1 {
            let ray = sampler.pinhole_ray(x, y);
            let meta = Meta::new(Vec2u::new(x, y), 1);

            let hdr_radiance = trace_radiance(&meta, &ray, &*scene, 4);
            let ldr_radiance = match scene.camera.tone_mapping {
                ToneMapping::Clamp => tone_map_clamp(hdr_radiance),
                ToneMapping::Reinhard => tone_map_reinhard(hdr_radiance),
                ToneMapping::Exposure(value) => tone_map_exposure(hdr_radiance, value),
            };

            let gamma_corrected = gamma_correction(ldr_radiance);
            let color = Pixel32::from_unit(gamma_corrected);

            backbuffer.add_pixel32_unsafe(x, y, color);
            backbuffer.assign_pixel8_unsafe(x, y, iso_factor);
        }
    }
}
