use std::sync::mpsc::{channel, Sender};
use std::thread::JoinHandle;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use std::thread;

use scene::HDRITexture;

enum AssetRequest {
    HDRITexture(String),
}

fn load_hdri_texture(path: String) -> HDRITexture {
    use stb_image::image;
    use stb_image::image::LoadResult;

    match image::load(&*path) {
        LoadResult::ImageU8(_) => panic!("Image \"{}\" is not an HDR image.", path),
        LoadResult::ImageF32(image) => {
            // @TODO: This should not be an assert.
            assert!(image.depth == 3);
            HDRITexture::new(image.data, image.width, image.height)
        },
        LoadResult::Error(message) => panic!("Could not load the image \"{}\"! Error: {}", path, message),
    }
}

pub enum AssetState<T> {
    Loading,
    Available(T),
}

pub struct Content {
    hdri_textures: Arc<Mutex<BTreeMap<String, AssetState<Arc<HDRITexture>>>>>,
    sender: Sender<AssetRequest>,
    #[allow(dead_code)]
    thread: JoinHandle<()>,
}

impl Content {
    pub fn new() -> Self {
        let (sender, receiver) = channel();

        let hdri_textures = Arc::new(Mutex::new(BTreeMap::new()));
        let hdri_textures2 = hdri_textures.clone();

        // @TODO: This thread should be joined.
        let thread = thread::spawn(move || {
            loop {
                match receiver.recv() {
                    Ok(asset_request) => {
                        match asset_request {
                            AssetRequest::HDRITexture(path) => {
                                let hdri_texture = load_hdri_texture(path.clone());
                                let mut hdri_textures = hdri_textures2.lock().expect("Could not aquire the lock for hdri_textures in Content.");
                                hdri_textures.insert(path, AssetState::Available(Arc::new(hdri_texture)));
                            },
                        }
                    },
                    Err(_) => (),
                }
            }
        });
        Content {
            hdri_textures: hdri_textures,
            thread: thread,
            sender: sender,
        }
    }

    pub fn get_hdri_texture(&self, path: &str) -> Option<Arc<HDRITexture>> {
        let mut hdri_textures = self.hdri_textures.lock().expect("Could not aquire the lock for hdri_textures in Content.");
        let (result, request_load) = match hdri_textures.get(path) {
            Some(AssetState::Loading) => (None, false),
            Some(AssetState::Available(hdri_texture)) => (Some(hdri_texture.clone()), false),
            None => (None, true),
        };
        if request_load {
            // @TODO: Handle the unwrap
            self.sender.send(AssetRequest::HDRITexture(String::from(path))).unwrap();
            hdri_textures.insert(String::from(path), AssetState::Loading);
        }
        result
    }
}