#[macro_use] extern crate derive_new;

extern crate glium;
extern crate glutin;
extern crate hmath;
extern crate rand;
extern crate time;
extern crate notify;
extern crate stb_image;

mod common;
mod window;
mod scene;
mod worker;
mod tracer;
mod parser;
mod content;

fn main() {
    window::start();
}