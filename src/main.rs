extern crate winit;
extern crate libc;

use winit::os::windows::WindowExt;

type COLORREF = libc::c_uint;
type HDC = *mut libc::c_void;
type HWND = *mut libc::c_void;

// Only windows is supported right now.
//#[cfg(all(target_os = "win32", target_arch = "x86"))]
#[link(name = "kernel32")]
#[allow(non_snake_case)]
extern "stdcall" {
    fn SetPixel(hdc: HDC, x: libc::c_int, y: libc::c_int, color: COLORREF) -> COLORREF;
    fn GetDC(hWnd: HWND) -> HDC;
}

struct Backbuffer {
    width: usize,
    height: usize,
    hdc: HDC,
}

fn main() {
    let mut events_loop = winit::EventsLoop::new();
    let window = winit::Window::new(&events_loop).unwrap();
    
    let h_wnd = window.get_hwnd();
    let hdc = unsafe { GetDC(h_wnd) };

    let size = window.get_inner_size().unwrap();

    let backbuffer = Backbuffer {
        width: size.width as usize, // from f64
        height: size.height as usize, // from f64
        hdc: hdc,
    };

    events_loop.run_forever(|event| {
        render(&backbuffer);

        match event {
            winit::Event::WindowEvent {
              event: winit::WindowEvent::CloseRequested,
              ..
            } => winit::ControlFlow::Break,
            _ => winit::ControlFlow::Continue,
        }
    });
}

fn render(backbuffer: &Backbuffer) {
    for y in 0..backbuffer.height as i32 {
        for x in 0..backbuffer.width as i32 {
            // 0x00bbggrr
            let color: COLORREF = 0x000000FF;

            unsafe {
                SetPixel(backbuffer.hdc, x, y, color);
            }
        }
    }
}