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

fn main() {
    let mut events_loop = winit::EventsLoop::new();
    let window = winit::Window::new(&events_loop).unwrap();
    
    let h_wnd = window.get_hwnd();
    let hdc = unsafe { GetDC(h_wnd) };

    // 0x00bbggrr
    let color: COLORREF = 0x000000FF;

    for i in 0..200 {
        unsafe {
            SetPixel(hdc, 25 + i, 25, color);
        }
    }

    events_loop.run_forever(|event| {
        match event {
            winit::Event::WindowEvent {
              event: winit::WindowEvent::CloseRequested,
              ..
            } => winit::ControlFlow::Break,
            _ => winit::ControlFlow::Continue,
        }
    });
}