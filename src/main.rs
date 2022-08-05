extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

//pub mod fluidstuff;
//pub mod fluidstuffbetter;
pub mod bestfluid;

//use fluidstuff::fluidobj;
//use fluidstuffbetter::*;
use bestfluid::*;


use glutin_window::GlutinWindow as Window;
use graphics::color::TRANSPARENT;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::{Button, PressEvent, MouseButton};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;

//crate::fluidstuff::*;

pub struct App {
    gl: GlGraphics,
    rotation: f64,
    fob: b_fluobj,
}

impl  App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;
        //use fluidstuffbetter::N;
        //use fluidstuff::N;
        //use fluidstuff::IX;
        

        const GREEN: [f32; 4] = [0.0, 1.0, 0.0, 1.0];
        const RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        // /* 
        let size: f64 = 40_f64 * 16_f64 / (N as f64);
        let square = rectangle::square(0_f64, 0_f64, size);
        let rotation = self.rotation;
        let (x, y) = (args.window_size[0] / 2.0, args.window_size[1] / 2.0);
        self.gl.draw(args.viewport(), |c, gl| {
            clear(GREEN, gl);

            /*

            let tranform = c
                .transform
                .trans(x, y)
                .rot_rad(rotation)
                .trans(-25.0, -25.0);

            // */
            for i in 1..=N {
                for j in 1..=N {
                    let tranform = c
                        .transform
                        .trans( size+ 0_f64 + (i as f64 * args.window_size[0] / (N+2) as f64) as f64,  size+ 0_f64 + (j as f64 * args.window_size[1] / (N+2) as f64) as f64)
                        .rot_rad(rotation)
                        .trans(-(size/2_f64), -(size/2_f64));

                    let transform2 = c
                        .transform
                        .trans(size + 1_f64 * i as f64 + (i as f64 * (size / 1_f64)), size + 1_f64 * j as f64 + (j as f64 * (size / 1_f64)) )
                        .trans(-(size/2_f64), -(size/2_f64));

                    let val = ((i*j) as f32) / (1 as f32);
                    let val2 = self.fob.x[IX(i,j)];
                    rectangle([1.0, 0.0, 0.0, (0.0 + val2 as f32) * 800 as f32], square, transform2, gl);
                }
            }
        });
        // */


    }
    fn update(&mut self, args: &UpdateArgs) {
        let dt: f64 = args.dt / 10_f64;
        //self.rotation += 2.0 * args.dt;
        //self.fob.clear_prev_vals();
        self.fob.vel_step_self(dt);
        self.fob.dens_step_self(dt);
        self.fob.clear_prev_vals();
        //print_matrix(self.fob.x.clone());
    }
}

fn main() {
    //use crate::fluidstuff::Nex;
    let opengl = OpenGL::V3_2;

    let mut window: Window = WindowSettings::new("spinning square", [1000, 1000])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();
    let mut app = App {
        gl: GlGraphics::new(opengl),
        rotation: 0.0,
        fob: b_fluobj::new(N),
    };

    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
        if let Some(Button::Mouse(button)) = e.press_args() {
            if button == MouseButton::Left {
                app.fob.x[150] += 1_f64;
            }
            else if button == MouseButton::Right {
                app.fob.x[150] -= 1_f64;
            }
            
        }
    }

}
