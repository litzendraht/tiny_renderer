use std::time;
use std::fs::File;
use std::io::BufReader;

use show_image::{WindowOptions, ImageView, ImageInfo, event, create_window};
use obj::{load_obj, Obj};

use crate::util::{Coord, Vector2, Vector3};
use crate::image::{Image, Color};

const WHITE: Color = Color { r: 255, g: 255, b: 255, };
const BLACK: Color = Color { r: 0,   g: 0,   b: 0,   };
const RED: Color =   Color { r: 255, g: 0,   b: 0,   };
const GREEN: Color = Color { r: 0,   g: 255, b: 0,   };
const BLUE: Color =  Color { r: 0,   g: 0,   b: 255, };

pub struct Context {
    pub width: u32,
    pub height: u32,
    pub print_fps: bool,
    pub filename: String,
}

/// Helper, defining exit event to be an Escape key press.
fn is_exit_event(window_event: event::WindowEvent) -> bool {
    if let event::WindowEvent::KeyboardInput(event) = window_event {
        // println!("{:#?}", event);
        if event.input.key_code == Some(event::VirtualKeyCode::Escape) && event.input.state.is_released() {
            return true;
        }
    }

    return false;
}

/// Transformtaion of a float in \[-1.0, 1.0\] to the pixel image coordinate
/// in range \[0, length - 1\].
fn to_rasterized_coord(float_coord: f32, scale: u32) -> i32 {
    return ((float_coord + 1.0) * ((scale - 1) as f32) / 2.0) as i32;
}

/// Actualy launches the window, showing images.
/// Takes struct, defining execution context.
pub fn run(context: Context) -> Result<(), Box<dyn std::error::Error>>{
    let mut image = Image::new(context.width, context.height);
    
    let model: Obj = load_obj(BufReader::new(File::open(context.filename)?))?;
    println!("Number of vertices - {}", model.vertices.len());
    println!("Number of indices  - {}", model.indices.len());

    let window_options: WindowOptions = WindowOptions {
        size: Some([context.width, context.height]),
        ..Default::default()
    };
    let window = create_window("output", window_options)?;
    let event_channel = window.event_channel()?;

    let mut exit = false;
    let time_begin = time::Instant::now();
    let mut frame_counter_time_begin = time::Instant::now();
    let mut frame_counter: u32 = 0;
    while !exit {
        let passed_time = time::Instant::now()
        .duration_since(time_begin)
        .as_secs_f32();

        // Drawing all faces of the model.
        for i in 0..model.indices.len() / 3 {
            // @TRASH can't quite believe that I need to explicitly deref and cast to usize to just
            // get an element of a vec.
            let coord_a_index = *model.indices.get(3 * i + 0).unwrap() as usize;
            let coord_b_index = *model.indices.get(3 * i + 1).unwrap() as usize;
            let coord_c_index = *model.indices.get(3 * i + 2).unwrap() as usize;
            let a = Vector3 {
                x: model.vertices.get(coord_a_index).unwrap().position[0],
                y: model.vertices.get(coord_a_index).unwrap().position[1],
                z: model.vertices.get(coord_a_index).unwrap().position[2],
            };
            let b = Vector3 {
                x: model.vertices.get(coord_b_index).unwrap().position[0],
                y: model.vertices.get(coord_b_index).unwrap().position[1],
                z: model.vertices.get(coord_b_index).unwrap().position[2],
            };
            let c = Vector3 {
                x: model.vertices.get(coord_c_index).unwrap().position[0],
                y: model.vertices.get(coord_c_index).unwrap().position[1],
                z: model.vertices.get(coord_c_index).unwrap().position[2],
            };            
            let coord_a = Coord { 
                x: to_rasterized_coord(a.x, context.width),
                y: to_rasterized_coord(a.y, context.height),
            };
            let coord_b = Coord { 
                x: to_rasterized_coord(b.x, context.width),
                y: to_rasterized_coord(b.y, context.height),
            };
            let coord_c = Coord { 
                x: to_rasterized_coord(c.x, context.width),
                y: to_rasterized_coord(c.y, context.height),
            };

            // Calculating normal projection on the face.
            let light_dir = Vector3 { x: 0.0, y: 0.0, z: 1.0 };  // Directed out of the screen(?).
            let face_normal = Vector3::cross(b - a, c - a);
            let mut normal_correction_coef = Vector3::dot(light_dir, face_normal);            
            // Backface culling.
            if normal_correction_coef > 0.0 {
                normal_correction_coef /= face_normal.norm();
                let normal_corrected_color = Color::blend(WHITE, BLACK, normal_correction_coef);
                image.draw_triangle(coord_a, coord_b, coord_c, normal_corrected_color);
            }
        }

        let image_data = ImageView::new(ImageInfo::rgb8(context.width, context.height), image.as_pixel_data());
        window.set_image("image", image_data)?;

        // Unloading all the garbage from event channel, that has piled up, looking for exit event.
        let exit_poll_result = event_channel.try_iter()
        .map(|window_event| is_exit_event(window_event))
        .reduce(|was_exit_event, is_exit_event| was_exit_event || is_exit_event);

        // If any event is Escape key press, then exiting.
        exit = match exit_poll_result {
            Some(value) => value,
            None => false,
        };
        
        if context.print_fps {
            // Counting frames to printout stats every seconds.
            frame_counter += 1;
            if time::Instant::now()
            .duration_since(frame_counter_time_begin)
            .as_secs_f32() > 1.0 {
                println!("FPS --- {}", frame_counter);
                frame_counter_time_begin = time::Instant::now();
                frame_counter = 0;
            }
        }
    }

    return Ok(());
}