use std::time;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use image;

use show_image::{WindowOptions, ImageView, ImageInfo, event, create_window};
use obj::raw::object::Polygon;
use obj::raw::{parse_obj};

use crate::util::{Vector2i, Vector2f, Vector3f};
use crate::scene::{Scene, Color};

// @TRASH double definition
const WHITE: Color = Color { r: 255, g: 255, b: 255, };
const BLACK: Color = Color { r: 0,   g: 0,   b: 0,   };
const RED: Color =   Color { r: 255, g: 0,   b: 0,   };
const GREEN: Color = Color { r: 0,   g: 255, b: 0,   };
const BLUE: Color =  Color { r: 0,   g: 0,   b: 255, };

// @TODO redo asset_path to be an actual Path object somehow
pub struct Context {
    pub width: u32,
    pub height: u32,
    pub print_fps: bool,
    pub asset_path: String,
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

/// Actualy launches the window, showing images.
/// Takes struct, defining execution context.
pub fn run(context: Context) -> Result<(), Box<dyn std::error::Error>>{
    let mut scene = Scene::new(context.width, context.height);
    
    let model_path = context.asset_path.clone() + "/model.obj";
    let texture_path = context.asset_path.clone() + "/texture.tga";

    println!("Loading model from: {}", model_path);
    let model = parse_obj(BufReader::new(File::open(model_path)?))?;
    println!("Number of vertices in a model: {}", model.positions.len());
    println!("Number of polygons in a model: {}", model.polygons.len());

    println!("Loading texture from: {}", texture_path);
    let texture = image::open(texture_path)?.into_rgb8();
    println!("Dimensions of loaded texture are: {} x {}", texture.width(), texture.height());

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
        // Clearing z-buffer and resetting rendered data to (0, 0, 0).
        scene.clear();

        let passed_time = time::Instant::now()
        .duration_since(time_begin)
        .as_secs_f32();

        // Drawing all polygons of the model.
        for polygon in model.polygons.iter() {
            let indices: &Vec<(usize, usize, usize)> = match polygon {
                Polygon::PTN(indices) => indices,
                _ => panic!("Encountered some garbage, while looking through polygons."),
            };
            let a = Vector3f {
                x: model.positions.get(indices[0].0).unwrap().0,
                y: model.positions.get(indices[0].0).unwrap().1,
                z: model.positions.get(indices[0].0).unwrap().2,
            };
            let b = Vector3f {
                x: model.positions.get(indices[1].0).unwrap().0,
                y: model.positions.get(indices[1].0).unwrap().1,
                z: model.positions.get(indices[1].0).unwrap().2,
            };
            let c = Vector3f {
                x: model.positions.get(indices[2].0).unwrap().0,
                y: model.positions.get(indices[2].0).unwrap().1,
                z: model.positions.get(indices[2].0).unwrap().2,
            };

            // Calculating normal projection on the face.
            let light_dir = Vector3f { x: 0.0, y: 0.0, z: 1.0 };  // Directed to us from the screen.
            let face_normal = Vector3f::cross(b - a, c - a);
            let mut normal_correction_coef = Vector3f::dot(light_dir, face_normal);            
            // Backface culling.
            if normal_correction_coef > 0.0 {
                normal_correction_coef /= face_normal.norm();
                // @TODO figure out, whu texture is upside down, lol? Why do I need to do 1 - y?
                let uv_a = Vector2f {
                    x: model.tex_coords.get(indices[0].1).unwrap().0,
                    y: 1.0 - model.tex_coords.get(indices[0].1).unwrap().1,
                };
                let uv_b = Vector2f {
                    x: model.tex_coords.get(indices[1].1).unwrap().0,
                    y: 1.0 - model.tex_coords.get(indices[1].1).unwrap().1,
                };
                let uv_c = Vector2f {
                    x: model.tex_coords.get(indices[2].1).unwrap().0,
                    y: 1.0 - model.tex_coords.get(indices[2].1).unwrap().1,
                };
                // scene.draw_triangle(a, b, c, WHITE, normal_correction_coef);
                scene.draw_triangle_textured(a, b, c, uv_a, uv_b, uv_c, &texture, normal_correction_coef);
            }
        }

        let image_data = ImageView::new(ImageInfo::rgb8(context.width, context.height), scene.as_render_data());
        // let image_data = ImageView::new(ImageInfo::rgb8(context.width, context.height), scene.as_depth_data());
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