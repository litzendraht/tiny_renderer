use std::time;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use show_image::{WindowOptions, ImageView, ImageInfo, event, create_window};
use obj::raw::{parse_obj};
use nalgebra as na;
use na::vector;

use crate::scene::Scene;

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
    let model_path = context.asset_path.clone() + "/model.obj";
    let texture_path = context.asset_path.clone() + "/texture.tga";

    println!("Loading model from: {}", model_path);
    let model = parse_obj(BufReader::new(File::open(model_path)?))?;
    println!("Number of vertices in a model: {}", model.positions.len());
    println!("Number of polygons in a model: {}", model.polygons.len());

    println!("Loading texture from: {}", texture_path);
    let texture = image::open(texture_path)?.into_rgb8();
    println!("Dimensions of loaded texture are: {} x {}", texture.width(), texture.height());

    let mut scene = Scene::new(context.width, context.height, model, texture);

    let window_options: WindowOptions = WindowOptions {
        size: Some([context.width, context.height]),
        ..Default::default()
    };
    let window = create_window("output", window_options)?;
    let event_channel = window.event_channel()?;

    // Stats.
    let mut exit = false;
    let time_begin = time::Instant::now();
    let mut frame_counter_time_begin = time::Instant::now();
    let mut frame_counter: u32 = 0;
    while !exit {
        let passed_time = time::Instant::now()
        .duration_since(time_begin)
        .as_secs_f32();

        // Clearing z-buffer and resetting rendered data to (0, 0, 0).
        scene.clear();        

        // Setting up camera position and direction.
        let look_from = vector![1.0 * passed_time.sin(), 0., 1.0 * passed_time.cos()];
        let look_at = vector![0., 0., 0.];
        let up = vector![0., 1., 0.];

        scene.prepare_camera(look_from, look_at, up);
        scene.render();

        let data = scene.render_data();
        let image_view = ImageView::new(ImageInfo::rgb8(context.width, context.height), data.as_raw());
        window.set_image("image", image_view)?;

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