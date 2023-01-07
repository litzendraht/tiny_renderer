use std::{time, collections::HashMap};
use std::fs::File;
use std::io::BufReader;

use show_image::{WindowOptions, ImageView, ImageInfo, event, create_window};
use obj::raw::parse_obj;
use nalgebra as na;
use na::vector;

use crate::scene::Scene;

const CAMERA_SPEED:       f32 = 3.0;
const LIGHT_SOURCE_SPEED: f32 = 3.0;

#[derive(PartialEq, Eq, Hash)]
enum Action {
    CameraLeft,
    CameraRight,
    LightLeft,
    LightRight,
    ExitApp,
}

pub struct Params {
    pub width:                u32,
    pub height:               u32,
    pub print_fps:            bool,
    pub asset_path:           String,
    pub shader_pipeline_name: String,
}

struct FrameActionBuffer {
    pub actions: HashMap<Action, bool>,
}

impl FrameActionBuffer {
    fn new() -> Self {
        return Self {
            actions: HashMap::from([
                (Action::CameraLeft,  false),
                (Action::CameraRight, false),
                (Action::LightLeft,   false),
                (Action::LightRight,  false),
                (Action::ExitApp,     false),
            ]),
        }
    }

    fn reset(&mut self) {
        for (_, value) in self.actions.iter_mut() {
            *value = false;
        }
    }

    fn process_window_event(&mut self, window_event: event::WindowEvent) {
        if let event::WindowEvent::KeyboardInput(event) = window_event {
            match (
                event.input.key_code, 
                event.input.state.is_pressed(), 
                event.input.state.is_released()
            ) {
                (Some(event::VirtualKeyCode::A),      true, _   ) => {
                    *self.actions.entry(Action::CameraLeft).or_insert(true) = true; 
                },
                (Some(event::VirtualKeyCode::D),      true, _   ) => {
                    *self.actions.entry(Action::CameraRight).or_insert(true) = true; 
                },
                (Some(event::VirtualKeyCode::Q),      true, _   ) => {
                    *self.actions.entry(Action::LightLeft).or_insert(true) = true; 
                },
                (Some(event::VirtualKeyCode::E),      true, _   ) => {
                    *self.actions.entry(Action::LightRight).or_insert(true) = true; 
                },
                (Some(event::VirtualKeyCode::Escape), _,    true) => {
                    *self.actions.entry(Action::ExitApp).or_insert(true) = true; 
                },
                _ => ()
            }
        }
    }
}

/// Actualy launches the window, showing images.
/// Takes struct, defining execution params.
pub fn run(params: Params) -> Result<(), Box<dyn std::error::Error>> {   
    let obj_path = params.asset_path.clone() + "/model.obj";
    let texture_path = params.asset_path.clone() + "/texture.tga";
    let normal_map_path = params.asset_path.clone() + "/normal_map.tga";
    let normal_map_tangent_path = params.asset_path.clone() + "/normal_map_tangent.tga";
    let specular_map_path = params.asset_path.clone() + "/specular_map.tga";

    println!("Loading model from: {}", obj_path);
    let obj = parse_obj(BufReader::new(File::open(obj_path)?))?;
    println!("Number of vertices in a model: {}", obj.positions.len());
    println!("Number of polygons in a model: {}", obj.polygons.len());

    println!("Loading texture from: {}", texture_path);
    let texture = image::open(texture_path)?.into_rgb8();
    println!("Dimensions of loaded texture are: {} x {}", texture.width(), texture.height());

    println!("Loading normal map from: {}", normal_map_path);
    let normal_map = image::open(normal_map_path)?.into_rgb8();
    println!("Dimensions of loaded normal map are: {} x {}", normal_map.width(), normal_map.height());

    println!("Loading normal map in tangent coordinates from: {}", normal_map_tangent_path);
    let normal_map_tangent = image::open(normal_map_tangent_path)?.into_rgb8();
    println!("Dimensions of loaded normal map in tangent coordinates are: {} x {}", normal_map.width(), normal_map.height());

    println!("Loading specular map from: {}", specular_map_path);
    let specular_map = image::open(specular_map_path)?.into_rgb8();
    println!("Dimensions of loaded specular map are: {} x {}", specular_map.width(), specular_map.height());

    println!("Cooking up a scene with '{}' shader pipeline", params.shader_pipeline_name);
    let mut scene = Scene::new(
        params.width, 
        params.height, 
        obj, 
        texture, 
        normal_map,
        normal_map_tangent, 
        specular_map,
        params.shader_pipeline_name
    );

    let window_options: WindowOptions = WindowOptions {
        size: Some([params.width, params.height]),
        ..Default::default()
    };
    let window = create_window("output", window_options)?;
    let event_channel = window.event_channel()?;

    // Buffer for tracking actionable window events.
    let mut frame_action_buffer = FrameActionBuffer::new();
    // Variables for convenience.
    let mut camera_angle: f32 = 0.0;
    let mut light_direction_angle: f32 = 0.0;
    // Stats.
    let mut exit = false;
    let mut frame_counter_time_begin = time::Instant::now();
    let mut frame_counter: u32 = 0;
    let mut frame_begin_time;
    let mut frame_time = 0.0;
    while !exit {
        frame_begin_time = time::Instant::now();

        // Clearing z-buffer and resetting rendered data to (0, 0, 0).
        scene.clear();        

        // Movement speed is proportional to previous frame dt for a smoother experience.
        if *frame_action_buffer.actions.get(&Action::CameraRight).unwrap() {
            camera_angle += CAMERA_SPEED * frame_time;
        }
        if *frame_action_buffer.actions.get(&Action::CameraLeft).unwrap() {
            camera_angle -= CAMERA_SPEED * frame_time;
        }
        // Direction is FROM surface TO source, so negative of true direction.
        // This simplifies math inside shaders somewhat by removing the need to place minus at some critical spots.
        // Easier to think of this as light source position on a unit sphere.
        if *frame_action_buffer.actions.get(&Action::LightRight).unwrap() {
            light_direction_angle += LIGHT_SOURCE_SPEED * frame_time;
        }
        if *frame_action_buffer.actions.get(&Action::LightLeft).unwrap() {
            light_direction_angle -= LIGHT_SOURCE_SPEED * frame_time;
        }
        let look_from       = vector![camera_angle.sin(), 0.0, camera_angle.cos()];
        let look_at         = vector![0.0, 0.0, 0.0];
        let up              = vector![0.0, 1.0, 0.0];        
        let light_direction = vector![light_direction_angle.sin(), 0.0, light_direction_angle.cos()];
        scene.set_light_direction(light_direction);
        scene.set_camera(look_from, look_at, up);
        scene.render();

        // Getting rendered data as a data slice and feeding it into window.
        let data = scene.get_frame_buffer();
        // let data = scene.get_z_buffer();
        // let data = scene.get_shadow_buffer();
        let image_view = ImageView::new(ImageInfo::rgb8(params.width, params.height), data.as_raw());
        window.set_image("image", image_view)?;

        // Unloading all the garbage from event channel, that has piled up, looking for actionable events.
        frame_action_buffer.reset();
        for window_event in event_channel.try_iter() {
            frame_action_buffer.process_window_event(window_event);
        }
        // If found Exit action, leaving the app in the beginning of the next frame.
        if *frame_action_buffer.actions.get(&Action::ExitApp).unwrap() {
            exit = true;
        }
        
        if params.print_fps {
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

        frame_time = time::Instant::now()
        .duration_since(frame_begin_time)
        .as_secs_f32();
    }

    return Ok(());
}