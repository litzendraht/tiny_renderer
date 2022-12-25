use std::time;
use show_image::{ImageView, ImageInfo, event, create_window};
use crate::image::{Image, Color, Coord};

const WHITE: Color = Color { r: 255, g: 255, b: 255, };
const BLACK: Color = Color { r: 0,   g: 0,   b: 0,   };
const RED: Color =   Color { r: 255, g: 0,   b: 0,   };
const GREEN: Color = Color { r: 0,   g: 255, b: 0,   };
const BLUE: Color =  Color { r: 0,   g: 0,   b: 255, };

pub struct Context {
    pub width: u32,
    pub height: u32,
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
    let mut image = Image::new(context.width, context.height);
    
    let window = create_window("output", Default::default())?;
    let event_channel = window.event_channel()?;

    let mut exit = false;
    let time_begin = time::Instant::now();
    while !exit {
        let passed_time = time::Instant::now()
        .duration_since(time_begin)
        .as_secs_f32();

        for i in 0..context.width {
            for j in 0..context.height {
                image.set_pixel(
                    Coord {x: i, y: j}, 
                    Color::blend(
                        WHITE, 
                        BLACK, 
                        (10.0 * ((i as f32 / context.width as f32) + passed_time)).sin()
                    )
                )?;
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
        }      
    }

    return Ok(());
}