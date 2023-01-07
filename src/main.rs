mod scene;
mod app;

use std::env;

const WIDTH: u32  = 800;
const HEIGHT: u32 = 800;

#[show_image::main]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Default values.
    let mut asset_path = String::from("assets/diablo");
    let mut shader_pipeline_name = String::from("default");

    let args: Vec<String> = env::args().collect();
    for i in 1..args.len() {
        match args[i].as_str() {
            "-p" => { asset_path = args[i + 1].clone(); }
            "-s" => { shader_pipeline_name = args[i + 1].clone(); }
            _ => ()
        }
    }

    let params = app::Params {
        width: WIDTH,
        height: HEIGHT,
        print_fps: false,
        asset_path,
        shader_pipeline_name,
    };

    app::run(params)?;

    return Ok(());
}
