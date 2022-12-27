mod util;
mod scene;
mod app;

const WIDTH: u32  = 800;
const HEIGHT: u32 = 800;

#[show_image::main]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // @TODO command line arguments to change context.
    let context = app::Context {
        width: WIDTH,
        height: HEIGHT,
        print_fps: true,
        asset_path: String::from("assets/african_head"),
    };

    app::run(context)?;

    return Ok(());
}
