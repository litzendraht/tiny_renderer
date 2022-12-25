mod app;
mod image;

const WIDTH: u32  = 800;
const HEIGHT: u32 = 600;

#[show_image::main]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    app::run(app::Context {
        width: WIDTH,
        height: HEIGHT,
    })?;

    return Ok(());
}
