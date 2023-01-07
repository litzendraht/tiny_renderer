mod scene;
mod app;

const WIDTH: u32  = 800;
const HEIGHT: u32 = 800;

#[show_image::main]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // @TODO command line arguments to change app parameters.
    let params = app::Params {
        width: WIDTH,
        height: HEIGHT,
        print_fps: true,
        // asset_path: String::from("assets/african_head"),
        asset_path: String::from("assets/diablo"),
        // shader_pipeline_name: "default",
        // shader_pipeline_name: "phong",
        // shader_pipeline_name: "normal_map",
        // shader_pipeline_name: "specular",
        // shader_pipeline_name: "darboux",
        // shader_pipeline_name: "shadow",
        shader_pipeline_name: "occlusion",
    };

    app::run(params)?;

    return Ok(());
}
