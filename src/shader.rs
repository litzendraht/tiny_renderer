use std::collections::HashMap;

use nalgebra as na;
use na::Vector3;

trait ShaderPipeline {
    fn vertex(&mut self, a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>, params: &mut HashMap<&str, f32>);
    fn fragment(&mut self, barycentric_coord: Vector3<f32>, params: &HashMap<&str, f32>) -> bool;
}

pub struct GouraudSP {
    // Global shader parameters.
    pub light_direction: Vector3<f32>,
    // Local buffer for passing values between vertex and fragment parts of the pipeline.
    face_intensity: f32, // Buffer between vertex and fragment shaders.
    pub fragment_color: Vector3<u8> // Access to color after application of fragment shader.
}

impl ShaderPipeline for GouraudSP {
    fn vertex(&mut self, a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>, params: &mut HashMap<&str, f32>) {
        self.face_intensity = 1.0;
    }
    fn fragment(&mut self, barycentric_coord: Vector3<f32>, params: &HashMap<&str, f32>) -> bool {
        return false;
    }
}