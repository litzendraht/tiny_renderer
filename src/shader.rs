use obj::raw::RawObj;
use image::RgbImage;
use nalgebra as na;
use na::{vector, Vector2, Vector3, Matrix4};

pub trait ShaderPipeline<'a> {
    fn new(model: &'a RawObj, texture: &'a RgbImage) -> Self;
    /// Getter for transformed coordinates obtained by aplying vertex part of the pipeline to the polygon.
    fn get_transformed_coords(&self) -> &[Vector2<i32>];
    /// Getter for z values of vertices obtained by aplying vertex part of the pipeline to the polygon.
    fn get_vertex_z_values(&self) -> &[f32];
    /// "Vertex" shader, taking 3 indices, defining the polygon of a model and returning bool,
    /// telling whether the polygon should be drawn or not. During processing saves transformed
    /// coordinates 
    fn vertex(&mut self, transform_matrix: Matrix4<f32>, indices: Vector3<usize>) -> bool;
    fn fragment(&mut self, bar_coord: Vector3<f32>) -> Option<Vector3<u8>>;
}

#[derive(Clone, Copy)]
pub struct DefaultSP<'a> {
    // Global shader parameters.
    model: &'a RawObj,
    texture: &'a RgbImage,
    // Global shader parameters.
    pub light_direction: Vector3<f32>,
    // Local buffer for passing values between vertex and fragment parts of the pipeline.    
    pub varying_intensity: Vector3<f32>, // Intensity in each vertex of a polygon.
    pub transformed_coords: [Vector2<i32>; 3], // Coordinates after all transformation, including viewport.
    pub vertex_z_values: [f32; 3], // Value used for comparison with existing z-buffer values.
    // Access to color after application of fragment shader.
    pub fragment_color: Vector3<u8> 
}

// #[derive(Clone, Copy)]
// pub struct GouraudSP<'a> {
//     // Global shader parameters.
//     model: &'a RawObj,
//     texture: &'a RgbImage,
//     pub light_direction: Vector3<f32>,
//     // Local buffer for passing values between vertex and fragment parts of the pipeline.

//     varying_intensity: Vector3<f32>, // Intensity in each vertex of a polygon.
//     // Access to color after application of fragment shader.
//     pub fragment_color: Vector3<u8> 
// }

impl<'a> ShaderPipeline<'a> for DefaultSP<'a> {
    fn new(model: &'a RawObj, texture: &'a RgbImage) -> Self {
        return DefaultSP {
            model,
            texture,
            light_direction: vector![0.0, 0.0, 1.0],
            varying_intensity: vector![0.0, 0.0, 0.0],
            transformed_coords: [vector![0, 0]; 3],
            vertex_z_values: [0.0; 3],
            fragment_color: vector![0, 0, 0],
        };
    }

    fn get_transformed_coords(&self) -> &[Vector2<i32>] {
        return &self.transformed_coords;
    }

    fn get_vertex_z_values(&self) -> &[f32] {
        return &self.vertex_z_values;
    }

    fn vertex(&mut self, transform_matrix: Matrix4<f32>, indices: Vector3<usize>) -> bool {
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertices = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertices[i] = vector![
                self.model.positions[indices[i]].0,
                self.model.positions[indices[i]].1,
                self.model.positions[indices[i]].2
            ];
        }

        // Calculating normal projection on the face.
        let face_normal = (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0]));
        let mut normal_correction_coef = self.light_direction.dot(&face_normal);
        // Backface culling.
        if normal_correction_coef < 0.0 {
            return false;
        }

        normal_correction_coef /= face_normal.norm();
        self.varying_intensity = vector![
            normal_correction_coef, 
            normal_correction_coef, 
            normal_correction_coef
        ];
        
        for i in 0..3 {
            let hom_v = vector![vertices[i].x, vertices[i].y, vertices[i].z, 1.0];
            let hom_transformed_v = transform_matrix * hom_v;
            let transformed_v = vector![
                hom_transformed_v.x / hom_transformed_v.w,
                hom_transformed_v.y / hom_transformed_v.w,
                hom_transformed_v.z / hom_transformed_v.w
            ];
            self.transformed_coords[0][0] = (hom_transformed_v.x / hom_transformed_v.w) as i32;
            self.transformed_coords[0][1] = (hom_transformed_v.y / hom_transformed_v.w) as i32;
            self.vertex_z_values[i] = transformed_v.z;
        }

        return true;
    }

    fn fragment(&mut self, bar_coord: Vector3<f32>) -> Option<Vector3<u8>> {
        return Some(vector![0, 0, 0]);
    }
}

// impl<'a> ShaderPipeline for GouraudSP<'a> {
//     fn new(model: &RawObj, texture: &RgbImage) -> Self {
//         return GouraudSP {
//             model,
//             texture,
//             light_direction: vector![0.0, 0.0, 1.0],
//             varying_intensity: vector![0.0, 0.0, 0.0],
//             fragment_color: vector![0, 0, 0],
//         };
//     }

//     fn vertex(&mut self, indices: Vector3<usize>) -> bool{
//         let a = vector![
//             self.model.positions[indices[0]].0,
//             self.model.positions[indices[0]].1,
//             self.model.positions[indices[0]].2
//         ];
//         let b = vector![
//             self.model.positions[indices[1]].0,
//             self.model.positions[indices[1]].1,
//             self.model.positions[indices[1]].2
//         ];
//         let c = vector![
//             self.model.positions[indices[2]].0,
//             self.model.positions[indices[2]].1,
//             self.model.positions[indices[2]].2
//         ];

//         // Calculating normal projection on the face.
//         let face_normal = (b - a).cross(&(c - a));
//         let mut normal_correction_coef = self.light_direction.dot(&face_normal);
//         // Backface culling.
//         if normal_correction_coef < 0.0 {
//             return false;
//         } else {
//             normal_correction_coef /= face_normal.norm();
//             self.varying_intensity = vector![
//                 normal_correction_coef, 
//                 normal_correction_coef, 
//                 normal_correction_coef
//                 ];
//             return true;
//         }
        
//     }

//     fn fragment(&mut self, bar_coord: Vector3<f32>) -> Option<Vector3<u8>> {
//         return None;
//     }
// }