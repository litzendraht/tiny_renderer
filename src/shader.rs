use obj::raw::RawObj;
use image::RgbImage;
use nalgebra as na;
use na::{vector, Vector2, Vector3, Matrix4};

/// Utility for getting convex combination of 2 Vector3<u8>'s
fn color_blend(color_1: Vector3<u8>, color_2: Vector3<u8>, t: f32) -> Vector3<u8>{
    return vector![
        (t * color_1.x as f32 + (1.0 - t) * color_2.x as f32) as u8,
        (t * color_1.y as f32 + (1.0 - t) * color_2.y as f32) as u8,
        (t * color_1.z as f32 + (1.0 - t) * color_2.z as f32) as u8
    ];
}

/// Pipeline, which defines "vertex" and "fragment" shaders as 2 methods and holds some values in order to pass them
/// between those methods, when applied in a scene. Methods require information about the model being rendered.
// @TODO figure out how to not have to pass model and texture to methods - approach with storing references to them
// in the struct produces complications with lifetimes and I was not able to overcome them immediately, which
// produced thi awkward unobviuous implementation.
pub trait ShaderPipeline {
    fn new() -> Self;
    fn get_buffer(&self) -> &Buffer;
    /// "Vertex" shader, taking model, 3 indices, defining the polygon of a model and transformation matrix,
    /// which produces vertex coordinates on a screen and z-buffer values. Basically is resposible for working
    /// on geometry.
    fn vertex(
        &mut self, 
        model: &RawObj,
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>,
        transform_matrix: Matrix4<f32>        
    ) -> bool;
    fn fragment(
        &mut self,
        texture: &RgbImage,
        bar_coord: Vector3<f32>
    ) -> bool;
}

#[derive(Clone, Copy)]
/// Buffer for passing values between different stages of a pipeline.
pub struct Buffer {
    // Local buffer for passing values between vertex and fragment parts of the pipeline.    
    pub vertex_intensities: Vector3<f32>, // Light intensity in each vertex of a polygon.
    pub vertex_transformed_coords: [Vector2<i32>; 3], // Coordinates after all transformation, including viewport.
    pub vertex_z_values: [f32; 3], // Value used for comparison with existing z-buffer values.
    pub vertex_uvs: [Vector2<f32>; 3], // UV coordiantes, defining where to look for a color of a vertex.
    // Access to color and vertex z-value after application of fragment shader.
    pub fragment_color: Vector3<u8> 
}

#[derive(Clone, Copy)]
pub struct DefaultSP {
    // Global shader parameters.
    light_direction: Vector3<f32>,
    buffer: Buffer
}

#[derive(Clone, Copy)]
pub struct GouraudSP {
    // Global shader parameters.
    light_direction: Vector3<f32>,
    buffer: Buffer
}

impl Buffer {
    fn new() -> Self {
        return Buffer {
            vertex_intensities:        vector![0.0, 0.0, 0.0],
            vertex_transformed_coords: [vector![0, 0]; 3],
            vertex_z_values:           [0.0; 3],
            vertex_uvs:                [vector![0.0, 0.0]; 3],
            fragment_color:            vector![0, 0, 0],
        };
    }
}

impl ShaderPipeline for DefaultSP {
    fn new() -> Self {
        return DefaultSP {
            light_direction: vector![0.0, 0.0, 1.0],
            buffer: Buffer::new(),
        };
    }

    fn get_buffer(&self) -> &Buffer {
        return &self.buffer;
    }

    fn vertex(
        &mut self, 
        model: &RawObj,         
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>,
        transform_matrix: Matrix4<f32>,
    ) -> bool {
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertices = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertices[i] = vector![
                model.positions[pos_indices[i]].0,
                model.positions[pos_indices[i]].1,
                model.positions[pos_indices[i]].2
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
        self.buffer.vertex_intensities = vector![
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
            self.buffer.vertex_transformed_coords[i][0] = (hom_transformed_v.x / hom_transformed_v.w) as i32;
            self.buffer.vertex_transformed_coords[i][1] = (hom_transformed_v.y / hom_transformed_v.w) as i32;
            self.buffer.vertex_z_values[i] = transformed_v.z;
        }

        // @TODO figure out, why texture is upside down, lol? Why do I need to do 1 - y to get correct vertex_uvs?
        for i in 0..3 {
            self.buffer.vertex_uvs[i] = vector![
                model.tex_coords[tex_indices[i]].0,
                1.0 - model.tex_coords[tex_indices[i]].1
            ];
        }

        return true;
    }

    fn fragment(
        &mut self, 
        texture: &RgbImage,
        bar_coord: Vector3<f32>
    ) -> bool {
        // Finding texture uv and coordinate with the help of calculated barycentric coordinates.
        let fragment_uv = bar_coord.x * self.buffer.vertex_uvs[0] + 
                          bar_coord.y * self.buffer.vertex_uvs[1] + 
                          bar_coord.z * self.buffer.vertex_uvs[2]; 

        // Converting uv into explicit tex_coords.
        let texture_coord = vector![
            (fragment_uv.x * texture.width() as f32) as u32,
            (fragment_uv.y * texture.height() as f32) as u32
        ];

        // Getting color at texture_coord.
        let texture_color = vector![
            texture.get_pixel(texture_coord.x, texture_coord.y).0[0],
            texture.get_pixel(texture_coord.x, texture_coord.y).0[1],
            texture.get_pixel(texture_coord.x, texture_coord.y).0[2]
        ];

        // Calculating intensity at the fragment by combining barycentric coordinates and 
        // calculated vertex intensities.
        let normal_correction_coef = bar_coord.x * self.buffer.vertex_intensities[0] + 
                                     bar_coord.y * self.buffer.vertex_intensities[1] + 
                                     bar_coord.z * self.buffer.vertex_intensities[2];
        self.buffer.fragment_color = color_blend(texture_color, vector![0, 0, 0], normal_correction_coef);

        return true;
    }
}

impl ShaderPipeline for GouraudSP {
    fn new() -> Self {
        return GouraudSP {
            light_direction: vector![0.0, 0.0, 1.0],
            buffer: Buffer::new(),
        };
    }

    fn get_buffer(&self) -> &Buffer {
        return &self.buffer;
    }

    fn vertex(
        &mut self, 
        model: &RawObj,         
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>,
        transform_matrix: Matrix4<f32>,
    ) -> bool {
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertices = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertices[i] = vector![
                model.positions[pos_indices[i]].0,
                model.positions[pos_indices[i]].1,
                model.positions[pos_indices[i]].2
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
        self.buffer.vertex_intensities = vector![
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
            self.buffer.vertex_transformed_coords[i][0] = (hom_transformed_v.x / hom_transformed_v.w) as i32;
            self.buffer.vertex_transformed_coords[i][1] = (hom_transformed_v.y / hom_transformed_v.w) as i32;
            self.buffer.vertex_z_values[i] = transformed_v.z;
        }

        // @TODO figure out, why texture is upside down, lol? Why do I need to do 1 - y to get correct vertex_uvs?
        for i in 0..3 {
            self.buffer.vertex_uvs[i] = vector![
                model.tex_coords[tex_indices[i]].0,
                1.0 - model.tex_coords[tex_indices[i]].1
            ];
        }

        return true;
    }

    fn fragment(
        &mut self, 
        texture: &RgbImage,
        bar_coord: Vector3<f32>
    ) -> bool {
        // Finding texture uv and coordinate with the help of calculated barycentric coordinates.
        let fragment_uv = bar_coord.x * self.buffer.vertex_uvs[0] + 
                          bar_coord.y * self.buffer.vertex_uvs[1] + 
                          bar_coord.z * self.buffer.vertex_uvs[2]; 

        // Converting uv into explicit tex_coords.
        let texture_coord = vector![
            (fragment_uv.x * texture.width() as f32) as u32,
            (fragment_uv.y * texture.height() as f32) as u32
        ];

        // Getting color at texture_coord.
        let texture_color = vector![
            texture.get_pixel(texture_coord.x, texture_coord.y).0[0],
            texture.get_pixel(texture_coord.x, texture_coord.y).0[1],
            texture.get_pixel(texture_coord.x, texture_coord.y).0[2]
        ];

        // Calculating intensity at the fragment by combining barycentric coordinates and 
        // calculated vertex intensities.
        let normal_correction_coef = bar_coord.x * self.buffer.vertex_intensities[0] + 
                                     bar_coord.y * self.buffer.vertex_intensities[1] + 
                                     bar_coord.z * self.buffer.vertex_intensities[2];
        self.buffer.fragment_color = color_blend(texture_color, vector![0, 0, 0], normal_correction_coef);

        return true;
    }
}
