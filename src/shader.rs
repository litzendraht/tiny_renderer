use nalgebra as na;
use na::{vector, Vector3, Vector4, Matrix4, Matrix2x3};

use crate::scene::Model;

/// Utility for getting convex combination of 2 Vector3<u8>'s
fn color_blend(color_1: Vector3<u8>, color_2: Vector3<u8>, t: f32) -> Vector3<u8>{
    return vector![
        (t * color_1.x as f32 + (1.0 - t) * color_2.x as f32) as u8,
        (t * color_1.y as f32 + (1.0 - t) * color_2.y as f32) as u8,
        (t * color_1.z as f32 + (1.0 - t) * color_2.z as f32) as u8
    ];
}

/// Transformation of a point to homogenous coordinates.
fn to_hom_point(v: Vector3<f32>) -> Vector4<f32> {
    return vector![v.x, v.y, v.z, 1.0];
}

/// Transformation of a vector to homogenous coordinates.
fn to_hom_vector(v: Vector3<f32>) -> Vector4<f32> {
    return vector![v.x, v.y, v.z, 0.0];
}

/// Transformation of a point from homogenous coordinates.
fn from_hom_point(v: Vector4<f32>) -> Vector3<f32> {
    return vector![v.x / v.w, v.y / v.w, v.z / v.w];
}

/// Transformation of a vector from homogenous coordinates.
fn from_hom_vector(v: Vector4<f32>) -> Vector3<f32> {
    return vector![v.x, v.y, v.z];
}

/// Pipeline, which defines "vertex" and "fragment" shaders as 2 methods and holds some values in order to pass them
/// between those methods, when applied in a scene. Methods require information about the model being rendered.
// @TODO figure out how to not have to pass model and texture to methods - approach with storing references to them
// in the struct produces complications with lifetimes and I was not able to overcome them immediately, which
// produced thi awkward unobviuous implementation.
pub trait ShaderPipeline {
    fn new() -> Self;
    fn get_buffer(&self) -> &Buffer;
    fn get_buffer_mut(&mut self) -> &mut Buffer;
    /// "Vertex" shader, taking model, 3 indices, defining the polygon of a model and transformation matrix,
    /// which produces vertex coordinates on a screen and z-buffer values. Basically is resposible for working
    /// on geometry.
    fn vertex(
        &mut self, 
        model: &Model,
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>      
    ) -> bool;
    /// "Fragment" shader, taking a texture and barycendric coordinates of a fragment.
    fn fragment(
        &mut self,
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool;
}

#[derive(Default, Clone, Copy)]
/// Buffer for passing values between different stages of a pipeline and setting up frame constants
/// like light direction and transform matrices.
pub struct Buffer {
    pub camera_direction:              Vector3<f32>,
    pub light_direction:               Vector3<f32>,
    pub vertex_transform_matrix:       Matrix4<f32>,   // Applied to vertcies to get final screen coodrdinates.
    pub direction_transform_matrix:    Matrix4<f32>,   // Applied to light direction.
    pub it_direction_transform_matrix: Matrix4<f32>,   // Applied to model normals.
    // Local buffer for passing values between vertex and fragment parts of the pipeline.    
    pub vertex_intensities:            Vector3<f32>,   // Light intensity in each vertex of a polygon.
    pub vertex_transformed_coords:     Matrix2x3<i32>, // Coordinates after all transformation, including viewport.
    pub vertex_z_values:               Vector3<f32>,   // Value used for comparison with existing z-buffer values.
    pub vertex_uvs:                    Matrix2x3<f32>, // UV coordinates, defining where to look for a color of a vertex.
    // Access to color and vertex z-value after application of fragment shader.
    pub fragment_color:                Vector3<u8> 
}

impl Buffer {
    fn new() -> Self {
        return Buffer {
            ..Default::default()
        };
    }
}

#[derive(Clone, Copy)]
pub struct DefaultSP {
    buffer: Buffer
}

#[derive(Clone, Copy)]
pub struct GouraudSP {
    buffer: Buffer
}

impl ShaderPipeline for DefaultSP {
    fn new() -> Self {
        return DefaultSP {
            buffer: Buffer::new(),
        };
    }

    fn get_buffer(&self) -> &Buffer {
        return &self.buffer;
    }

    fn get_buffer_mut(&mut self) -> &mut Buffer {
        return &mut self.buffer;
    }

    fn vertex(
        &mut self, 
        model: &Model,         
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>
    ) -> bool {
        let transformed_light_direction = from_hom_vector(
            self.buffer.direction_transform_matrix * to_hom_vector(self.buffer.light_direction)
        ).normalize();

        // Getting actual vertex coordinates from the model using given indices.
        let mut vertices = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertices[i] = vector![
                model.obj.positions[pos_indices[i]].0,
                model.obj.positions[pos_indices[i]].1,
                model.obj.positions[pos_indices[i]].2
            ];
        }

        // Calculating normal projection on the face.
        let face_normal = (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0]));
        let transformed_face_normal = from_hom_vector(
            self.buffer.it_direction_transform_matrix * to_hom_vector(face_normal)
        ).normalize();
        let normal_correction_coef = transformed_light_direction.dot(&transformed_face_normal);

        // @TODO no culling, because I don't pas camera position to shader, worth fixing.
        // // Backface culling.
        // if normal_correction_coef < 0.0 {
        //     return false;
        // }

        self.buffer.vertex_intensities = vector![
            normal_correction_coef, 
            normal_correction_coef, 
            normal_correction_coef
        ];
        
        for i in 0..3 {
            let transformed_v = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertices[i])
            );
            self.buffer.vertex_transformed_coords.set_column(
                i,
                &vector![
                    transformed_v.x as i32,
                    transformed_v.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = transformed_v.z;
        }

        // @TODO figure out, why texture is upside down, lol? Why do I need to do 1 - y to get correct vertex_uvs?
        for i in 0..3 {
            self.buffer.vertex_uvs.set_column(
                i, 
                &vector![
                    model.obj.tex_coords[tex_indices[i]].0,
                    1.0 - model.obj.tex_coords[tex_indices[i]].1
                ]
            );
        }

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        // Finding texture uv and coordinate with the help of calculated barycentric coordinates.
        let fragment_uv = self.buffer.vertex_uvs * bar_coord; 

        // Converting uv into explicit tex_coords.
        let texture_coord = vector![
            (fragment_uv.x * model.texture.width() as f32) as u32,
            (fragment_uv.y * model.texture.height() as f32) as u32
        ];

        // Getting color at texture_coord.
        let texture_color = Vector3::<u8>::from_row_slice(
            &model.texture.get_pixel(texture_coord.x, texture_coord.y).0[0..3]
        );

        // Calculating intensity at the fragment by combining barycentric coordinates and 
        // calculated vertex intensities.
        let normal_correction_coef = bar_coord.dot(&self.buffer.vertex_intensities);
        self.buffer.fragment_color = color_blend(texture_color, vector![0, 0, 0], normal_correction_coef);

        return true;
    }
}

impl ShaderPipeline for GouraudSP {
    fn new() -> Self {
        return GouraudSP {
            buffer: Buffer::new(),
        };
    }

    fn get_buffer(&self) -> &Buffer {
        return &self.buffer;
    }

    fn get_buffer_mut(&mut self) -> &mut Buffer {
        return &mut self.buffer;
    }

    fn vertex(
        &mut self, 
        model: &Model,         
        pos_indices: Vector3<usize>,
        tex_indices: Vector3<usize>,
        normal_indices: Vector3<usize>
    ) -> bool {
        // let transformed_light_direction = from_hom(to_hom(self.buffer.light_direction));
        let transformed_light_direction = from_hom_vector(
            self.buffer.direction_transform_matrix * to_hom_vector(self.buffer.light_direction)
        ).normalize();

        // Getting actual vertex coordinates from the model using given indices.
        let mut vertices = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertices[i] = vector![
                model.obj.positions[pos_indices[i]].0,
                model.obj.positions[pos_indices[i]].1,
                model.obj.positions[pos_indices[i]].2
            ];
        }

        // Calculating normal from the face.
        let face_normal = (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0])).normalize();

        // // Backface culling.
        // if transformed_light_direction.dot(&face_normal) < 0.0 {
        //     return false;
        // }

        // Calculating light intensities at each vertex to then interpolate them in fragment shader.
        for i in 0..3 {
            let vertex_normal = vector![
                model.obj.normals[normal_indices[i]].0,
                model.obj.normals[normal_indices[i]].1,
                model.obj.normals[normal_indices[i]].2
            ];
            let transformed_vertex_normal = from_hom_vector(
                self.buffer.it_direction_transform_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            self.buffer.vertex_intensities[i] = transformed_light_direction.dot(&transformed_vertex_normal);
        }
        
        for i in 0..3 {
            let transformed_v = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertices[i])
            );
            self.buffer.vertex_transformed_coords.set_column(
                i,
                &vector![
                    transformed_v.x as i32,
                    transformed_v.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = transformed_v.z;
        }

        // @TODO figure out, why texture is upside down, lol? Why do I need to do 1 - y to get correct vertex_uvs?
        for i in 0..3 {
            self.buffer.vertex_uvs.set_column(
                i, 
                &vector![
                    model.obj.tex_coords[tex_indices[i]].0,
                    1.0 - model.obj.tex_coords[tex_indices[i]].1
                ]
            );
        }

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        // Finding texture uv and coordinate with the help of calculated barycentric coordinates.
        let fragment_uv = self.buffer.vertex_uvs * bar_coord; 

        // Converting uv into explicit tex_coords.
        let texture_coord = vector![
            (fragment_uv.x * model.texture.width() as f32) as u32,
            (fragment_uv.y * model.texture.height() as f32) as u32
        ];

        // Getting color at texture_coord.
        let texture_color = Vector3::<u8>::from_row_slice(
            &model.texture.get_pixel(texture_coord.x, texture_coord.y).0[0..3]
        );

        // Calculating intensity at the fragment by combining barycentric coordinates and 
        // calculated vertex intensities.
        let normal_correction_coef = bar_coord.dot(&self.buffer.vertex_intensities);
        self.buffer.fragment_color = color_blend(texture_color, vector![0, 0, 0], normal_correction_coef);

        return true;
    }
}
