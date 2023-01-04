// @TODO no culling at the momement, drawing every face, maybe worth fixing.

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
    pub t_light_direction:             Vector3<f32>,   // Light direction with model and view transformations applied.
    pub vertex_transform_matrix:       Matrix4<f32>,   // Applied to vertcies to get final screen coodrdinates.
    pub direction_transform_matrix:    Matrix4<f32>,   // Applied to light direction.
    pub it_direction_transform_matrix: Matrix4<f32>,   // Applied to model normals.
    // Local buffer for passing values between vertex and fragment parts of the pipeline.    
    pub vertex_intensities:            Vector3<f32>,   // Light intensity in each vertex of a polygon.
    pub vertex_t_coords:               Matrix2x3<i32>, // Coordinates after all transformation, including viewport.
    pub vertex_z_values:               Vector3<f32>,   // Value used for comparison with existing z-buffer values.
    pub vertex_uvs:                    Matrix2x3<f32>, // UV coordinates, defining where to look for a color of a vertex.
    // Access to color after application of fragment shader.
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

#[derive(Clone, Copy)]
pub struct TrueNormalSP {
    buffer: Buffer
}

#[derive(Clone, Copy)]
pub struct SpecularSP {
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
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }

        // Calculating normal projection on the face.
        let face_normal = (vertex_positions[1] - vertex_positions[0]).cross(
            &(vertex_positions[2] - vertex_positions[0])
        );
        let t_face_normal = from_hom_vector(
            self.buffer.it_direction_transform_matrix * to_hom_vector(face_normal)
        ).normalize();
        let diff_coef = self.buffer.t_light_direction.dot(&t_face_normal);

        self.buffer.vertex_intensities = vector![
            diff_coef, 
            diff_coef, 
            diff_coef
        ];
        
        for i in 0..3 {
            let t_vertex_position = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertex_positions[i])
            );
            self.buffer.vertex_t_coords.set_column(
                i,
                &vector![
                    t_vertex_position.x as i32,
                    t_vertex_position.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = t_vertex_position.z;
        }

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
        // Correcting texture color by intensity, but in this pipeline normal is just a face normal, so
        // intensities are all equal and we can just take one at index 0 in vertex_intensities vector.
        let uv = self.buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = self.buffer.vertex_intensities[0];
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

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
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }

        // Calculating light intensities at each vertex to then interpolate them in fragment shader.
        for i in 0..3 {
            let vertex_normal = vector![
                model.obj.normals[normal_indices[i]].0,
                model.obj.normals[normal_indices[i]].1,
                model.obj.normals[normal_indices[i]].2
            ];
            let t_vertex_normal = from_hom_vector(
                self.buffer.it_direction_transform_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            self.buffer.vertex_intensities[i] = self.buffer.t_light_direction.dot(
                &t_vertex_normal
            );
        }
        
        for i in 0..3 {
            let t_vertex_position = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertex_positions[i])
            );
            self.buffer.vertex_t_coords.set_column(
                i,
                &vector![
                    t_vertex_position.x as i32,
                    t_vertex_position.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = t_vertex_position.z;
        }

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
        // Correcting texture color by intensity, where intensity is based on the interpolated normal
        // between the vertex_positions, which is produced via barycentric coordinates of the fragment.
        let uv = self.buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = bar_coord.dot(&self.buffer.vertex_intensities);
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}

impl ShaderPipeline for TrueNormalSP {
    fn new() -> Self {
        return TrueNormalSP {
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
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }
        
        for i in 0..3 {
            let t_vertex_position = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertex_positions[i])
            );
            self.buffer.vertex_t_coords.set_column(
                i,
                &vector![
                    t_vertex_position.x as i32,
                    t_vertex_position.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = t_vertex_position.z;
        }

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
        // Correcting texture color by intensity, where intensity is based on the normal, produced
        // by the full normal map of the model.
        let uv = self.buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        // @TODO why whis does not work correctly?
        let t_fragment_normal = from_hom_vector(
            self.buffer.it_direction_transform_matrix * to_hom_vector(fragment_normal)
        ).normalize();
        let diff_coef = self.buffer.t_light_direction.dot(&t_fragment_normal);
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}

impl ShaderPipeline for SpecularSP {
    fn new() -> Self {
        return SpecularSP {
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
        // Getting actual vertex coordinates from the model using given indices.
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }
        
        for i in 0..3 {
            let t_vertex_position = from_hom_point(
                self.buffer.vertex_transform_matrix * to_hom_point(vertex_positions[i])
            );
            self.buffer.vertex_t_coords.set_column(
                i,
                &vector![
                    t_vertex_position.x as i32,
                    t_vertex_position.y as i32
                ]
            );
            self.buffer.vertex_z_values[i] = t_vertex_position.z;
        }

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
        // Correcting texture color by intensity, where intensity is based on the normal, produced
        // by the full normal map of the model.
        let uv = self.buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        // @TODO why whis does not work correctly?
        let t_fragment_normal = from_hom_vector(
            self.buffer.it_direction_transform_matrix * to_hom_vector(fragment_normal)
        ).normalize();
        // Important minus here - light direction is from source to 
        let reflected_light_direction = (
            2.0 * (
                t_fragment_normal * 
                self.buffer.t_light_direction.dot(&t_fragment_normal)
            ) - self.buffer.t_light_direction
        ).normalize();
        let diff_coef = self.buffer.t_light_direction.dot(&t_fragment_normal);
        let spec_coef = 0.6 * reflected_light_direction.z.max(0.0).powf(model.get_specular_value_at_uv(uv));
        let corrected_color = vector![
            ((diff_coef + spec_coef) * color[0] as f32).min(255.0) as u8,
            ((diff_coef + spec_coef) * color[1] as f32).min(255.0) as u8,
            ((diff_coef + spec_coef) * color[2] as f32).min(255.0) as u8
        ];
        self.buffer.fragment_color = corrected_color;

        return true;
    }
}
