// @TODO no culling at the momement, drawing every face, maybe worth fixing.
use super::util::{to_hom_point, from_hom_point, to_hom_vector, from_hom_vector, color_blend};

use nalgebra as na;
use na::{vector, Vector3, Matrix3, Matrix4, Matrix2x3};

use crate::scene::Model;

/// Pipeline, which defines "vertex" and "fragment" shaders as 2 methods and holds some values in a buffer 
/// in order to pass them between those methods, when applied in a scene. 
/// Methods require information about the model being rendered.
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
    /// "Fragment" shader, taking model struct and barycendric coordinates of a fragment.
    fn fragment(
        &mut self,
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool;
}

#[derive(Default, Clone, Copy)]
/// Buffer for passing values between different stages of a pipeline and setting up frame constants
/// like light direction and transform matrices.
// @TODO not all buffer is used in every shader, but for simplicity it just has all options. Maybe
// worth thinking about how to get a more fine-grained buffer structure for each shader.
pub struct Buffer {
    pub camera_direction:  Vector3<f32>,
    pub t_light_direction: Vector3<f32>,   // Light direction with model and view transformations applied.
    pub vpmv_matrix:       Matrix4<f32>,   // Applied to vertices to get final screen coodrdinates.
    pub mv_matrix:         Matrix4<f32>,   // Applied to light direction.
    pub it_mv_matrix:      Matrix4<f32>,   // Applied to model normals.
    // Local buffer for passing values between vertex and fragment parts of the pipeline.    
    vertex_intensities:    Vector3<f32>,   // Light intensity in each vertex of a polygon.
    vertex_t_positions:    Matrix3<f32>,   // Transformed vertex positions as columns.
    vertex_t_normals:      Matrix3<f32>,   // Transformed vertex normals at each vertex as columns.
    vertex_uvs:            Matrix2x3<f32>, // UV coordinates, defining where to look for a color of a vertex as columns.
    pub vertex_t_coords:   Matrix2x3<i32>, // Coordinates after all transformation, including viewport as columns.
    pub vertex_z_values:   Vector3<f32>,   // Value used for comparison with existing z-buffer values.
    // Access to color after application of fragment shader.
    pub fragment_color:    Vector3<u8> 
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
pub struct PhongSP {
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

#[derive(Clone, Copy)]
pub struct DarbouxSP {
    buffer: Buffer
}

/// Boilerplate for moving uvs to a 2x3 matrix buffer.
fn store_vertex_uvs(
    uvs_buffer: &mut Matrix2x3<f32>, 
    tex_coords: &Vec<(f32, f32, f32)>,
    indices: Vector3<usize>
) {
    for i in 0..3 {
        uvs_buffer.set_column(
            i, 
            &vector![
                tex_coords[indices[i]].0,
                1.0 - tex_coords[indices[i]].1
            ]
        );
    }
}

/// Boilerplate for transforming and moving vertex information, namely screen coords and z-values into buffers.
fn store_vertex_transformation_results(
    vertex_positions: [Vector3<f32>; 3],
    vpmv_matrix: Matrix4<f32>,
    t_coords_buffer: &mut Matrix2x3<i32>, 
    z_values_buffer: &mut Vector3<f32>
) {
    for i in 0..3 {
        let vertex_t_position = from_hom_point(
            vpmv_matrix * to_hom_point(vertex_positions[i])
        );
        t_coords_buffer.set_column(
            i,
            &vector![
                vertex_t_position.x as i32,
                vertex_t_position.y as i32
            ]
        );
        z_values_buffer[i] = vertex_t_position.z;
    }
}

/// Pipeline calculates normal for a face and corrects diffuse texture color for a fragment by using
/// intensity based on the dot product between light direction and face normal.
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
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }

        // Calculating normal projection on the face.
        let face_normal = (vertex_positions[1] - vertex_positions[0]).cross(
            &(vertex_positions[2] - vertex_positions[0])
        );
        let t_face_normal = from_hom_vector(
            self.buffer.it_mv_matrix * to_hom_vector(face_normal)
        ).normalize();
        let diff_coef = self.buffer.t_light_direction.dot(&t_face_normal);

        self.buffer.vertex_intensities = vector![diff_coef, diff_coef, diff_coef];

        store_vertex_transformation_results(
            vertex_positions,
            self.buffer.vpmv_matrix, 
            &mut self.buffer.vertex_t_coords, 
            &mut self.buffer.vertex_z_values
        );
        store_vertex_uvs(&mut self.buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        let uv = self.buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = self.buffer.vertex_intensities[0];
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}

/// Pipeline interpolates fragment normal by using barycentric coordinates and vertex normals
/// and corrects fragment color using resulting light intensity.
impl ShaderPipeline for PhongSP {
    fn new() -> Self {
        return PhongSP {
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
            let vertex_t_normal = from_hom_vector(
                self.buffer.it_mv_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            self.buffer.vertex_intensities[i] = self.buffer.t_light_direction.dot(
                &vertex_t_normal
            );
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            self.buffer.vpmv_matrix, 
            &mut self.buffer.vertex_t_coords, 
            &mut self.buffer.vertex_z_values
        );
        store_vertex_uvs(&mut self.buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        let uv = self.buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = bar_coord.dot(&self.buffer.vertex_intensities);
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}

/// Pipeline corrects fragemnt color by using the normal, provided by the normal map of the model.
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
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            self.buffer.vpmv_matrix, 
            &mut self.buffer.vertex_t_coords, 
            &mut self.buffer.vertex_z_values
        );
        store_vertex_uvs(&mut self.buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        let uv = self.buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        let t_fragment_normal = from_hom_vector(
            self.buffer.it_mv_matrix * to_hom_vector(fragment_normal)
        ).normalize();
        let diff_coef = self.buffer.t_light_direction.dot(&t_fragment_normal);
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}

/// Pipeline calculates fragment color using a model accounting for a reflected specular light component.
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
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            self.buffer.vpmv_matrix, 
            &mut self.buffer.vertex_t_coords, 
            &mut self.buffer.vertex_z_values
        );
        store_vertex_uvs(&mut self.buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        let uv = self.buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        let t_fragment_normal = from_hom_vector(
            self.buffer.it_mv_matrix * to_hom_vector(fragment_normal)
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

/// Pipeline using normals stored in Darboux frame for color intensity calculation.
impl ShaderPipeline for DarbouxSP {
    fn new() -> Self {
        return DarbouxSP {
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
        let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
        for i in 0..3 {
            vertex_positions[i] = model.get_vertex_position_at_index(pos_indices[i]);
            // Collecting transformed vertex positions in a buffer to use in local basis calculation.
            self.buffer.vertex_t_positions.set_column(i, &from_hom_point(
                    self.buffer.mv_matrix * to_hom_point(vertex_positions[i])
                )
            );
        }

        // Collecting transformed normals at each vertex into a single matrix for subsequent interpolation
        // in a fragment shader.
        for i in 0..3 {
            let vertex_normal = vector![
                model.obj.normals[normal_indices[i]].0,
                model.obj.normals[normal_indices[i]].1,
                model.obj.normals[normal_indices[i]].2
            ];
            let vertex_t_normal = from_hom_vector(
                self.buffer.it_mv_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            self.buffer.vertex_t_normals.set_column(i, &vertex_t_normal);
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            self.buffer.vpmv_matrix, 
            &mut self.buffer.vertex_t_coords, 
            &mut self.buffer.vertex_z_values
        );
        store_vertex_uvs(&mut self.buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    }

    fn fragment(
        &mut self, 
        model: &Model,
        bar_coord: Vector3<f32>
    ) -> bool {
        let uv = self.buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal_tangent = model.get_normal_tangent_at_uv(uv);
        // Calculating the matrix, giving required transformation from Darboux basis to the global one.
        // A is a matrix with rows representing a local basis at each fragment.
        let mut local_basis_matrix: Matrix3<f32> = Default::default();
        let local_z = self.buffer.vertex_t_normals * bar_coord;
        local_basis_matrix.set_row(
            0, &(self.buffer.vertex_t_positions * vector![-1.0, 1.0, 0.0]).normalize().transpose()
        );
        local_basis_matrix.set_row(
            1, &(self.buffer.vertex_t_positions * vector![-1.0, 0.0, 1.0]).normalize().transpose()
        );
        local_basis_matrix.set_row(
            2, &(self.buffer.vertex_t_normals * bar_coord).normalize().transpose()
        );

        let i_local_basis_matrix = local_basis_matrix.try_inverse().unwrap();
        let local_x = i_local_basis_matrix * vector![
            self.buffer.vertex_uvs.m12 - self.buffer.vertex_uvs.m11,
            self.buffer.vertex_uvs.m13 - self.buffer.vertex_uvs.m11,
            0.0
        ];
        let local_y = i_local_basis_matrix * vector![
            self.buffer.vertex_uvs.m22 - self.buffer.vertex_uvs.m21,
            self.buffer.vertex_uvs.m23 - self.buffer.vertex_uvs.m21,
            0.0
        ];

        let mut local_transform_matrix: Matrix3<f32> = Default::default();
        local_transform_matrix.set_column(0, &local_x.normalize());
        local_transform_matrix.set_column(1, &local_y.normalize());
        local_transform_matrix.set_column(2, &local_z.normalize());
        let t_fragment_normal = (local_transform_matrix * fragment_normal_tangent).normalize();
        
        let diff_coef = self.buffer.t_light_direction.dot(&t_fragment_normal);
        self.buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    }
}
