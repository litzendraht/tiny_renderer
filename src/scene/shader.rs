// @TODO no culling at the momement, drawing every face, maybe worth fixing.
use super::util::{Model, to_hom_point, from_hom_point, to_hom_vector, from_hom_vector, color_blend};

use nalgebra as na;
use na::{vector, Vector3, matrix, Matrix3, Matrix4, Matrix2x3};

/// Buffer for passing values between different stages of a pipeline and setting up frame constants
/// like light direction and transform matrices.
// @TODO not all buffer is used in every shader, but for simplicity it just has all options. Maybe
// worth thinking about how to get a more fine-grained buffer structure for each shader.
#[derive(Default)]
pub struct Buffer {
    // Pointers to some fat buffers like z-buffer and shadow buffer.    
    pub z_buffer:          Vec<f32>,
    pub shadow_buffer:     Vec<f32>,
    // Small stuff.
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
    fn new(width: u32, height: u32) -> Self {
        let frame_buffer_size = (width * height) as usize;
        return Self {
            z_buffer:      vec![0.0; frame_buffer_size],
            shadow_buffer: vec![0.0; frame_buffer_size],
            ..Default::default()
        }
    }
}

/// Type representing a funtion, which is called in order to prepare pipeline for application of vertex
/// and fragment shaders.
type Prepare = dyn Fn(
    &mut Buffer,    // Buffer.
    u32,            // Screen width.
    u32,            // Screen height.
    Vector3<f32>,   // Light direction.   
    Vector3<f32>,   // Camera placement.
    Vector3<f32>,   // Camera direction.
    Vector3<f32>    // Camera up direction.
);

/// Type representing vertex shader.
type VertexShader = dyn Fn(
    &mut Buffer,    // Buffer.
    &Model,         // Model info.
    Vector3<usize>, // Position indices.
    Vector3<usize>, // Diffuse texture indices.
    Vector3<usize>  // Vertex normal indices.
) -> bool;

/// Type representing fragment shader.
type FragmentShader = dyn Fn(
    &mut Buffer,    // Buffer
    &Model,         // Model info.
    Vector3<f32>,   // Barycentric coordinates.
    usize           // Index of the fragment in the framebuffer.
) -> bool;

/// Representation of one pass in the shader pipeline storing closures, representing a 3 steps - 
/// preparation of the pipeline buffer, vertex shader application and fragment shader application.
pub struct ShaderPass {
    pub prepare:  Box<Prepare>,
    pub vertex:   Box<VertexShader>,
    pub fragment: Box<FragmentShader>,
}

/// Simple struct to organize several passes of the pipeline and provide a reference to the buffer.
pub struct ShaderPipeline {
    pub buffer: Buffer,
    pub passes: Vec<ShaderPass>
}

/// Simple backface culling.
fn should_cull_face(vertex_positions: [Vector3<f32>; 3], camera_direction: Vector3<f32>) -> bool {
    let face_normal = (vertex_positions[1] - vertex_positions[0]).cross(
        &(vertex_positions[2] - vertex_positions[0])
    );
    if camera_direction.dot(&face_normal) <= 0.0 {
        return true;
    } else {
        return false;
    }
}

/// Boilerplate for accessing vertex positions from model vertex list.
fn get_vertex_positions(model: &Model, indices: Vector3<usize>) -> [Vector3<f32>; 3] {
    let mut vertex_positions = [vector![0.0, 0.0, 0.0]; 3];
    for i in 0..3 {
        vertex_positions[i] = model.get_vertex_position_at_index(indices[i]);
    }
    return vertex_positions;
}

/// Boilerplate for moving uvs to a 2x3 matrix buffer.
fn store_vertex_uvs(
    uvs_buffer: &mut Matrix2x3<f32>, 
    tex_coords: &Vec<(f32, f32, f32)>,
    indices:    Vector3<usize>
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
    vpmv_matrix:      Matrix4<f32>,
    t_coords_buffer:  &mut Matrix2x3<i32>, 
    z_values_buffer:  &mut Vector3<f32>
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

/// Boilerplate for checking z-value of the fragment against the z-buffer.
/// Returns false if there is no need to update the framebuffer.
fn process_z_value(buffer: &mut Buffer, bar_coord: Vector3<f32>, index: usize) -> bool {
    // Checking fragment z-value in the pipeline buffer and comparing it to the value in the 
    // buffer, on failure returning false, signifying that no further fragment processing
    // should be done
    let z_to_screen = bar_coord.dot(&buffer.vertex_z_values); 
    if z_to_screen <= buffer.z_buffer[index] {
        return false;
    }
    buffer.z_buffer[index] = z_to_screen;
    return true;
}

impl ShaderPipeline {
    pub fn new(pipeline_name: &str, width: u32, height: u32) -> Self {
        let buffer =        Buffer::new(width, height);
        let passes:         Vec<ShaderPass>;
        match pipeline_name {
            "default"      => passes = get_default_pipeline_passes(),
            "phong"        => passes = get_phong_pipeline_passes(),
            "normal_map"   => passes = get_normal_map_pipeline_passes(),
            "specular"     => passes = get_specular_pipeline_passes(),
            "darboux"      => passes = get_darboux_pipeline_passes(),
            "shadow"       => passes = get_shadow_pipeline_passes(),
            _ => panic!("Provided pipeline name is not supported!"),
        }

        return Self {
            buffer,
            passes,
        };
    }
}

/// Standard setup which prepares transforms to the basis relative to the camera.
fn get_default_prepare() -> Box<Prepare> {
    let prepare = |
        buffer:          &mut Buffer,
        width:           u32,
        height:          u32,
        light_direction: Vector3<f32>,
        look_from:       Vector3<f32>,
        look_at:         Vector3<f32>,
        up:              Vector3<f32>        
    | {
        // New coordinate system a, b, c around camera position.
        let new_z = (look_from - look_at).normalize();
        let new_y = (up - new_z.dot(&up) * new_z).normalize();
        let new_x = new_y.cross(&new_z).normalize();        
        let model_matrix = matrix![new_x.x, new_x.y, new_x.z, 0.0;
                                   new_y.x, new_y.y, new_y.z, 0.0;
                                   new_z.x, new_z.y, new_z.z, 0.0;
                                   0.0,     0.0,     0.0,     1.0];
        let view_matrix = matrix![1.0, 0.0, 0.0, -look_from.x;
                                  0.0, 1.0, 0.0, -look_from.y;
                                  0.0, 0.0, 1.0, -look_from.z;
                                  0.0, 0.0, 0.0, 1.0];
        // @TODO figure out, how this actually works, OpenGL tutorials immediately give more complicated camera
        // with fov, near/far clipping planes, etc. For now I just know, that dividing by 5.0 works ok.
        let coef = -1.0 / 5.0;
        let projection_matrix = matrix![1.0, 0.0, 0.0,  0.0;
                                        0.0, 1.0, 0.0,  0.0;
                                        0.0, 0.0, 1.0,  0.0;
                                        0.0, 0.0, coef, 1.0];
        // Viewport matrix depends only on constants.
        // Setting z-buffer resolution to 255.
        // Redef for convenience.
        let w = (width - 1) as f32;
        let h = (height - 1) as f32;
        let d = 255.;
        let viewport_matrix = matrix![w / 2.0, 0.0,     0.0,     w / 2.0;
                                      0.0,     h / 2.0, 0.0,     h / 2.0;
                                      0.0,     0.0,     d / 2.0, d / 2.0;
                                      0.0,     0.0,     0.0,     1.0];

        // Preparing shader pipeline for the render pass.
        buffer.vpmv_matrix = viewport_matrix *
                             projection_matrix * 
                             model_matrix * 
                             view_matrix;
        // Not interested in projection and rasterization, when transformaing light direction and normals.
        buffer.mv_matrix = model_matrix * view_matrix;
        buffer.it_mv_matrix = (model_matrix * view_matrix).transpose().try_inverse().unwrap();
        buffer.camera_direction = new_z;   
        buffer.t_light_direction = from_hom_vector(
            buffer.mv_matrix * to_hom_vector(light_direction)
        ).normalize();
    };

    return Box::new(prepare);
}

/// Pipeline preparation for the pass, where we want to get depth information as if our camera was placed at
/// the light source.
fn get_shadow_pass_prepare() -> Box<Prepare> {
    let prepare = |
        buffer:          &mut Buffer,
        width:           u32,
        height:          u32,
        light_direction: Vector3<f32>,
        look_from:       Vector3<f32>,
        look_at:         Vector3<f32>,
        up:              Vector3<f32>        
    | {
        // Light direction is inverted so we don't have to place annoying minuses in shaders, so we need
        // to invert it back, if it is to serve as a camera for the pass.
        let new_z = -light_direction.normalize();
        let new_y = (up - new_z.dot(&up) * new_z).normalize();
        let new_x = new_y.cross(&new_z).normalize();        
        let model_matrix = matrix![new_x.x, new_x.y, new_x.z, 0.0;
                                   new_y.x, new_y.y, new_y.z, 0.0;
                                   new_z.x, new_z.y, new_z.z, 0.0;
                                   0.0,     0.0,     0.0,     1.0];
        let view_matrix = matrix![1.0, 0.0, 0.0, light_direction.x;
                                  0.0, 1.0, 0.0, light_direction.y;
                                  0.0, 0.0, 1.0, light_direction.z;
                                  0.0, 0.0, 0.0, 1.0];
        // @TODO figure out, how this actually works, OpenGL tutorials immediately give more complicated camera
        // with fov, near/far clipping planes, etc. For now I just know, that dividing by 5.0 works ok.
        let coef = -1.0 / 5.0;
        let projection_matrix = matrix![1.0, 0.0, 0.0,  0.0;
                                        0.0, 1.0, 0.0,  0.0;
                                        0.0, 0.0, 1.0,  0.0;
                                        0.0, 0.0, coef, 1.0];
        // Viewport matrix depends only on constants.
        // Setting z-buffer resolution to 255.
        // Redef for convenience.
        let w = (width - 1) as f32;
        let h = (height - 1) as f32;
        let d = 255.;
        let viewport_matrix = matrix![w / 2.0, 0.0,     0.0,     w / 2.0;
                                      0.0,     h / 2.0, 0.0,     h / 2.0;
                                      0.0,     0.0,     d / 2.0, d / 2.0;
                                      0.0,     0.0,     0.0,     1.0];

        // Preparing shader pipeline for the render pass.
        buffer.vpmv_matrix = viewport_matrix *
                             projection_matrix * 
                             model_matrix * 
                             view_matrix;
        // Not interested in projection and rasterization, when transformaing light direction and normals.
        buffer.mv_matrix = model_matrix * view_matrix;
        buffer.it_mv_matrix = (model_matrix * view_matrix).transpose().try_inverse().unwrap();
        buffer.camera_direction = new_z;   
        buffer.t_light_direction = from_hom_vector(
            buffer.mv_matrix * to_hom_vector(light_direction)
        ).normalize();
    };

    return Box::new(prepare);
}

/// Calculating diffuse coefficient based on the face normal and light direction.
fn get_default_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        let vertex_positions = get_vertex_positions(model, pos_indices);
        if should_cull_face(vertex_positions, buffer.camera_direction) { return false; }

        // Calculating normal projection on the face.
        let face_normal = (vertex_positions[1] - vertex_positions[0]).cross(
            &(vertex_positions[2] - vertex_positions[0])
        );
        let t_face_normal = from_hom_vector(
            buffer.it_mv_matrix * to_hom_vector(face_normal)
        ).normalize();
        let diff_coef = buffer.t_light_direction.dot(&t_face_normal);
        buffer.vertex_intensities = vector![diff_coef, diff_coef, diff_coef];

        store_vertex_transformation_results(
            vertex_positions,
            buffer.vpmv_matrix, 
            &mut buffer.vertex_t_coords, 
            &mut buffer.vertex_z_values
        );
        store_vertex_uvs(&mut buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        let uv = buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = buffer.vertex_intensities[0];
        buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    };

    passes.push(ShaderPass {
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });

    return passes;
}

/// Calculating diffuse coefficient based on the interpolation of vertex normals in a particular fragment
/// via barycentric coordinates.
fn get_phong_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        let vertex_positions = get_vertex_positions(model, pos_indices);
        if should_cull_face(vertex_positions, buffer.camera_direction) { return false; }

        // Calculating light intensities at each vertex to then interpolate them in fragment shader.
        for i in 0..3 {
            let vertex_normal = vector![
                model.obj.normals[normal_indices[i]].0,
                model.obj.normals[normal_indices[i]].1,
                model.obj.normals[normal_indices[i]].2
            ];
            let vertex_t_normal = from_hom_vector(
                buffer.it_mv_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            buffer.vertex_intensities[i] = buffer.t_light_direction.dot(
                &vertex_t_normal
            );
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            buffer.vpmv_matrix, 
            &mut buffer.vertex_t_coords, 
            &mut buffer.vertex_z_values
        );
        store_vertex_uvs(&mut buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        let uv = buffer.vertex_uvs * bar_coord;
        let color = model.get_color_at_uv(uv);
        let diff_coef = bar_coord.dot(&buffer.vertex_intensities);
        buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    };

    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });

    return passes;
}

/// Calculating diffuse coefficient by looking up normal from normal map at the fragment uv.
fn get_normal_map_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        let vertex_positions = get_vertex_positions(model, pos_indices);
        if should_cull_face(vertex_positions, buffer.camera_direction) { return false; }
        
        store_vertex_transformation_results(
            vertex_positions,
            buffer.vpmv_matrix, 
            &mut buffer.vertex_t_coords, 
            &mut buffer.vertex_z_values
        );
        store_vertex_uvs(&mut buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        let uv = buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        let t_fragment_normal = from_hom_vector(
            buffer.it_mv_matrix * to_hom_vector(fragment_normal)
        ).normalize();
        let diff_coef = buffer.t_light_direction.dot(&t_fragment_normal);
        buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    };

    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });

    return passes;
}

/// Calculating diffuse coefficient as in the Phong pipeline, but also accounting for a specular component
/// by using specular map texture.
fn get_specular_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        let vertex_positions = get_vertex_positions(model, pos_indices);
        if should_cull_face(vertex_positions, buffer.camera_direction) { return false; }
        
        store_vertex_transformation_results(
            vertex_positions,
            buffer.vpmv_matrix, 
            &mut buffer.vertex_t_coords, 
            &mut buffer.vertex_z_values
        );
        store_vertex_uvs(&mut buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        let uv = buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal = model.get_normal_at_uv(uv);
        let t_fragment_normal = from_hom_vector(
            buffer.it_mv_matrix * to_hom_vector(fragment_normal)
        ).normalize();
        // Calculated reflection direction, immediately in a new camera frame.
        let reflected_t_light_direction = (
            2.0 * (
                t_fragment_normal * 
                buffer.t_light_direction.dot(&t_fragment_normal)
            ) - buffer.t_light_direction
        ).normalize();
        let diff_coef = buffer.t_light_direction.dot(&t_fragment_normal);
        // Accesing only .z, since in the new frame camera direction is always [0.0, 0.0, -1.0].
        let spec_coef = 0.6 * reflected_t_light_direction.z.max(0.0).powf(model.get_specular_value_at_uv(uv));
        let corrected_color = vector![
            ((diff_coef + spec_coef) * color[0] as f32).min(255.0) as u8,
            ((diff_coef + spec_coef) * color[1] as f32).min(255.0) as u8,
            ((diff_coef + spec_coef) * color[2] as f32).min(255.0) as u8
        ];
        buffer.fragment_color = corrected_color;

        return true;
    };

    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });

    return passes;
}

/// Calculating normal at a fragment by transforming normal from local Darboux basis to a global one.
fn get_darboux_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        let vertex_positions = get_vertex_positions(model, pos_indices);
        if should_cull_face(vertex_positions, buffer.camera_direction) { return false; }
        
        // Collecting transformed vertex positions in a buffer to use in local basis calculation.
        for i in 0..3 {
            buffer.vertex_t_positions.set_column(i, &from_hom_point(
                    buffer.mv_matrix * to_hom_point(vertex_positions[i])
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
                buffer.it_mv_matrix * to_hom_vector(vertex_normal)
            ).normalize();
            buffer.vertex_t_normals.set_column(i, &vertex_t_normal);
        }
        
        store_vertex_transformation_results(
            vertex_positions,
            buffer.vpmv_matrix, 
            &mut buffer.vertex_t_coords, 
            &mut buffer.vertex_z_values
        );
        store_vertex_uvs(&mut buffer.vertex_uvs, &model.obj.tex_coords, tex_indices);

        return true;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        let uv = buffer.vertex_uvs * bar_coord; 
        let color = model.get_color_at_uv(uv);
        let fragment_normal_tangent = model.get_normal_tangent_at_uv(uv);
        // Calculating the matrix, giving required transformation from Darboux basis to the global one.
        // A is a matrix with rows representing a local basis at each fragment.
        let mut local_basis_matrix: Matrix3<f32> = Default::default();
        let local_z = buffer.vertex_t_normals * bar_coord;
        local_basis_matrix.set_row(
            0, &(buffer.vertex_t_positions * vector![-1.0, 1.0, 0.0]).normalize().transpose()
        );
        local_basis_matrix.set_row(
            1, &(buffer.vertex_t_positions * vector![-1.0, 0.0, 1.0]).normalize().transpose()
        );
        local_basis_matrix.set_row(
            2, &(buffer.vertex_t_normals * bar_coord).normalize().transpose()
        );

        let i_local_basis_matrix = local_basis_matrix.try_inverse().unwrap();
        let local_x = i_local_basis_matrix * vector![
            buffer.vertex_uvs.m12 - buffer.vertex_uvs.m11,
            buffer.vertex_uvs.m13 - buffer.vertex_uvs.m11,
            0.0
        ];
        let local_y = i_local_basis_matrix * vector![
            buffer.vertex_uvs.m22 - buffer.vertex_uvs.m21,
            buffer.vertex_uvs.m23 - buffer.vertex_uvs.m21,
            0.0
        ];

        let mut local_transform_matrix: Matrix3<f32> = Default::default();
        local_transform_matrix.set_column(0, &local_x.normalize());
        local_transform_matrix.set_column(1, &local_y.normalize());
        local_transform_matrix.set_column(2, &local_z.normalize());
        let t_fragment_normal = (local_transform_matrix * fragment_normal_tangent).normalize();
        
        let diff_coef = buffer.t_light_direction.dot(&t_fragment_normal);
        buffer.fragment_color = color_blend(color, vector![0, 0, 0], diff_coef);

        return true;
    };

    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });

    return passes;
}

/// Two pass pipeline, doing render, placing camera at the light position and then using obtained z-buffer
/// to augment render result from the camera position.
fn get_shadow_pipeline_passes() -> Vec::<ShaderPass> {
    let mut passes = Vec::<ShaderPass>::new();

    let vertex_pass_1 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        return false;
    };

    let fragment_pass_1 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        return false;
    };

    let vertex_pass_2 = |
        buffer:         &mut Buffer,
        model:          &Model,
        pos_indices:    Vector3<usize>,
        tex_indices:    Vector3<usize>,
        normal_indices: Vector3<usize>
    | -> bool {
        return false;
    };

    let fragment_pass_2 = |
        buffer:    &mut Buffer,
        model:     &Model,
        bar_coord: Vector3<f32>,
        index:     usize
    | -> bool {
        if !process_z_value(buffer, bar_coord, index) { return false; }
        return false;
    };

    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_1), 
        fragment: Box::new(fragment_pass_1),
    });
    passes.push(ShaderPass { 
        prepare:  get_default_prepare(),
        vertex:   Box::new(vertex_pass_2), 
        fragment: Box::new(fragment_pass_2),
    });

    return passes;
}
