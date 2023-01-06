/// Buffer for small things, like passing values between different stages of a shader pipeline, 
/// setting up frame constants like light direction and transform matrices.
// @TODO not all buffer is used in every shader, but for simplicity it just has all options. Maybe
// worth thinking about how to get a more fine-grained buffer structure for each shader.
#[derive(Default)]
pub struct Buffer {
    pub camera_direction:  Vector3<f32>,
    pub t_light_direction: Vector3<f32>,   // Light direction with model and view transformations applied.
    pub vpmv_matrix:       Matrix4<f32>,   // Applied to vertices to get final screen coodrdinates.
    pub mv_matrix:         Matrix4<f32>,   // Applied to light direction.
    pub it_mv_matrix:      Matrix4<f32>,   // Applied to model normals.
    // Big ol' fat z-buffer
    pub z_buffer:          Vec<f32>,
    // Big ol' fat shadow buffer
    pub shadow_buffer:     Vec<f32>,
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
    fn new(buffer_size: usize) -> Self {
        return Self {
            z_buffer:      vec![0.0; buffer_size],
            shadow_buffer: vec![0.0; buffer_size],
            ..Default::default()
        }
    }
}