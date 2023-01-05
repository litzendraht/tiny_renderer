mod util;
use util::{to_hom_vector, from_hom_vector};
pub mod shader;
pub use self::shader::{ShaderPipeline, DefaultSP, PhongSP, TrueNormalSP, SpecularSP, DarbouxSP};

use std::cmp::{min, max};

use obj::raw::RawObj;
use obj::raw::object::Polygon;
use image::{ImageBuffer, Rgb, RgbImage};
use nalgebra as na;
use na::{vector, Vector2, Vector3, matrix, Matrix2x3};

/// Struct, holding all information about the model, including geometry, texture and normal and specular maps.
pub struct Model {
    pub obj: RawObj,
    pub texture: RgbImage,
    pub normal_map: RgbImage,
    pub normal_map_tangent: RgbImage,
    pub specular_map: RgbImage,
}

impl Model {
    pub fn get_vertex_position_at_index(&self, index: usize) -> Vector3<f32> {
        return vector![
            self.obj.positions[index].0,
            self.obj.positions[index].1,
            self.obj.positions[index].2
        ];
    }

    /// Returns texture color from texture file at uv.
    pub fn get_color_at_uv(&self, uv: Vector2<f32>) -> Vector3<u8> {
        let coord = vector![
            (uv.x * self.texture.width() as f32) as u32,
            (uv.y * self.texture.height() as f32) as u32
        ];

        return Vector3::<u8>::from_row_slice(
            &self.texture.get_pixel(coord.x, coord.y).0[0..3]
        );
    }

    /// Returns normalized normal from normal map at uv.
    pub fn get_normal_at_uv(&self, uv: Vector2<f32>) -> Vector3<f32> {
        let coord = vector![
            (uv.x * self.normal_map.width() as f32) as u32,
            (uv.y * self.normal_map.height() as f32) as u32
        ];

        // Subtracting 0.5 to get from [0, 255] to [-0.5, 0.5]
        return vector![
            (self.normal_map.get_pixel(coord.x, coord.y).0[0]) as f32 / 255.0 - 0.5,
            (self.normal_map.get_pixel(coord.x, coord.y).0[1]) as f32 / 255.0 - 0.5,
            (self.normal_map.get_pixel(coord.x, coord.y).0[2]) as f32 / 255.0 - 0.5
        ].normalize();
    }

    /// Returns normalized normal from normal map in tangent coordinates at uv.
    pub fn get_normal_tangent_at_uv(&self, uv: Vector2<f32>) -> Vector3<f32> {
        let coord = vector![
            (uv.x * self.normal_map.width() as f32) as u32,
            (uv.y * self.normal_map.height() as f32) as u32
        ];

        // Subtracting 0.5 to get from [0, 255] to [-0.5, 0.5]
        return vector![
            (self.normal_map_tangent.get_pixel(coord.x, coord.y).0[0]) as f32 / 255.0 - 0.5,
            (self.normal_map_tangent.get_pixel(coord.x, coord.y).0[1]) as f32 / 255.0 - 0.5,
            (self.normal_map_tangent.get_pixel(coord.x, coord.y).0[2]) as f32 / 255.0 - 0.5
        ].normalize();
    }

    /// Returns a [0, 1] from specular mpa at uv.
    pub fn get_specular_value_at_uv(&self, uv: Vector2<f32>) -> f32 {
        let coord = vector![
            (uv.x * self.specular_map.width() as f32) as u32,
            (uv.y * self.specular_map.height() as f32) as u32
        ];

        return self.specular_map.get_pixel(coord.x, coord.y).0[0] as f32;
    }
}

/// Scene, holding its width, height and private flat array(vec) of pixel data,
/// showing the rendered image.
/// (0, 0) is the bottom left coordinate.
pub struct Scene<T> where 
    T: ShaderPipeline 
{
    width: u32,
    height: u32,
    model: Model,
    // Pipeline, specifying vertex and fragment shaders
    shader_pipeline: T,
    // Scene settings.
    light_direction: Vector3<f32>,
    // z-buffer, which continuously fills out after clear() call with every new primitive drawn.
    // Takes f32 values in [0, 255] for ease of debugging.
    z_buffer: Vec<f32>,
    // u8 vesrion of z-buffer directly passed to image_show.
    depth_data: Vec<u8>,
    // Storing flat array.
    render_data: Vec<u8>,
}

impl<T> Scene<T> where
    T: ShaderPipeline
{
    /// Generates new Scene struct with specified width and height.
    /// Pixel data format is assumed to be rgb8.
    pub fn new(
        width: u32, 
        height: u32,
        obj: RawObj,
        texture: RgbImage,
        normal_map: RgbImage,
        normal_map_tangent: RgbImage,
        specular_map: RgbImage
    ) -> Self {
        let model = Model {
            obj,
            texture,
            normal_map,
            normal_map_tangent,
            specular_map,
        };
        let shader_pipeline = T::new();
        let light_direction = vector![0.0, 0.0, -1.0];  // Directed from us to the screen.
        let n_pixels = (width * height) as usize;
        let z_buffer: Vec<f32>   = vec![f32::MIN; n_pixels];
        let depth_data: Vec<u8>  = vec![0; 3 * n_pixels];
        let render_data: Vec<u8> = vec![0; 3 * n_pixels];
        return Scene {
            width,
            height,
            model,
            shader_pipeline,
            light_direction,
            z_buffer,
            depth_data,
            render_data,
        }
    }
    
    /// Get rendered scene as a slice of color values of size 3 * (number of pixels).
    /// Flips the image, so (0, 0) is the lower left corner.
    pub fn get_render_data(&self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        let mut buffer: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_vec(
            self.width, self.height, self.render_data.clone()
        ).unwrap();
        image::imageops::flip_vertical_in_place(&mut buffer);
        return buffer;
    }

    /// Get image, representing z-buffer values.
    /// Lazy in a sense, that color data for the image is calculated only if this call is made.
    // @TODO figure out, why it doesn't want to work.
    pub fn get_depth_data(&mut self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        for i in 0..self.z_buffer.len() {
            self.depth_data[3 * i + 0] = self.z_buffer[i] as u8;
            self.depth_data[3 * i + 1] = self.z_buffer[i] as u8;
            self.depth_data[3 * i + 2] = self.z_buffer[i] as u8;
        }
        let mut buffer: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_vec(
            self.width, self.height, self.depth_data.clone()
        ).unwrap();
        image::imageops::flip_vertical_in_place(&mut buffer);
        return buffer;
    }

    /// Setter for light direction.
    pub fn set_light_direction(&mut self, light_direction: Vector3<f32>) {
        self.light_direction = light_direction;
    }

    /// Sets all rendered pixels data to (0, 0, 0) and clears z-buffer.
    pub fn clear(&mut self) {
        let n_pixels = (self.width * self.height) as usize;
        for i in 0..n_pixels {
            self.z_buffer[i] = f32::MIN;
            self.render_data[3 * i + 0] = 0;
            self.render_data[3 * i + 1] = 0;
            self.render_data[3 * i + 2] = 0;
        }
    }

    /// Recalculating model, view, projection and viewport matrices to be used when drawing primitives.
    /// look_from is camera position, look_at defines directon.
    pub fn prepare_render(&mut self, look_from: Vector3<f32>, look_at: Vector3<f32>, up: Vector3<f32>) {
        // New coordinate system a, b, c around camera position.
        // @TODO figure out consistent naming for basis vectors.
        let c = (look_from - look_at).normalize();
        let b = (up - c.dot(&up) * c).normalize();
        let a = b.cross(&c).normalize();        
        let model_matrix = matrix![a.x, a.y, a.z, 0.0;
                                   b.x, b.y, b.z, 0.0;
                                   c.x, c.y, c.z, 0.0;
                                   0.0, 0.0, 0.0, 1.0];
        let view_matrix = matrix![1.0, 0.0, 0.0, -look_from.x;
                                  0.0, 1.0, 0.0, -look_from.y;
                                  0.0, 0.0, 1.0, -look_from.z;
                                  0.0, 0.0, 0.0, 1.0];
        // @TODO figure out, how this actually works, OpenGL tutorials immediately give more complicated camera
        // with fov, near/far clipping planes, etc. For now I just know, that dividing by 5.0 work ok.
        let coef = -1.0 / 5.0;
        let projection_matrix = matrix![1.0, 0.0, 0.0,  0.0;
                                        0.0, 1.0, 0.0,  0.0;
                                        0.0, 0.0, 1.0,  0.0;
                                        0.0, 0.0, coef, 1.0];
        // Viewport matrix depends only on constants.
        // Setting z-buffer resolution to 255.
        // Redef for convenience.
        let w = (self.width - 1) as f32;
        let h = (self.height - 1) as f32;
        let d = 255.;
        let viewport_matrix = matrix![w / 2.0, 0.0,     0.0,     w / 2.0;
                                      0.0,     h / 2.0, 0.0,     h / 2.0;
                                      0.0,     0.0,     d / 2.0, d / 2.0;
                                      0.0,     0.0,     0.0,     1.0];

        // Preparing shader pipeline for the render pass.
        let buffer = self.shader_pipeline.get_buffer_mut();
        buffer.vpmv_matrix = viewport_matrix *
                                         projection_matrix * 
                                         model_matrix * 
                                         view_matrix;
        buffer.mv_matrix = model_matrix * view_matrix;
        // Not interested in projection and rasterization, when transformaing light direction and normals.
        buffer.it_mv_matrix = (model_matrix * view_matrix).transpose().try_inverse().unwrap();       
        buffer.t_light_direction = from_hom_vector(
            buffer.mv_matrix * to_hom_vector(self.light_direction)
        ).normalize();
    }

    pub fn render(&mut self) {
        // Simple local bounding box struct for convenience.
        #[derive(Debug)]
        struct BoundingBox {
            ll: Vector2<i32>, // lower left corner
            ur: Vector2<i32>, // upper right corner
        }

        // Helper used to find bounding box of a triangle.
        fn get_triangle_bounding_box(coords: Matrix2x3<i32>) -> BoundingBox {
            return BoundingBox {
                ll: vector![
                    min(min(coords.m11, coords.m12), coords.m13),
                    min(min(coords.m21, coords.m22), coords.m23)
                ],
                ur: vector![
                    max(max(coords.m11, coords.m12), coords.m13),
                    max(max(coords.m21, coords.m22), coords.m23)
                ]
            };
        }

        // Getting barycentric coordinates for a point in relation to a rasterized triangle coordinates.
        fn to_barycentric_coord(
            internal_point: Vector2<i32>, 
            coords: Matrix2x3<i32>
        ) -> Vector3<f32> {
            let raw_cross = vector![
                    (coords.m12 - coords.m11) as f32,
                    (coords.m13 - coords.m11) as f32,
                    (coords.m11 - internal_point.x) as f32
                ].cross(&vector![
                    (coords.m22 - coords.m21) as f32,
                    (coords.m23 - coords.m21) as f32,
                    (coords.m21 - internal_point.y) as f32
                ]);
            if raw_cross.z.abs() < 1.0 {
                // Degenerate triangle, returning something with negative coordinate.
                return vector![-1.0, 1.0, 1.0];
            }
            return vector![
                1.0 - (raw_cross.x + raw_cross.y) / raw_cross.z,
                raw_cross.x / raw_cross.z,
                raw_cross.y / raw_cross.z
            ];
        }

        // Drawing all polygons of the model.
        for polygon in &self.model.obj.polygons {
            // Indices are &Vec((usize, usize, usize)), where first item corresponds to indices for
            // positions, second to indices for texture uv coords and third to indices for normals
            // which results in a bloated call to vertex shader.
            // @TODO make it not bloated.
            let indices: &Vec<(usize, usize, usize)> = match polygon {
                Polygon::PTN(indices) => &indices,
                _ => panic!("Encountered some garbage, while looking through polygons."),
            };

            if !self.shader_pipeline.vertex(
                &self.model,
                vector![indices[0].0, indices[1].0, indices[2].0],
                vector![indices[0].1, indices[1].1, indices[2].1],
                vector![indices[0].2, indices[1].2, indices[2].2]               
            ) {
                // Vertex shader decided, that whole polygon shouldn't be rendered.
                continue;
            }

            let vertex_t_coords = self.shader_pipeline.get_buffer().vertex_t_coords;
            let bbox = get_triangle_bounding_box(vertex_t_coords);

            // Accounting for possibility that bbox can reach outside of the screen.
            for i in max(0, bbox.ll.x)..=min(bbox.ur.x, (self.width - 1) as i32) {
                for j in max(0, bbox.ll.y)..=min(bbox.ur.y, (self.height - 1) as i32) {
                    let pixel_index = (i + j * self.width as i32) as usize;
                    let bar_coord = to_barycentric_coord(
                        vector![i, j], 
                        vertex_t_coords
                    );

                    // If any of the coordinates are negative, point is not in the triangle, so skipping it.
                    if bar_coord.x < 0.0 || bar_coord.y < 0.0 || bar_coord.z < 0.0 {
                        continue;
                    }

                    // z-buffer section - checking fragment z-value in the pipeline buffer and comparing it
                    // to the value in the buffer, on failure skipping the pixel, else updating buffer and
                    // invoking the fragment part of the pipeline.
                    let vertex_z_values = self.shader_pipeline.get_buffer().vertex_z_values;
                    let z_to_screen = bar_coord.dot(&vertex_z_values); 
                    if z_to_screen <= self.z_buffer[pixel_index] {
                        continue;
                    }
                    self.z_buffer[pixel_index] = z_to_screen;

                    // If fragment shader returns true, getting color from the pipeline and coloring the
                    // pixel.
                    if !self.shader_pipeline.fragment(&self.model, bar_coord) {
                        continue;
                    }
                    let fragment_color = self.shader_pipeline.get_buffer().fragment_color;
                    self.render_data[3 * pixel_index + 0] = fragment_color.x;
                    self.render_data[3 * pixel_index + 1] = fragment_color.y;
                    self.render_data[3 * pixel_index + 2] = fragment_color.z;
                }
            }
        }
    }
}
