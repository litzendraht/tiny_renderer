use std::cmp::{min, max};

use obj::raw::RawObj;
use obj::raw::object::Polygon;
use image::{ImageBuffer, Rgb, RgbImage};
use nalgebra as na;
use na::{vector, Vector2, Vector3, matrix, Matrix2x3};

use crate::shader::ShaderPipeline;

/// Struct, holding all information about the model, including geometry, texture and normal and specular maps.
pub struct Model {
    pub obj: RawObj,
    pub texture: RgbImage,
    pub normal_map: RgbImage,
    pub spec_map: RgbImage,
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
        spec_map: RgbImage
    ) -> Self {
        let model = Model {
            obj,
            texture,
            normal_map,
            spec_map,
        };
        let shader_pipeline = T::new();
        let light_direction = vector![0.0, 0.0, 1.0];  // Directed to us from the screen.
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
        // New coordinate system around camera position, with z
        let y = up.normalize();
        let z = (look_from - look_at).normalize();
        let x = up.cross(&z).normalize();        
        let model_matrix = matrix![x.x, x.y, x.z, 0.0;
                                   y.x, y.y, y.z, 0.0;
                                   z.x, z.y, z.z, 0.0;
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
        buffer.light_direction = self.light_direction;
        buffer.vertex_transform_matrix = viewport_matrix *
                                         projection_matrix * 
                                         model_matrix * 
                                         view_matrix;
        buffer.direction_transform_matrix = model_matrix * view_matrix;
        // @TODO something is not quite right here, light direction seems to move in relation 
        // to the model, when it shouldn't. 
        buffer.it_direction_transform_matrix = ( 
            model_matrix * 
            view_matrix
        ).transpose().try_inverse().unwrap();
    }

    /// Sets Scene pixel to a color at specifed coordinate.
    /// 
    /// Assumes, that pixel data is rgb8.
    pub fn set_pixel(&mut self, p: Vector2<i32>, color: Vector3<u8>) {
        // Pixel data is rgb8, so we find the starting pixel_index of a 3-tuple and do 3 assignments.
        // @OPTI maybe I can set 3 values simultaneously somehow?
        let pixel_index = (3 * (p.x + p.y * self.width as i32)) as usize;
        self.render_data[pixel_index + 0] = color.x;
        self.render_data[pixel_index + 1] = color.y;
        self.render_data[pixel_index + 2] = color.z;
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
                // Vertex shader decided, that whole polygon wouldn't be rendered.
                continue;
            }

            let vertex_transformed_coords = self.shader_pipeline.get_buffer().vertex_transformed_coords;
            let bbox = get_triangle_bounding_box(vertex_transformed_coords);

            // Accounting for possibility that bbox can reach outside of the screen.
            for i in max(0, bbox.ll.x)..=min(bbox.ur.x, (self.width - 1) as i32) {
                for j in max(0, bbox.ll.y)..=min(bbox.ur.y, (self.height - 1) as i32) {
                    let pixel_index = (i + j * self.width as i32) as usize;
                    let bar_coord = to_barycentric_coord(
                        vector![i, j], 
                        vertex_transformed_coords
                    );

                    // If any of the coordinates are negative, point is not in the triangle, so skipping it.
                    if bar_coord.x < 0.0 || bar_coord.y < 0.0 || bar_coord.z < 0.0 {
                        continue;
                    }

                    let vertex_z_values = self.shader_pipeline.get_buffer().vertex_z_values;
                    let z_to_screen = bar_coord.dot(&vertex_z_values);                    

                    // Checking z-buffer, on failure skipping the pixel.
                    if z_to_screen <= self.z_buffer[pixel_index] {
                        continue;
                    }

                    // z-buffer check is passed, so setting the pixel and new z-buffer value.
                    self.z_buffer[pixel_index] = z_to_screen;

                    // Fragment shader has an ability to skip a pixel.
                    if !self.shader_pipeline.fragment(&self.model, bar_coord) {
                        continue;
                    }

                    // If fragment shader returned true, taking fragment color form the buffer and coloring
                    // in the appropriate pixel.
                    let fragment_color = self.shader_pipeline.get_buffer().fragment_color;
                    self.render_data[3 * pixel_index + 0] = fragment_color.x;
                    self.render_data[3 * pixel_index + 1] = fragment_color.y;
                    self.render_data[3 * pixel_index + 2] = fragment_color.z;
                }
            }
        }
    }
}
