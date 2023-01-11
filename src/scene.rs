// @TODO similarly to shader.rs crate, this crate, which is closely coupled to it is also hot garbage,
// requiring some refactoring.

mod util;
mod shader;

use util::Model;
use self::shader::ShaderPipeline;

use std::{cmp::{min, max}, thread::available_parallelism};

use obj::raw::RawObj;
use obj::raw::object::Polygon;
use image::{ImageBuffer, Rgb, RgbImage};
use nalgebra as na;
use na::{vector, Vector2, Vector3, Matrix2x3};
use threadpool::ThreadPool;

/// Scene, holding its width, height and private flat array(vec) of pixel data,
/// showing the rendered image.
/// (0, 0) is the bottom left coordinate.
pub struct Scene
{
    width:           u32,
    height:          u32,
    model:           Model,
    // Pipeline, specifying vertex and fragment shaders
    shader_pipeline: ShaderPipeline,
    // Lighting and camera settings.
    light_direction: Vector3<f32>,
    look_from:       Vector3<f32>,
    look_at:         Vector3<f32>,
    up:              Vector3<f32>,
    // u8 version of z-buffer directly passed to image_show.
    depth_data:      Vec<u8>,
    // Storing flat array.
    frame_buffer:    Vec<u8>,
    // Threadpool for multi-threaded fragment shader execution.
    thread_pool:     ThreadPool,
    
}

impl Scene
{
    /// Generates new Scene struct with specified width and height.
    /// Pixel data format is assumed to be rgb8.
    pub fn new(
        width:                u32, 
        height:               u32,
        obj:                  RawObj,
        texture:              RgbImage,
        normal_map:           RgbImage,
        normal_map_tangent:   RgbImage,
        specular_map:         RgbImage,
        shader_pipeline_name: String
    ) -> Self {
        let model = Model {
            obj,
            texture,
            normal_map,
            normal_map_tangent,
            specular_map,
        };
        let frame_buffer_size = (width * height) as usize;
        let shader_pipeline = ShaderPipeline::new(shader_pipeline_name, width, height);
        let light_direction = vector![0.0, 0.0, -1.0];   
        let look_from       = vector![0.0, 0.0, 1.0];
        let look_at         = vector![0.0, 0.0, 0.0];
        let up              = vector![0.0, 1.0, 0.0];
        let depth_data:   Vec<u8> = vec![0; 3 * frame_buffer_size];
        let frame_buffer: Vec<u8> = vec![0; 3 * frame_buffer_size];
        let n_threads = available_parallelism().unwrap().get();
        println!("scene is creating thread pool with {} threads", n_threads);
        let thread_pool     = ThreadPool::new(n_threads);
        return Scene {
            width,
            height,
            model,
            shader_pipeline,
            light_direction,
            look_from,
            look_at,
            up,
            depth_data,
            frame_buffer,
            thread_pool,
        }
    }
    
    /// Get rendered scene as a slice of color values of size 3 * (number of pixels).
    /// Flips the image, so (0, 0) is the lower left corner.
    pub fn get_frame_buffer(&self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        let mut buffer: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_vec(
            self.width, self.height, self.frame_buffer.clone()
        ).unwrap();
        image::imageops::flip_vertical_in_place(&mut buffer);
        return buffer;
    }

    /// Get image, representing z-buffer values.
    /// Lazy in a sense, that color data for the image is calculated only if this call is made.
    pub fn get_z_buffer(&mut self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        for i in 0..self.shader_pipeline.buffer.z_buffer.len() {
            self.depth_data[3 * i + 0] = self.shader_pipeline.buffer.z_buffer[i] as u8;
            self.depth_data[3 * i + 1] = self.shader_pipeline.buffer.z_buffer[i] as u8;
            self.depth_data[3 * i + 2] = self.shader_pipeline.buffer.z_buffer[i] as u8;
        }
        let mut buffer: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_vec(
            self.width, self.height, self.depth_data.clone()
        ).unwrap();
        image::imageops::flip_vertical_in_place(&mut buffer);
        return buffer;
    }

    /// Get image, representing shadow-buffer values.
    /// Lazy in a sense, that color data for the image is calculated only if this call is made.
    pub fn get_shadow_buffer(&mut self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        for i in 0..self.shader_pipeline.buffer.shadow_buffer.len() {
            self.depth_data[3 * i + 0] = self.shader_pipeline.buffer.shadow_buffer[i] as u8;
            self.depth_data[3 * i + 1] = self.shader_pipeline.buffer.shadow_buffer[i] as u8;
            self.depth_data[3 * i + 2] = self.shader_pipeline.buffer.shadow_buffer[i] as u8;
        }
        let mut buffer: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_vec(
            self.width, self.height, self.depth_data.clone()
        ).unwrap();
        image::imageops::flip_vertical_in_place(&mut buffer);
        return buffer;
    }

    /// Sets all rendered pixels data to (0, 0, 0) and clears z-buffer.
    pub fn clear(&mut self) {
        let frame_buffer_size = (self.width * self.height) as usize;
        for i in 0..frame_buffer_size {
            self.shader_pipeline.buffer.z_buffer[i] = f32::MIN;
            self.shader_pipeline.buffer.shadow_buffer[i] = f32::MIN;
            self.frame_buffer[3 * i + 0] = 0;
            self.frame_buffer[3 * i + 1] = 0;
            self.frame_buffer[3 * i + 2] = 0;
        }
    }

    /// Settign light parameters for the scene.
    pub fn set_light_direction(&mut self, light_direction: Vector3<f32>) {
        self.light_direction = light_direction;
    }

    /// Setting camera parameters for the scene,
    pub fn set_camera(&mut self, look_from: Vector3<f32>, look_at: Vector3<f32>, up: Vector3<f32>) {
        self.look_from = look_from;
        self.look_at   = look_at;
        self.up        = up;
    }

    pub fn render(&mut self) {
        // Simple local bounding box struct for convenience.
        #[derive(Debug)]
        struct BoundingBox {
            ll: Vector2<i32>, // lower left corner
            ur: Vector2<i32>, // upper right corner
        }

        // Helper used to find bounding box of a triangle. Can reach outside of the screen.
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

        // Applying all passes of the shader pipeline.
        for pass in &self.shader_pipeline.passes {
            // Filling the buffer with pass constants.
            (pass.prepare)(
                &mut self.shader_pipeline.buffer,
                self.width,
                self.height,
                self.light_direction,
                self.look_from,
                self.look_at, 
                self.up                
            );
            // Drawing all polygons of the model.
            for polygon in &self.model.obj.polygons {
                // Indices are &Vec((usize, usize, usize)), where first item corresponds to indices for
                // positions, second to indices for texture uv coords and third to indices for normals
                // which results in a bloated call to vertex shader.
                let indices: &Vec<(usize, usize, usize)> = match polygon {
                    Polygon::PTN(indices) => &indices,
                    _ => panic!("Encountered some garbage, while looking through polygons."),
                };

                if !(pass.vertex)(
                    &mut self.shader_pipeline.buffer,
                    &self.model,
                    vector![indices[0].0, indices[1].0, indices[2].0],
                    vector![indices[0].1, indices[1].1, indices[2].1],
                    vector![indices[0].2, indices[1].2, indices[2].2]               
                ) {
                    // Vertex shader decided, that whole polygon shouldn't be rendered.
                    continue;
                }

                let vertex_t_raster = self.shader_pipeline.buffer.vertex_t_raster;
                let bbox = get_triangle_bounding_box(vertex_t_raster);

                // Accounting for possibility that bbox can reach outside of the screen.
                let x_min = max(0, bbox.ll.x);
                let x_max = min(bbox.ur.x, (self.width - 1) as i32);
                let y_min = max(0, bbox.ll.y);
                let y_max = min(bbox.ur.y, (self.height - 1) as i32);
                for i in x_min..=x_max {
                    for j in y_min..=y_max {
                        let bar_coord = to_barycentric_coord(
                            vector![i, j], 
                            vertex_t_raster
                        );

                        // If any of the coordinates are negative, point is not in the triangle, so skipping it.
                        if bar_coord.x < 0.0 || bar_coord.y < 0.0 || bar_coord.z < 0.0 {
                            continue;
                        }                        

                        // If fragment shader returns true, getting color from the pipeline and coloring the
                        // pixel, else skipping the pixel.
                        if !(pass.fragment)(
                            &mut self.shader_pipeline.buffer, 
                            &self.model,
                            vector![i as u32, j as u32],
                            bar_coord
                        ) {
                            continue;
                        }
                        let fragment_color = self.shader_pipeline.buffer.fragment_color;
                        let pixel_index = (i + j * self.width as i32) as usize;
                        self.frame_buffer[3 * pixel_index + 0] = fragment_color.x;
                        self.frame_buffer[3 * pixel_index + 1] = fragment_color.y;
                        self.frame_buffer[3 * pixel_index + 2] = fragment_color.z;
                    }
                }
            }
        }
    }
}
