use std::ops::Bound;
use std::cmp::{min, max};
use image::{ImageBuffer, Rgb};

use crate::util::{Vector2i, Vector2f, Vector3f};

/// Struct, representing raw rgb8 pixel data.
#[derive(Clone, Copy)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

// @TRASH double definition
const WHITE: Color = Color { r: 255, g: 255, b: 255, };
const BLACK: Color = Color { r: 0,   g: 0,   b: 0,   };
const RED:   Color = Color { r: 255, g: 0,   b: 0,   };
const GREEN: Color = Color { r: 0,   g: 255, b: 0,   };
const BLUE:  Color = Color { r: 0,   g: 0,   b: 255, };

impl Color {
    /// Get convex combination of two colors: t * c_1 + (1 - t) * c_2.
    /// t is unrestricted.
    pub fn blend(color_1: Color, color_2: Color, t: f32) -> Color {
        return Color {
            r: (t * color_1.r as f32 + (1.0 - t) * color_2.r as f32) as u8,
            g: (t * color_1.g as f32 + (1.0 - t) * color_2.g as f32) as u8,
            b: (t * color_1.b as f32 + (1.0 - t) * color_2.b as f32) as u8,
        }
    }
}

/// Point in a scene - x, y give pixel index and z gives distance to the camera.
#[derive(Debug, Clone, Copy)]
struct ScenePoint {
    x: i32,
    y: i32,
    z: f32,
}

/// Scene, holding its width, height and private flat array(vec) of pixel data,
/// showing the rendered image.
/// (0, 0) is the bottom left coordinate.
pub struct Scene {
    pub width: u32,
    pub height: u32,
    z_buffer: Vec<f32>,    // z-buffer, which continuously fills out after clear() call with every new primitive drawn.
    depth_data: Vec<u8>,   // Normalized values of the z-buffer for visualization.
    render_data: Vec<u8>,  // Storing flat array.
}

impl Scene {
    /// Generates new Scene struct with specified width and height.
    /// Pixel data format is assumed to be rgb8.
    pub fn new(width: u32, height: u32) -> Scene {
        let n_pixels = (width * height) as usize;
        let z_buffer: Vec<f32>   = vec![f32::MIN; n_pixels];
        let depth_data: Vec<u8>  = vec![0; 3 * n_pixels];
        let render_data: Vec<u8> = vec![0; 3 * n_pixels];
        return Scene {
            width,
            height,
            z_buffer,
            depth_data,
            render_data,
        }
    }
    
    /// Get rendered scene as a slice of color values of size 3 * (number of pixels).
    pub fn as_render_data(&self) -> &[u8] {
        return &self.render_data[..];
    }

    /// Get image, representing z-buffer values.
    /// Lazy in a sense, that color data for the image is calculated only if this call is made.
    pub fn as_depth_data(&mut self) -> &[u8] {
        let z_max = self.z_buffer.iter().fold(f32::MIN, |max_value, value| value.max(max_value));
        let z_min = self.z_buffer.iter().fold(f32::MAX, |min_value, value| value.min(min_value));
        let scale = z_max - z_min;
        for i in 0..self.z_buffer.len() {
            let scaled_z = ((self.z_buffer[i] + z_min) / scale) * 255.0;
            self.depth_data[3 * i + 0] = scaled_z as u8;
            self.depth_data[3 * i + 1] = scaled_z as u8;
            self.depth_data[3 * i + 2] = scaled_z as u8;
        }
        return &self.depth_data[..];
    }

    /// Sets all rendered pixels data to (0, 0, 0) and clears z-buffer.
    pub fn clear(&mut self) {
        let capacity = (self.width * self.height) as usize;
        for i in 0..capacity {
            self.z_buffer[i] = f32::MIN;
            self.render_data[3 * i + 0] = 0;
            self.render_data[3 * i + 1] = 0;
            self.render_data[3 * i + 2] = 0;
        }
    }

    /// Checking if coordinate is in Scene bounds.
    fn coord_in_bounds(&self, v: Vector2i) -> bool {
        if v.x > self.width as i32 {
            return false;
        }
        if v.y > self.height as i32 {
            return false;
        }
        return true;
    }

    /// Tranfromation of Vector3 with x, y in [-1.0, 1.0] to ScenePoint.
    fn to_scene_point(&self, v: Vector3f) -> ScenePoint {
        // Transformtaion of a float in \[-1.0, 1.0\] to the pixel image coordinate
        // in range \[0, scale - 1\].
        fn to_rasterized_coord(float_coord: f32, scale: u32) -> i32 {
            return ((float_coord + 1.0) * ((scale - 1) as f32) / 2.0) as i32;
}
        return ScenePoint {
            x: to_rasterized_coord(v.x, self.width),
            y: to_rasterized_coord(v.y, self.height),
            z: v.z,
        }
    }

    /// Sets Scene pixel to a color at specifed coordinate.
    /// 
    /// Assumes, that pixel data is rgb8.
    pub fn set_pixel(&mut self, p: Vector2i, color: Color) {
        // Pixel data is rgb8, so we find the starting index of a 3-tuple and do 3 assignments.
        // @OPTI maybe I can set 3 values simultaneously somehow?
        // @TODO do something with y inversion, seems very hacky here, maybe need to explicitly
        // invert the image, when providing a slice in as_render_data.
        let index = (3 * (p.x + (self.height as i32 - 1 - p.y) * self.width as i32)) as usize;
        self.render_data[index + 0] = color.r;
        self.render_data[index + 1] = color.g;
        self.render_data[index + 2] = color.b;
    }

    /// Draws a line between v_1 and v_2 coordinates with specified color
    /// via Bresenham's algorithm as presented in https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    /// Draws line over anything, but consequent primitive draw calls can disrespct it, since z-buffer doesn't change.
    // @CLEAN strange function, maybe need to standardize all draw functions to have z_buffer_ignore param.
    pub fn draw_line_z_ignore(&mut self, a: Vector3f, b: Vector3f, color: Color) {
        let scene_point_a = self.to_scene_point(a);
        let scene_point_b = self.to_scene_point(b);
        let mut x_0 = scene_point_a.x;
        let x_1 = scene_point_b.x;
        let mut y_0 = scene_point_a.y;
        let y_1 = scene_point_b.y;
        let dx: i32 = (x_1 - x_0).abs();
        let sx: i8 = match x_0 < x_1 {
            true  =>  1,
            false => -1,
        };
        let dy: i32 = -(y_1 - y_0).abs();
        let sy: i8 = match y_0 < y_1 {
            true  =>  1,
            false => -1,
        };
        let mut error: i32 = dx + dy;

        let mut e2;
        loop {
            self.set_pixel(Vector2i { x: x_0, y: y_0 }, color);
            if (x_0 == x_1) && (y_0 == y_1) {
                break;
            }
            e2 = 2 * error;
            if e2 >= dy {
                if x_0 == x_1 {
                    break;
                }
                error += dy;
                x_0 += sx as i32;
            }
            if e2 <= dx {
                if y_1 == y_0 {
                    break;
                }
                error += dx;
                y_0 += sy as i32;
            }
        }
    }

    pub fn draw_point(&mut self, a: Vector3f, color: Color) {
        let scene_point = self.to_scene_point(a);
        if scene_point.z > self.z_buffer[(scene_point.x + scene_point.y * self.width as i32) as usize] {
            // Drawing over the point that is further away.
            self.set_pixel(Vector2i { x: scene_point.x, y: scene_point.y }, color);
        }
    }

    /// Draw a triangle with vertices at specified coordiantes and filled with specified color.
    pub fn draw_triangle(&mut self, a: Vector3f, b: Vector3f, c: Vector3f, color: Color, normal_correction_coef: f32) {
        // Simple local bounding box struct for convenience.
        #[derive(Debug)]
        struct BoundingBox {
            ll: Vector2i, // lower left corner
            ur: Vector2i, // upper right corner
        }

        // Helper used to find bounding box of a triangle.
        fn get_triangle_bounding_box(coord_a: Vector2i, coord_b: Vector2i, coord_c: Vector2i) -> BoundingBox {
            return BoundingBox {
                ll: Vector2i {
                    x: min(min(coord_a.x, coord_b.x), coord_c.x),
                    y: min(min(coord_a.y, coord_b.y), coord_c.y),
                },
                ur: Vector2i {
                    x: max(max(coord_a.x, coord_b.x), coord_c.x),
                    y: max(max(coord_a.y, coord_b.y), coord_c.y),
                }
            }
        }

        fn to_barycentric_coord(coord_point: Vector2i, coord_a: Vector2i, coord_b: Vector2i, coord_c: Vector2i) -> Vector3f {
            let raw_cross = Vector3f::cross(
                Vector3f {
                    x: (coord_b.x - coord_a.x) as f32,
                    y: (coord_c.x - coord_a.x) as f32,
                    z: (coord_a.x - coord_point.x) as f32,
                },
                Vector3f {
                    x: (coord_b.y - coord_a.y) as f32,
                    y: (coord_c.y - coord_a.y) as f32,
                    z: (coord_a.y - coord_point.y) as f32,
                }
            );
            if raw_cross.z.abs() < 1.0 {
                // Degenerate triangle, returning something with negative coordinate.
                return Vector3f { x: -1.0, y: 1.0, z: 1.0};
            }
            return Vector3f {
                x: 1.0 - (raw_cross.x + raw_cross.y) / raw_cross.z,
                y: raw_cross.x / raw_cross.z,
                z: raw_cross.y / raw_cross.z,
            }
        }

        let scene_point_a = self.to_scene_point(a);
        let scene_point_b = self.to_scene_point(b);
        let scene_point_c = self.to_scene_point(c);
        let coord_a = Vector2i { x: scene_point_a.x, y: scene_point_a.y };
        let coord_b = Vector2i { x: scene_point_b.x, y: scene_point_b.y };
        let coord_c = Vector2i { x: scene_point_c.x, y: scene_point_c.y };
        let bbox = get_triangle_bounding_box(coord_a, coord_b, coord_c);
        for i in bbox.ll.x..=bbox.ur.x {
            for j in bbox.ll.y..=bbox.ur.y {
                let coord_point = Vector2i { x:i, y:j };
                let barycentric_coord = to_barycentric_coord(coord_point, coord_a, coord_b, coord_c);
                if barycentric_coord.x < 0.0 || barycentric_coord.y < 0.0 || barycentric_coord.z < 0.0 {
                    // If any of the coordinates are negative, point is not in the triangle, so skipping it.
                    continue;
                } else {
                    let z_to_screen = scene_point_a.z * barycentric_coord.x + 
                                      scene_point_b.z * barycentric_coord.y +
                                      scene_point_c.z * barycentric_coord.z;
                    let index = (i + j * self.width as i32) as usize;
                    // Checking z-buffer, on success redrawing the pixel and updating the buffer.
                    if z_to_screen > self.z_buffer[index] {
                        self.z_buffer[index] = z_to_screen;
                        self.set_pixel(coord_point, Color::blend(color, BLACK, normal_correction_coef));
                    }
                }
            }
        }
    }

    /// Draw a triangle with vertices at specified coordiantes and filled with specified color.
    pub fn draw_triangle_textured(
        &mut self, 
        a: Vector3f, b: Vector3f, c: Vector3f, 
        uv_a: Vector2f, uv_b: Vector2f, uv_c: Vector2f,
        texture: &ImageBuffer<Rgb<u8>, Vec<u8>>,
        normal_correction_coef: f32
    ) {
        // Simple local bounding box struct for convenience.
        #[derive(Debug)]
        struct BoundingBox {
            ll: Vector2i, // lower left corner
            ur: Vector2i, // upper right corner
        }

        // Helper used to find bounding box of a triangle.
        fn get_triangle_bounding_box(coord_a: Vector2i, coord_b: Vector2i, coord_c: Vector2i) -> BoundingBox {
            return BoundingBox {
                ll: Vector2i {
                    x: min(min(coord_a.x, coord_b.x), coord_c.x),
                    y: min(min(coord_a.y, coord_b.y), coord_c.y),
                },
                ur: Vector2i {
                    x: max(max(coord_a.x, coord_b.x), coord_c.x),
                    y: max(max(coord_a.y, coord_b.y), coord_c.y),
                }
            }
        }

        fn to_barycentric_coord(coord_point: Vector2i, coord_a: Vector2i, coord_b: Vector2i, coord_c: Vector2i) -> Vector3f {
            let raw_cross = Vector3f::cross(
                Vector3f {
                    x: (coord_b.x - coord_a.x) as f32,
                    y: (coord_c.x - coord_a.x) as f32,
                    z: (coord_a.x - coord_point.x) as f32,
                },
                Vector3f {
                    x: (coord_b.y - coord_a.y) as f32,
                    y: (coord_c.y - coord_a.y) as f32,
                    z: (coord_a.y - coord_point.y) as f32,
                }
            );
            if raw_cross.z.abs() < 1.0 {
                // Degenerate triangle, returning something with negative coordinate.
                return Vector3f { x: -1.0, y: 1.0, z: 1.0};
            }
            return Vector3f {
                x: 1.0 - (raw_cross.x + raw_cross.y) / raw_cross.z,
                y: raw_cross.x / raw_cross.z,
                z: raw_cross.y / raw_cross.z,
            }
        }

        let scene_point_a = self.to_scene_point(a);
        let scene_point_b = self.to_scene_point(b);
        let scene_point_c = self.to_scene_point(c);
        let coord_a = Vector2i { x: scene_point_a.x, y: scene_point_a.y };
        let coord_b = Vector2i { x: scene_point_b.x, y: scene_point_b.y };
        let coord_c = Vector2i { x: scene_point_c.x, y: scene_point_c.y };
        let bbox = get_triangle_bounding_box(coord_a, coord_b, coord_c);
        for i in bbox.ll.x..=bbox.ur.x {
            for j in bbox.ll.y..=bbox.ur.y {
                let coord_point = Vector2i { x:i, y:j };
                let barycentric_coord = to_barycentric_coord(coord_point, coord_a, coord_b, coord_c);
                if barycentric_coord.x < 0.0 || barycentric_coord.y < 0.0 || barycentric_coord.z < 0.0 {
                    // If any of the coordinates are negative, point is not in the triangle, so skipping it.
                    continue;
                } else {
                    let z_to_screen = barycentric_coord.x * scene_point_a.z + 
                                      barycentric_coord.y * scene_point_b.z +
                                      barycentric_coord.z * scene_point_c.z;
                    let index = (i + j * self.width as i32) as usize;
                    // Checking z-buffer, on success redrawing the pixel and updating the buffer.
                    if z_to_screen > self.z_buffer[index] {
                        self.z_buffer[index] = z_to_screen;
                        let point_uv = barycentric_coord.x * uv_a + 
                                       barycentric_coord.y * uv_b + 
                                       barycentric_coord.z * uv_c; 
                        let texture_coord = Vector2i {
                            x: (point_uv.x * texture.width() as f32) as i32,
                            y: (point_uv.y * texture.height() as f32) as i32,
                        };
                        let texture_color = Color {
                            r: texture.get_pixel(texture_coord.x as u32, texture_coord.y as u32).0[0],
                            g: texture.get_pixel(texture_coord.x as u32, texture_coord.y as u32).0[1],
                            b: texture.get_pixel(texture_coord.x as u32, texture_coord.y as u32).0[2],
                        };
                        self.set_pixel(coord_point, Color::blend(texture_color, BLACK, normal_correction_coef));
                    }
                }
            }
        }
    }
}
