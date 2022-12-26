use std::ops::Bound;
use std::cmp::{min, max};

use crate::util::{Coord, Vector2, Vector3};

/// Struct, representing raw rgb8 pixel data.
#[derive(Clone, Copy)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

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

/// Image, holding its width, height and private flat array(vec) of pixel data.
/// /// (0, 0) is the bottom left coord.
pub struct Image {
    pub width: u32,
    pub height: u32,
    pixel_data: Vec<u8>,  // Storing flat array.
}

impl Image {
    /// Generates new image struct with specified width and height.
    /// Pixel data format is assumed to be rgb8.
    pub fn new(width: u32, height: u32) -> Image {
        let capacity = (3 * width * height) as usize;
        let pixel_data: Vec<u8> = vec![0; capacity];
        return Image {
            width,
            height,
            pixel_data,
        }
    }

    /// Checking if coord is in image bounds.
    fn coord_in_bounds(&self, coord: Coord) -> bool {
        if coord.x > self.width as i32 {
            return false;
        }
        if coord.y > self.height as i32 {
            return false;
        }
        return true;
    }

    
    pub fn as_pixel_data(&self) -> &[u8]{
        return &self.pixel_data[..];
    }

    /// Sets all pixels data to (0, 0, 0).
    pub fn flush(&mut self) {
        for elem in &mut self.pixel_data {
            *elem = 0;
        }
    }

    /// Sets image pixel to a color at specifed coordinate.
    /// 
    /// Assumes, that pixel data is rgb8.
    pub fn set_pixel(&mut self, coord: Coord, color: Color) {
        // Pixel data is rgb8, so we find the starting index of a 3-tuple and do 3 assignments.
        // @OPTI maybe I can set 3 values simultaneously somehow?
        // Forcing (0, 0) to be in the bottom left here by inverting coord.y.
        let index = (3 * (coord.x + ((self.height as i32 - 1) - coord.y) * self.width as i32)) as usize;
        self.pixel_data[index + 0] = color.r;
        self.pixel_data[index + 1] = color.g;
        self.pixel_data[index + 2] = color.b;
    }

    /// Draws a line between v_1 and v_2 coordinates with specified color
    /// via Bresenham's algorithm as presented in https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    pub fn draw_line(&mut self, coord_a: Coord, coord_b: Coord, color: Color) {
        let mut x_0 = coord_a.x;
        let x_1 = coord_b.x;
        let mut y_0 = coord_a.y;
        let y_1 = coord_b.y;
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
            self.set_pixel(Coord { x: x_0, y: y_0}, color);
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

    /// Draw a triangle with vertices at specified cooediantes and filled with specified color.
    pub fn draw_triangle(&mut self, coord_a: Coord, coord_b: Coord, coord_c: Coord, color: Color) {
        // Simple local bounding box struct for convenience.
        #[derive(Debug)]
        struct BoundingBox {
            ll: Coord, // lower left corner
            ur: Coord, // upper right corner
}

        // Helper used to find bounding box of a triangle.
        fn get_triangle_bounding_box(coord_a: Coord, coord_b: Coord, coord_c: Coord) -> BoundingBox {
            return BoundingBox {
                ll: Coord {
                    x: min(min(coord_a.x, coord_b.x), coord_c.x),
                    y: min(min(coord_a.y, coord_b.y), coord_c.y),
                },
                ur: Coord {
                    x: max(max(coord_a.x, coord_b.x), coord_c.x),
                    y: max(max(coord_a.y, coord_b.y), coord_c.y),
                }
            }
        }

        // Function to check if position of a point is inside a triangle.
        // Calculates barycentric coordinates and checks that they are all positive.
        fn point_is_inside_triangle(coord_point: Coord, coord_a: Coord, coord_b: Coord, coord_c: Coord) -> bool {
            fn to_barycentric_coord(coord_point: Coord, coord_a: Coord, coord_b: Coord, coord_c: Coord) -> Vector3 {
                let raw_cross = Vector3::cross(
                    Vector3 {
                        x: (coord_b.x - coord_a.x) as f32,
                        y: (coord_c.x - coord_a.x) as f32,
                        z: (coord_a.x - coord_point.x) as f32,
                    },
                    Vector3 {
                        x: (coord_b.y - coord_a.y) as f32,
                        y: (coord_c.y - coord_a.y) as f32,
                        z: (coord_a.y - coord_point.y) as f32,
                    }
                );
                if raw_cross.z.abs() < 1.0 {
                    // Degenerate triangle.
                    return Vector3 { x: -1.0, y: 1.0, z: 1.0};
                }
                return Vector3 {
                    x: 1.0 - (raw_cross.x + raw_cross.y) / raw_cross.z,
                    y: raw_cross.x / raw_cross.z,
                    z: raw_cross.y / raw_cross.z,
                }
            }

            let coord_barycentric = to_barycentric_coord(coord_point, coord_a, coord_b, coord_c);
            if coord_barycentric.x < 0.0 || coord_barycentric.y < 0.0 || coord_barycentric.z < 0.0
            {
                return false;
            } else {
                return true;
            }
        }

        let bbox = get_triangle_bounding_box(coord_a, coord_b, coord_c);
        for i in bbox.ll.x..=bbox.ur.x {
            for j in bbox.ll.y..=bbox.ur.y {
                let coord_point = Coord { x:i, y:j };
                if point_is_inside_triangle(coord_point, coord_a, coord_b, coord_c) {
                    self.set_pixel(Coord { x:i, y:j }, color);
                }
            }
        }
    }
}
