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

/// 2D coordinate of an image pixel. 
#[derive(Clone, Copy)]
pub struct Coord {
    pub x: u32,
    pub y: u32,
}

/// Image, holding its width, height and private flat array(vec) of pixel data.
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
        if coord.x > self.width {
            return false;
        }
        if coord.y > self.height {
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
    pub fn set_pixel(&mut self, coord: Coord, color: Color) -> Result<(), &'static str>{
        if !self.coord_in_bounds(coord) {
            return Err("Given coordinate is not in image bounds");
        }
        // Pixel data is rgba, so we find the starting index of a 3-tuple and do 3 assignments.
        // @OPTI maybe I can set 4 values simultaneously somehow? 
        let index = (3 * (coord.x + coord.y * self.width)) as usize;
        self.pixel_data[index + 0] = color.r;
        self.pixel_data[index + 1] = color.g;
        self.pixel_data[index + 2] = color.b;
        return Ok(());
    }

    /// Draws a line between start and end coordinates with specified color.
    pub fn draw_line(&mut self, start: Coord, end: Coord, color: Color) -> Result<(), &'static str>{
        if !self.coord_in_bounds(start) {
            return Err("Start coordinate is not in image bounds.");
        }
        if !self.coord_in_bounds(end) {
            return Err("End coordinate is not in image bounds.");
        }

        return Ok(());
    }
}
