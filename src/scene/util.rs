use obj::raw::RawObj;
use image::RgbImage;
use nalgebra as na;
use na::{vector, Vector2, Vector3, Vector4};

// @TODO remove double definition.
/// Transformation of a point to homogenous coordinates.
pub fn to_hom_point(v: Vector3<f32>) -> Vector4<f32> {
    return vector![v.x, v.y, v.z, 1.0];
}

/// Transformation of a vector to homogenous coordinates.
pub fn to_hom_vector(v: Vector3<f32>) -> Vector4<f32> {
    return vector![v.x, v.y, v.z, 0.0];
}

/// Transformation of a point from homogenous coordinates.
pub fn from_hom_point(v: Vector4<f32>) -> Vector3<f32> {
    return vector![v.x / v.w, v.y / v.w, v.z / v.w];
}

/// Transformation of a vector from homogenous coordinates.
pub fn from_hom_vector(v: Vector4<f32>) -> Vector3<f32> {
    return vector![v.x, v.y, v.z];
}

/// Utility for getting convex combination of 2 Vector3<u8>'s
pub fn color_blend(color_1: Vector3<u8>, color_2: Vector3<u8>, t: f32) -> Vector3<u8>{
    return vector![
        (t * color_1.x as f32 + (1.0 - t) * color_2.x as f32) as u8,
        (t * color_1.y as f32 + (1.0 - t) * color_2.y as f32) as u8,
        (t * color_1.z as f32 + (1.0 - t) * color_2.z as f32) as u8
    ];
}

/// Struct, holding all information about the model, including geometry, texture and normal and specular maps.
pub struct Model {
    pub obj:                RawObj,
    pub texture:            RgbImage,
    pub normal_map:         RgbImage,
    pub normal_map_tangent: RgbImage,
    pub specular_map:       RgbImage,
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