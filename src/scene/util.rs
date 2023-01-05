use nalgebra as na;
use na::{vector, Vector3, Vector4};

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