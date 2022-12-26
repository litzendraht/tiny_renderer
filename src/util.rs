use std::ops;

/// 2D coordinate of an image pixel.
/// i32 to allow coordiantes outside of the image and to remove a lot of casts.
#[derive(Debug, Clone, Copy)]
pub struct Coord {
    pub x: i32,
    pub y: i32,
}

/// Vector2 storing 2 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector2 {
    pub x: f32,
    pub y: f32,
}

/// Vector3 storing 3 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector3 {
    /// Dot product of 2 Vector3's.
    pub fn dot(a: Vector3, b: Vector3) -> f32 {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /// Norm of a Vector3.
    pub fn norm(&self) -> f32 {
        return Vector3::dot(*self, *self).sqrt();
    }

    /// Cross product of 2 Vector3's.
    pub fn cross(a: Vector3, b: Vector3) -> Vector3 {
        return Vector3 { 
            x: a.y * b.z - a.z * b.y, 
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x,
        };
    }
}

impl std::ops::Add<Vector3> for Vector3 {
    type Output = Vector3;

    fn add(self, _rhs: Vector3) -> Vector3 {
        return Vector3 {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
            z: self.z + _rhs.z,
        };
    }
}

impl std::ops::Sub<Vector3> for Vector3 {
    type Output = Vector3;

    fn sub(self, _rhs: Vector3) -> Vector3 {
        return Vector3 {
            x: self.x - _rhs.x,
            y: self.y - _rhs.y,
            z: self.z - _rhs.z,
        };
    }
}