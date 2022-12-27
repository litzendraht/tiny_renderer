/// Vector2 storing f32
#[derive(Debug, Clone, Copy)]
pub struct Vector2f {
    pub x: f32,
    pub y: f32,
}

/// Vector2 storing 2 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector2i {
    pub x: i32,
    pub y: i32,
}

/// Vector3 storing 3 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector3f {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector3f {
    /// Dot product of 2 Vector3's.
    pub fn dot(a: Vector3f, b: Vector3f) -> f32 {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /// Norm of a Vector3.
    pub fn norm(&self) -> f32 {
        return Vector3f::dot(*self, *self).sqrt();
    }

    /// Cross product of 2 Vector3's.
    pub fn cross(a: Vector3f, b: Vector3f) -> Vector3f {
        return Vector3f { 
            x: a.y * b.z - a.z * b.y, 
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x,
        };
    }
}

impl std::ops::Add<Vector3f> for Vector3f {
    type Output = Vector3f;

    fn add(self, _rhs: Vector3f) -> Vector3f {
        return Vector3f {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
            z: self.z + _rhs.z,
        };
    }
}

impl std::ops::Sub<Vector3f> for Vector3f {
    type Output = Vector3f;

    fn sub(self, _rhs: Vector3f) -> Vector3f {
        return Vector3f {
            x: self.x - _rhs.x,
            y: self.y - _rhs.y,
            z: self.z - _rhs.z,
        };
    }
}

impl std::ops::Mul<Vector3f> for f32 {
    type Output = Vector3f;

    fn mul(self, _rhs: Vector3f) -> Vector3f {
        return Vector3f {
            x: self * _rhs.x,
            y: self * _rhs.y,
            z: self * _rhs.z,
        };
    }
}