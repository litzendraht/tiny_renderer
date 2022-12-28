#[macro_export]
macro_rules! matrix4f {
    ( $( $x:expr ),* ) => {
        {
            let mut matrix = Vec::new();
            $(
                temp_vec.push($x);
            )*
            temp_vec
        }
    };
}

/// Vector2 storing 2 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector2i {
    pub x: i32,
    pub y: i32,
}

/// Vector2 storing f32
#[derive(Debug, Clone, Copy)]
pub struct Vector2f {
    pub x: f32,
    pub y: f32,
}

/// Vector3 storing 3 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector3f {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

/// Vector4 storing 4 f32
#[derive(Debug, Clone, Copy)]
pub struct Vector4f {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

/// Matrix4x4 storing 16 f32
#[derive(Clone, Copy)]
pub struct Matrix4f {
    // @TRASH this is suffering.
    pub xx: f32,
    pub xy: f32,
    pub xz: f32,
    pub xw: f32,
    pub yx: f32,
    pub yy: f32,
    pub yz: f32,
    pub yw: f32,
    pub zx: f32,
    pub zy: f32,
    pub zz: f32,
    pub zw: f32,
    pub wx: f32,
    pub wy: f32,
    pub wz: f32,
    pub ww: f32,
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

impl std::ops::Add<Vector2f> for Vector2f {
    type Output = Vector2f;

    fn add(self, _rhs: Vector2f) -> Vector2f {
        return Vector2f {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
        };
    }
}

impl std::ops::Mul<Vector2f> for f32 {
    type Output = Vector2f;

    fn mul(self, _rhs: Vector2f) -> Vector2f {
        return Vector2f {
            x: self * _rhs.x,
            y: self * _rhs.y,
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
