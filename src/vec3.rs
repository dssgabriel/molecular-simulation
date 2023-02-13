use std::fmt::Display;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub fn splat(val: f64) -> Self {
        Self::new(val, val, val)
    }

    #[inline]
    pub fn unit_x() -> Self {
        Self {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }

    #[inline]
    pub fn unit_y() -> Self {
        Self {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    }

    #[inline]
    pub fn unit_z() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    }

    #[inline]
    pub fn dot(&self, rhs: Self) -> f64 {
        (self.x * rhs.x) + (self.y * rhs.y) + (self.z * rhs.z)
    }

    /// The wedge (aka exterior) product of two vectors.
    ///
    /// This operation results in a bivector, which represents
    /// the plane parallel to the two vectors, and which has a
    /// 'oriented area' equal to the parallelogram created by extending
    /// the two vectors, oriented such that the positive direction is the
    /// one which would move `self` closer to `rhs`.
    #[inline]
    pub fn wedge(&self, rhs: Self) -> Self {
        Self::new(
            (self.x * rhs.y) - (self.y * rhs.x),
            (self.x * rhs.z) - (self.z * rhs.x),
            (self.y * rhs.z) - (self.z * rhs.y),
        )
    }

    #[inline]
    pub fn cross(&self, rhs: Self) -> Self {
        Self::new(
            (self.y * rhs.z) + (-self.z * rhs.y),
            (self.z * rhs.x) + (-self.x * rhs.z),
            (self.x * rhs.y) + (-self.y * rhs.x),
        )
    }

    #[inline]
    pub fn distance_square(&self, rhs: Self) -> f64 {
        (rhs.x - self.x) * (rhs.x - self.x)
            + (rhs.y - self.y) * (rhs.y - self.y)
            + (rhs.z - self.z) * (rhs.z - self.z)
    }

    #[inline]
    pub fn distance(&self, rhs: Self) -> f64 {
        self.distance_square(rhs).sqrt()
    }

    #[inline]
    pub fn mag_sq(&self) -> f64 {
        (self.x * self.x) + (self.y * self.y) + (self.z * self.z)
    }

    #[inline]
    pub fn mag(&self) -> f64 {
        self.mag_sq().sqrt()
    }

    #[inline]
    pub fn normalize(&mut self) {
        let r_mag = 1.0 / self.mag();
        self.x *= r_mag;
        self.y *= r_mag;
        self.z *= r_mag;
    }

    #[inline]
    pub fn normalized(&self) -> Self {
        let mut r = self.clone();
        r.normalize();
        r
    }

    #[inline]
    pub fn add(&self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }

    #[inline]
    pub fn sub(&self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }

    #[inline]
    pub fn mul(&self, rhs: Self) -> Self {
        Self::new(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
    }

    #[inline]
    pub fn mul_add(&self, mul: Self, add: Self) -> Self {
        Self::new(
            self.x.mul_add(mul.x, add.x),
            self.y.mul_add(mul.y, add.y),
            self.z.mul_add(mul.z, add.z),
        )
    }

    #[inline]
    pub fn abs(&self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    #[inline]
    pub fn clamp(&mut self, min: Self, max: Self) {
        self.x = self.x.max(min.x).min(max.x);
        self.y = self.y.max(min.y).min(max.y);
        self.z = self.z.max(min.z).min(max.z);
    }

    #[inline]
    pub fn clamped(mut self, min: Self, max: Self) -> Self {
        self.clamp(min, max);
        self
    }

    #[inline]
    pub fn map<F>(&self, mut f: F) -> Self
    where
        F: FnMut(f64) -> f64,
    {
        Self::new(f(self.x), f(self.y), f(self.z))
    }

    #[inline]
    pub fn apply<F>(&mut self, mut f: F)
    where
        F: FnMut(f64) -> f64,
    {
        self.x = f(self.x);
        self.y = f(self.y);
        self.z = f(self.z);
    }

    #[inline]
    pub fn max_by_component(mut self, rhs: Self) -> Self {
        self.x = self.x.max(rhs.x);
        self.y = self.y.max(rhs.y);
        self.z = self.z.max(rhs.z);
        self
    }

    #[inline]
    pub fn min_by_component(mut self, rhs: Self) -> Self {
        self.x = self.x.min(rhs.x);
        self.y = self.y.min(rhs.y);
        self.z = self.z.min(rhs.z);
        self
    }

    #[inline]
    pub fn component_max(&self) -> f64 {
        self.x.max(self.y).max(self.z)
    }

    #[inline]
    pub fn component_min(&self) -> f64 {
        self.x.min(self.y).min(self.z)
    }

    #[inline]
    pub fn zero() -> Self {
        Self::splat(0.0)
    }

    #[inline]
    pub fn one() -> Self {
        Self::splat(1.0)
    }
}

impl Add for Vec3 {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl Mul for Vec3 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl MulAssign for Vec3 {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        self.x *= rhs.x;
        self.y *= rhs.y;
        self.z *= rhs.z;
    }
}

impl Neg for Vec3 {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl Display for Vec3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{ x: {:.5e}, y: {:.5e}, z: {:.5e} }}",
            self.x, self.y, self.z
        )
    }
}
