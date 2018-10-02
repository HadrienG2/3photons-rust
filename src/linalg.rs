//! This module is a helper around nalgebra which smooths out some rough edges
//! of that library (e.g. lack of vectors of size 8) and implements some domain-
//! specific 4-momentum handling logic.


/// Re-export of nalgebra's type-level integers
pub mod dimension {
    pub use nalgebra::dimension::*;
}

/// Re-export of nalgebra's matrix and vector types + some extra definitions
pub mod vecmat {
    use nalgebra::{
        dimension::*,
        MatrixMN,
        VectorN,
    };

    // Re-export of various matrix and vector types from nalgebra
    pub use nalgebra::{
        Vector2,
        Vector3,
        Vector5,
        Matrix2x3,
        Matrix2x4,
        Matrix2x5,
        Matrix3,
        Matrix3x2,
        Matrix3x4,
        Matrix4x3,
        Matrix5,
        Matrix5x4,
        MatrixSlice,
    };

    /// An 8-dimensional vector type
    pub type Vector8<T> = VectorN<T, U8>;

    /// A 5x8 matrix type
    pub type Matrix5x8<T> = MatrixMN<T, U5, U8>;

    /// Convenience shorthand for defining vector slices
    pub type VectorSlice<'a, T, SliceDim, ParentDim> =
        MatrixSlice<'a, T, SliceDim, U1, U1, ParentDim>;
}

/// Handling of relativistic 4-momenta
pub mod momentum {
    use ::numeric::Real;

    use nalgebra::{
        dimension::*,
        Vector4
    };

    /// Relativistic 4-momentum
    pub type Momentum = Vector4<Real>;

    /// Convenience const for accessing the X coordinate of a 4-vector
    pub const X: usize = 0;

    /// Convenience const for accessing the Y coordinate of a 4-vector
    pub const Y: usize = 1;

    /// Convenience const for accessing the Z coordinate of a 4-vector
    pub const Z: usize = 2;

    /// Convenience const for accessing the E coordinate of a 4-vector
    pub const E: usize = 3;

    /// Get a read-only view on the spatial part of a 4-momentum
    pub fn xyz(m: &Momentum) -> super::vecmat::VectorSlice<Real, U3, U4> {
        m.fixed_rows::<U3>(X)
    }
}