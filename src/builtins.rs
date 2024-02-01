// Builtin types

use num::{Signed, traits::NumAssignOps, Float, NumCast, ToPrimitive};
pub use num::{Complex, Num};

pub trait Rank1
where
    Self: Num + NumAssignOps + Copy + Clone + Sized + Signed + Float + PartialOrd + NumCast + ToPrimitive,
{
}

impl Rank1 for f32 {}
impl Rank1 for f64 {}

pub trait ComplexOp<T: Rank1> {
    fn dot(&self, other: &Complex<T>) -> T;
}

impl<T: Rank1> ComplexOp<T> for Complex<T>{
    fn dot(&self, other: &Complex<T>) -> T {
        self.re * other.re - self.im * other.im
    }
}

/// Storage type enum
/// Storage::R -> Row Major storage
/// Storage::C -> Row Major storage
#[derive(Clone, Copy)]
pub enum Storage {
    R, // Row Major storage
    C, // Column Major storage
}

/// Transpose enum
/// Transpose::N -> No Transpose
/// Transpose::T -> Transpose the matrix
/// Transpose::C -> Take Conjugate transpose
#[derive(Clone, Copy)]
pub enum Transpose {
    N, // No transpose
    T, // Transpose
    C, // Conjugate Transpose
}

/// Uplo enum
/// Uplo::U -> Access only Upper Triangular part
/// Uplo::L -> Access only Lower Triangular part
#[derive(Clone, Copy)]
pub enum Uplo {
    U, // Access the Upper Triangular part
    L, // Access the Lower Triangular part
}

/// `Matrix` struct is the most basic matrix storage
pub struct Matrix<T: Rank1> {
    pub data: Vec<T>,                  // Matrix Data Storage
    pub n_rows: usize,                 // Number of rows
    pub n_cols: usize,                 // Number of cols
    pub storage_type: Option<Storage>, // Storage type
}

/// `BandMatrix` struct for band matrix storage
/// supplied with additional metadata along with
/// core matrix as `Matrix` struct
pub struct BandMatrix<T: Rank1> {
    pub matrix: Matrix<T>, // Core matrix data
    pub ku: usize,         // # Super-diagonals
    pub kl: usize,         // # Sub-diagonals
    pub uplo: Uplo,        // Whether UTM or LTM for triangular band matrix
}

/// `Vector` struct for simple vector data
pub struct Vector<T: Rank1> {
    pub data: Vec<T>, // Vector data
    pub n: usize,     // # elements
    pub incx: i32,    // Increment
}

/// `FromRowVector` defines methods for different data types
/// that are implemented to load dense matrix inputs from a
/// row-major layout to raw data in target DataType.
pub trait FromRowVector<T: Rank1> {
    /// Fills data field in `T` in column major layout,
    /// from values in `data` which is dense matrix in row-major layout.
    fn as_column_major(&mut self, data: &[T]);

    /// Fills data field in `T` in row major layout,
    /// from values in `data` which is dense matrix in row-major layout.
    fn as_row_major(&mut self, data: &[T]);

    /// Fills data field in `T` in default or current layout,
    /// from values in `data` which is dense matrix in row-major layout.
    fn with_data(&mut self, data: &[T]);
}

pub struct VecIterator {
    n: usize,
    inc: usize,
    rev: bool,

    // Iterator State
    curr: usize,
    next: usize,
    data_size: usize
}

impl VecIterator{
    pub fn new(n: usize, inc: i64) -> Self {
        let mut itr = VecIterator { n, inc: inc.wrapping_abs() as usize, rev: (inc < 0), curr: 0, next: 0, data_size: 0 };
        itr.data_size = (itr.n - 1) * itr.inc;
        itr
    }
}

impl Iterator for VecIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let mut val: Option<usize> = None;

        if self.curr < self.n {
            if self.rev {
                val = Some(self.data_size - self.curr * self.inc - 1);
            } else {
                val = Some(self.curr * self.inc);
            }
        }

        self.curr = self.next;
        self.next += 1;

        val
    }
}
