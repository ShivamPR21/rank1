use num::{Complex, Zero};

// Level-1 BLAS subprogram implementations
use crate::builtins::{Rank1, ComplexOp, VecIterator};
use super::LV1;

impl<T: Rank1> LV1<Complex<T>, T> for Complex<T> {
    fn asum(n: usize, x: &[Complex<T>], incx: i64, result: &mut T) {
        let itrx = VecIterator::new(n, incx);

        for idx in itrx {    *result += x[idx].l1_norm();
        }
    }

    fn axpy(
        n: usize,
        alpha: Complex<T>,
        x: &[Complex<T>],
        incx: i64,
        y: &mut [Complex<T>],
        incy: i64,
    ) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {    y[idy] = alpha * x[idx] + y[idy];
        }
    }

    fn copy(n: usize, x: &[Complex<T>], incx: i64, y: &mut [Complex<T>], incy: i64) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            y[idy] = x[idx];
        }
    }

    fn dotc(n: usize, x: &[Complex<T>], incx: i64, y: &mut [Complex<T>], incy: i64, result: &mut T) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            *result += x[idx].conj().dot(&y[idy]);
        }
    }

    fn dotu(n: usize, x: &[Complex<T>], incx: i64, y: &mut [Complex<T>], incy: i64, result: &mut T) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            *result += x[idx].dot(&y[idy]);
        }
    }

    fn nrm2(n: usize, x: &[Complex<T>], incx: i64, result: &mut T) {
        result.set_zero();
        let itrx = VecIterator::new(n, incx);

        for idx in itrx {
            *result += x[idx].dot(&x[idx]);
        }

        *result = result.sqrt();
    }

    fn rot(n: usize, x: &mut [Complex<T>], incx: i64, y: &mut [Complex<T>], incy: i64, c: T, s: T) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            let x_tmp = x[idx] * c + y[idy] * s;
            let y_tmp = - x[idx] * s + y[idy] * c;
            (x[idx], y[idy]) = (x_tmp, y_tmp);
        }
    }
}

impl<T: Rank1> LV1<Complex<T>, Complex<T>> for Complex<T> {
    fn dot(n: usize, x: &[Complex<T>], incx: i64, y: &mut [Complex<T>], incy: i64, result: &mut Complex<T>) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            *result += x[idx] * y[idy];
        }
    }

    fn nrm2(n: usize, x: &[Complex<T>], incx: i64, result: &mut Complex<T>) {
        result.set_zero();
        let itrx = VecIterator::new(n, incx);

        for idx in itrx {
            *result += x[idx] * x[idx];
        }

        *result = result.sqrt();
    }
}
