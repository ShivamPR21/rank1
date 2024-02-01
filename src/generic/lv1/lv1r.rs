// Level-1 BLAS subprogram implementations
use crate::builtins::{Rank1, VecIterator};
use super::LV1;

impl<T: Rank1> LV1<T, T> for T {
    fn asum(n: usize, x: &[T], incx: i64, result: &mut T) {
        result.set_zero();
        let itrx = VecIterator::new(n, incx);

        for idx in itrx {
            *result += x[idx].abs();
        }
    }

    fn axpy(n: usize, alpha: T, x: &[T], incx: i64, y: &mut [T], incy: i64) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            y[idy] = alpha * x[idx] + y[idy];
        }
    }

    fn copy(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            y[idy] = x[idx];
        }
    }

    fn dot(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut T) {
        result.set_zero();
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            *result += x[idx] * y[idy];
        }
    }

    fn sdsdot(n: usize, sb: T, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut T) {
        result.set_zero();

        *result = sb;
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            *result += x[idx] * y[idy];
        }
    }

    fn nrm2(n: usize, x: &[T], incx: i64, result: &mut T) {
        result.set_zero();
        let itrx = VecIterator::new(n, incx);

        for idx in itrx {
            *result += x[idx] * x[idx];
        }

        *result = result.sqrt();
    }

    fn rot(n: usize, x: &mut [T], incx: i64, y: &mut [T], incy: i64, c: T, s: T) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        for (idx, idy) in itrx.zip(itry) {
            let x_tmp = c * x[idx] + s * y[idy];
            let y_tmp = -s * x[idx] + c * y[idy];
            (x[idx], y[idy]) = (x_tmp, y_tmp);
        }
    }

    fn rotg(a: &mut T, b: &mut T, c: &mut T, s: &mut T) {
        let r = (*a * *a + *b * *b).sqrt();
        *c = *a/r;
        *s = *b/r;

        if a.abs() > b.abs() {
            *b = *s;
        } else if *c != T::zero() {
            *b = T::one()/(*c);
        } else {
            *b = T::one();
        }

        *a = r;
    }

    fn rotm(n: usize, x: &mut [T], incx: i64, y: &mut [T], incy: i64, param: &[T;5]) {
        let itrx = VecIterator::new(n, incx);
        let itry = VecIterator::new(n, incy);

        let flag = param[0].to_f32().unwrap();
        let _fn = |_x: &T, _y: &T| {
            if flag == -1.0 {
                (
                    param[1] * *_x + param[3] * *_y,
                    param[2] * *_x + param[4] * *_y
                )
            } else if flag == 0.0 {
                (
                    *_x + param[3] * *_y,
                    param[2] * *_x + *_y
                )
            } else if flag == 1.0 {
                (
                    param[1] * *_x + *_y,
                    - *_x + param[4] * *_y
                )
            } else if flag == -2.0 {
                (
                    *_x,
                    *_y
                )
            } else {
                panic!();
            }
        };

        for (idx, idy) in itrx.zip(itry) {
            (x[idx], y[idy]) = _fn(&x[idx], &y[idy]);
        }
    }

    fn rotmg(d1: &mut T, d2: &mut T, x1: &mut T, y1: &T, param: &[T;5]) {

    }
}
