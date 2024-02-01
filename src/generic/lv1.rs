pub trait LV1<T, R> {
    fn asum(n: usize, x: &[T], incx: i64, result: &mut R) {
        unimplemented!();
    }

    fn axpy(n: usize, alpha: T, x: &[T], incx: i64, y: &mut [T], incy: i64) {
        unimplemented!();
    }

    fn copy(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64) {
        unimplemented!();
    }

    fn dot(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut R) {
        unimplemented!();
    }

    fn sdsdot(n: usize, sb: T, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut R){
        unimplemented!();
    }

    fn dotc(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut R){
        unimplemented!();
    }

    fn dotu(n: usize, x: &[T], incx: i64, y: &mut [T], incy: i64, result: &mut R){
        unimplemented!();
    }

    fn nrm2(n: usize, x: &[T], incx: i64, result: &mut R){
        unimplemented!();
    }

    fn rot(n: usize, x: &mut [T], incx: i64, y: &mut [T], incy: i64, c: R, s: R){
        unimplemented!();
    }

    fn rotg(a: &mut T, b: &mut T, c: &mut T, s: &mut T) {
        unimplemented!();
    }

    fn rotm(n: usize, x: &mut [T], incx: i64, y: &mut [T], incy: i64, param: &[T;5]) {
        unimplemented!();
    }

    fn rotmg(d1: &mut T, d2: &mut T, x1: &mut T, y1: &T, param: &[T;5]) {
        unimplemented!();
    }
}

mod lv1c;
mod lv1r;

pub use lv1c::*;
pub use lv1r::*;
