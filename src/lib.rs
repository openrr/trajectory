#![warn(missing_debug_implementations, rust_2018_idioms)]

mod cubic_spline;
mod linear;
mod traits;
mod utils;

pub use self::cubic_spline::*;
pub use self::linear::*;
pub use self::traits::*;
