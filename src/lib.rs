extern crate num_traits;

mod traits;
mod utils;
mod linear;
mod cubic_spline;

pub use self::traits::*;
pub use self::linear::*;
pub use self::cubic_spline::*;
