use super::traits::*;
use super::utils::*;
use num_traits::Float;

fn convert<T>(val: f64) -> T
where
    T: Float,
{
    T::from(val).unwrap()
}

#[derive(Debug)]
pub struct CubicSpline<T>
where
    T: Float,
{
    times: Vec<T>,
    points: Vec<Vec<T>>,
    a: Vec<Vec<T>>,
    b: Vec<Vec<T>>,
    c: Vec<Vec<T>>,
    d: Vec<Vec<T>>,
}

// from https://en.wikipedia.org/wiki/Spline_(mathematics)
impl<T> CubicSpline<T>
where
    T: Float,
{
    #[allow(clippy::many_single_char_names, clippy::just_underscores_and_digits)]
    pub fn new(times: Vec<T>, points: Vec<Vec<T>>) -> Option<Self> {
        if !is_inputs_valid(&times, &points) {
            return None;
        }
        let dim = points[0].len();
        let _0: T = T::zero();
        let _1: T = T::one();
        let _2: T = convert(2.0);
        let _3: T = convert(3.0);

        // x_i = times[i]
        // y_i = points[i][j]
        let n = points.len() - 1;

        // size of a is n+1
        let a = points.clone();
        // size of h is n
        let mut h = vec![_0; n];
        for i in 0..n {
            h[i] = times[i + 1] - times[i];
        }
        let mut alpha = vec![vec![_0; dim]; n];
        for i in 1..n {
            for j in 0..dim {
                alpha[i][j] = (_3 / h[i] * (a[i + 1][j] - a[i][j]))
                    - (_3 / h[i - 1] * (a[i][j] - a[i - 1][j]));
            }
        }
        let mut b = vec![vec![_0; dim]; n];
        let mut c = vec![vec![_0; dim]; n + 1];
        let mut d = vec![vec![_0; dim]; n];
        let mut l = vec![_1; n + 1];
        let mut m = vec![_0; n + 1];
        let mut z = vec![vec![_0; dim]; n + 1];
        for i in 1..n {
            l[i] = _2 * (times[i + 1] - times[i - 1]) - h[i - 1] * m[i - 1];
            m[i] = h[i] / l[i];
            for j in 0..dim {
                z[i][j] = (alpha[i][j] - h[i - 1] * z[i - 1][j]) / l[i];
            }
        }
        l[n] = _1;
        for j_ in 1..n + 1 {
            let j = n - j_;
            for k in 0..dim {
                c[j][k] = z[j][k] - m[j] * c[j + 1][k];
            }
            for k in 0..dim {
                b[j][k] =
                    (a[j + 1][k] - a[j][k]) / h[j] - (h[j] * (c[j + 1][k] + _2 * c[j][k]) / _3);
            }
            for k in 0..dim {
                d[j][k] = (c[j + 1][k] - c[j][k]) / (_3 * h[j]);
            }
        }
        Some(Self {
            times,
            points,
            a,
            b,
            c,
            d,
        })
    }
}

impl<T> Trajectory for CubicSpline<T>
where
    T: Float,
{
    type Point = Vec<T>;
    type Time = T;

    fn position(&self, t: T) -> Option<Vec<T>> {
        if t < self.times[0] {
            return None;
        }
        let dim = self.points[0].len();
        for i in 0..(self.points.len() - 1) {
            if t >= self.times[i] && t <= self.times[i + 1] {
                let mut pt = vec![T::zero(); dim];
                let dx = t - self.times[i];
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = self.a[i][j]
                        + self.b[i][j] * dx
                        + self.c[i][j] * dx * dx
                        + self.d[i][j] * dx * dx * dx;
                }
                return Some(pt);
            }
        }
        None
    }

    fn velocity(&self, t: T) -> Option<Vec<T>> {
        if t < self.times[0] {
            return None;
        }
        let dim = self.points[0].len();
        for i in 0..(self.points.len() - 1) {
            if t >= self.times[i] && t <= self.times[i + 1] {
                let mut pt = vec![T::zero(); dim];
                let dx = t - self.times[i];
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = self.b[i][j]
                        + convert::<T>(2.0) * self.c[i][j] * dx
                        + convert::<T>(3.0) * self.d[i][j] * dx * dx;
                }
                return Some(pt);
            }
        }
        None
    }

    fn acceleration(&self, t: T) -> Option<Vec<T>> {
        if t < self.times[0] {
            return None;
        }
        let dim = self.points[0].len();
        for i in 0..(self.points.len() - 1) {
            if t >= self.times[i] && t <= self.times[i + 1] {
                let mut pt = vec![T::zero(); dim];
                let dx = t - self.times[i];
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter =
                        convert::<T>(2.0) * self.c[i][j] + convert::<T>(6.0) * self.d[i][j] * dx;
                }
                return Some(pt);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_spline() {
        let times = vec![0.0, 1.0, 3.0, 4.0];
        let points = vec![
            vec![0.0, -1.0],
            vec![2.0, -3.0],
            vec![3.0, 3.0],
            vec![1.0, 5.0],
        ];
        let ip = CubicSpline::new(times, points).unwrap();
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.position(t).unwrap();
            println!("{t} {} {}", p[0], p[1]);
        }
    }

    #[test]
    fn test_velocity() {
        let times = vec![0.0, 1.0, 3.0, 4.0];
        let points = vec![
            vec![0.0, -1.0],
            vec![2.0, -3.0],
            vec![3.0, 3.0],
            vec![1.0, 5.0],
        ];
        let ip = CubicSpline::new(times, points).unwrap();
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.velocity(t).unwrap();
            println!("{t} {} {}", p[0], p[1]);
        }
    }

    #[test]
    fn test_acceleration() {
        let times = vec![0.0, 1.0, 3.0, 4.0];
        let points = vec![
            vec![0.0, -1.0],
            vec![2.0, -3.0],
            vec![3.0, 3.0],
            vec![1.0, 5.0],
        ];
        let ip = CubicSpline::new(times, points).unwrap();
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.acceleration(t).unwrap();
            println!("{t} {} {}", p[0], p[1]);
        }
    }
}
