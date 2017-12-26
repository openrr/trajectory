extern crate num_traits;
use num_traits::Float;

pub trait Trajectory {
    type Point;
    type Time;
    fn position(&self, t: Self::Time) -> Option<Self::Point>;
    fn velocity(&self, t: Self::Time) -> Option<Self::Point>;
    fn acceleration(&self, t: Self::Time) -> Option<Self::Point>;
}

pub struct TimedPoint<K, T> {
    pub time: K,
    pub point: T,
}

impl<K, T> TimedPoint<K, T> {
    pub fn new(time: K, point: T) -> Self {
        Self { time, point }
    }
}

pub struct LinearVectorInterpolator<T>
where
    T: Float,
{
    points: Vec<TimedPoint<T, Vec<T>>>,
    dim: usize,
}

impl<T> LinearVectorInterpolator<T>
where
    T: Float,
{
    pub fn new(points: Vec<TimedPoint<T, Vec<T>>>) -> Self {
        Self {
            dim: points[0].point.len(),
            points,
        }
    }
}

impl<T> Trajectory for LinearVectorInterpolator<T>
where
    T: Float,
{
    type Point = Vec<T>;
    type Time = T;
    fn position(&self, t: T) -> Option<Vec<T>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = self.points[i].point.clone();
                let p0 = &self.points[i].point;
                let t0 = self.points[i].time;
                let p1 = &self.points[i + 1].point;
                let t1 = self.points[i + 1].time;
                for j in 0..p0.len() {
                    pt[j] = pt[j] + (p1[j] - p0[j]) / (t1 - t0) * (t - t0);
                }
                return Some(pt);
            }
        }
        None
    }
    fn velocity(&self, t: T) -> Option<Vec<T>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = vec![T::zero(); self.dim];
                let p0 = &self.points[i].point;
                let t0 = self.points[i].time;
                let p1 = &self.points[i + 1].point;
                let t1 = self.points[i + 1].time;
                for j in 0..p0.len() {
                    pt[j] = (p1[j] - p0[j]) / (t1 - t0);
                }
                return Some(pt);
            }
        }
        None
    }
    fn acceleration(&self, _t: T) -> Option<Vec<T>> {
        None
    }
}


pub struct CSplineVectorInterpolator<T>
where
    T: Float,
{
    points: Vec<TimedPoint<T, Vec<T>>>,
    dim: usize,
    a: Vec<Vec<T>>,
    b: Vec<Vec<T>>,
    c: Vec<Vec<T>>,
    d: Vec<Vec<T>>,
}

fn convert<T>(val: f64) -> T
where
    T: Float,
{
    T::from(val).unwrap()
}

// from https://en.wikipedia.org/wiki/Spline_(mathematics)
impl<T> CSplineVectorInterpolator<T>
where
    T: Float,
{
    pub fn new(points: Vec<TimedPoint<T, Vec<T>>>) -> Self {
        let _0: T = T::zero();
        let _1: T = T::one();
        let _2: T = convert(2.0);
        let _3: T = convert(3.0);

        // x_i = points[i].time
        // y_i = points[i].point[j]
        let n = points.len() - 1;
        let dim = points[0].point.len();
        // size of a is n+1
        let a = points.iter().map(|p| p.point.clone()).collect::<Vec<_>>();
        // size of h is n
        let mut h: Vec<T> = vec![_0; n];
        for i in 0..n {
            h[i] = points[i + 1].time - points[i].time;
        }
        let mut alpha = Vec::<Vec<T>>::new();
        alpha.reserve(n);
        for i in 0..n {
            // alpha[0] is never used
            let mut alpha_j = Vec::<T>::new();
            for j in 0..dim {
                let second_elm: T = if i == 0 {
                    _0
                } else {
                    _3 / h[i - 1] * (a[i][j] - a[i - 1][j])
                };
                alpha_j.push((_3 / h[i] * (a[i + 1][j] - a[i][j])) - second_elm);
            }
            alpha.push(alpha_j);
        }
        let mut c = Vec::<Vec<T>>::new();
        let mut b = Vec::<Vec<T>>::new();
        let mut d = Vec::<Vec<T>>::new();
        b.resize(n, vec![_0; dim]);
        c.resize(n + 1, vec![_0; dim]);
        d.resize(n, vec![_0; dim]);
        let mut l = vec![_1; n + 1];
        let mut m = vec![_0; n + 1];
        let mut z = Vec::<Vec<T>>::new();
        z.push(vec![_0; dim]);
        for i in 1..n {
            l[i] = _2 * (points[i + 1].time - points[i - 1].time) - h[i - 1] * m[i - 1];
            m[i] = h[i] / l[i];
            let mut z_j = Vec::<T>::new();
            z_j.reserve(dim);
            for j in 0..dim {
                z_j.push((alpha[i][j] - h[i - 1] * z[i - 1][j]) / l[i])
            }
            z.push(z_j);
        }
        l[n] = _1;
        z.push(vec![_0; dim]);
        for j_ in 1..n + 1 {
            let j = n - j_;
            for k in 0..dim {
                c[j][k] = z[j][k] - m[j] * c[j + 1][k];
            }
            for k in 0..dim {
                b[j][k] = (a[j + 1][k] - a[j][k]) / h[j] -
                    (h[j] * (c[j + 1][k] + _2 * c[j][k]) / _3);
            }
            for k in 0..dim {
                d[j][k] = (c[j + 1][k] - c[j][k]) / (_3 * h[j]);
            }
        }
        Self {
            points,
            dim,
            a,
            b,
            c,
            d,
        }
    }
}

impl<T> Trajectory for CSplineVectorInterpolator<T>
where
    T: Float,
{
    type Point = Vec<T>;
    type Time = T;
    fn position(&self, t: T) -> Option<Vec<T>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<T>::new();
                pt.resize(self.dim, T::zero());
                let dx = t - self.points[i].time;
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = self.a[i][j] + self.b[i][j] * dx + self.c[i][j] * dx * dx +
                        self.d[i][j] * dx * dx * dx;
                }
                return Some(pt);
            }
        }
        None
    }
    fn velocity(&self, t: T) -> Option<Vec<T>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<T>::new();
                pt.resize(self.dim, T::zero());
                let dx = t - self.points[i].time;
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = self.b[i][j] + convert::<T>(2.0) * self.c[i][j] * dx +
                        convert::<T>(3.0) * self.d[i][j] * dx * dx;
                }
                return Some(pt);
            }
        }
        None
    }
    fn acceleration(&self, t: T) -> Option<Vec<T>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<T>::new();
                pt.resize(self.dim, T::zero());
                let dx = t - self.points[i].time;
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = convert::<T>(2.0) * self.c[i][j] +
                        convert::<T>(6.0) * self.d[i][j] * dx;
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
    type TimedVec = TimedPoint<f64, Vec<f64>>;
    #[test]
    fn test_linear() {
        let points = vec![
            TimedVec::new(0.0, vec![0.0, -1.0]),
            TimedVec::new(1.0, vec![2.0, -3.0]),
            TimedVec::new(3.0, vec![3.0, 3.0]),
            TimedVec::new(4.0, vec![1.0, 5.0]),
        ];
        let ip = LinearVectorInterpolator::new(points);
        let p0 = ip.position(-0.5);
        assert!(p0.is_none());
        let p1 = ip.position(0.5).unwrap();
        assert_eq!(p1[0], 1.0);
        assert_eq!(p1[1], -2.0);
        let p2 = ip.position(1.0).unwrap();
        assert_eq!(p2[0], 2.0);
        assert_eq!(p2[1], -3.0);
        let p3 = ip.position(3.0).unwrap();
        assert_eq!(p3[0], 3.0);
        assert_eq!(p3[1], 3.0);
        let p4 = ip.position(3.5).unwrap();
        assert_eq!(p4[0], 2.0);
        assert_eq!(p4[1], 4.0);
    }

    #[test]
    fn test_spline() {
        let points = vec![
            TimedVec::new(0.0, vec![0.0, -1.0]),
            TimedVec::new(1.0, vec![2.0, -3.0]),
            TimedVec::new(3.0, vec![3.0, 3.0]),
            TimedVec::new(4.0, vec![1.0, 5.0]),
        ];
        let ip = CSplineVectorInterpolator::new(points);
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.position(t).unwrap();
            println!("{} {} {}", t, p[0], p[1]);
        }
    }

    #[test]
    fn test_velocity() {
        let points = vec![
            TimedVec::new(0.0, vec![0.0, -1.0]),
            TimedVec::new(1.0, vec![2.0, -3.0]),
            TimedVec::new(3.0, vec![3.0, 3.0]),
            TimedVec::new(4.0, vec![1.0, 5.0]),
        ];
        let ip = CSplineVectorInterpolator::new(points);
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.velocity(t).unwrap();
            println!("{} {} {}", t, p[0], p[1]);
        }
    }

    #[test]
    fn test_acceleration() {
        let points = vec![
            TimedVec::new(0.0, vec![0.0, -1.0]),
            TimedVec::new(1.0, vec![2.0, -3.0]),
            TimedVec::new(3.0, vec![3.0, 3.0]),
            TimedVec::new(4.0, vec![1.0, 5.0]),
        ];
        let ip = CSplineVectorInterpolator::new(points);
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.acceleration(t).unwrap();
            println!("{} {} {}", t, p[0], p[1]);
        }
    }

}
