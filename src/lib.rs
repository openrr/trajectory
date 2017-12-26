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
    pub fn new(times: Vec<T>, points: Vec<Vec<T>>) -> Option<Self> {
        let timed_points = from_times_and_points(times, points)?;
        Self::from_points(timed_points)
    }
    pub fn from_points(points: Vec<TimedPoint<T, Vec<T>>>) -> Option<Self> {
        if !is_points_valid(&points) {
            return None;
        }
        Some(Self {
            dim: points[0].point.len(),
            points,
        })
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

fn is_points_valid<T>(points: &Vec<TimedPoint<T, Vec<T>>>) -> bool {
    if points.is_empty() {
        return false;
    }
    if points[0].point.is_empty() {
        return false;
    }
    let dim = points[0].point.len();
    for p in points {
        if p.point.len() != dim {
            return false;
        }
    }
    true
}

fn from_times_and_points<T>(
    times: Vec<T>,
    points: Vec<Vec<T>>,
) -> Option<Vec<TimedPoint<T, Vec<T>>>>
where
    T: Float,
{
    if times.len() != points.len() {
        return None;
    }
    let mut timed_points: Vec<TimedPoint<T, Vec<T>>> = Vec::new();
    for (i, point) in points.into_iter().enumerate() {
        timed_points.push(TimedPoint::new(times[i], point));
    }
    Some(timed_points)
}

// from https://en.wikipedia.org/wiki/Spline_(mathematics)
impl<T> CSplineVectorInterpolator<T>
where
    T: Float,
{
    pub fn new(times: Vec<T>, points: Vec<Vec<T>>) -> Option<Self> {
        let timed_points = from_times_and_points(times, points)?;
        Self::from_points(timed_points)
    }

    pub fn from_points(points: Vec<TimedPoint<T, Vec<T>>>) -> Option<Self> {
        if !is_points_valid(&points) {
            return None;
        }
        let dim = points[0].point.len();
        let _0: T = T::zero();
        let _1: T = T::one();
        let _2: T = convert(2.0);
        let _3: T = convert(3.0);

        // x_i = points[i].time
        // y_i = points[i].point[j]
        let n = points.len() - 1;

        // size of a is n+1
        let a = points.iter().map(|p| p.point.clone()).collect::<Vec<_>>();
        // size of h is n
        let mut h = vec![_0; n];
        for i in 0..n {
            h[i] = points[i + 1].time - points[i].time;
        }
        let mut alpha = vec![vec![_0; dim]; n];
        for i in 1..n {
            for j in 0..dim {
                alpha[i][j] = (_3 / h[i] * (a[i + 1][j] - a[i][j])) -
                    (_3 / h[i - 1] * (a[i][j] - a[i - 1][j]));
            }
        }
        let mut b = vec![vec![_0; dim]; n];
        let mut c = vec![vec![_0; dim]; n + 1];
        let mut d = vec![vec![_0; dim]; n];
        let mut l = vec![_1; n + 1];
        let mut m = vec![_0; n + 1];
        let mut z = vec![vec![_0; dim]; n + 1];
        for i in 1..n {
            l[i] = _2 * (points[i + 1].time - points[i - 1].time) - h[i - 1] * m[i - 1];
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
                b[j][k] = (a[j + 1][k] - a[j][k]) / h[j] -
                    (h[j] * (c[j + 1][k] + _2 * c[j][k]) / _3);
            }
            for k in 0..dim {
                d[j][k] = (c[j + 1][k] - c[j][k]) / (_3 * h[j]);
            }
        }
        Some(Self {
            points,
            dim,
            a,
            b,
            c,
            d,
        })
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
                let mut pt = vec![T::zero(); self.dim];
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
                let mut pt = vec![T::zero(); self.dim];
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
                let mut pt = vec![T::zero(); self.dim];
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
        let ip = LinearVectorInterpolator::from_points(points).unwrap();
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
        let ip = CSplineVectorInterpolator::from_points(points).unwrap();
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
        let ip = CSplineVectorInterpolator::from_points(points).unwrap();
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.velocity(t).unwrap();
            println!("{} {} {}", t, p[0], p[1]);
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
        let ip = CSplineVectorInterpolator::new(times, points).unwrap();
        for i in 0..400 {
            let t = i as f64 * 0.01f64;
            let p = ip.acceleration(t).unwrap();
            println!("{} {} {}", t, p[0], p[1]);
        }
    }

}
