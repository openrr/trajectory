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

type TimedVec = TimedPoint<f64, Vec<f64>>;

pub struct LinearVectorInterpolator {
    points: Vec<TimedVec>,
    dim: usize,
}

impl LinearVectorInterpolator {
    pub fn new(points: Vec<TimedVec>) -> Self {
        Self {
            dim: points[0].point.len(),
            points,
        }
    }
}

impl Trajectory for LinearVectorInterpolator {
    type Point = Vec<f64>;
    type Time = f64;
    fn position(&self, t: f64) -> Option<Vec<f64>> {
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
                    pt[j] += (p1[j] - p0[j]) / (t1 - t0) * (t - t0);
                }
                return Some(pt);
            }
        }
        None
    }
    fn velocity(&self, t: f64) -> Option<Vec<f64>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = vec![0.0; self.dim];

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
    fn acceleration(&self, _t: f64) -> Option<Vec<f64>> {
        None
    }
}


pub struct CSplineVectorInterpolator {
    points: Vec<TimedVec>,
    dim: usize,
    a: Vec<Vec<f64>>,
    b: Vec<Vec<f64>>,
    c: Vec<Vec<f64>>,
    d: Vec<Vec<f64>>,
}

// from https://en.wikipedia.org/wiki/Spline_(mathematics)
impl CSplineVectorInterpolator {
    pub fn new(points: Vec<TimedVec>) -> Self {
        // x_i = points[i].time
        // y_i = points[i].point[j]
        let n = points.len() - 1;
        let dim = points[0].point.len();
        // size of a is n+1
        let a = points.iter().map(|p| p.point.clone()).collect::<Vec<_>>();
        // size of h is n
        let mut h: Vec<f64> = vec![0.0; n];
        for i in 0..n {
            h[i] = points[i + 1].time - points[i].time;
        }
        let mut alpha = Vec::<Vec<f64>>::new();
        alpha.reserve(n);
        for i in 0..n {
            // alpha[0] is never used
            let mut alpha_j = Vec::<f64>::new();
            for j in 0..dim {
                let second_elm = if i == 0 {
                    0.0
                } else {
                    3.0 / h[i - 1] * (a[i][j] - a[i - 1][j])
                };
                alpha_j.push((3.0 / h[i] * (a[i + 1][j] - a[i][j])) - second_elm);
            }
            alpha.push(alpha_j);
        }
        let mut c = Vec::<Vec<f64>>::new();
        let mut b = Vec::<Vec<f64>>::new();
        let mut d = Vec::<Vec<f64>>::new();
        b.resize(n, vec![0.0; dim]);
        c.resize(n + 1, vec![0.0; dim]);
        d.resize(n, vec![0.0; dim]);
        let mut l = vec![1.0; n + 1];
        let mut m = vec![0.0; n + 1];
        let mut z = Vec::<Vec<f64>>::new();
        z.push(vec![0.0; dim]);
        for i in 1..n {
            l[i] = 2.0 * (points[i + 1].time - points[i - 1].time) - h[i - 1] * m[i - 1];
            m[i] = h[i] / l[i];
            let mut z_j = Vec::<f64>::new();
            z_j.reserve(dim);
            for j in 0..dim {
                z_j.push((alpha[i][j] - h[i - 1] * z[i - 1][j]) / l[i])
            }
            z.push(z_j);
        }
        l[n] = 1.0;
        z.push(vec![0.0; dim]);
        for j_ in 1..n + 1 {
            let j = n - j_;
            for k in 0..dim {
                c[j][k] = z[j][k] - m[j] * c[j + 1][k];
            }
            for k in 0..dim {
                b[j][k] = (a[j + 1][k] - a[j][k]) / h[j] -
                    (h[j] * (c[j + 1][k] + 2.0 * c[j][k]) / 3.0);
            }
            for k in 0..dim {
                d[j][k] = (c[j + 1][k] - c[j][k]) / (3.0 * h[j]);
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

impl Trajectory for CSplineVectorInterpolator {
    type Point = Vec<f64>;
    type Time = f64;
    fn position(&self, t: f64) -> Option<Vec<f64>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<f64>::new();
                pt.resize(self.dim, 0.0);
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
    fn velocity(&self, t: f64) -> Option<Vec<f64>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<f64>::new();
                pt.resize(self.dim, 0.0);
                let dx = t - self.points[i].time;
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = self.b[i][j] + 2.0 * self.c[i][j] * dx +
                        3.0 * self.d[i][j] * dx * dx;
                }
                return Some(pt);
            }
        }
        None
    }
    fn acceleration(&self, t: f64) -> Option<Vec<f64>> {
        if t < self.points[0].time {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            if t >= self.points[i].time && t <= self.points[i + 1].time {
                let mut pt = Vec::<f64>::new();
                pt.resize(self.dim, 0.0);
                let dx = t - self.points[i].time;
                for (j, point_iter) in pt.iter_mut().enumerate() {
                    *point_iter = 2.0 * self.c[i][j] + 6.0 * self.d[i][j] * dx;
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
