extern crate nalgebra as na;

use na::DVector;

pub trait Interpolator {
    fn interpolate(&self, t: f64) -> Option<Vec<f64>>;
}

pub struct LinearInterpolator {
    times: Vec<f64>,
    points: Vec<Vec<f64>>,
}

impl LinearInterpolator {
    pub fn new(times: Vec<f64>, points: Vec<Vec<f64>>) -> Self {
        Self { times, points }
    }
}

impl Interpolator for LinearInterpolator {
    fn interpolate(&self, t: f64) -> Option<Vec<f64>> {
        if t < self.times[0] || t > *self.times.last().unwrap() {
            return None;
        }
        for i in 0..(self.times.len() - 1) {
            if t >= self.times[i] && t <= self.times[i + 1] {
                let dim = self.points[0].len();
                let p0 = DVector::from_column_slice(dim, &self.points[i]);
                let p1 = DVector::from_column_slice(dim, &self.points[i + 1]);
                let diff = (p1 - p0.clone()) / (self.times[i + 1] - self.times[i]) *
                    (t - self.times[i]);
                let pt = p0 + diff;
                return Some(pt.as_slice().to_vec());
            }
        }
        return None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        let times = vec![0.0, 1.0, 3.0, 4.0];
        let points = vec![
            vec![0.0, -1.0],
            vec![2.0, -3.0],
            vec![3.0, 3.0],
            vec![1.0, 5.0],
        ];
        let ip = LinearInterpolator::new(times, points);
        let p0 = ip.interpolate(-0.5);
        assert!(p0.is_none());
        let p1 = ip.interpolate(0.5).unwrap();
        assert_eq!(p1[0], 1.0);
        assert_eq!(p1[1], -2.0);
        let p2 = ip.interpolate(1.0).unwrap();
        assert_eq!(p2[0], 2.0);
        assert_eq!(p2[1], -3.0);
        let p3 = ip.interpolate(3.0).unwrap();
        assert_eq!(p3[0], 3.0);
        assert_eq!(p3[1], 3.0);
        let p4 = ip.interpolate(3.5).unwrap();
        assert_eq!(p4[0], 2.0);
        assert_eq!(p4[1], 4.0);
    }
}
