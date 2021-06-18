use super::traits::*;
use super::utils::*;
use num_traits::Float;

#[derive(Debug)]
pub struct Linear<T>
where
    T: Float,
{
    times: Vec<T>,
    points: Vec<Vec<T>>,
}

impl<T> Linear<T>
where
    T: Float,
{
    pub fn new(times: Vec<T>, points: Vec<Vec<T>>) -> Option<Self> {
        if !is_inputs_valid(&times, &points) {
            return None;
        }
        Some(Self { times, points })
    }
}

impl<T> Trajectory for Linear<T>
where
    T: Float,
{
    type Point = Vec<T>;
    type Time = T;

    fn position(&self, t: T) -> Option<Vec<T>> {
        if t < self.times[0] {
            return None;
        }
        for i in 0..(self.points.len() - 1) {
            let t0 = self.times[i];
            let t1 = self.times[i + 1];
            if t >= t0 && t <= t1 {
                let mut pt = self.points[i].clone();
                let p0 = &self.points[i];
                let p1 = &self.points[i + 1];
                for j in 0..p0.len() {
                    pt[j] = p0[j] + (p1[j] - p0[j]) / (t1 - t0) * (t - t0);
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
            let t0 = self.times[i];
            let t1 = self.times[i + 1];
            if t >= t0 && t <= t1 {
                let mut pt = vec![T::zero(); dim];
                let p0 = &self.points[i];
                let p1 = &self.points[i + 1];
                for j in 0..p0.len() {
                    pt[j] = (p1[j] - p0[j]) / (t1 - t0);
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
        for tm in &self.times {
            if t == *tm {
                return None;
            }
        }
        Some(vec![T::zero(); dim])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_linear() {
        let times = vec![0.0, 1.0, 3.0, 4.0];
        let points = vec![
            vec![0.0, -1.0],
            vec![2.0, -3.0],
            vec![3.0, 3.0],
            vec![1.0, 5.0],
        ];
        let ip = Linear::new(times, points).unwrap();
        let p0 = ip.position(-0.5);
        assert!(p0.is_none());
        let p1 = ip.position(0.5).unwrap();
        assert_approx_eq!(p1[0], 1.0);
        assert_approx_eq!(p1[1], -2.0);
        let p2 = ip.position(1.0).unwrap();
        assert_approx_eq!(p2[0], 2.0);
        assert_approx_eq!(p2[1], -3.0);
        let p3 = ip.position(3.0).unwrap();
        assert_approx_eq!(p3[0], 3.0);
        assert_approx_eq!(p3[1], 3.0);
        let p4 = ip.position(3.5).unwrap();
        assert_approx_eq!(p4[0], 2.0);
        assert_approx_eq!(p4[1], 4.0);
    }
}
