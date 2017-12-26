extern crate gnuplot;
extern crate trajectory;

use gnuplot::{Figure, Caption, Color};
use trajectory::CSplineVectorInterpolator;
use trajectory::Trajectory;
type TimedVec = trajectory::TimedPoint<f64, Vec<f64>>;

fn main() {
    let points = vec![
        TimedVec::new(0.0, vec![0.0, -1.0]),
        TimedVec::new(1.0, vec![2.0, -3.0]),
        TimedVec::new(3.0, vec![3.0, 3.0]),
        TimedVec::new(4.0, vec![1.0, 5.0]),
    ];
    let ip = CSplineVectorInterpolator::new(points);

    let mut times = Vec::new();
    let mut positions0 = Vec::new();
    let mut positions1 = Vec::new();
    let mut velocities0 = Vec::new();
    let mut velocities1 = Vec::new();
    let mut accs0 = Vec::new();
    let mut accs1 = Vec::new();
    for i in 0..400 {
        let t = i as f64 * 0.01f64;
        let p = ip.position(t).unwrap();
        let v = ip.velocity(t).unwrap();
        let a = ip.acceleration(t).unwrap();
        times.push(t);
        positions0.push(p[0]);
        positions1.push(p[1]);
        velocities0.push(v[0]);
        velocities1.push(v[1]);
        accs0.push(a[0]);
        accs1.push(a[1]);
    }
    let mut fg = Figure::new();
    fg.axes2d()
        .lines(&times, &positions0, &[Caption("Position"), Color("red")])
        .lines(&times, &velocities0, &[Caption("Velocity"), Color("green")])
        .lines(&times, &accs0, &[Caption("Acceleration"), Color("blue")]);
    fg.show();
}