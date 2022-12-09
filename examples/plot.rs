use gnuplot::{Caption, Color, Figure};
use trajectory::CubicSpline;
use trajectory::Trajectory;

fn main() {
    let times = vec![0.0, 1.0, 3.0, 4.0];
    let points = vec![
        vec![0.0, -1.0],
        vec![2.0, -3.0],
        vec![3.0, 3.0],
        vec![1.0, 5.0],
    ];
    let ip = CubicSpline::new(times, points).unwrap();

    let mut times = Vec::new();
    let mut positions0 = Vec::new();
    let mut positions1 = Vec::new();
    let mut velocities0 = Vec::new();
    let mut velocities1 = Vec::new();
    let mut accelerations0 = Vec::new();
    let mut accelerations1 = Vec::new();
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
        accelerations0.push(a[0]);
        accelerations1.push(a[1]);
    }
    let mut fg = Figure::new();
    fg.axes2d()
        .lines(&times, &positions0, &[Caption("Position"), Color("red")])
        .lines(&times, &velocities0, &[Caption("Velocity"), Color("green")])
        .lines(
            &times,
            &accelerations0,
            &[Caption("Acceleration"), Color("blue")],
        );
    let _ = fg.show();
    let mut fg = Figure::new();
    fg.axes2d()
        .lines(&times, &positions1, &[Caption("Position"), Color("red")])
        .lines(&times, &velocities1, &[Caption("Velocity"), Color("green")])
        .lines(
            &times,
            &accelerations1,
            &[Caption("Acceleration"), Color("blue")],
        );
    let _ = fg.show();
}
