use trajectory::CubicSpline;
use trajectory::Trajectory;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rec = rerun::RecordingStreamBuilder::new("trajectory").spawn()?;

    let times = vec![0.0, 1.0, 3.0, 4.0];
    let points = vec![
        vec![0.0, -1.0],
        vec![2.0, -3.0],
        vec![3.0, 3.0],
        vec![1.0, 5.0],
    ];
    let ip = CubicSpline::new(times, points).unwrap();

    for i in 0..400 {
        let t = i as f64 * 0.01;
        let p = ip.position(t).unwrap();
        let v = ip.velocity(t).unwrap();
        let a = ip.acceleration(t).unwrap();

        rec.set_time_seconds("time", t);
        rec.log("x/trajectory/position", &rerun::Scalar::new(p[0]))?;
        rec.log("y/trajectory/position", &rerun::Scalar::new(p[1]))?;
        rec.log("x/trajectory/velocity", &rerun::Scalar::new(v[0]))?;
        rec.log("y/trajectory/velocity", &rerun::Scalar::new(v[1]))?;
        rec.log("x/trajectory/acceleration", &rerun::Scalar::new(a[0]))?;
        rec.log("y/trajectory/acceleration", &rerun::Scalar::new(a[1]))?;
    }

    Ok(())
}
