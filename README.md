# trajectory

[![Build Status](https://img.shields.io/github/actions/workflow/status/openrr/trajectory/ci.yml?branch=main&logo=github)](https://github.com/openrr/trajectory/actions) [![crates.io](https://img.shields.io/crates/v/trajectory.svg?logo=rust)](https://crates.io/crates/trajectory) [![docs](https://docs.rs/trajectory/badge.svg)](https://docs.rs/trajectory) [![discord](https://dcbadge.vercel.app/api/server/8DAFFKc88B?style=flat)](https://discord.gg/8DAFFKc88B)

Trajectory interpolator for Rust.

## Code example

```rust
use trajectory::{CubicSpline, Trajectory};

let times = vec![0.0_f64, 1.0, 3.0, 4.0];
let points = vec![
    vec![0.0, -1.0],
    vec![2.0, -3.0],
    vec![3.0, 3.0],
    vec![1.0, 5.0],
];
let ip = CubicSpline::new(times, points).unwrap();
for i in 0..400 {
    let t = i as f64 * 0.01_f64;
    let p = ip.position(t).unwrap();
    let v = ip.velocity(t).unwrap();
    let a = ip.acceleration(t).unwrap();
}
```

## Run example

It requires `gnuplot`.

```bash
cargo run --example plot
```

![plot1](https://github.com/openrr/trajectory/raw/main/img/plot1.png)
![plot2](https://github.com/openrr/trajectory/raw/main/img/plot2.png)

## `OpenRR` Community

[Here](https://discord.gg/8DAFFKc88B) is a discord server for `OpenRR` users and developers.
