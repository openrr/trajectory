use num_traits::Float;

pub fn is_inputs_valid<T>(times: &Vec<T>, points: &Vec<Vec<T>>) -> bool
where
    T: Float,
{
    if times.len() != points.len() {
        return false;
    }
    for i in 0..times.len() - 1 {
        // not sorted time
        if times[i] >= times[i + 1] {
            return false;
        }
    }
    if points.is_empty() {
        return false;
    }
    if points[0].is_empty() {
        return false;
    }
    let dim = points[0].len();
    for p in points {
        if p.len() != dim {
            return false;
        }
    }
    true
}
