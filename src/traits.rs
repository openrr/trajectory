pub trait Trajectory {
    type Point;
    type Time;
    fn position(&self, t: Self::Time) -> Option<Self::Point>;
    fn velocity(&self, t: Self::Time) -> Option<Self::Point>;
    fn acceleration(&self, t: Self::Time) -> Option<Self::Point>;
}
