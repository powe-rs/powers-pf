pub trait ProgressMonitor {
    fn update(&self, i: usize, norm_f: f64);
}
