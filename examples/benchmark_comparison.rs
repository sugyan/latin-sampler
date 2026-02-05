//! Benchmark comparison between `sample()` and `Sampler` iterator.
//!
//! Usage: cargo run --release --example benchmark_comparison -- [n]
//!
//! This example measures the performance difference between:
//! - `sample()`: One-shot function that performs burn-in (n³) each call
//! - `Sampler`: Iterator that performs burn-in once, then thinning (3n²) between samples
//!
//! For generating multiple samples, the iterator is significantly faster due to
//! amortized burn-in cost.

use latin_sampler::{Sampler, SamplerParams, sample};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::env;
use std::time::{Duration, Instant};

const TRIALS: usize = 5;

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(10);

    println!("Benchmark: sample() vs Sampler iterator (n={})", n);
    println!("=========================================================");
    println!();

    for k in [1, 10, 100, 1_000, 10_000] {
        benchmark_k_samples(n, k);
    }

    println!("Theoretical Analysis (default params: burn_in=n³, thinning=3n²):");
    println!("  n={}, burn_in={}, thinning={}", n, n * n * n, 3 * n * n);
    println!();
    for k in [1, 10, 100, 1_000, 10_000] {
        let sample_steps = k * n * n * n;
        let iter_steps = n * n * n + (k - 1) * 3 * n * n;
        let speedup = sample_steps as f64 / iter_steps as f64;
        println!(
            "  k={:>5}: {:.2}x expected (sample: {:>9} steps, iter: {:>9} steps)",
            k, speedup, sample_steps, iter_steps
        );
    }
}

fn benchmark_k_samples(n: usize, k: usize) {
    println!("--- Generating {} sample(s) ---", k);

    let params = SamplerParams::default();

    // Benchmark sample() (one-shot)
    let sample_times = run_trials(TRIALS, || {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        for _ in 0..k {
            let sq = sample(n, &mut rng, &params);
            std::hint::black_box(sq);
        }
    });

    // Benchmark Sampler iterator
    let iter_times = run_trials(TRIALS, || {
        let rng = ChaCha20Rng::seed_from_u64(42);
        let sampler = Sampler::new(n, rng, params.clone());
        for sq in sampler.take(k) {
            std::hint::black_box(sq);
        }
    });

    let sample_stats = Stats::from(&sample_times);
    let iter_stats = Stats::from(&iter_times);

    println!(
        "  sample()   mean={:>10.2}ms  min={:>10.2}ms  max={:>10.2}ms",
        sample_stats.mean.as_secs_f64() * 1000.0,
        sample_stats.min.as_secs_f64() * 1000.0,
        sample_stats.max.as_secs_f64() * 1000.0
    );
    println!(
        "  Sampler    mean={:>10.2}ms  min={:>10.2}ms  max={:>10.2}ms",
        iter_stats.mean.as_secs_f64() * 1000.0,
        iter_stats.min.as_secs_f64() * 1000.0,
        iter_stats.max.as_secs_f64() * 1000.0
    );

    let speedup = sample_stats.mean.as_secs_f64() / iter_stats.mean.as_secs_f64();
    println!("  Speedup (Sampler vs sample()): {:.2}x", speedup);
    println!();
}

fn run_trials<F: FnMut()>(trials: usize, mut f: F) -> Vec<Duration> {
    (0..trials)
        .map(|_| {
            let start = Instant::now();
            f();
            start.elapsed()
        })
        .collect()
}

struct Stats {
    mean: Duration,
    min: Duration,
    max: Duration,
}

impl Stats {
    fn from(times: &[Duration]) -> Self {
        let sum: Duration = times.iter().sum();
        let mean = sum / times.len() as u32;
        let min = *times.iter().min().unwrap();
        let max = *times.iter().max().unwrap();
        Self { mean, min, max }
    }
}
