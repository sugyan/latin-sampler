//! Benchmarks for Latin square sampling.
//!
//! This file benchmarks individual operations for regression testing:
//! - `sample()`: One-shot sampling (includes burn-in each call)
//! - `Sampler` first `next()`: Initial sample with burn-in
//! - `Sampler` subsequent `next()`: Samples with thinning only
//!
//! For comparing `sample()` vs `Sampler` for generating k samples, see:
//!   cargo run --release --example benchmark_comparison

#![feature(test)]

extern crate test;

use latin_sampler::{Sampler, SamplerParams, sample};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use test::Bencher;

#[bench]
fn bench_sample_n4(b: &mut Bencher) {
    let params = SamplerParams::default();
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    b.iter(|| {
        let sq = sample(4, &mut rng, &params);
        test::black_box(sq)
    });
}

#[bench]
fn bench_sample_n7(b: &mut Bencher) {
    let params = SamplerParams::default();
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    b.iter(|| {
        let sq = sample(7, &mut rng, &params);
        test::black_box(sq)
    });
}

#[bench]
fn bench_sample_n10(b: &mut Bencher) {
    let params = SamplerParams::default();
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    b.iter(|| {
        let sq = sample(10, &mut rng, &params);
        test::black_box(sq)
    });
}

#[bench]
fn bench_sampler_first_n10(b: &mut Bencher) {
    let params = SamplerParams::default();

    b.iter(|| {
        let rng = ChaCha20Rng::seed_from_u64(42);
        let mut sampler = Sampler::new(10, rng, params.clone());
        let sq = sampler.next().unwrap();
        test::black_box(sq)
    });
}

#[bench]
fn bench_sampler_subsequent_n10(b: &mut Bencher) {
    let params = SamplerParams::default();
    let rng = ChaCha20Rng::seed_from_u64(42);
    let mut sampler = Sampler::new(10, rng, params);

    // Consume first sample (with burn-in) outside the benchmark loop
    let _ = sampler.next();

    b.iter(|| {
        let sq = sampler.next().unwrap();
        test::black_box(sq)
    });
}
