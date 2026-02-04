# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
cargo build                    # Build the library
cargo test                     # Run all tests
cargo test <test_name>         # Run a specific test
cargo bench                    # Run benchmarks (requires nightly)
```

### Examples

```bash
cargo run --release --example generate -- <n> [seed]     # Generate a Latin square
cargo run --release --example uniformity_test -- <n>     # Chi-square uniformity test
cargo run --release --example independence_check         # Autocorrelation analysis
cargo run --release --example benchmark_comparison -- [n] # sample() vs Sampler comparison
```

## Architecture

This is an MCMC sampler for generating approximately uniform Latin squares using the Jacobson-Matthews algorithm.

### Core Modules

- **`jacobson_matthews.rs`**: The Jacobson-Matthews algorithm implementation using a 3D {-1, 0, 1} array representation (`JMState`). Can pass through "improper" states where one position has value -1, enabling ergodicity.

- **`sampler.rs`**: Public API with two sampling modes:
  - `sample()`: One-shot function with fresh burn-in each call
  - `Sampler`: Iterator that performs burn-in once, then uses thinning between samples

- **`square.rs`**: `LatinSquare` type representing an n×n array where each row and column is a permutation of {0..n-1}.

### Key Concepts

- **Burn-in** (default n³): Steps to reach equilibrium from initial cyclic state
- **Thinning** (default 3n²): Steps between consecutive samples for independence
- **p_do_nothing** (default 0.01): Probability of null move for aperiodicity

### Constraints

- Order n must be in range 2..=255 (stored as u8)
- Requires Rust 1.85+ (edition 2024)
