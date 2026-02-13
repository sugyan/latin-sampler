#![doc = include_str!("../README.md")]

mod jacobson_matthews;
mod sampler;
mod square;
#[cfg(target_arch = "wasm32")]
mod wasm;

pub use sampler::{Sampler, SamplerParams, sample};
pub use square::LatinSquare;
