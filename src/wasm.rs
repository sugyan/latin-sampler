use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use wasm_bindgen::prelude::*;

use crate::{LatinSquare, Sampler, SamplerParams, sample};

/// Convert a LatinSquare to a JsValue (2D array of u8).
fn latin_square_to_js(sq: &LatinSquare) -> Result<JsValue, JsError> {
    let rows: Vec<Vec<u8>> = (0..sq.n())
        .map(|r| (0..sq.n()).map(|c| sq.get(r, c)).collect())
        .collect();
    serde_wasm_bindgen::to_value(&rows).map_err(|e| JsError::new(&e.to_string()))
}

/// Generate a Latin square of order n with the given seed.
/// Returns a 2D array directly usable in JavaScript.
#[wasm_bindgen]
pub fn generate(n: u8, seed: u64) -> Result<JsValue, JsError> {
    if n < 2 {
        return Err(JsError::new("n must be at least 2"));
    }
    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let params = SamplerParams::default();
    let sq = sample(n as usize, &mut rng, &params);

    latin_square_to_js(&sq)
}

/// A stateful sampler that produces approximately uniform Latin squares.
///
/// Uses MCMC sampling with the Jacobson-Matthews algorithm.
/// Burn-in is performed automatically on the first call to `next()`,
/// and thinning steps are applied between subsequent samples.
#[wasm_bindgen]
pub struct WasmSampler {
    sampler: Sampler<ChaCha20Rng>,
}

#[wasm_bindgen]
impl WasmSampler {
    /// Create a new sampler for Latin squares of order `n`.
    ///
    /// `n` must be in range 2..=255. The `seed` determines the random sequence.
    #[wasm_bindgen(constructor)]
    pub fn new(n: u8, seed: u64) -> Result<WasmSampler, JsError> {
        if n < 2 {
            return Err(JsError::new("n must be at least 2"));
        }
        let rng = ChaCha20Rng::seed_from_u64(seed);
        let params = SamplerParams::default();
        let sampler = Sampler::new(n as usize, rng, params);
        Ok(WasmSampler { sampler })
    }

    /// Get the next Latin square sample.
    ///
    /// The first call performs burn-in; subsequent calls apply thinning steps.
    /// Returns a 2D array of u8 values.
    pub fn next(&mut self) -> Result<JsValue, JsError> {
        let sq = self
            .sampler
            .next()
            .expect("Sampler is an infinite iterator");
        latin_square_to_js(&sq)
    }
}
