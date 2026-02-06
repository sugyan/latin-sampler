use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use wasm_bindgen::prelude::*;

use crate::{SamplerParams, sample};

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

    // Vec<Vec<u8>>に変換してJsValueへ
    let rows: Vec<Vec<u8>> = (0..sq.n())
        .map(|r| (0..sq.n()).map(|c| sq.get(r, c)).collect())
        .collect();

    serde_wasm_bindgen::to_value(&rows).map_err(|e| JsError::new(&e.to_string()))
}
