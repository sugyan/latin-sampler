import { useEffect, useRef, useState } from "react";
import init, { generate, WasmSampler } from "latin-sampler";

type Mode = "oneshot" | "sampler";

function cellColor(value: number, n: number): string {
  const hue = (value / n) * 360;
  return `hsl(${hue}, 65%, 50%)`;
}

function App() {
  const [ready, setReady] = useState(false);
  const [mode, setMode] = useState<Mode>("oneshot");
  const [n, setN] = useState(5);
  const [seed, setSeed] = useState(42);
  const [grid, setGrid] = useState<number[][] | null>(null);
  const [error, setError] = useState<string | null>(null);
  const samplerRef = useRef<WasmSampler | null>(null);

  useEffect(() => {
    init().then(() => setReady(true));
  }, []);

  useEffect(() => {
    return () => {
      samplerRef.current?.free();
      samplerRef.current = null;
    };
  }, []);

  const handleGenerate = () => {
    setError(null);
    try {
      const result = generate(n, BigInt(seed));
      setGrid(result);
    } catch (e) {
      setError(String(e));
    }
  };

  const handleCreateSampler = () => {
    setError(null);
    try {
      samplerRef.current?.free();
      samplerRef.current = new WasmSampler(n, BigInt(seed));
      const result = samplerRef.current.next();
      setGrid(result);
    } catch (e) {
      setError(String(e));
    }
  };

  const handleNext = () => {
    setError(null);
    try {
      if (!samplerRef.current) return;
      const result = samplerRef.current.next();
      setGrid(result);
    } catch (e) {
      setError(String(e));
    }
  };

  if (!ready) {
    return <p>Loading wasm...</p>;
  }

  return (
    <>
      <h1>Latin Sampler Demo</h1>

      <div className="mode-toggle">
        <button
          className={mode === "oneshot" ? "active" : ""}
          onClick={() => setMode("oneshot")}
        >
          One-shot
        </button>
        <button
          className={mode === "sampler" ? "active" : ""}
          onClick={() => setMode("sampler")}
        >
          Sampler
        </button>
      </div>

      <div className="controls">
        <label>
          n (2-255)
          <input
            type="number"
            min={2}
            max={255}
            value={n}
            onChange={(e) => setN(Number(e.target.value))}
          />
        </label>
        <label>
          seed
          <input
            type="number"
            min={0}
            value={seed}
            onChange={(e) => setSeed(Number(e.target.value))}
          />
        </label>
        {mode === "oneshot" ? (
          <button onClick={handleGenerate}>Generate</button>
        ) : (
          <>
            <button onClick={handleCreateSampler}>Create Sampler</button>
            <button onClick={handleNext} disabled={!samplerRef.current}>
              Next
            </button>
          </>
        )}
      </div>

      {error && <p className="error">{error}</p>}

      {grid && (
        <div
          className="grid"
          style={{ gridTemplateColumns: `repeat(${grid.length}, 40px)` }}
        >
          {grid.flatMap((row, i) =>
            row.map((val, j) => (
              <div
                key={`${i}-${j}`}
                className="cell"
                style={{ backgroundColor: cellColor(val, grid.length) }}
              >
                {val}
              </div>
            ))
          )}
        </div>
      )}
    </>
  );
}

export default App;
