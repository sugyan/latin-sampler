import { useEffect, useRef, useState } from "react";
import init, { WasmSampler } from "latin-sampler";

function cellColor(value: number, n: number): string {
  const hue = (value / n) * 360;
  return `hsl(${hue}, 65%, 50%)`;
}

function App() {
  const [ready, setReady] = useState(false);
  const [n, setN] = useState(5);
  const [seed, setSeed] = useState(42);
  const [grid, setGrid] = useState<number[][] | null>(null);
  const [sampleNumber, setSampleNumber] = useState(0);
  const [dirty, setDirty] = useState(false);
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

  const needsReset = !samplerRef.current || dirty;

  const handleClick = () => {
    setError(null);
    try {
      if (needsReset) {
        samplerRef.current?.free();
        samplerRef.current = new WasmSampler(n, BigInt(seed));
        setDirty(false);
        setSampleNumber(0);
      }
      const result = samplerRef.current!.next();
      setGrid(result);
      setSampleNumber((prev) => prev + 1);
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

      <div className="controls">
        <label>
          n (2-20)
          <input
            type="number"
            min={2}
            max={20}
            value={n}
            onChange={(e) => {
              setN(Number(e.target.value));
              setDirty(true);
            }}
          />
        </label>
        <label>
          seed
          <input
            type="number"
            min={0}
            value={seed}
            onChange={(e) => {
              setSeed(Number(e.target.value));
              setDirty(true);
            }}
          />
        </label>
        <button onClick={handleClick}>{needsReset ? "Generate" : "Next"}</button>
      </div>

      {error && <p className="error">{error}</p>}

      {grid && (
        <>
          <p className="sample-number">#{sampleNumber}</p>
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
        </>
      )}
    </>
  );
}

export default App;
