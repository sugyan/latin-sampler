import { useEffect, useRef, useState } from "react";
import init, { WasmSampler } from "latin-sampler";

function cellColor(value: number, n: number): string {
  const hue = (value / n) * 360;
  return `hsl(${hue}, 75%, 55%)`;
}

function cellSize(n: number): number {
  if (n <= 5) return 56;
  if (n <= 10) return 44;
  if (n <= 15) return 36;
  return 28;
}

function App() {
  const [ready, setReady] = useState(false);
  const [n, setN] = useState(5);
  const [seed, setSeed] = useState(42);
  const [grid, setGrid] = useState<number[][] | null>(null);
  const [sampleNumber, setSampleNumber] = useState(0);
  const [dirty, setDirty] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [hasSampler, setHasSampler] = useState(false);
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

  const needsReset = !hasSampler || dirty;

  const handleClick = () => {
    setError(null);
    try {
      if (needsReset) {
        samplerRef.current?.free();
        samplerRef.current = new WasmSampler(n, BigInt(seed));
        setHasSampler(true);
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
    return (
      <div className="min-h-screen flex items-center justify-center">
        <p className="text-neutral-400 animate-pulse text-lg">
          Loading wasm...
        </p>
      </div>
    );
  }

  const size = grid ? cellSize(grid.length) : 0;

  return (
    <div className="min-h-screen flex items-center justify-center p-4">
      <div className="bg-surface rounded-2xl border border-border p-8 shadow-2xl w-full max-w-2xl animate-fade-in">
        <h1 className="text-3xl font-bold bg-gradient-to-r from-indigo-400 to-purple-400 bg-clip-text text-transparent mb-1">
          Latin Sampler
        </h1>
        <p className="text-sm text-neutral-500 mb-6">
          Jacobson-Matthews MCMC Algorithm
        </p>

        <div className="flex gap-3 items-end flex-wrap mb-6">
          <label className="flex flex-col gap-1.5 text-sm font-semibold text-neutral-300">
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
              className="w-28 px-3 py-2 bg-neutral-800 border border-neutral-700 text-white rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent transition-colors"
            />
          </label>
          <label className="flex flex-col gap-1.5 text-sm font-semibold text-neutral-300">
            seed
            <input
              type="number"
              min={0}
              value={seed}
              onChange={(e) => {
                setSeed(Number(e.target.value));
                setDirty(true);
              }}
              className="w-28 px-3 py-2 bg-neutral-800 border border-neutral-700 text-white rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-accent focus:border-transparent transition-colors"
            />
          </label>
          <button
            onClick={handleClick}
            className="px-5 py-2 bg-accent hover:bg-accent-hover text-white font-semibold rounded-lg text-sm transition-colors cursor-pointer disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {needsReset ? "Generate" : "Next"}
          </button>
        </div>

        {error && (
          <div className="bg-red-950/50 border border-red-800 text-red-300 rounded-lg px-4 py-3 text-sm mb-4">
            {error}
          </div>
        )}

        {grid && (
          <div className="animate-fade-in">
            <p className="mb-3">
              <span className="bg-neutral-800 text-neutral-300 rounded-full px-3 py-1 font-mono text-sm">
                #{sampleNumber}
              </span>
            </p>
            <div
              className="inline-grid gap-0.5"
              style={{
                gridTemplateColumns: `repeat(${grid.length}, ${size}px)`,
              }}
            >
              {grid.flatMap((row, i) =>
                row.map((val, j) => (
                  <div
                    key={`${i}-${j}`}
                    className="flex items-center justify-center rounded-lg font-semibold text-white animate-scale-in"
                    style={{
                      width: size,
                      height: size,
                      fontSize: size <= 28 ? "0.7rem" : "0.85rem",
                      backgroundColor: cellColor(val, grid.length),
                      textShadow: "0 1px 2px rgba(0,0,0,0.4)",
                      transition: "background-color 0.3s ease",
                    }}
                  >
                    {val}
                  </div>
                ))
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default App;
