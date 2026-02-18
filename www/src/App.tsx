import { useEffect, useRef, useState } from "react";
import init, { WasmSampler } from "latin-sampler";

function cellColor(value: number, n: number): string {
  const hue = (value / n) * 360;
  return `oklch(0.72 0.16 ${hue})`;
}

function cellSize(n: number): number {
  return Math.min(56, Math.max(20, Math.round(400 / n)));
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
  const [hoverValue, setHoverValue] = useState<number | null>(null);
  const samplerRef = useRef<WasmSampler | null>(null);

  useEffect(() => {
    init().then(() => setReady(true));
  }, []);

  useEffect(() => {
    return () => {
      samplerRef.current?.free();
      samplerRef.current = null;
      setHasSampler(false);
    };
  }, []);

  const needsReset = !hasSampler || dirty;

  const handleClick = () => {
    setError(null);
    if (n < 2 || n > 20) {
      setError("n must be between 2 and 20");
      return;
    }
    try {
      if (needsReset) {
        samplerRef.current?.free();
        samplerRef.current = null;
        samplerRef.current = new WasmSampler(n, BigInt(seed));
        setHasSampler(true);
        setDirty(false);
        setSampleNumber(0);
      }
      const result = samplerRef.current!.next();
      setGrid(result);
      setSampleNumber((prev) => prev + 1);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
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
            disabled={!ready}
            className="px-5 py-2 bg-accent hover:bg-accent-hover text-white font-semibold rounded-lg text-sm transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
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
              className="inline-grid rounded-xl"
              style={{
                gridTemplateColumns: `repeat(${grid.length}, ${size}px)`,
                gap: grid.length <= 10 ? "3px" : "2px",
                background: "oklch(0.18 0.01 260)",
                border: "1px solid oklch(0.3 0.02 260)",
                padding: grid.length <= 10 ? "3px" : "2px",
              }}
              onMouseLeave={() => setHoverValue(null)}
            >
              {grid.flatMap((row, i) =>
                row.map((val, j) => {
                  const dimmed =
                    hoverValue !== null && hoverValue !== val;
                  const addDepth = grid.length <= 15;
                  return (
                    <div
                      key={`${i}-${j}`}
                      className="flex items-center justify-center font-semibold text-white animate-scale-in"
                      style={{
                        width: size,
                        height: size,
                        borderRadius: Math.round(size * 0.1),
                        fontSize: `${Math.max(10, Math.round(size * 0.35))}px`,
                        backgroundColor: cellColor(val, grid.length),
                        textShadow: "0 1px 2px rgba(0,0,0,0.4)",
                        opacity: dimmed ? 0.25 : 1,
                        transition:
                          "opacity 0.15s ease, background-color 0.3s ease",
                        border: addDepth
                          ? "1px solid rgba(255,255,255,0.1)"
                          : "none",
                        boxShadow: addDepth
                          ? "inset 0 1px 0 rgba(255,255,255,0.06), 0 2px 4px rgba(0,0,0,0.3)"
                          : "0 1px 3px rgba(0,0,0,0.2)",
                        animationDelay: `${i * 40}ms`,
                        animationFillMode: "backwards",
                        cursor: "default",
                      }}
                      onMouseEnter={() => setHoverValue(val)}
                    >
                      {val}
                    </div>
                  );
                })
              )}
            </div>
          </div>
        )}
        <footer className="mt-6 pt-4 border-t border-border flex items-center gap-4 text-xs text-neutral-500">
          <a
            href="https://github.com/sugyan/latin-sampler"
            target="_blank"
            rel="noopener noreferrer"
            className="flex items-center gap-1.5 hover:text-neutral-300 transition-colors"
            aria-label="GitHub repository"
          >
            <svg viewBox="0 0 16 16" width="14" height="14" fill="currentColor" aria-hidden="true">
              <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z" />
            </svg>
            sugyan/latin-sampler
          </a>
          <a
            href="https://crates.io/crates/latin-sampler"
            target="_blank"
            rel="noopener noreferrer"
            className="flex items-center gap-1.5 hover:text-neutral-300 transition-colors"
            aria-label="crates.io"
          >
            <svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" aria-hidden="true">
              <path d="M21 8.25l-9-5.25-9 5.25v7.5l9 5.25 9-5.25v-7.5zm-9 3.75l-6-3.5 6-3.5 6 3.5-6 3.5zm-7 1.75l6 3.5v5l-6-3.5v-5zm8 8.5v-5l6-3.5v5l-6 3.5z" />
            </svg>
            crates.io
          </a>
        </footer>
      </div>
    </div>
  );
}

export default App;
