// Futoshiki: fill an n×n grid with 1..n so each row and column holds every
// number once, while respecting the inequality signs between neighbours.
// `givens` is a flat array (0 = empty). `h[r*(n-1)+c]` relates cell (r,c) to
// (r,c+1) and `v[r*n+c]` relates (r,c) to (r+1,c): 0 = no sign,
// LT = first < second, GT = first > second.

import { shuffle, range, popcount } from './rng.js';

export const LT = 1;
export const GT = 2;

export const meta = {
  id: 'futoshiki',
  name: 'Futoshiki',
  tagline: 'A Latin square with greater-than signs.',
  rules: [
    'Fill each row and each column with the numbers 1 to N, each exactly once.',
    'Every inequality sign between two cells must hold (the sign points at the smaller number).',
  ],
  options: [
    { id: '4', label: '4×4' },
    { id: '5', label: '5×5' },
    { id: '6', label: '6×6' },
    { id: '7', label: '7×7' },
  ],
  defaultOption: '5',
};

// Constraints as a list per cell: [otherCell, relation] meaning cell REL other.
function relations(n, h, v) {
  const rel = range(n * n).map(() => []);
  for (let r = 0; r < n; r++) {
    for (let c = 0; c < n - 1; c++) {
      const s = h[r * (n - 1) + c];
      if (!s) continue;
      const a = r * n + c, b = a + 1;
      rel[a].push([b, s]);
      rel[b].push([a, s === LT ? GT : LT]);
    }
  }
  for (let r = 0; r < n - 1; r++) {
    for (let c = 0; c < n; c++) {
      const s = v[r * n + c];
      if (!s) continue;
      const a = r * n + c, b = a + n;
      rel[a].push([b, s]);
      rel[b].push([a, s === LT ? GT : LT]);
    }
  }
  return rel;
}

const lowBit = (m) => 31 - Math.clz32(m & -m);
const highBit = (m) => 31 - Math.clz32(m);

// Domains are bitmasks per cell (bit d = value d+1). Propagation repeats until
// nothing changes: placed values leave their row and column, a value with a
// single possible place in a line goes there, and every sign trims the bounds
// of both of its cells.
function propagate(dom, n, rel, lines) {
  for (let changed = true; changed;) {
    changed = false;
    for (let i = 0; i < n * n; i++) {
      let m = dom[i];
      if (!m) return false;
      for (const [j, s] of rel[i]) {
        const o = dom[j];
        if (!o) return false;
        m &= s === LT ? (1 << highBit(o)) - 1 : ~((2 << lowBit(o)) - 1);
      }
      if (m !== dom[i]) {
        if (!m) return false;
        dom[i] = m; changed = true;
      }
    }
    for (const line of lines) {
      for (const i of line) {
        const m = dom[i];
        if (m & (m - 1)) continue;
        for (const j of line) {
          if (j !== i && dom[j] & m) {
            dom[j] &= ~m;
            if (!dom[j]) return false;
            changed = true;
          }
        }
      }
      for (let d = 0; d < n; d++) {
        const bit = 1 << d;
        let place = -1, count = 0;
        for (const i of line) if (dom[i] & bit) { place = i; count++; }
        if (count === 0) return false;
        if (count === 1 && dom[place] !== bit) { dom[place] = bit; changed = true; }
      }
    }
  }
  return true;
}

export function solve(puzzle, limit = 2, rng = null) {
  const { size: n, givens } = puzzle;
  const h = puzzle.h || new Array(n * (n - 1)).fill(0);
  const v = puzzle.v || new Array(n * n).fill(0);
  const rel = relations(n, h, v);
  const full = (1 << n) - 1;
  const lines = [
    ...range(n).map((r) => range(n).map((c) => r * n + c)),
    ...range(n).map((c) => range(n).map((r) => r * n + c)),
  ];
  let count = 0;
  let solution = null;

  function rec(dom) {
    if (!propagate(dom, n, rel, lines)) return;
    let best = -1, bestCount = 99;
    for (let i = 0; i < n * n; i++) {
      const k = popcount(dom[i]);
      if (k > 1 && k < bestCount) { best = i; bestCount = k; }
    }
    if (best === -1) {
      count++;
      if (!solution) solution = Array.from(dom, (m) => lowBit(m) + 1);
      return;
    }
    const vals = [];
    for (let d = 0; d < n; d++) if (dom[best] & (1 << d)) vals.push(d);
    if (rng) shuffle(vals, rng);
    for (const d of vals) {
      const next = dom.slice();
      next[best] = 1 << d;
      rec(next);
      if (count >= limit) return;
    }
  }

  rec(Int32Array.from(givens, (g) => (g ? 1 << (g - 1) : full)));
  return { count, solution };
}

export function countSolutions(puzzle, limit = 2) {
  return solve(puzzle, limit).count;
}

export function generate(option, rng) {
  const n = Number(option) || 5;
  const solution = solve({ size: n, givens: new Array(n * n).fill(0) }, 1, rng).solution;
  const rel = (a, b) => (solution[a] < solution[b] ? LT : GT);

  // Start from every number and every sign, then strip clues while the
  // solution stays unique.
  const puzzle = {
    size: n,
    givens: solution.slice(),
    h: range(n * (n - 1)).map((k) => {
      const r = Math.floor(k / (n - 1)), c = k % (n - 1);
      return rel(r * n + c, r * n + c + 1);
    }),
    v: range(n * n).map((k) => (k < n * (n - 1) ? rel(k, k + n) : 0)),
  };
  // Numbers go first so that the finished puzzle leans on its signs.
  const clues = [
    ...shuffle(range(n * n).map((i) => ['givens', i]), rng),
    ...shuffle([
      ...range(n * (n - 1)).map((i) => ['h', i]),
      ...range(n * (n - 1)).map((i) => ['v', i]),
    ], rng),
  ];
  for (const [key, i] of clues) {
    const keep = puzzle[key][i];
    puzzle[key][i] = 0;
    if (countSolutions(puzzle) !== 1) puzzle[key][i] = keep;
  }
  return { ...puzzle, solution };
}

// Cells that clash with a row/column duplicate or a violated sign.
export function conflicts(values, puzzle) {
  const n = puzzle.size;
  const rel = relations(n, puzzle.h, puzzle.v);
  const bad = new Set();
  for (let i = 0; i < n * n; i++) {
    if (!values[i]) continue;
    for (let j = i + 1; j < n * n; j++) {
      if (values[i] === values[j] && (Math.floor(i / n) === Math.floor(j / n) || i % n === j % n)) {
        bad.add(i); bad.add(j);
      }
    }
    for (const [j, s] of rel[i]) {
      if (values[j] && (s === LT ? values[i] >= values[j] : values[i] <= values[j])) {
        bad.add(i); bad.add(j);
      }
    }
  }
  return bad;
}
