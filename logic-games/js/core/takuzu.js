// Takuzu (a.k.a. Binairo / Tango): fill an even-sized grid with two symbols so
// that every row and column holds equally many of each, no three identical
// symbols are adjacent in a line, and no two rows or two columns are identical.
// Grids are flat arrays with -1 for empty, 0 and 1 for the two symbols.

import { shuffle, range, popcount } from './rng.js';

export const meta = {
  id: 'takuzu',
  name: 'Takuzu',
  tagline: 'Balance suns and moons in every line.',
  rules: [
    'Fill every cell with a sun or a moon.',
    'Each row and each column holds the same number of suns and moons.',
    'No more than two identical symbols may be next to each other in a row or column.',
    'No two rows are identical, and no two columns are identical.',
    'Tap a cell to cycle sun → moon → empty.',
  ],
  options: [
    { id: '6', label: '6×6' },
    { id: '8', label: '8×8' },
    { id: '10', label: '10×10' },
    { id: '12', label: '12×12' },
  ],
  defaultOption: '8',
};

function lineCells(n, isRow, k) {
  return isRow ? range(n).map((c) => k * n + c) : range(n).map((r) => r * n + k);
}

// All valid lines of length n as bitmasks: balanced and without three in a row.
const patternCache = new Map();
function patterns(n) {
  if (!patternCache.has(n)) {
    const out = [];
    const tripleMask = 0b111;
    for (let p = 0; p < 1 << n; p++) {
      if (popcount(p) !== n / 2) continue;
      let ok = true;
      for (let s = 0; s + 3 <= n && ok; s++) {
        const t = (p >> s) & tripleMask;
        ok = t !== 0 && t !== tripleMask;
      }
      if (ok) out.push(p);
    }
    patternCache.set(n, out);
  }
  return patternCache.get(n);
}

// Lines 0..n-1 are rows, n..2n-1 are columns; bit j of a line is its j-th cell.
// Each line keeps a domain of patterns still compatible with the known cells.
function propagate(g, n, domains, cells) {
  for (let changed = true; changed;) {
    changed = false;
    for (let L = 0; L < 2 * n; L++) {
      const lc = cells[L];
      let kmask = 0, kval = 0;
      for (let j = 0; j < n; j++) {
        const v = g[lc[j]];
        if (v !== -1) { kmask |= 1 << j; if (v) kval |= 1 << j; }
      }
      const dom = domains[L];
      const filtered = dom.filter((p) => (p & kmask) === kval);
      if (!filtered.length) return false;
      if (filtered.length !== dom.length) { domains[L] = filtered; changed = true; }

      let and = -1, or = 0;
      for (const p of filtered) { and &= p; or |= p; }
      for (let j = 0; j < n; j++) {
        if (kmask & (1 << j)) continue;
        if (and & (1 << j)) { g[lc[j]] = 1; changed = true; }
        else if (!(or & (1 << j))) { g[lc[j]] = 0; changed = true; }
      }

      // A settled line removes its pattern from every parallel line.
      if (filtered.length === 1) {
        const p = filtered[0];
        const base = L < n ? 0 : n;
        for (let M = base; M < base + n; M++) {
          if (M === L || !domains[M].includes(p)) continue;
          domains[M] = domains[M].filter((q) => q !== p);
          if (!domains[M].length) return false;
          changed = true;
        }
      }
    }
  }
  return true;
}

export function solve(grid, n, limit = 2, rng = null) {
  const cells = [...range(n).map((k) => lineCells(n, true, k)), ...range(n).map((k) => lineCells(n, false, k))];
  let count = 0;
  let solution = null;

  function rec(g, domains) {
    if (!propagate(g, n, domains, cells)) return;
    let best = -1;
    for (let L = 0; L < 2 * n; L++) {
      if (domains[L].length > 1 && (best === -1 || domains[L].length < domains[best].length)) best = L;
    }
    if (best === -1) {
      count++;
      if (!solution) solution = g.slice();
      return;
    }
    const options = rng ? shuffle(domains[best].slice(), rng) : domains[best];
    for (const p of options) {
      const g2 = g.slice();
      cells[best].forEach((i, j) => { g2[i] = (p >> j) & 1; });
      const d2 = domains.slice();
      d2[best] = [p];
      rec(g2, d2);
      if (count >= limit) return;
    }
  }

  rec(grid.slice(), new Array(2 * n).fill(patterns(n)));
  return { count, solution };
}

export function countSolutions(grid, n, limit = 2) {
  return solve(grid, n, limit).count;
}

export function generate(option, rng) {
  const n = Number(option) || 8;
  const solution = solve(new Array(n * n).fill(-1), n, 1, rng).solution;
  const givens = solution.slice();
  for (const i of shuffle(range(n * n), rng)) {
    const keep = givens[i];
    givens[i] = -1;
    if (countSolutions(givens, n) !== 1) givens[i] = keep;
  }
  return { size: n, givens, solution };
}

// Cells involved in a visible rule violation.
export function conflicts(values, n) {
  const bad = new Set();
  const half = n / 2;
  for (const isRow of [true, false]) {
    const lines = range(n).map((k) => lineCells(n, isRow, k));
    lines.forEach((cells, k) => {
      for (let p = 2; p < n; p++) {
        const v = values[cells[p]];
        if (v !== -1 && v === values[cells[p - 1]] && v === values[cells[p - 2]]) {
          bad.add(cells[p]); bad.add(cells[p - 1]); bad.add(cells[p - 2]);
        }
      }
      for (const sym of [0, 1]) {
        if (cells.filter((i) => values[i] === sym).length > half) {
          cells.forEach((i) => { if (values[i] === sym) bad.add(i); });
        }
      }
      if (cells.every((i) => values[i] !== -1)) {
        for (let o = k + 1; o < n; o++) {
          const oc = lines[o];
          if (oc.every((i, p) => values[i] === values[cells[p]])) {
            cells.forEach((i) => bad.add(i)); oc.forEach((i) => bad.add(i));
          }
        }
      }
    });
  }
  return bad;
}
