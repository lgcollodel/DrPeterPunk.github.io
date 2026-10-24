// Queens: an N×N board split into N coloured regions. Place one queen in every
// row, column and region so that no two queens touch, not even diagonally.
// `regions` is a flat array of n*n region ids (0..n-1).

import { shuffle, range, randInt } from './rng.js';

export const meta = {
  id: 'queens',
  name: 'Queens',
  tagline: 'One queen per row, column and colour.',
  rules: [
    'Place exactly one queen in every row, every column and every coloured region.',
    'Queens may not touch each other, not even diagonally.',
    'Tap a cell once to mark it ✕, twice to place a queen.',
  ],
  options: [
    { id: '6', label: '6×6' },
    { id: '7', label: '7×7' },
    { id: '8', label: '8×8' },
    { id: '9', label: '9×9' },
    { id: '10', label: '10×10' },
  ],
  defaultOption: '8',
};

// Returns up to `limit` solutions, each an array col[row].
export function solveAll(regions, n, limit = 2) {
  const solutions = [];
  const cols = new Array(n).fill(-1);
  let usedCols = 0, usedRegions = 0;

  function rec(r) {
    if (r === n) {
      solutions.push(cols.slice());
      return;
    }
    const prev = r > 0 ? cols[r - 1] : -10;
    for (let c = 0; c < n; c++) {
      if (usedCols & (1 << c)) continue;
      if (Math.abs(c - prev) <= 1) continue;
      const reg = regions[r * n + c];
      if (usedRegions & (1 << reg)) continue;
      cols[r] = c;
      usedCols |= 1 << c; usedRegions |= 1 << reg;
      rec(r + 1);
      usedCols &= ~(1 << c); usedRegions &= ~(1 << reg);
      if (solutions.length >= limit) return;
    }
    cols[r] = -1;
  }

  rec(0);
  return solutions;
}

export function countSolutions(regions, n, limit = 2) {
  return solveAll(regions, n, limit).length;
}

function neighbours4(i, n) {
  const r = Math.floor(i / n), c = i % n, out = [];
  if (r > 0) out.push(i - n);
  if (r < n - 1) out.push(i + n);
  if (c > 0) out.push(i - 1);
  if (c < n - 1) out.push(i + 1);
  return out;
}

function randomQueens(n, rng) {
  const cols = [];
  const used = new Set();
  function rec(r) {
    if (r === n) return true;
    for (const c of shuffle(range(n), rng)) {
      if (used.has(c) || (r > 0 && Math.abs(c - cols[r - 1]) <= 1)) continue;
      cols.push(c); used.add(c);
      if (rec(r + 1)) return true;
      cols.pop(); used.delete(c);
    }
    return false;
  }
  rec(0);
  return cols;
}

// Grow one region from every queen by random flood fill.
function growRegions(queens, n, rng) {
  const regions = new Array(n * n).fill(-1);
  queens.forEach((c, r) => { regions[r * n + c] = r; });
  let left = n * n - n;
  while (left > 0) {
    const frontier = [];
    for (let i = 0; i < n * n; i++) {
      if (regions[i] !== -1) continue;
      for (const j of neighbours4(i, n)) if (regions[j] !== -1) frontier.push([i, regions[j]]);
    }
    // Half of the time favour the smallest growing region to keep sizes varied
    // but not lopsided.
    let pool = frontier;
    if (rng() < 0.5) {
      const sizes = new Array(n).fill(0);
      regions.forEach((g) => { if (g >= 0) sizes[g]++; });
      const min = Math.min(...frontier.map(([, g]) => sizes[g]));
      pool = frontier.filter(([, g]) => sizes[g] === min);
    }
    const [cell, reg] = pool[randInt(rng, pool.length)];
    regions[cell] = reg;
    left--;
  }
  return regions;
}

function stillConnectedWithout(regions, n, reg, removed) {
  const cells = [];
  for (let i = 0; i < n * n; i++) if (regions[i] === reg && i !== removed) cells.push(i);
  if (cells.length === 0) return false;
  const seen = new Set([cells[0]]);
  const stack = [cells[0]];
  while (stack.length) {
    const i = stack.pop();
    for (const j of neighbours4(i, n)) {
      if (j !== removed && regions[j] === reg && !seen.has(j)) {
        seen.add(j); stack.push(j);
      }
    }
  }
  return seen.size === cells.length;
}

// Relabel regions in reading order so ids carry no hint about the solution.
function normalise(regions) {
  const map = new Map();
  return regions.map((r) => {
    if (!map.has(r)) map.set(r, map.size);
    return map.get(r);
  });
}

export function generate(option, rng) {
  const n = Number(option) || 8;
  for (;;) {
    const queens = randomQueens(n, rng);
    const seeds = new Set(queens.map((c, r) => r * n + c));
    const regions = growRegions(queens, n, rng);

    // Hill-climb: while another solution exists, move one of its queen cells into
    // a neighbouring region. The intended queens stay alone in their regions,
    // so the intended solution always remains valid.
    for (let iter = 0; iter < 400; iter++) {
      const sols = solveAll(regions, n, 2);
      if (sols.length === 1) {
        return { size: n, regions: normalise(regions), solution: queens.slice() };
      }
      const alt = sols.find((s) => s.some((c, r) => c !== queens[r]));
      const candidates = shuffle(
        alt.map((c, r) => r * n + c).filter((i) => !seeds.has(i)),
        rng,
      );
      let moved = false;
      for (const cell of candidates) {
        const from = regions[cell];
        const targets = shuffle(
          [...new Set(neighbours4(cell, n).map((j) => regions[j]))].filter((g) => g !== from),
          rng,
        );
        if (!targets.length || !stillConnectedWithout(regions, n, from, cell)) continue;
        // Usually grow the smallest neighbouring region; sometimes pick at
        // random so the search does not get stuck cycling.
        const size = (g) => regions.filter((x) => x === g).length;
        regions[cell] = rng() < 0.7 ? targets.reduce((a, b) => (size(b) < size(a) ? b : a)) : targets[0];
        moved = true;
        break;
      }
      if (!moved) break;
    }
  }
}

// Indices of queens that break a rule given the current placements.
export function conflicts(queenCells, regions, n) {
  const bad = new Set();
  const list = [...queenCells];
  for (let a = 0; a < list.length; a++) {
    for (let b = a + 1; b < list.length; b++) {
      const i = list[a], j = list[b];
      const ri = Math.floor(i / n), ci = i % n, rj = Math.floor(j / n), cj = j % n;
      if (ri === rj || ci === cj || regions[i] === regions[j] ||
          (Math.abs(ri - rj) <= 1 && Math.abs(ci - cj) <= 1)) {
        bad.add(i); bad.add(j);
      }
    }
  }
  return bad;
}
