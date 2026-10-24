// Sudoku: fill a 9x9 grid so every row, column and 3x3 box holds 1-9 once.
// Grids are flat arrays of 81 numbers, 0 meaning empty.

import { shuffle, range, popcount } from './rng.js';

const N = 9;
const ALL = 0x1ff;
const boxOf = range(81).map((i) => Math.floor(Math.floor(i / 9) / 3) * 3 + Math.floor((i % 9) / 3));

export const meta = {
  id: 'sudoku',
  name: 'Sudoku',
  tagline: 'Fill 1–9 in every row, column and box.',
  rules: [
    'Fill every empty cell with a digit from 1 to 9.',
    'Each row, each column and each 3×3 box contains every digit exactly once.',
  ],
  options: [
    { id: 'easy', label: 'Easy' },
    { id: 'medium', label: 'Medium' },
    { id: 'hard', label: 'Hard' },
  ],
  defaultOption: 'medium',
};

const TARGET_CLUES = { easy: 38, medium: 30, hard: 0 };

// Counts solutions up to `limit`. With an rng the candidate order is randomised,
// which is how full solution grids are generated.
export function solve(grid, limit = 2, rng = null) {
  const g = grid.slice();
  const rows = new Int32Array(N), cols = new Int32Array(N), boxes = new Int32Array(N);
  for (let i = 0; i < 81; i++) {
    const v = g[i];
    if (!v) continue;
    const bit = 1 << (v - 1), r = Math.floor(i / 9), c = i % 9, b = boxOf[i];
    if ((rows[r] | cols[c] | boxes[b]) & bit) return { count: 0, solution: null };
    rows[r] |= bit; cols[c] |= bit; boxes[b] |= bit;
  }
  let count = 0;
  let solution = null;

  function rec() {
    let best = -1, bestMask = 0, bestCount = 10;
    for (let i = 0; i < 81; i++) {
      if (g[i]) continue;
      const mask = ~(rows[Math.floor(i / 9)] | cols[i % 9] | boxes[boxOf[i]]) & ALL;
      const c = popcount(mask);
      if (c === 0) return;
      if (c < bestCount) {
        best = i; bestMask = mask; bestCount = c;
        if (c === 1) break;
      }
    }
    if (best === -1) {
      count++;
      if (!solution) solution = g.slice();
      return;
    }
    const digits = [];
    for (let d = 0; d < 9; d++) if (bestMask & (1 << d)) digits.push(d);
    if (rng) shuffle(digits, rng);
    const r = Math.floor(best / 9), c = best % 9, b = boxOf[best];
    for (const d of digits) {
      const bit = 1 << d;
      g[best] = d + 1;
      rows[r] |= bit; cols[c] |= bit; boxes[b] |= bit;
      rec();
      rows[r] &= ~bit; cols[c] &= ~bit; boxes[b] &= ~bit;
      g[best] = 0;
      if (count >= limit) return;
    }
  }

  rec();
  return { count, solution };
}

export function countSolutions(grid, limit = 2) {
  return solve(grid, limit).count;
}

export function generate(option, rng) {
  const solution = solve(new Array(81).fill(0), 1, rng).solution;
  const puzzle = solution.slice();

  // Remove clues in symmetric pairs as long as the solution stays unique.
  const order = shuffle(range(41), rng);
  for (const i of order) {
    const j = 80 - i;
    const keepI = puzzle[i], keepJ = puzzle[j];
    puzzle[i] = 0; puzzle[j] = 0;
    if (countSolutions(puzzle) !== 1) {
      puzzle[i] = keepI; puzzle[j] = keepJ;
    }
  }

  // Easier levels add clues back; adding clues never breaks uniqueness.
  const target = TARGET_CLUES[option] ?? TARGET_CLUES.medium;
  const empty = shuffle(range(41).filter((i) => !puzzle[i]), rng);
  let clues = puzzle.filter(Boolean).length;
  for (const i of empty) {
    if (clues >= target) break;
    puzzle[i] = solution[i];
    puzzle[80 - i] = solution[80 - i];
    clues += i === 40 ? 1 : 2;
  }

  return { size: N, givens: puzzle, solution };
}

// Indices of filled cells that clash with another cell in their row, column or box.
export function conflicts(values) {
  const bad = new Set();
  for (let i = 0; i < 81; i++) {
    if (!values[i]) continue;
    for (let j = i + 1; j < 81; j++) {
      if (values[i] !== values[j]) continue;
      if (Math.floor(i / 9) === Math.floor(j / 9) || i % 9 === j % 9 || boxOf[i] === boxOf[j]) {
        bad.add(i); bad.add(j);
      }
    }
  }
  return bad;
}

export { boxOf };
