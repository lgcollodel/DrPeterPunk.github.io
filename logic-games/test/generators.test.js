// Every generated puzzle must have exactly one solution, and that solution must
// obey the rules. Rule checks here are written independently of the solvers.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { makeRng } from '../js/core/rng.js';
import * as sudoku from '../js/core/sudoku.js';
import * as queens from '../js/core/queens.js';
import * as takuzu from '../js/core/takuzu.js';
import * as futoshiki from '../js/core/futoshiki.js';

const SEEDS = [1, 2, 3, 42, 1234];

function isPermutation(values, n) {
  return [...values].sort((a, b) => a - b).every((v, i) => v === i + 1) && values.length === n;
}

test('sudoku puzzles are valid and unique', () => {
  for (const option of ['easy', 'medium', 'hard']) {
    for (const seed of SEEDS) {
      const p = sudoku.generate(option, makeRng(seed));
      const s = p.solution;
      for (let k = 0; k < 9; k++) {
        assert.ok(isPermutation(s.slice(k * 9, k * 9 + 9), 9));
        assert.ok(isPermutation(s.filter((_, i) => i % 9 === k), 9));
        assert.ok(isPermutation(s.filter((_, i) => sudoku.boxOf[i] === k), 9));
      }
      p.givens.forEach((v, i) => v && assert.equal(v, s[i]));
      assert.equal(sudoku.countSolutions(p.givens, 3), 1);
    }
  }
});

test('sudoku difficulty changes clue count', () => {
  const clues = (o) => sudoku.generate(o, makeRng(7)).givens.filter(Boolean).length;
  assert.ok(clues('easy') > clues('medium'));
  assert.ok(clues('medium') >= clues('hard'));
});

// Brute force over all column permutations, independent of the solver.
function bruteQueens(regions, n) {
  let count = 0;
  const cols = [];
  const used = new Array(n).fill(false);
  (function rec(r) {
    if (r === n) {
      const regs = new Set(cols.map((c, row) => regions[row * n + c]));
      const touching = cols.some((c, row) => row > 0 && Math.abs(c - cols[row - 1]) <= 1);
      if (regs.size === n && !touching) count++;
      return;
    }
    for (let c = 0; c < n; c++) {
      if (used[c]) continue;
      used[c] = true; cols.push(c); rec(r + 1); cols.pop(); used[c] = false;
    }
  })(0);
  return count;
}

function regionsConnected(regions, n) {
  for (let g = 0; g < n; g++) {
    const cells = regions.map((r, i) => (r === g ? i : -1)).filter((i) => i >= 0);
    if (!cells.length) return false;
    const seen = new Set([cells[0]]);
    const stack = [cells[0]];
    while (stack.length) {
      const i = stack.pop();
      const r = Math.floor(i / n), c = i % n;
      for (const [dr, dc] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const rr = r + dr, cc = c + dc, j = rr * n + cc;
        if (rr >= 0 && rr < n && cc >= 0 && cc < n && regions[j] === g && !seen.has(j)) {
          seen.add(j); stack.push(j);
        }
      }
    }
    if (seen.size !== cells.length) return false;
  }
  return true;
}

test('queens puzzles are valid and unique', () => {
  for (const option of ['6', '7', '8']) {
    for (const seed of SEEDS) {
      const p = queens.generate(option, makeRng(seed));
      const n = p.size;
      assert.equal(bruteQueens(p.regions, n), 1, `n=${n} seed=${seed}`);
      assert.ok(regionsConnected(p.regions, n));
      assert.equal(new Set(p.solution.map((c, r) => p.regions[r * n + c])).size, n);
    }
  }
  for (const option of ['9', '10']) {
    const p = queens.generate(option, makeRng(5));
    assert.equal(queens.countSolutions(p.regions, p.size, 3), 1);
    assert.ok(regionsConnected(p.regions, p.size));
  }
});

function takuzuValid(s, n) {
  const lines = [];
  for (let k = 0; k < n; k++) {
    lines.push(s.slice(k * n, k * n + n).join(''));
    lines.push(s.filter((_, i) => i % n === k).join(''));
  }
  const rows = lines.filter((_, i) => i % 2 === 0), cols = lines.filter((_, i) => i % 2 === 1);
  return lines.every((l) => !/000|111/.test(l) && l.split('0').length - 1 === n / 2) &&
    new Set(rows).size === n && new Set(cols).size === n;
}

test('takuzu puzzles are valid and unique', () => {
  for (const option of ['6', '8', '10', '12']) {
    for (const seed of SEEDS.slice(0, option === '12' ? 2 : 5)) {
      const p = takuzu.generate(option, makeRng(seed));
      assert.ok(takuzuValid(p.solution, p.size));
      p.givens.forEach((v, i) => v !== -1 && assert.equal(v, p.solution[i]));
      assert.equal(takuzu.countSolutions(p.givens, p.size, 3), 1);
    }
  }
});

test('futoshiki puzzles are valid and unique', () => {
  for (const option of ['4', '5', '6', '7']) {
    for (const seed of SEEDS) {
      const p = futoshiki.generate(option, makeRng(seed));
      const n = p.size, s = p.solution;
      for (let k = 0; k < n; k++) {
        assert.ok(isPermutation(s.slice(k * n, k * n + n), n));
        assert.ok(isPermutation(s.filter((_, i) => i % n === k), n));
      }
      p.h.forEach((sign, k) => {
        const a = Math.floor(k / (n - 1)) * n + (k % (n - 1));
        if (sign) assert.equal(sign === futoshiki.LT, s[a] < s[a + 1]);
      });
      p.v.forEach((sign, k) => {
        if (sign) assert.equal(sign === futoshiki.LT, s[k] < s[k + n]);
      });
      assert.equal(futoshiki.countSolutions(p, 3), 1);
    }
  }
});

test('same seed gives the same puzzle', () => {
  assert.deepEqual(sudoku.generate('hard', makeRng(9)), sudoku.generate('hard', makeRng(9)));
  assert.deepEqual(queens.generate('8', makeRng(9)), queens.generate('8', makeRng(9)));
});
