import { conflicts } from '../core/queens.js';
import { el, History, flash } from './common.js';

const EMPTY = 0, MARK = 1, QUEEN = 2;

export function mount(root, puzzle, { onChange, onWin }) {
  const { size: n, regions, solution } = puzzle;
  const target = new Set(solution.map((c, r) => r * n + c));
  const history = new History();
  let state = new Array(n * n).fill(EMPTY);
  let done = false;

  const board = el('div', 'board queens');
  board.style.setProperty('--n', n);
  const cells = regions.map((reg, i) => {
    const cell = el('button', `cell region-${reg}`);
    cell.type = 'button';
    cell.setAttribute('aria-label', `row ${Math.floor(i / n) + 1} column ${(i % n) + 1}`);
    const r = Math.floor(i / n), c = i % n;
    if (c < n - 1 && regions[i + 1] !== reg) cell.classList.add('wall-r');
    if (r < n - 1 && regions[i + n] !== reg) cell.classList.add('wall-b');
    cell.addEventListener('click', () => cycle(i, 1));
    cell.addEventListener('contextmenu', (e) => { e.preventDefault(); cycle(i, 2); });
    board.append(cell);
    return cell;
  });
  root.append(board, el('p', 'hint-text', 'Tap: ✕ → ♛ → empty. Right-click goes backwards.'));

  function cycle(i, step) {
    if (done) return;
    history.push(state);
    state[i] = (state[i] + step) % 3;
    changed();
  }

  function queenCells() {
    return state.map((s, i) => (s === QUEEN ? i : -1)).filter((i) => i >= 0);
  }

  function changed() {
    render();
    onChange();
    const q = queenCells();
    if (q.length === n && q.every((i) => target.has(i))) {
      done = true;
      render();
      onWin();
    }
  }

  function render() {
    const bad = conflicts(queenCells(), regions, n);
    cells.forEach((cell, i) => {
      cell.textContent = state[i] === QUEEN ? '♛' : state[i] === MARK ? '✕' : '';
      cell.classList.toggle('queen', state[i] === QUEEN);
      cell.classList.toggle('mark', state[i] === MARK);
      cell.classList.toggle('conflict', bad.has(i));
      cell.classList.toggle('solved', done);
    });
  }
  render();

  return {
    undo() {
      const s = history.pop();
      if (!s || done) return;
      state = s;
      render(); onChange();
    },
    hint() {
      if (done) return;
      history.push(state);
      // First remove a misplaced queen, otherwise place a missing one.
      const wrongQueen = queenCells().find((i) => !target.has(i));
      const missing = [...target].find((i) => state[i] !== QUEEN);
      const i = wrongQueen ?? missing;
      if (i === undefined) return;
      state[i] = wrongQueen !== undefined ? MARK : QUEEN;
      flash([cells[i]], 'hinted');
      changed();
    },
    check() {
      const wrong = cells.filter((_, i) => (state[i] === QUEEN && !target.has(i)) || (state[i] === MARK && target.has(i)));
      flash(wrong, 'wrong');
      return wrong.length;
    },
    reset() {
      history.push(state);
      state = new Array(n * n).fill(EMPTY); done = false;
      render(); onChange();
    },
    reveal() {
      state = state.map((_, i) => (target.has(i) ? QUEEN : EMPTY)); done = true;
      render();
    },
    destroy() {},
  };
}
