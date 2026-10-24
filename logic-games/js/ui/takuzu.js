import { conflicts } from '../core/takuzu.js';
import { el, History, flash } from './common.js';

const SYMBOLS = { 0: '☀', 1: '☾' };

export function mount(root, puzzle, { onChange, onWin }) {
  const { size: n, givens, solution } = puzzle;
  const history = new History();
  let values = givens.slice();
  let done = false;

  const board = el('div', 'board takuzu');
  board.style.setProperty('--n', n);
  const cells = givens.map((g, i) => {
    const cell = el('button', 'cell');
    cell.type = 'button';
    if (g !== -1) cell.classList.add('given');
    cell.addEventListener('click', () => cycle(i, [0, 1, -1]));
    cell.addEventListener('contextmenu', (e) => { e.preventDefault(); cycle(i, [1, -1, 0]); });
    board.append(cell);
    return cell;
  });
  root.append(board, el('p', 'hint-text', 'Tap: ☀ → ☾ → empty. Right-click goes backwards.'));

  // next[v + 1] is the value that follows v.
  function cycle(i, next) {
    if (done || givens[i] !== -1) return;
    history.push(values);
    values[i] = next[values[i] + 1];
    changed();
  }

  function changed() {
    render();
    onChange();
    if (values.every((v, i) => v === solution[i])) {
      done = true;
      render();
      onWin();
    }
  }

  function render() {
    const bad = conflicts(values, n);
    cells.forEach((cell, i) => {
      const v = values[i];
      cell.textContent = v === -1 ? '' : SYMBOLS[v];
      cell.classList.toggle('sun', v === 0);
      cell.classList.toggle('moon', v === 1);
      cell.classList.toggle('conflict', bad.has(i));
      cell.classList.toggle('solved', done);
    });
  }
  render();

  return {
    undo() {
      const s = history.pop();
      if (!s || done) return;
      values = s;
      render(); onChange();
    },
    hint() {
      if (done) return;
      const i = values.findIndex((v, k) => v !== solution[k]);
      if (i < 0) return;
      history.push(values);
      values[i] = solution[i];
      flash([cells[i]], 'hinted');
      changed();
    },
    check() {
      const wrong = cells.filter((_, i) => values[i] !== -1 && values[i] !== solution[i]);
      flash(wrong, 'wrong');
      return wrong.length;
    },
    reset() {
      history.push(values);
      values = givens.slice(); done = false;
      render(); onChange();
    },
    reveal() {
      values = solution.slice(); done = true;
      render();
    },
    destroy() {},
  };
}
