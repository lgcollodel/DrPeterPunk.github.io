import { conflicts, LT } from '../core/futoshiki.js';
import { el, History, numberPad, moveSelection, flash } from './common.js';

export function mount(root, puzzle, { onChange, onWin }) {
  const { size: n, givens, solution, h, v } = puzzle;
  const history = new History();
  let values = givens.slice();
  let selected = givens.indexOf(0);
  let done = false;

  // The board interleaves cells with sign slots: (2n-1) × (2n-1) tracks.
  const board = el('div', 'board futoshiki');
  board.style.setProperty('--n', n);
  const tracks = `repeat(${n - 1}, var(--track) calc(var(--track) * 0.4)) var(--track)`;
  board.style.gridTemplateColumns = tracks;
  board.style.gridTemplateRows = tracks;
  const cells = [];
  for (let R = 0; R < 2 * n - 1; R++) {
    for (let C = 0; C < 2 * n - 1; C++) {
      const r = R >> 1, c = C >> 1;
      if (R % 2 === 0 && C % 2 === 0) {
        const i = r * n + c;
        const cell = el('button', 'cell');
        cell.type = 'button';
        if (givens[i]) cell.classList.add('given');
        cell.addEventListener('click', () => { selected = i; render(); });
        cells[i] = cell;
        board.append(cell);
      } else if (R % 2 === 0) {
        const s = h[r * (n - 1) + c];
        board.append(el('span', 'sign', s ? (s === LT ? '<' : '>') : ''));
      } else if (C % 2 === 0) {
        const s = v[r * n + c];
        board.append(el('span', 'sign', s ? (s === LT ? '∧' : '∨') : ''));
      } else {
        board.append(el('span', 'gap'));
      }
    }
  }
  const pad = numberPad(n, (d) => enter(d), [{ label: '⌫', title: 'Erase (Backspace)', onClick: () => enter(0) }]);
  root.append(board, pad);

  function enter(d) {
    if (done || selected < 0 || givens[selected]) return;
    history.push(values);
    values[selected] = values[selected] === d ? 0 : d;
    changed();
  }

  function changed() {
    render();
    onChange();
    if (values.every((x, i) => x === solution[i])) {
      done = true;
      render();
      onWin();
    }
  }

  function render() {
    const bad = conflicts(values, puzzle);
    cells.forEach((cell, i) => {
      cell.textContent = values[i] ? String(values[i]) : '';
      cell.classList.toggle('selected', i === selected);
      cell.classList.toggle('conflict', bad.has(i) && !givens[i]);
      cell.classList.toggle('solved', done);
    });
  }

  function onKey(e) {
    if (e.target.closest && e.target.closest('select, input')) return;
    const d = Number(e.key);
    if (d >= 1 && d <= n) { enter(d); e.preventDefault(); return; }
    if (e.key === 'Backspace' || e.key === 'Delete' || e.key === '0') { enter(0); e.preventDefault(); return; }
    const next = moveSelection(e.key, Math.max(selected, 0), n);
    if (next !== null) { selected = next; render(); e.preventDefault(); }
  }
  document.addEventListener('keydown', onKey);
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
      const wrong = (i) => values[i] !== solution[i];
      const i = selected >= 0 && wrong(selected) ? selected : values.findIndex((_, k) => wrong(k));
      if (i < 0) return;
      history.push(values);
      values[i] = solution[i];
      selected = i;
      flash([cells[i]], 'hinted');
      changed();
    },
    check() {
      const wrong = cells.filter((_, i) => values[i] && values[i] !== solution[i]);
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
    destroy() { document.removeEventListener('keydown', onKey); },
  };
}
