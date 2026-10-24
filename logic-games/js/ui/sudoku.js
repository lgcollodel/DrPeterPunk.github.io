import { conflicts, boxOf } from '../core/sudoku.js';
import { el, History, numberPad, moveSelection, flash } from './common.js';

export function mount(root, puzzle, { onChange, onWin }) {
  const { givens, solution } = puzzle;
  const history = new History();
  let values = givens.slice();
  let notes = new Array(81).fill(0);
  let selected = givens.indexOf(0);
  let notesMode = false;
  let done = false;

  const board = el('div', 'board sudoku');
  const cells = givens.map((g, i) => {
    const cell = el('button', 'cell');
    cell.type = 'button';
    const r = Math.floor(i / 9), c = i % 9;
    if (c % 3 === 2 && c < 8) cell.classList.add('edge-r');
    if (r % 3 === 2 && r < 8) cell.classList.add('edge-b');
    if (g) cell.classList.add('given');
    cell.addEventListener('click', () => { selected = i; render(); });
    board.append(cell);
    return cell;
  });

  const notesBtn = { label: '✎', title: 'Pencil marks (N)', className: 'toggle', onClick: (b) => {
    notesMode = !notesMode;
    b.classList.toggle('on', notesMode);
  } };
  const pad = numberPad(9, (d) => enter(d), [
    { label: '⌫', title: 'Erase (Backspace)', onClick: () => enter(0) },
    notesBtn,
  ]);
  root.append(board, pad);

  function snapshot() { return { values, notes }; }

  function enter(d) {
    if (done || selected < 0 || givens[selected]) return;
    history.push(snapshot());
    if (notesMode && d) {
      values[selected] = 0;
      notes[selected] ^= 1 << (d - 1);
    } else {
      values[selected] = values[selected] === d ? 0 : d;
      notes[selected] = 0;
      if (d) clearNotes(selected, d);
    }
    changed();
  }

  // Placing a digit removes it from pencil marks of the cell's peers.
  function clearNotes(i, d) {
    const r = Math.floor(i / 9), c = i % 9;
    for (let j = 0; j < 81; j++) {
      if (Math.floor(j / 9) === r || j % 9 === c || boxOf[j] === boxOf[i]) notes[j] &= ~(1 << (d - 1));
    }
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
    const bad = conflicts(values);
    const selVal = selected >= 0 ? values[selected] : 0;
    cells.forEach((cell, i) => {
      const v = values[i];
      cell.classList.toggle('selected', i === selected);
      cell.classList.toggle('peer', selected >= 0 && i !== selected && (
        Math.floor(i / 9) === Math.floor(selected / 9) || i % 9 === selected % 9 || boxOf[i] === boxOf[selected]));
      cell.classList.toggle('same', !!v && v === selVal && i !== selected);
      cell.classList.toggle('conflict', bad.has(i) && !givens[i]);
      cell.classList.toggle('solved', done);
      if (v) {
        cell.textContent = String(v);
      } else if (notes[i]) {
        cell.textContent = '';
        const grid = el('span', 'notes');
        for (let d = 1; d <= 9; d++) grid.append(el('span', '', notes[i] & (1 << (d - 1)) ? String(d) : ''));
        cell.append(grid);
      } else {
        cell.textContent = '';
      }
    });
  }

  function onKey(e) {
    if (e.target.closest && e.target.closest('select, input')) return;
    if (/^[1-9]$/.test(e.key)) { enter(Number(e.key)); e.preventDefault(); return; }
    if (e.key === 'Backspace' || e.key === 'Delete' || e.key === '0') { enter(0); e.preventDefault(); return; }
    if (e.key === 'n' || e.key === 'N') { pad.querySelector('.toggle').click(); return; }
    const next = moveSelection(e.key, Math.max(selected, 0), 9);
    if (next !== null) { selected = next; render(); e.preventDefault(); }
  }
  document.addEventListener('keydown', onKey);
  render();

  return {
    undo() {
      const s = history.pop();
      if (!s || done) return;
      ({ values, notes } = s);
      render(); onChange();
    },
    hint() {
      if (done) return;
      const wrong = (i) => !givens[i] && values[i] !== solution[i];
      const target = selected >= 0 && wrong(selected) ? selected : values.findIndex((_, i) => wrong(i));
      if (target < 0) return;
      history.push(snapshot());
      values[target] = solution[target];
      notes[target] = 0;
      clearNotes(target, solution[target]);
      selected = target;
      flash([cells[target]], 'hinted');
      changed();
    },
    check() {
      const wrong = cells.filter((_, i) => values[i] && values[i] !== solution[i]);
      flash(wrong, 'wrong');
      return wrong.length;
    },
    reset() {
      history.push(snapshot());
      values = givens.slice(); notes = new Array(81).fill(0); done = false;
      render(); onChange();
    },
    reveal() {
      values = solution.slice(); notes = new Array(81).fill(0); done = true;
      render();
    },
    destroy() { document.removeEventListener('keydown', onKey); },
  };
}
