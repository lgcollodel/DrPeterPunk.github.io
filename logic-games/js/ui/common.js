// Small DOM helpers shared by the game views.

export function el(tag, className = '', text = '') {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text) node.textContent = text;
  return node;
}

// Undo stack of snapshots taken before each move.
export class History {
  constructor() { this.stack = []; }
  push(snapshot) { this.stack.push(structuredClone(snapshot)); }
  pop() { return this.stack.pop(); }
  clear() { this.stack = []; }
}

export function numberPad(count, onDigit, extras = []) {
  const pad = el('div', 'numpad');
  for (let d = 1; d <= count; d++) {
    const b = el('button', 'numkey', String(d));
    b.type = 'button';
    b.addEventListener('click', () => onDigit(d));
    pad.append(b);
  }
  for (const { label, title, onClick, className } of extras) {
    const b = el('button', `numkey ${className || ''}`, label);
    b.type = 'button';
    if (title) b.title = title;
    b.addEventListener('click', () => onClick(b));
    pad.append(b);
  }
  return pad;
}

// Arrow-key movement inside an n×n grid; returns the new index or null.
export function moveSelection(key, index, n) {
  const r = Math.floor(index / n), c = index % n;
  switch (key) {
    case 'ArrowUp': return ((r + n - 1) % n) * n + c;
    case 'ArrowDown': return ((r + 1) % n) * n + c;
    case 'ArrowLeft': return r * n + ((c + n - 1) % n);
    case 'ArrowRight': return r * n + ((c + 1) % n);
    default: return null;
  }
}

export function flash(cells, className, ms = 1500) {
  cells.forEach((c) => c.classList.add(className));
  setTimeout(() => cells.forEach((c) => c.classList.remove(className)), ms);
}
