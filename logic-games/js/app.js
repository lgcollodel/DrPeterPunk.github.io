import { makeRng, randomSeed } from './core/rng.js';
import * as sudoku from './core/sudoku.js';
import * as queens from './core/queens.js';
import * as takuzu from './core/takuzu.js';
import * as futoshiki from './core/futoshiki.js';
import * as sudokuView from './ui/sudoku.js';
import * as queensView from './ui/queens.js';
import * as takuzuView from './ui/takuzu.js';
import * as futoshikiView from './ui/futoshiki.js';
import { el } from './ui/common.js';

const GAMES = [
  { core: queens, view: queensView, icon: '♛' },
  { core: sudoku, view: sudokuView, icon: '9' },
  { core: takuzu, view: takuzuView, icon: '☀' },
  { core: futoshiki, view: futoshikiView, icon: '>' },
];
const byId = Object.fromEntries(GAMES.map((g) => [g.core.meta.id, g]));

const app = document.getElementById('app');
let current = null;
let timer = null;
let renderToken = 0;

function storage(key, value) {
  try {
    if (value === undefined) return JSON.parse(localStorage.getItem(key));
    localStorage.setItem(key, JSON.stringify(value));
  } catch {
    return null;
  }
  return null;
}

function formatTime(ms) {
  const s = Math.floor(ms / 1000);
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
}

function teardown() {
  renderToken++;
  if (current) current.destroy();
  current = null;
  clearInterval(timer);
}

function renderHome() {
  teardown();
  document.title = 'Logic Games';
  app.replaceChildren();
  const intro = el('p', 'intro',
    'Every puzzle is generated on the spot and checked by a solver to have exactly one solution.');
  const grid = el('div', 'cards');
  for (const { core: { meta }, icon } of GAMES) {
    const card = el('a', 'card');
    card.href = `#/${meta.id}/${storage(`lg-last-${meta.id}`) || meta.defaultOption}`;
    card.append(el('span', `card-icon icon-${meta.id}`, icon), el('h2', '', meta.name), el('p', '', meta.tagline));
    grid.append(card);
  }
  app.append(intro, grid);
}

function renderGame(game, option, seed) {
  teardown();
  const { meta } = game.core;
  document.title = `${meta.name} · Logic Games`;
  storage(`lg-last-${meta.id}`, option);

  const toolbar = el('div', 'toolbar');
  const back = el('a', 'back', '← Games');
  back.href = '#/';
  const title = el('h2', 'game-title', meta.name);
  const select = el('select');
  select.setAttribute('aria-label', 'Difficulty');
  for (const o of meta.options) {
    const opt = el('option', '', o.label);
    opt.value = o.id;
    opt.selected = o.id === option;
    select.append(opt);
  }
  select.addEventListener('change', () => { location.hash = `#/${meta.id}/${select.value}`; });
  const newBtn = el('button', 'primary', 'New puzzle');
  newBtn.type = 'button';
  newBtn.addEventListener('click', () => { location.hash = `#/${meta.id}/${option}/${randomSeed()}`; });
  toolbar.append(back, title, select, newBtn);

  const status = el('div', 'status');
  const clock = el('span', 'clock', '0:00');
  const best = el('span', 'best');
  const message = el('span', 'message');
  status.append(clock, message, best);

  const stage = el('div', 'stage');
  const loading = el('p', 'loading', 'Generating a puzzle with a unique solution…');
  stage.append(loading);

  const actions = el('div', 'actions');
  const rules = el('details', 'rules');
  rules.append(el('summary', '', 'How to play'));
  const list = el('ul');
  meta.rules.forEach((r) => list.append(el('li', '', r)));
  rules.append(list);
  const seedInfo = el('p', 'seed');

  app.replaceChildren(toolbar, status, stage, actions, rules, seedInfo);

  const bestKey = `lg-best-${meta.id}-${option}`;
  const showBest = () => {
    const b = storage(bestKey);
    best.textContent = b ? `Best ${formatTime(b)}` : '';
  };
  showBest();

  // Let the "Generating" message paint before the synchronous work starts.
  const token = renderToken;
  setTimeout(() => {
    if (token !== renderToken) return;
    const puzzle = game.core.generate(option, makeRng(seed));
    stage.replaceChildren();
    const started = Date.now();
    let finished = false;
    let assisted = false;
    const tick = () => { clock.textContent = formatTime(Date.now() - started); };
    timer = setInterval(tick, 1000);

    const say = (text, cls = '') => { message.textContent = text; message.className = `message ${cls}`; };
    const stop = () => { finished = true; clearInterval(timer); tick(); };

    current = game.view.mount(stage, puzzle, {
      onChange: () => { if (!finished) say(''); },
      onWin: () => {
        const elapsed = Date.now() - started;
        stop();
        if (assisted) {
          say(`Solved with hints in ${formatTime(elapsed)}`, 'win');
          return;
        }
        const prev = storage(bestKey);
        if (!prev || elapsed < prev) storage(bestKey, elapsed);
        showBest();
        say(`Solved in ${formatTime(elapsed)}! 🎉`, 'win');
      },
    });

    const action = (label, fn, cls = '') => {
      const b = el('button', cls, label);
      b.type = 'button';
      b.addEventListener('click', fn);
      actions.append(b);
    };
    action('Undo', () => current.undo());
    action('Hint', () => { if (!finished) assisted = true; current.hint(); });
    action('Check', () => {
      if (finished) return;
      const wrong = current.check();
      say(wrong ? `${wrong} mistake${wrong > 1 ? 's' : ''} highlighted` : 'No mistakes so far', wrong ? 'bad' : 'good');
    });
    action('Restart', () => current.reset());
    action('Show solution', () => {
      if (finished || !confirm('Reveal the solution?')) return;
      stop();
      current.reveal();
      say('Solution revealed');
    }, 'subtle');

    seedInfo.textContent = `Puzzle #${seed} · share this page's link to share the puzzle`;
  }, 30);
}

function route() {
  const [, id, option, seedText] = location.hash.split('/');
  const game = byId[id];
  if (!game) return renderHome();
  const { meta } = game.core;
  const opt = meta.options.some((o) => o.id === option) ? option : meta.defaultOption;
  const seed = Number(seedText);
  if (!Number.isInteger(seed) || seed <= 0 || opt !== option) {
    history.replaceState(null, '', `#/${id}/${opt}/${randomSeed()}`);
    return route();
  }
  renderGame(game, opt, seed);
}

window.addEventListener('hashchange', route);
route();
