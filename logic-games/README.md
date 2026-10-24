# Logic Games

Queens, Sudoku, Takuzu and Futoshiki in the browser. Every puzzle is generated
on the fly from a seed and is guaranteed to have exactly one solution.

Open `index.html` through any static server (ES modules don't load from
`file://`), e.g. `python3 -m http.server` in the repo root, then visit
`/logic-games/`. On GitHub Pages it is served at `/logic-games/`.

## How generation works

Each game has a solver that counts solutions up to a limit of 2. The
generator only returns a puzzle when the count is exactly 1.

| Game | Generation |
| --- | --- |
| Sudoku | Random full grid, then remove clues in symmetric pairs while the solution stays unique. Easier levels add clues back. |
| Queens | Random non-touching queen placement, regions grown from each queen, then hill-climbing on region borders until no other solution exists. |
| Takuzu | Random valid grid built from line patterns, then remove cells while the solution stays unique. The solver propagates over line-pattern domains. |
| Futoshiki | Random Latin square with every number and sign, then remove numbers first and signs second while the solution stays unique. |

Puzzles are reproducible: the URL `#/<game>/<option>/<seed>` always gives the
same puzzle, so links can be shared.

## Code layout

- `js/core/` – generators and solvers, no DOM (`rng.js` is the seeded PRNG)
- `js/ui/` – one view per game
- `js/app.js` – routing, toolbar, timer, best times
- `test/` – checks that generated puzzles are valid and unique

Run the tests with `npm test` (Node 18+, no dependencies).
