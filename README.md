# CaRaSh spherical shell experiment

Local implementation of the Case/Rajan/Shende idea: a 3-D cubic lattice in
which a wavefront expands from a point source at (approximately) one cell per
tick, producing integer-radius spherical shells of thickness 1.

## What is tested

- `tick` equals the integer Euclidean radius `r = floor(sqrt(dx^2+dy^2+dz^2))`.
- No `sqrt`, `isqrt` or general multiplication in the propagation rule.
- Only additions and shifts are used: `r2_new = r2 + 2*|coord| + 1` and
  `next_sq = next_sq + 2*r + 1`.
- Because cells on the same shell can activate each other tangentially,
  each main tick is iterated to a fixed point (virtual sub-layers).

## Build and run

```bash
make          # L=17 by default
make L=31     # larger lattice
```

## Output

- Number of active cells per tick.
- Cross-section showing the concentric integer-radius shells.
- Verification count against a `sqrt`-based reference.

## Notes

- This is a *forward-only* expansion. Pulsation (expansion + contraction) is
  not implemented here; it is the next step for the `automaton` model.
- The algorithm uses each cell's signed axial offset from the source
  (`ax, ay, az`) to decide the outward directions. In a full CA these would
  be propagated along with the wave, just like the radius.
