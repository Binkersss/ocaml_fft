# OCaml Array-Based FFT

> This implementation is based on `njwilkins18/nas-mpl` Parallel ML FT benchmark
> You can find this implementation on GitHub

- Grid data: Two `float array`s per grid (`re` and `im` in the `c3` record)
- Scratch blocks and tables: The FFT work buffers (`y0r`, `y0i`, `y1r`, `y1i`),
  roots-of-unity tabls (`ur`, `ui`), `ex` component table, and RNG buffer `tmp`
  are all plain float arrays
- Index map: `int array`

This means that each array is standard, not allocating a heap block for each
element.

`sums` is `(float * float) array`
`ref_values` is `float list`, converted to arrays in `verify`.

## To run

Assuming you have dune, OCaml, etc.

```
dune build --release ./bin/ft.exe
./_build/default/bin/ft.exe
```
