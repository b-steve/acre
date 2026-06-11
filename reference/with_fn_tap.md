# Wrap an objective to emit buffered trace rows on each evaluation

Creates an idempotent wrapper around an objective function
`fn(par, ...)` that forwards the call unchanged, then notifies a tap
object produced by
[`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md)
via `tap$on_fn_evaluation(par_vec, fval, par_names)`. Use this with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (or similar
optimizers) to stream live, correctly paired trace rows when combined
with
[`with_gr_tap()`](https://b-steve.github.io/acre/reference/with_gr_tap.md).

## Usage

``` r
with_fn_tap(fn, tap, par_names, G_matrix_list)
```

## Arguments

- fn:

  Objective function with signature `function(par, ...) -> scalar`.

- tap:

  Tap/list returned by
  [`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md),
  expected to provide a method
  `on_fn_evaluation(par_vec, fval, par_names)`.

- par_names:

  Character vector of names corresponding to `par`. These are used by
  the tap to map positions in `par` onto the `trace_cols` it prints.

- G_matrix_list:

  List of back-transformation matrices used to back transform parameter
  values in the case they have been scaled and centered for model
  fitting.

## Value

A function with the same signature and return value as `fn`, which also
triggers the tap side-effect after each evaluation.

## Details

The wrapper is **idempotent**: if it detects the function was already
wrapped (attribute `.__fn_tapped`), it returns `fn` unchanged. Arguments
are [`force()`](https://rdrr.io/r/base/force.html)d so the wrapper
closes over the current objects.

## See also

[`with_gr_tap()`](https://b-steve.github.io/acre/reference/with_gr_tap.md),
[`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md),
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)
