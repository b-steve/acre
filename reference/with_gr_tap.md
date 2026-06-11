# Wrap a gradient to emit buffered trace rows (mgc) on each evaluation

Creates an idempotent wrapper around a gradient function `gr(par, ...)`
that forwards the call unchanged, then notifies a tap object produced by
[`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md)
via `tap$on_gr_evaluation(par_vec, grad_vec, par_names)`. When used
together with
[`with_fn_tap()`](https://b-steve.github.io/acre/reference/with_fn_tap.md),
the tap can pair `fn` and `gr` evaluations at the same parameter vector
and print a trace row that includes `mgc = max(abs(grad))`.

## Usage

``` r
with_gr_tap(gr, tap, par_names, G_matrix_list)
```

## Arguments

- gr:

  Gradient function with signature `function(par, ...) -> numeric` (same
  length/order as `par`), or `NULL` if no gradient is supplied.

- tap:

  Tap/list returned by
  [`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md),
  expected to provide a method
  `on_gr_evaluation(par_vec, grad_vec, par_names)`.

- par_names:

  Character vector of names corresponding to `par`. These are used by
  the tap to map gradient components to parameter columns.

- G_matrix_list:

  List of back-transformation matrices used to back transform parameter
  values in the case they have been scaled and centered for model
  fitting.

## Value

A function with the same signature and return value as `gr`, which also
triggers the tap side-effect after each evaluation; or `NULL` if `gr`
was `NULL`.

## Details

If `gr` is `NULL`, this helper returns `NULL` (no wrapping). The wrapper
is **idempotent**: if it detects the function was already wrapped
(attribute `.__gr_tapped`), it returns `gr` unchanged. Arguments are
[`force()`](https://rdrr.io/r/base/force.html)d so the wrapper closes
over the current objects.

## See also

[`with_fn_tap()`](https://b-steve.github.io/acre/reference/with_fn_tap.md),
[`make_buffer_printer()`](https://b-steve.github.io/acre/reference/make_buffer_printer.md),
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)
