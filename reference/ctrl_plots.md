# Plot methods

All plots are built with 'ggplot2' and returned as `ggplot` objects, so
they can be modified or saved and work in noninteractive sessions.
Discrete time signals are drawn as steps (sample and hold); continuous
signals as lines. Axis labels state the time domain and units (rad per
time unit for frequencies).

## Usage

``` r
# S3 method for class 'ctrl_kalman_filter'
plot(x, truth = NULL, ...)

# S3 method for class 'ctrl_idcompare'
plot(x, ...)

# S3 method for class 'ctrl_resid'
plot(x, ...)

# S3 method for class 'ctrl_mpc_simulation'
plot(x, ...)

# S3 method for class 'ctrl_nl_response'
plot(x, ...)

# S3 method for class 'ctrl_pid_simulation'
plot(x, ...)

# S3 method for class 'ctrl_response'
plot(x, what = c("outputs", "states", "inputs", "all"), ...)

# S3 method for class 'ctrl_bode'
plot(x, ...)

# S3 method for class 'ctrl_nyquist'
plot(x, ...)

# S3 method for class 'ctrl_nichols'
plot(x, ...)

# S3 method for class 'ctrl_pzmap'
plot(x, ...)

# S3 method for class 'ctrl_rlocus'
plot(x, ...)
```

## Arguments

- x:

  An object returned by a `ctrl_*` analysis function.

- truth:

  Optional matrix of true states to overlay.

- ...:

  Unused.

- what:

  For responses: `"outputs"`, `"states"`, `"inputs"` or `"all"`.

## Value

A `ggplot` object.

## Examples

``` r
G <- ctrl_tf(4, c(1, 1.2, 4))
plot(ctrl_step(G))

plot(ctrl_bode(G))

plot(ctrl_nyquist(G))

plot(ctrl_pzmap(G))

plot(ctrl_rlocus(ctrl_tf(1, c(1, 3, 2, 0))))
```
