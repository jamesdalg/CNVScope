## Submission

This fixes the examples ERROR reported for 3.7.7 on
r-devel-linux-x86_64-debian-gcc and r-devel-linux-x86_64-fedora-gcc:

```
> rebinGenomicInteractions(gint=gint_small_chr1, ...)
Error in as.vector(x) : no method for coercing this S4 class to a vector
Calls: rebinGenomicInteractions ... as.data.frame.Vector -> .as.data.frame.Vector -> as.vector
```

S4Vectors 0.51 (Bioconductor devel) changed the S3 method
`as.data.frame.Vector()` so that it no longer re-dispatches to the S4
`as.data.frame()` generic. `rebinGenomicInteractions()` called base
`as.data.frame()` on a GenomicInteractions object, which therefore no longer
reached the GInteractions method from InteractionSet. The function now calls
the S4 generic, `S4Vectors::as.data.frame()`, directly. S4Vectors was already
in Suggests, so the dependencies are unchanged, and the function's output is
identical to before.

## Test environments

* Ubuntu Linux, R-devel with Bioconductor devel (S4Vectors 0.51.10,
  InteractionSet 1.41.0, GenomicRanges 1.65.4), R-release, R-oldrel
  (GitHub Actions)
* macOS and Windows, R-release (GitHub Actions)
* Local Ubuntu Linux, R 4.6.2, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 0 notes on all GitHub Actions environments.

The local `--as-cran` check gave 1 NOTE:

```
Examples with CPU (user + system) or elapsed time > 5s
                     user system elapsed
importBreakpointBed 6.621  0.538   7.161
```

This time is taken up by loading the GenomicInteractions and InteractionSet
namespaces (Suggests), which this is the first example to use. The example
itself reads a six-line BED file. The same timing appeared in the checks of
3.7.7.
