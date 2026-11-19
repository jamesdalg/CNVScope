## Resubmission

This is a resubmission of 3.7.9, whose incoming pretest gave one NOTE on
r-devel-linux-x86_64-debian-gcc:

```
Examples with CPU (user + system) or elapsed time > 5s
                     user system elapsed
importBreakpointBed 6.024  0.208   6.236
```

Most of that time was loading the GenomicInteractions namespace (Suggests),
which pulls in about 50 further Bioconductor namespaces and takes about 5
seconds on its own. The rest was a closeAllConnections() call inside
importBreakpointBed(), which should not have been there. In this version:

* importBreakpointBed() no longer calls closeAllConnections(). The connection
  that rtracklayer::import.bed() opens is already closed by rtracklayer.
* The examples of importBreakpointBed() and rebinGenomicInteractions(), the
  only two that load GenomicInteractions, are wrapped in \donttest{}.

The 3.7.9 fix for the r-devel ERROR is unchanged. It is described below.

## Original submission (3.7.9)

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

* Local Ubuntu Linux, R 4.6.2, `R CMD check --as-cran` (3.7.10)
* CRAN incoming pretest of 3.7.9: r-devel-windows-x86_64 OK;
  r-devel-linux-x86_64-debian-gcc 1 NOTE (the example timing addressed above)
* Ubuntu Linux, R-devel with Bioconductor devel (S4Vectors 0.51.10,
  InteractionSet 1.41.0, GenomicRanges 1.65.4), R-release, R-oldrel; macOS and
  Windows, R-release (GitHub Actions, 3.7.9)

## R CMD check results

0 errors | 0 warnings | 0 notes
