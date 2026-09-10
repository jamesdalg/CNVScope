## Submission

This submission fixes the NOTE reported for version 3.7.6 at
<https://cran.r-project.org/web/checks/check_results_CNVScope.html> on
r-devel-linux-x86_64-debian-clang and r-devel-linux-x86_64-debian-gcc:

```
Rd files without \usage:
  'CNVScopeserver.Rd' 'calcCNVKernelProbDist.Rd'
  'downsample_genomic_matrix.Rd' 'formSampleMatrixFromRawGDCData.Rd'
  'getBlockAverageMatrixFromBreakpoints.Rd'
  'getInterchromosomalInteractivePlot.Rd' 'importBreakpointBed.Rd'
  'rebinGenomicInteractions.Rd'
\arguments should not be documented without \usage.
```

In each of the eight affected sources, a top-level `globalVariables()` (or
`dontCheck()`) call sat between the roxygen block and the function it
documents. roxygen2 therefore attached the block to that call rather than to
the function and emitted no `\usage` section, leaving `\arguments` documented
without `\usage`. Those calls now precede their roxygen blocks, and the
documentation has been regenerated; all eight topics have `\usage` again.

Restoring `\usage` surfaced two arguments of `formSampleMatrixFromRawGDCData()`
(`parallel` and `cnlabel`) that had never been documented; both are now
documented. Stray top-level test code that ran at package build time was also
removed from `R/downsample_genomic_matrix.R`.

There are no user-visible changes to any function's behaviour.

## Test environments

* Local: Ubuntu Linux, R 4.6.1 -- 0 errors, 0 warnings, 1 note (see below)
* win-builder: R-devel and R-release -- 0 errors, 0 warnings, 0 notes

## R CMD check results

0 errors | 0 warnings | 1 note (local only)

The note is the elapsed-time one, and only on the local machine:

```
* checking examples ... [52s/52s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
                     user system elapsed
importBreakpointBed 6.076  0.362    6.44
```

Nearly all of that is loading the `GenomicInteractions` namespace, not the
package's own work: the example's input is a 6-line, 306-byte BED file shipped
in `inst/extdata`, and the `importBreakpointBed()` call itself takes about 1.5s
once the namespace is attached. The example uses no external resources.

None of CRAN's own flavors reported this note for 3.7.6 -- the only note there
was the `Rd files without \usage` one this submission fixes -- so the example is
unchanged from the version currently on CRAN.
