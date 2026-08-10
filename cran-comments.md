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

* Local: Ubuntu Linux, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

`R CMD check --as-cran` is clean locally. CRAN's incoming checks may add the
usual "Days since last update" NOTE; this submission is solely to clear the
r-devel Rd NOTE above.
