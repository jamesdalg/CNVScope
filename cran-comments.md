## Submission

This supersedes the 3.7.7 submission currently in the queue -- please discard
it. The package is unchanged apart from the version number; only these comments
differ.

This fixes the NOTE reported for 3.7.6 on r-devel-linux-x86_64-debian-clang and
r-devel-linux-x86_64-debian-gcc:

```
Rd files without \usage:
  'CNVScopeserver.Rd' 'calcCNVKernelProbDist.Rd'
  'downsample_genomic_matrix.Rd' 'formSampleMatrixFromRawGDCData.Rd'
  'getBlockAverageMatrixFromBreakpoints.Rd'
  'getInterchromosomalInteractivePlot.Rd' 'importBreakpointBed.Rd'
  'rebinGenomicInteractions.Rd'
\arguments should not be documented without \usage.
```

In each of the eight sources a top-level `globalVariables()` or `dontCheck()`
call sat between the roxygen block and the function it documents, so roxygen2
attached the block to that call and emitted no `\usage`. Those calls now precede
their roxygen blocks and the documentation has been regenerated; all eight
topics have `\usage` again.

This also surfaced two arguments of `formSampleMatrixFromRawGDCData()`
(`parallel`, `cnlabel`) that had never been documented; both are now documented.
No function's behaviour changes.

## Test environments

* Ubuntu Linux, R-devel, R-release, R-oldrel (GitHub Actions)
* Windows, R-devel and R-release (win-builder)
* macOS, R-release (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes
