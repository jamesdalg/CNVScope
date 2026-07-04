## Resubmission

This submission fixes the Windows installation failure reported for version
3.7.5 at <https://cran.r-project.org/web/checks/check_results_CNVScope.html>:

```
Error in loadNamespace(i, ...) : there is no package called 'RSQLite'
ERROR: lazy loading failed for package 'CNVScope'
```

The RSQLite dependency was pulled in transitively (biomaRt -> BiocFileCache ->
RSQLite, and via GenomicInteractions), so CNVScope's own lazy loading failed
whenever that chain was unavailable on a build machine. This is the same class
of failure previously seen as "there is no package called 'dbplyr'".

To fix it durably, biomaRt and GenomicInteractions have been moved from
Imports to Suggests and are now used behind `requireNamespace()` guards at
their points of use. This removes the entire BiocFileCache chain (RSQLite,
dbplyr) from the package's strong dependencies, so installation no longer
depends on those packages being present.

## Test environments

* Local: Ubuntu Linux, R 4.6.1
* GitHub Actions:
  * Ubuntu Linux, R-release, R-devel, R-oldrel-1
  * macOS, R-release
  * Windows, R-release
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 2 notes

`R CMD check --as-cran` reports two NOTEs, both benign:

* "Days since last update: 1". This is a fast resubmission because 3.7.5 fails
  to install on r-oldrel-windows (the RSQLite/lazy-loading error above). This
  submission is solely to fix that installation failure.
* One example (`importBreakpointBed`) occasionally exceeds 5s elapsed on a
  loaded machine (~6.5s locally). It is a small, self-contained example with
  no external resources.
