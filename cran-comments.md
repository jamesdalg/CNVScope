## Resubmission

This submission updates the maintainer email address to a working one, as
requested by CRAN (email to the previous address bounced). The version number
has been bumped accordingly.

It also addresses the issues reported for the previous version (3.7.2) at
<https://cran.r-project.org/web/checks/check_results_CNVScope.html>:

* Windows installation now succeeds. The earlier "lazy loading failed ...
  there is no package called 'dbplyr'" error came from a transitive
  dependency load and no longer reproduces; the package installs and checks
  cleanly on Windows (verified on win-builder and on the Windows R-release
  GitHub Actions runner).
* The "Namespace in Imports field not imported from: 'Hmisc'" NOTE is
  resolved; 'Hmisc' is no longer declared in Imports.
* DESCRIPTION now provides an `Authors@R` field (the previous version had
  none).

## Test environments

* Local: Ubuntu Linux, R 4.6.1
* GitHub Actions:
  * Ubuntu Linux, R-release, R-devel, R-oldrel-1
  * macOS, R-release
  * Windows, R-release
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* The note reports the new maintainer address:

  ```
  New maintainer:
    James Dalgleish <jamesdalg@gmail.com>
  Old maintainer(s):
    James Dalgleish <james.dalgleish@nih.gov>
  ```

  This is the intended change and the reason for this resubmission.
