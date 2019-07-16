CNVScope v2.4.6 (Release date: 2019-07-15)
==============
*Minor fix to freadGDCfile documentation

CNVScope v2.4.6 (Release date: 2019-07-15)
==============
*More information about BLCA, SKCM, and AML in additional examples vignette
**Their processed files have been placed within the package
**These are available only on the github site as the CRAN size limitations forbid it.
*Added the ability to specify chrM and chrY in formSampleMatrixFromRawGDCData function.

CNVScope v2.4.4 (Release date: 2019-07-09)
==============

Changes:
* Fixed vignettes not properly linking to the package.
CNVScope v2.4.4 (Release date: 2019-07-09)
==============

Changes:

* Added a power analysis vignette
* Fixed issues with formSampleMatrixFromRawGDCData
** These issues arose from changing types with underlying Bioconductor dependencies from character to factor Rle.
* Updated the build for R-base v3.6.1
