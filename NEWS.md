CNVScope v3.3.0 (Release Date 2020-11-28)
==============
*We have now added smoothie::kernel2d smooth
as an alternative smoother to spatialfil (now disabled) 
and now our package should once again be CRAN ready.

CNVScope v3.2.9 (Release Date 2020-11-27)
==============
*With spatialfil removed from CRAN, we have disabled the spatialfil requirement
to run the additional examples.
*This algorithm is neat and we'll look at other ways to calculate interaction
probability, but for the moment,
we'll just disable things to keep the package on CRAN.

CNVScope v3.2.6 (Release Date 2020-8-14)
==============
*Fixed the gene search in the shiny application
*Added gene listing for row and column genes to the left of the plot.
**These gene data table outputs are searchable and the size
is adjustable.
*Disabled automatic switching to the sample data tab on click so that
users can properly use the new gene column feature.
*CRAN fixes to keep the package publicly available



CNVScope v3.1.8 (Release Date 2020-7-2)
==============
*Reduced imports to 20 and moved the remainder of packages into suggests.
*Fixed whole genome map.
*Optimized load time to plot a whole genome map (should be roughly a 50% reduction, according to profvis).

CNVScope v2.9.4 (Release Date 2019-12-7)
==============
*Removed Biostrings dependency to clear NOTE 
(although packages that CNVScope depends on require it). The prior Biostrings issue has been
fixed.

CNVScope v2.9.2 (Release Date 2019-12-5)
==============
*Added an explicit dependency, Biostrings to keep CRAN check errors on debian from occurring.
*Made a small fix to rebinGenomicInteractions, fixing a check at the end to handle a vector input
**This will make the package compatible with R 4.0.

CNVScope v2.8.2 (Release Date 2019-10-23)
==============
*Removed all blockseg references in code. It was a nice feature to have a third algorithm, but blockseg will not be on CRAN shortly.

CNVScope v2.8.1 (Release Date 2019-10-23)
==============
*Removed blockseg as a strong dependency, per CRAN.


CNVScope v2.7.7 (Release Date 2019-10-18)
==============
*Updated documentation (some new images have been added to show the effect of bin sizes and resolution).
*Fixed a broken URL that was causing a CRAN note.

CNVScope v2.7.5 (Release Date 2019-10-12)
==============
*Updated documentation for finely binned matrices (see the input matrix vignette).
*Several sizes of NBL whole genome matrices are available within the package (1e5,2.5e5,1e6,1e7,1e8).
**These are for package use only. The HD matrices are not speed-optimal for app use.
**Functionality of basic package functions remains the same as 2.7.2.5.
*passes CRAN checks again with flying colors (no notes, warnings, or errors).

CNVScope v2.7.2.5 (Release Date 2019-09-21)
==============
*Added P-value filter
*Added correlatoin options
*Added additional app visualization options.
*Added material on how to make higher resolution maps.

CNVScope v2.5.4 (Release date: 2019-07-19)
==============
*Minor changes to package dependency CNVScopePublicData.

CNVScope v2.5.2 (Release date: 2019-07-18)
==============
*Another minor cran fix to CNVScopeserver.

CNVScope v2.5.1 (Release date: 2019-07-18)
==============
*Small changes to CNVScopeserver to fix a CRAN check.

CNVScope v2.5.0 (Release date: 2019-07-18)
==============
*Small changes in 2.4.9 to the README, for publicity.
*Major changes to many functions to pass CRAN checks (global variable issues).
CNVScope v2.4.8 (Release date: 2019-07-17)
==============
*Small changes to dependencies (including pwr package) and vignettes.

CNVScope v2.4.7 (Release date: 2019-07-16)
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
*prior to this, this package now works with CNVScopePublicData, a large package containing
*all that is requisite to run the server locally.
