## Test environments
* Local machine
  * Running under: Ubuntu 22.04.2 LTS
  * Platform: x86_64-pc-linux-gnu (64-bit)
  * R version 4.3.1 (2023-06-16)
* Github Actions:
  * Windows-latest (R release)
  * MacOS-latest (R release)
  * Ubuntu-latest (R release, devel, oldrel-1)
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 1 note

* This is an update.

NOTE comments:
* win-builder found three invalid URLs (error 403) in vignette d_other_sens. However after checking it, URLs are ok.

## Downstream dependencies
One downstream dependency, package apisensr, not affected by this update.
