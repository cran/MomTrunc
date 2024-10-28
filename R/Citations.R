.onAttach <- function(libname, pkgname) {
  methodology_citation <- "Galarza, C. E., Matos, L. A., Castro, L. M., & Lachos, V. H. (2022). Moments of the doubly truncated selection elliptical distributions with emphasis on the unified multivariate skew-t distribution. Journal of Multivariate Analysis, 189, 104944. doi:10.1016/j.jmva.2021.104944."
  package_citation <- "Galarza, C. E., Matos, L. A., Castro, L. M., & Lachos, V. H. (2024). MomTrunc: Moments of Truncated Multivariate Distributions. R package version 6.1. R Foundation for Statistical Computing. https://CRAN.R-project.org/package=MomTrunc"
  packageStartupMessage("Thank you for using MomTrunc!")
  packageStartupMessage("To acknowledge our work, please cite the following:")
  packageStartupMessage("  ")
  packageStartupMessage("1. Methodology behind the package:")
  packageStartupMessage(methodology_citation)
  packageStartupMessage("  ")
  packageStartupMessage("2. The R package itself:")
  packageStartupMessage(package_citation)
}