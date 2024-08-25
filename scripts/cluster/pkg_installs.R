pkgs <- c('dplyr', 'tidyr', 'readr', 'tibble', 'fs', 'pls', 'spls')
package.check <- lapply(
  pkgs,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x,
                       repos = 'https://ftp.osuosl.org/pub/cran/',
                       dependencies = T)
      library(x, character.only = TRUE)
    }
  }
)
