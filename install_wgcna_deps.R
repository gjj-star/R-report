lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))
BiocManager::install(c("impute", "preprocessCore"), lib = lib, ask = FALSE, update = FALSE)
cat("WGCNA dependencies installed.\n")
