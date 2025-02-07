list.files("E:/RaMBat")
setwd("E:/RaMBat")
getwd()

usethis::use_data(all_13datasets, overwrite = TRUE, compress = "xz")
file.exists("data/all_13datasets.rda")
usethis::use_data(GSE85217, overwrite = TRUE, compress = "xz")
file.exists("data/GSE85217.rda")
usethis::use_data(samp_13, overwrite = TRUE, compress = "xz")
usethis::use_data(sampAnnot_GSE85217, overwrite = TRUE, compress = "xz")

usethis::use_data(MB_RANK_GP, overwrite = TRUE, compress = "xz")

usethis::use_data(all_rank_t_genes, overwrite = TRUE, compress = "xz")

usethis::use_data(all_reversed_gp_genes, overwrite = TRUE, compress = "xz")




usethis::use_data(all_13datasets, GSE85217, samp_13, sampAnnot_GSE85217, overwrite = TRUE, compress = "xz")

usethis::create_package("E:/R-Studio/RStudio/new")



devtools::install("E:/RaMBat", force = TRUE)
