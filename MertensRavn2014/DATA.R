pklist <- c("tidyverse", "data.table", "gdata",  "curl")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")
file <- "jme2014_data.xls"
url <- "https://karelmertenscom.files.wordpress.com/2017/09/"
curl_download(paste(url, file, sep = ""), file, quiet = FALSE)
DATA <- read.xls(file, sheet = 1, skip = 1)
unlink(file)