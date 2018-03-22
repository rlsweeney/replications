pklist <- c("dtplyr", "curl", "foreign")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")

url.Saiz2010 <- "http://web.archive.org/web/20100619052721/http://real.wharton.upenn.edu/~saiz/"
filename.zip <- "SUPPLYDATA.zip"
filename.dta <- "HOUSING_SUPPLY.dta"

curl_download(paste(url.Saiz2010, filename.zip, sep = ""), destfile = filename.zip, quiet = FALSE)
unzip(filename.zip)
housing.supply <- read.dta(filename.dta)
unlink(filename.dta)
unlink(filename.zip)
rm(filename.dta, filename.zip, url.Saiz2010)
