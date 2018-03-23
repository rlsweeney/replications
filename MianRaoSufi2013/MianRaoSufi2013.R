
pklist <- c("dtplyr", "foreign")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")

url.amir.sufi.data <- "http://faculty.chicagobooth.edu/amir.sufi/data-and-appendices/"

filename <- "MianRaoSufiQJE_housingnetworthshock.dta"
filename2 <- "CYd2i2006nyf_frommiansufi_20120621.dta"

curl_download(paste(url.amir.sufi.data, filename, sep = ""), filename, quiet = FALSE)
curl_download(paste(url.amir.sufi.data, filename2, sep = ""), filename2, quiet = FALSE)

housing.networth.shock <- read.dta(filename)
mian.sufi.2013 <- read.dta(filename2)

unlink(filename)
unlink(filename2)

rm(filename, filename2)
