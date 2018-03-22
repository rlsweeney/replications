pklist <- c("dtplyr", "curl", "foreign")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")

url.Saiz2010 <- "http://real.wharton.upenn.edu/~gyourko/WRLURI/"
filename <- "WHARTON LAND REGULATION DATA_1_24_2008"
filename.html <- gsub(" ", "%20", filename)
filename.zip <- paste(filename, ".zip", sep = "")
filename.dta <- paste(filename, ".dta", sep = "")

curl_download(paste(url.Saiz2010, filename.html, ".zip", sep = ""), destfile = filename.zip, quiet = FALSE)

unzip(filename.zip)
wharton.land.regulation <- read.dta(filename.dta)
unlink(filename.zip)
unlink(filename.dta)
rm(url.Saiz2010, filename, filename.html, filename.zip, filename.dta)
