rm(list = ls())
pklist <- c("tidyverse", "data.table", "gdata",  "curl", "pracma", "nleqslv", "matrixStats")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")

source("https://raw.githubusercontent.com/fgeerolf/econ-221/master/course2/plotirf.R")
load(url("https://github.com/fgeerolf/replications/raw/master/MertensRavn2014/DATA.RData"))
load(url("https://github.com/fgeerolf/datasets/raw/master/oecd/EO.RData"))


# fntarget ----------------------

fntarget <- function(x, Sigma){
  thetaG <- x[1] 
  sigmaG <- x[2]
  sigmaT <- x[3]
  sigmaY <- x[4]
  zetaT <- x[5]
  zetaG <- x[6]
  A <- matrix(c(1, 0, -thetaY, 0, 1, -gammaY, -zetaT, -zetaG, 1), ncol = 3, byrow = TRUE)
  B <- matrix(c(sigmaT, thetaG*sigmaG, 0, gammaT*sigmaT, sigmaG, 0, 0, 0, sigmaY), 
              ncol = 3, byrow = TRUE)
  D <- mldivide(A, B, pinv = TRUE)
  DD <- D%*%t(D)
  return(list("objective" = DD[lower.tri(DD, diag = TRUE)] - Sigma[lower.tri(Sigma, diag = TRUE)], 
              "A" = A, "B" = B, "D" = D))
}

# calculateirf ----------------------

calculateirf <- function(dates, VARS){
  T <- dim(VARS)[1]
  n <- dim(VARS)[2]
  
  VARSLAGS1 <- xts::lag.xts(VARS, 1)
  colnames(VARSLAGS1) <- paste(colnames(VARSLAGS1), "L", as.character(1), sep = "")
  VARSLAGS2 <- xts::lag.xts(VARS, 2)
  colnames(VARSLAGS2) <- paste(colnames(VARSLAGS2), "L", as.character(2), sep = "")
  VARSLAGS3 <- xts::lag.xts(VARS, 3)
  colnames(VARSLAGS3) <- paste(colnames(VARSLAGS3), "L", as.character(3), sep = "")
  VARSLAGS4 <- xts::lag.xts(VARS, 4)
  colnames(VARSLAGS4) <- paste(colnames(VARSLAGS4), "L", as.character(4), sep = "")
  
  VARSLAGS1 <- VARSLAGS1[(p+1):T,]
  VARSLAGS2 <- VARSLAGS2[(p+1):T,]
  VARSLAGS3 <- VARSLAGS3[(p+1):T,]
  VARSLAGS4 <- VARSLAGS4[(p+1):T,]
  VARS <- VARS[(p+1):T,]
  VAR.T <- dim(VARS)[1]
  VAR.n <- dim(VARS)[2]
  VAR.DET <- matrix(c(rep.int(1, T), 1:T, (1:T)^2, (dates==1975.25)), nrow = T, ncol = n+1, byrow = FALSE, 
                    dimnames = list(NULL,c("constant", "trend", "trend2", "dummy")))
  VAR.DET <- VAR.DET[(p+1):T,]
  VAR.bet <- mldivide(cbind(VARSLAGS1, VARSLAGS2, VARSLAGS3, VARSLAGS4, VAR.DET), VARS, pinv = FALSE)
  VAR.bet2 <- mldivide(cbind(VARSLAGS1, VARSLAGS2, VARSLAGS3, VARSLAGS4, VAR.DET), VARS, pinv = TRUE)
  ALLDATA <- as.data.table(cbind(VARS, VARSLAGS1, VARSLAGS2, VARSLAGS3, VARSLAGS4, VAR.DET))
  summary(lm(Tax.Revenues ~ .  - Govt.Spending - Output - 1, data = ALLDATA))
  
  res <- VARS - cbind(VARSLAGS1, VARSLAGS2, VARSLAGS3, VARSLAGS4, VAR.DET)%*%VAR.bet
  Sigma <- (t(res)%*%res)/(T - n*p - 1)
  
  #results <- nleqslv(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), fntarget, control = list(btol = 0.01))
  results <- nleqslv(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                     function(x) fntarget(x, Sigma)$objective, 
                     method = "Newton", 
                     control = list(xtol = 1e-200))
  x <- results$x
  thetaG <- x[1] 
  sigmaG <- x[2]
  sigmaT <- x[3]
  sigmaY <- x[4]
  zetaT <- x[5]
  zetaG <- x[6]
  A <- matrix(c(1, 0, -thetaY, 0, 1, -gammaY, -zetaT, -zetaG, 1), ncol = 3, byrow = TRUE)
  B <- matrix(c(sigmaT, thetaG*sigmaG, 0, gammaT*sigmaT, sigmaG, 0, 0, 0, sigmaY), 
              ncol = 3, byrow = TRUE)
  D <- mldivide(A, B, pinv = FALSE)
  DD <- D%*%t(D)
  
  
  irs <- matrix(0, nrow = p + irfhorizon, ncol = n)
  irs[p+1,] <- -100*D[,1]/D[1,1]*VAR.tshocksize
  for (j in 2:(irfhorizon)){
    lvars <- t(irs[seq(p+j-1, j, -1),])
    irs[p+j,] <- t(as.vector(lvars))%*%VAR.bet[1:(p*n),]
  }
  irs <- irs[(p+1):nrow(irs),]
  irsg <- matrix(0, nrow = p + irfhorizon, ncol = n)
  irsg[p+1,] <- 100*D[,2]/D[2,2]*VAR.gshocksize
  for (j in 2:(irfhorizon)){
    lvarsg <- t(irsg[seq(p+j-1, j, -1),])
    irsg[p+j,] <- t(as.vector(lvarsg))%*%VAR.bet[1:(p*n),]
  }
  irsg <- irsg[(p+1):nrow(irsg),]
  VAR.e <- t(mldivide(D, t(res), pinv = FALSE)) # structural shocks
  return(list("irfT" = irs, "irfG" = irsg, "shocks" = VAR.e, "D" = D, "VAR.bet" = VAR.bet))
}


# Blanchard-Perotti (2002) ---------------------------------

thetaY <- 2.08
gammaT <- 0
gammaY <- 0
p <- 4
irfhorizon <- 20
DAT.TRY <- 0.1746  # Average ratio of federal tax revenues to GDP
DAT.GY <- 0.0981  # Average ratio of federal expenditures to GDP
VAR.tshocksize <- 0.01 / DAT.TRY
VAR.gshocksize <- 0.01 / DAT.GY 

# Import Data ---------------------------

#DATA <- read.xlsx('JME2014_data.xlsx', sheetIndex = 1, startRow = 2)


DATA2 <- EO %>% 
  filter(VARIABLE %in% c("GDPV", "GDP", "CGV", "IGV") & FREQUENCY=="Q") %>% 
  mutate(countryname = paste(Country),
         Date = as.numeric(substr(TIME,1,4))+(as.numeric(substr(TIME,7,7))-1)/4,
         VARIABLE = paste(VARIABLE),
         Value = log(Value)) %>% 
  filter(countryname == "United States") %>% 
  select(Date, variable = VARIABLE, value = Value) %>% 
  dcast(Date ~ variable, value.var = "value") %>% 
  mutate(Govt.Spending = log(exp(CGV) + exp(IGV)),
         Output = GDPV) %>% 
  select(Date, Output, Govt.Spending) %>% 
  as.data.table() %>% 
  setkey(Date)


dates <- as.matrix(DATA[,1])
VARS <- as.matrix(DATA[,2:4])
VARS

irs <- calculateirf(dates, VARS)
irf <- as.vector(irs$irfT)


# Confidence Intervals by Recursive WILD bootstrap --------------

nboot <- 100
p1 <- 0.68
p2 <- 0.90 

# Using a bootstrap : multiply by 1 and -1 (adding noise)
jj <- 1
irfci <- matrix(NA, nrow = nboot, ncol = size(irf)[2])
while (jj <= nboot){
  VARS <- DATA[,2:4]
  T <- dim(VARS)[1] - p
  n <- dim(VARS)[2]
  rr <- 1 - 2*(runif(T, min = 0, max = 1) > 0.5)
  eb <- (irs$shocks)*(rr%*%matrix(1, nrow = 1, ncol = 3))
  resb <- (irs$D)%*%t(eb)
  
  varsb <- matrix(0, nrow = p+T, ncol = n)
  varsb[1:p,] <- as.matrix(VARS[1:p,])
  VAR.DET <- matrix(c(rep.int(1, T+p), 1:(T+p), (1:(T+p))^2, (dates==1975.25)), 
                    nrow = T+p, ncol = n+1, byrow = FALSE, 
                    dimnames = list(NULL, c("constant", "trend", "trend2", "dummy")))
  VAR.DET <- VAR.DET[(p+1):(T+p),]
  VAR.bet <- irs$VAR.bet
  for (j in ((p+1):(p+T))){
    lvars <- t(varsb[seq(j-1, j-p, -1),])
    varsb[j,] <- t(as.vector(lvars))%*%VAR.bet[1:(p*n),] + VAR.DET[j-p,]%*%VAR.bet[(p*n+1):nrow(VAR.bet),] + t(resb[,j-p])
  }
  VARS <- varsb
  colnames(VARS) <- colnames(VAR.bet)
  irfci[jj,] <- as.vector(calculateirf(dates, VARS)$irfT)
  jj <- jj+1
}

irfcip1low <- matrix(colQuantiles(irfci, probs = (1-p1)/2), nrow = 20, byrow = FALSE)
irfcip1high <- matrix(colQuantiles(irfci, probs = 1-(1-p1)/2), nrow = 20, byrow = FALSE)
irfcip2low <- matrix(colQuantiles(irfci, probs = (1-p2)/2), nrow = 20, byrow = FALSE)
irfcip2high <- matrix(colQuantiles(irfci, probs = 1-(1-p2)/2), nrow = 20, byrow = FALSE)

bindData <- cbind(1:20 ,irfcip2low[,3], irfcip1low[,3], irs$irfT[,3], irfcip1high[,3], irfcip2high[,3])
colnames(bindData) <- c("lag", 
                        paste("IRF",as.character(p2*100),sep=""),
                        paste("IRF",as.character(p1*100),sep=""),
                        "IRF",
                        paste("IRF",as.character((1-p1)*100),sep=""),
                        paste("IRF",as.character((1-p2)*100),sep=""))

bindData <- as.data.table(bindData)

bindData <- rbind(list(0, 0, 0, 0, 0, 0), bindData)
bindData[,2:6] <- - bindData[,2:6]                                 

bindDatamelted <- bindData %>% 
  as.data.table %>% 
  melt(id.vars = "lag")


# plotirf -------------------------------------------------------


Q <- 20
Q <- Q + 1

maxtickers <- 40
bindDatamelted

plotirf(bindDatamelted, "Impulse Response from Blanchard - Perotti (2002)", Q)

Q <- Q - 1


