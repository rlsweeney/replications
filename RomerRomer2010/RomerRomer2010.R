rm(list = ls())
pklist <- c("tidyverse", "dynlm", "data.table", "matrixStats", "boot")
source("https://raw.githubusercontent.com/fgeerolf/R/master/load-packages.R")
load(url("https://github.com/fgeerolf/data/raw/master/oecd/EO.RData"))

p1 <- 0.95
p2 <- 0.68

problist <- c((1-p1)/2, (1-p2)/2, 1-(1-p2)/2, 1-(1-p1)/2)
problist

make.ts <- function(DataNew){
  DataNew.ts <- ts(DataNew[,-1], 
                   start = c((min(DataNew[,1])%/%1), (min(DataNew[,1])%%1)*4+1), 
                   end = c((max(DataNew[,1])%/%1), (max(DataNew[,1])%%1)*4+1), frequency = 4) 
  return(DataNew.ts)
}

# Replicating Romer-Romer (2010) ----------------------

DATA <- JME2014_data %>% 
  as.data.table %>% 
  setkey(Date) %>% 
  select(Date, Tax.Revenues, Govt.Spending, Output) %>% 
  rename(yearqtr = Date) %>% 
  as.data.table %>% 
  setkey(yearqtr)

DATA <- RRTAXSHOCKS %>% 
  dplyr::select(Date = V1, RRshocks = V2)  %>% 
  rename(yearqtr = Date) %>% 
  as.data.table %>% 
  setkey(yearqtr) %>% 
  merge(DATA, all.x = TRUE) %>% 
  filter(!is.na(Output))

DATA.ts <- ts(DATA[,-1], 
              start = c((min(DATA[,1])%/%1), (min(DATA[,1])%%1)*4+1), 
              end = c((max(DATA[,1])%/%1), (max(DATA[,1])%%1)*4+1), frequency=4) 
DATA.ts

lm.1 <- dynlm(d(Output) ~ L(RRshocks, 1:20), DATA.ts)
summary(lm.1)
irf1 <- cbind(as.matrix(unname(lm.1$coefficients[-1])*100), as.matrix(unname(cumsum(lm.1$coefficients[-1])*100)))
irf1

# International Implications ----------------------

listvariables <- c("GDPV", "GDP", "CGV", "CPV", "ITV","IGV", "MGSV", "XGSV", "EXCHER", "EXCHEB", "IHV")
listvariablesLN <- paste("LN", listvariables, sep = "")

dataRomer <- EO %>% 
  filter(VARIABLE %in% listvariables & FREQUENCY=="Q") %>% 
  mutate(countryname = paste(Country),
         Date = as.numeric(substr(TIME,1,4)) + (as.numeric(substr(TIME,7,7))-1)/4,
         VARIABLE = paste("LN", paste(VARIABLE), sep = ""),
         Value = log(Value)) %>% 
  filter(countryname == "United States") %>% 
  dplyr::select(Date, variable = VARIABLE, value = Value) %>% 
  dcast(Date ~ variable, value.var = "value") %>% 
  mutate(Govt.Spending = log(exp(LNCGV) + exp(LNIGV))) %>% 
  as.data.table %>% 
  setkey(Date)  %>% 
  rename(yearqtr = Date) %>% 
  as.data.table() %>% 
  setkey(yearqtr) %>% 
  merge(DATA, all = TRUE)

dataRomer <- make.ts(dataRomer)
dataRomer

colnames(dataRomer)


Q <- 12
resultsRomer <- array(0, c(Q, 2, length(listvariablesLN)), dimnames = list(paste("t=",1:Q,sep=""), c("Estimates", "Cumsum"), listvariablesLN))
resultsRomer

for (var in listvariablesLN){
  cat("\nCurrently showing IRF:", var,"\n ")
  lm <- dynlm(d(get(var)) ~ L(RRshocks, 1:12), dataRomer)
  summary(lm)
  resultsRomer[,,var] <- cbind(as.matrix(unname(lm$coefficients[-1])*100), as.matrix(unname(cumsum(lm$coefficients[-1])*100)))
  print(resultsRomer[,,var])
}

resultsRomer


# Now, calculate the standard errors -- by Bootstrapping those standard errors ----------

# Bootstrap 1 ------------

Q <- 12
IRFs.Romer <- array(0, c(Q+1, 6, length(listvariablesLN)), dimnames = list(paste("t=", 0:Q, sep=""), c("lag", "5%", "32%", "Pt", "68%", "95%"), listvariablesLN))



for (var in listvariablesLN){
  cat("Currently Dealing with: ", var, "\n")
  
  dataRomer2 <- dataRomer %>%
    as.data.table %>%
    mutate(D1var = get(var) - lag(get(var),1)) %>%
    do(data.frame(., setNames(shift(.$RRshocks, 1:Q), c(paste("RRshocksL", 1:Q, sep = ""))))) %>%
    select(D1var, paste("RRshocksL", 1:12, sep="")) %>%
    filter(!is.na(RRshocksL12)) %>%
    filter(!is.na(D1var))
  
  nrep <- 1000
  matRep <- array(0,c(nrep, Q))
  
  for (r in 1:nrep){
    dataRomer_sample <- cbind(D1var = lm(D1var ~ ., dataRomer2)$fitted.values + 
                                 sample(lm(D1var ~ ., dataRomer2)$residuals, nrow(dataRomer2), replace = TRUE), dataRomer2[,-1])
    matRep[r,] <- cumsum(lm(D1var ~ ., dataRomer_sample)$coefficients[-1])
  }
 
  IRFs.Romer[2:(Q+1) , c("lag"), var] <- 1:Q
  IRFs.Romer[2:(Q+1) , c("5%", "32%", "68%", "95%"), var] <- colQuantiles(matRep, probs = problist)
  IRFs.Romer[2:(Q+1) , "Pt", var] <- cumsum(lm(D1var ~ ., dataRomer2)$coefficients[-1])
}


IRFs.Romer



plotirf <- function(plotirfmelted, plottitle, Q){
  ggplot(plotirfmelted, aes(x = lag, y = value), group = variable) + 
    geom_line(aes(colour = factor(variable))) + 
    scale_color_manual(values = c("#D3D3D3", "#AAAAAA", "#000000", "#AAAAAA", "#D3D3D3")) +
    ggtitle(plottitle) +
    geom_polygon(data = plotirfmelted[c(1:Q, (2*Q):(Q+1)), ], aes(x = lag, y = value), alpha = 0.25)  +
    geom_polygon(data = plotirfmelted[c((Q+1):(2*Q), (3*Q):(2*Q+1)), ], aes(x = lag, y = value), alpha = 0.5) +
    geom_polygon(data = plotirfmelted[c((2*Q+1):(3*Q), (4*Q):(3*Q+1)),], aes(x = lag, y = value), alpha = 0.5) +
    geom_polygon(data = plotirfmelted[c((3*Q+1):(4*Q), (5*Q):(4*Q+1)),], aes(x = lag, y = value), alpha = 0.25) + 
    theme_bw() +
    theme(legend.position="none",
          axis.text.x = element_text(face = "bold", color="#000000", size = 14, angle = 0),
          axis.text.y = element_text(face = "bold", color="#000000", size = 14, angle = 0),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x =  element_text(face = "bold", color = "#000000", size = 14, angle = 0),
          axis.title.y =  element_text(face = "bold", color = "#000000", size = 14, angle = 90)) +
    xlab("Lags") + ylab("% Deviation from Trend") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_x_discrete(limits = c(seq(0, Q, floor(Q/5))))
}



for (var in listvariablesLN){
  IRF.melted <- as.data.frame(IRFs.Romer[,,var]) %>%
    as.data.table %>%
    melt(id.vars = "lag")
  plotirf(IRF.melted, paste("Response of: ", var, " in Romer/Romer (2010)", sep = ""), Q+1)
  ggsave(paste("graphs/", var, ".pdf", sep = ""), width = 9.57, height = 6.35)
}

# From Romer-Romer 2010 -----------


system("pwd")
system("ssconvert -S DataSet/RomerRomerFiscalData.xlsx DataSet/RomerRomerFiscalData.csv")

RomerRomerFiscalData1 <- read.csv('DataSet/RomerRomerFiscalData1.csv', skip = 11) %>%
  rename(yearqtr = X) %>%
  select(-starts_with("X.")) %>%
  mutate(yearqtr = (row_number()-1)*0.25 + 1945) %>%
  select(yearqtr, everything()) %>%
  filter(yearqtr != 2008)

RomerRomerFiscalData1.extract <- RomerRomerFiscalData1 %>%
  select(yearqtr, EXOGENR)

RomerRomerFiscalData2 <- read.csv('DataSet/RomerRomerFiscalData2.csv', skip = 18) %>%
  rename(yearqtr = X) %>%
  select(-starts_with("X.")) %>%
  mutate(yearqtr = (row_number()-1)*0.25 + 1945) %>%
  select(yearqtr, everything()) %>%
  filter(yearqtr != 2008)

Q <- 20
method <- "D1"

coeff.to.irf <- function(coeff){ # export IRFs
  if (method == "D1") return(cumsum(coeff[2:(Q+1), 1])) 
}


bootfun <- function(formula, data, indices) {
  d <- data[indices,]
  lm.1 <- lm(formula, data = d)
  return(coeff.to.irf(summary(lm.1)$coefficients))
}

plotirf <- function(irfdata, plottitle, Q){
  ggplot(irfdata, aes(x = lag, y = value), group = variable)  +  theme_bw() + 
    geom_line(aes(colour = factor(variable))) + 
    ggtitle(plottitle) + scale_color_manual(values = c("#D3D3D3", "#AAAAAA", "#000000", "#AAAAAA", "#D3D3D3")) + 
    geom_polygon(data = irfdata[c(1:Q, (2*Q):(Q+1)), ], aes(x = lag, y = value), alpha = 0.15)  +
    geom_polygon(data = irfdata[c((Q+1):(2*Q), (3*Q):(2*Q+1)), ], aes(x = lag, y = value), alpha = 0.3) +
    geom_polygon(data = irfdata[c((2*Q+1):(3*Q), (4*Q):(3*Q+1)),], aes(x = lag, y = value), alpha = 0.3) +
    geom_polygon(data = irfdata[c((3*Q+1):(4*Q), (5*Q):(4*Q+1)),], aes(x = lag, y = value), alpha = 0.15) + 
    theme(legend.position = "none",
          axis.text.x = element_text(face = "bold", color="#000000", size = 14, angle = 0),
          axis.text.y = element_text(face = "bold", color="#000000", size = 14, angle = 0),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x =  element_text(face = "bold", color = "#000000", size = 14, angle = 0),
          axis.title.y =  element_text(face = "bold", color = "#000000", size = 14, angle = 90)) +
    xlab("Lags") + ylab("% Deviation from Trend") + geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(limits = c(seq(0, Q, floor(Q/5))))
}

dir.create("graphs/romer-romer")

romerromer2010 <- function(variable){
  formul <- as.formula(paste(variable, " ~ . - yearqtr", sep = ""))
  
  RomerRomerFiscalData2.extract <- RomerRomerFiscalData2 %>%
    select(yearqtr, variable, NOMGDP) %>%
    merge(RomerRomerFiscalData1.extract, by = "yearqtr", all = TRUE) %>%
    mutate(EXOGENR_GDP = EXOGENR / NOMGDP) %>%
    mutate_at(vars(variable), funs(log(.) - lag(log(.),1))) %>%
    do(data.frame(., setNames(shift(.$EXOGENR_GDP, 1:Q), paste("EXOGENR_GDPL", 1:Q, sep = "")))) %>%
    select(yearqtr, variable, starts_with("EXOGENR_GDP")) %>%
    select(-EXOGENR_GDP) %>%
    na.omit
  
  nboot <- 1000
  irf <- matrix(0, Q+1, 6)
  irf[, 1] <- 0:Q
  irf[2:(Q+1),4] <- coeff.to.irf(summary(lm(formul, data = RomerRomerFiscalData2.extract))$coefficients)
  results <- boot(data = RomerRomerFiscalData2.extract, statistic = bootfun, R = nboot, formula = formul)
  irf[2:(Q+1),c(2:3, 5:6)] <- colQuantiles(results$t, probs = problist)
  
  irf %<>% as.data.table %>%
    rename(lag = V1, p5 = V2, p32 = V3, point = V4, p68 = V5, p95 = V6) %>%
    melt(id.vars = "lag")
  
  return(plotirf(irf, paste("Romer - Romer (2010): ", variable, sep = ""), Q+1))
}

ggplotly(romerromer2010("GDP"))

romerromer2010("EX")
romerromer2010("IM")
romerromer2010("PCE")
romerromer2010("DUR")
romerromer2010("NONDUR")
romerromer2010("PGDP")
romerromer2010("NOMGDP")


