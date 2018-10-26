## Setup
check.packages <- function(package){
  new.package <- package[!(package %in% installed.packages()[, "Package"])]
  if (length(new.package)) 
    install.packages(new.package, dependencies = TRUE)
  sapply(package, require, character.only = TRUE)
}

packages<-c("tidyverse", "readxl", "car", "MASS", "lme4", "arm","here")

check.packages(packages)

raw.inputs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 2)

CH4.inputs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 1)

LoQs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 3)
LOQs <- structure(as.numeric(LoQs), names=colnames(LoQs))

codings <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 5)

## Data Munging

#clean up raw measurements
clean_data <- raw.inputs %>%
  dplyr::filter(!is.na(PFBA)) %>%
  dplyr::select(-`Sample Date`, -`Date Analyzed`)

for (i in 1:length(LOQs)) {
  n = i+2
  clean_data[which(clean_data[,(n)] == "<LOQ"),(n)] <- LOQs[i]/sqrt(2) 
}

clean_data[clean_data == "ND" | is.na(clean_data)] <- 0

clean_data2 <- clean_data %>%
  mutate_at(colnames(clean_data[2:72]), as.numeric)

# calculate methane phase levels
CH4.phase <- CH4.inputs %>%
  group_by(Reactor) %>%
  summarize(max.CH4 = max(CH4)) %>%
  left_join(CH4.inputs) %>%
  mutate(peak.day = ifelse((CH4 == max.CH4), Days.Elapsed,NA)) %>%
  dplyr::filter(!is.na(peak.day)) %>%
  dplyr::select(peak.day,Reactor)%>%
  full_join(CH4.inputs) %>%
  mutate(phase = ifelse(CH4 < 0.2 & Days.Elapsed < peak.day, "lag", "peak"),
         phase = ifelse(CH4 < 0.2 & Days.Elapsed > peak.day, "decline", phase)) %>%
  dplyr::select(Days.Elapsed,Reactor,phase) %>%
  dplyr::filter(!is.na(phase))

phase_factoring <- function(df,reactor.list,lag.knot,decline.knot) {
  df %>% mutate(phase = case_when(Reactor %in% reactor.list & Days.Elapsed <= lag.knot ~ "lag",
                                  Reactor %in% reactor.list & Days.Elapsed > lag.knot & Days.Elapsed < decline.knot ~ "peak",
                                  Reactor %in% reactor.list & Days.Elapsed >= decline.knot ~ "decline",
                                  TRUE ~ as.character(phase)))
                
}

final_data <- clean_data2 %>% full_join(codings) %>%
  mutate(phase = "waiting") %>%
  phase_factoring(c("P2"), 74, 220)%>%
  phase_factoring(c("P3"), 311, 458)%>%
  phase_factoring(c("P20"), 75, 335)%>%
  phase_factoring(c("P21"), 67, 257)%>%
  phase_factoring(c("P24"), 20, 140)%>%
  phase_factoring(c("P25"), 33, 208)%>%
  phase_factoring(c("P25"), 33, 208)%>%
  phase_factoring(c("P7"), 20, 132)%>%
  phase_factoring(c("P6"), 54, 193)%>%
  phase_factoring(c("P4","P5"), 74, 220) %>%
  phase_factoring(c("P22","P23"), 109, 296)%>%
  phase_factoring(c("P26","P27"),26, 174) %>%
  phase_factoring(c("P8","P9"), 37, 163) %>%
  mutate_at(c("Reactor", "Season", "biotype","phase"), as.factor)

PFASs <- colnames(final_data)[3:72]

PFASs.raw <- dplyr::select(final_data, Days.Elapsed, Reactor, phase, biotype, Season, dil, one_of(PFASs)) %>%
  gather(one_of(PFASs), key = cmp, value = raw) %>%
  group_by(cmp,Reactor) %>%
  arrange(Days.Elapsed) %>%
  mutate(conc = raw*dil,
         lag.conc = conc - dplyr::lag(conc, n=1, default = NA, order_by=Days.Elapsed),
         lag.days = Days.Elapsed - dplyr::lag(Days.Elapsed, n=1, default = NA, order_by=Days.Elapsed)) %>%
  mutate(rate = lag.conc/lag.days) %>%
  #dplyr::filter(!is.na(rate)) %>%
  mutate(rel.rate = rate/conc) %>%
  group_by(Reactor) %>%
  mutate(tt = factor(seq_along(Days.Elapsed))) %>%ungroup()

PFASs.raw$phase <- factor(PFASs.raw$phase, levels = c("lag","peak","decline"))

chosen <- "PFHxA"

stats.raw <- dplyr::filter(PFASs.raw, cmp == chosen)

BC.transform <- boxcox(conc ~ phase + biotype + Season + Reactor,
                       plotit = TRUE,
                       data = stats.raw)

lambda <- BC.transform$x[which(BC.transform$y == max(BC.transform$y))]

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}

bc2 <- boxcoxfit(stats.raw$conc, lambda2 = TRUE)

lambda1 <- bc2$lambda[1]
lambda2 <- bc2$lambda[2]

response <- powerTransform(stats.raw$conc,lambda)

stats.raw <- mutate(stats.raw, response = response)

ggplot(stats.raw) +
  theme_bw()+
  geom_line(aes(x = Days.Elapsed, y = conc, linetype=biotype, group = Reactor))+
  geom_point(aes(x = Days.Elapsed, y= conc, color = phase))+
  #scale_linetype_manual( values = rep(1,length(levels(PFOA.raw$Reactor))))+
  #guides(linetype = FALSE)+
  facet_wrap(~Season)


fm.null <- lmer(conc ~ (1|Reactor),
                data = stats.raw)

scale.data <- mutate(stats.raw,
                     response_sc = scale(stats.raw$conc, scale = TRUE, center = TRUE))

fm.null2 <- lmer(conc ~ (phase|Reactor),
                 data = stats.raw)

anova(fm.null,fm.null2)

fm.lme1 <- lmer(log10(conc) ~ biotype + Season + phase + (phase|Reactor),
                data = stats.raw)



anova(fm.lme1,fm.null2)
Anova(fm.lme1)

mod <- emmeans(fm.lme1, c("biotype"), type = "response")

CLD(mod)

mod2 <- emmeans(fm.lme1, c("phase"), type = "response")

CLD(mod2)

mod3 <- emmeans(fm.lme1, ~ biotype + phase, type = "response")

CLD(mod3)

summary(fm.lme1)


xyplot(fitted(fm.lme1) ~ biotype|phase,
       groups = Reactor,
       data = stats.raw)

bwplot(fitted(fm.lme1) ~ biotype|phase,
       groups = Season,
       data = stats.raw)

ggplot(stats.raw) +
  theme_bw()+
  geom_line(aes(x = Days.Elapsed, y = response, linetype=biotype, group = Reactor))+
  geom_point(aes(x = Days.Elapsed, y= response, color = phase))+
  #scale_linetype_manual( values = rep(1,length(levels(PFOA.raw$Reactor))))+
  #guides(linetype = FALSE)+
  facet_wrap(~Season)
  



## test boxcox instead of log for later PFASs

BC.transform <- boxcox(PFOA ~ phase + biotype + Season + Reactor,
       plotit = FALSE,
       data = PFOA.raw)

lambda <- BC.transform$x[which(BC.transform$y == max(BC.transform$y))]

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}

mypar <- c(1,2)

qqnorm(powerTransform(PFOA.raw$PFOA, lambda))
qqline(powerTransform(PFOA.raw$PFOA, lambda))

qqnorm(log10(PFOA.raw$PFOA))
qqline(log10(PFOA.raw$PFOA))

library("geoR")
bc2 <- boxcoxfit(PFOA.raw$PFOA, lambda2 = TRUE)

lambda1 <- bc2$lambda[1]
lambda2 <- bc2$lambda[2]

qqnorm(powerTransform(PFOA.raw$PFOA, lambda1, lambda2))
qqline(powerTransform(PFOA.raw$PFOA, lambda1, lambda2))
