## Setup
check.packages <- function(package){
  new.package <- package[!(package %in% installed.packages()[, "Package"])]
  if (length(new.package)) 
    install.packages(new.package, dependencies = TRUE)
  sapply(package, require, character.only = TRUE)
}

packages<-c("tidyverse", "readxl", "car", "MASS", "lme4", "arm", "here")

check.packages(packages)

raw.inputs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 2)

CH4.inputs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 1)

LoQs <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 3)
LOQs <- structure(as.numeric(LoQs), names=colnames(LoQs))

codings <- read_excel('Reactor Raws/Johnsie Clean Data 2.xlsx', sheet = 5) %>%
  mutate(crossed_factors = paste0(Season,"_",biotype)) %>%
  group_by(crossed_factors)%>%
  mutate(dummy_rep = seq_along(crossed_factors))


## Data Munging

#clean up raw measurements
clean_data <- raw.inputs %>%
  dplyr::filter(!is.na(PFBA)) %>%
  dplyr::select(-`Sample Date`, -`Date Analyzed`)

for (i in 1:length(LOQs)) {
  n = i+2
  clean_data[which(clean_data[,(n)] == "<LOQ"),(n)] <- LOQs[i]/sqrt(2) 
  clean_data[which(clean_data[,(n)] == "LOQ"),(n)] <- LOQs[i]/sqrt(2) 
  clean_data[which(clean_data[,(n)] == "ND"),(n)] <- LOQs[i]/10
  clean_data[which(clean_data[,(n)] == "0"),(n)] <- LOQs[i]/10
}

#clean_data[clean_data == "ND" | is.na(clean_data)] <- 0

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

final_raw <- clean_data2 %>% full_join(codings) %>%
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

write_csv(final_raw,"raw_dataframe.csv")
