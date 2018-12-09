source("rawDataGeneration.R")

library(fitdistrplus)

##### PFAS Measurement Summaries and subsetting ####

PFASs <- colnames(final_raw)[3:72]

means <- final_raw[PFASs] %>%
  summarize_all(.funs = funs(mean))
median <- final_raw[PFASs] %>%
  summarize_all(.funs = funs(median))
sd <- final_raw[PFASs] %>%
  summarize_all(.funs = funs(sd))
min <- final_raw[PFASs] %>%
  summarize_all(.funs = funs(min))
max <- final_raw[PFASs] %>%
  summarize_all(.funs = funs(max))

summaries <- bind_rows(means,median,sd,min,max) %>% as.data.frame()
rownames(summaries) <- c("mean","median","sd","min","max")

PFASs_subset <- c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA",
  "FHEA","FHUEA","FPePA","FOEA","FHpPA",
  "PFBS","PFHxS","PFOS",
  "MeFBSAA","FTS.6.2","diPAP.6.2","diPAP.6.2.8.2","diPAP.8.2")

subset2 <- PFASs_subset[which(summaries["max",PFASs_subset] > 10)]

working_raws <- final_raw[subset2] %>%
  cbind(final_raw[c("Days.Elapsed","dil","Reactor","Season","biotype","phase")])

summary_plot <- working_raws %>%
  gather(one_of(subset2), key = "PFAS", value = "raw_measure") %>%
  mutate(conc = raw_measure * dil)

ggplot(summary_plot) +
  theme_minimal()+
  geom_point(aes(x = Days.Elapsed, y = conc, group = Reactor), size = 0.5) +
  geom_line(aes(x = Days.Elapsed, y = conc, group = Reactor, color = biotype)) +
  scale_y_log10() +
  facet_wrap(~PFAS)

##### LEGACY COMPOUNDS #####

legacy_summary <- filter(summary_plot, PFAS %in% c("PFBA","PFDA","PFOA","PFHpA","PFHxA","PFNA","PFOA","PFOS","PFPeA", "PFHxS"))

ggplot(filter(legacy_summary, phase != "lag")) +
  theme_minimal()+
  geom_point(aes(x = Days.Elapsed, y = conc, group = Reactor), size = 0.5) +
  geom_line(aes(x = Days.Elapsed, y = conc, group = Reactor, color = biotype)) +
  scale_y_log10() +
  facet_wrap(~PFAS)


##### F COMPOUNDS #####
telomer_summary <- filter(summary_plot, PFAS %in% c("FHEA","FHUEA","FPePA","FHpPA","FTS.6.2"))

ggplot(filter(telomer_summary, phase != "lag")) +
  theme_minimal()+
  geom_point(aes(x = Days.Elapsed, y = conc, group = Reactor), size = 0.5) +
  geom_line(aes(x = Days.Elapsed, y = conc, group = Reactor, color = biotype)) +
  scale_y_log10() +
  facet_wrap(~PFAS)

###### OTHER COMPOUNDS #####

other_summary <- filter(summary_plot, ! PFAS %in% telomer_summary$PFAS & ! PFAS %in% legacy_summary$PFAS)

ggplot(filter(other_summary, phase != "lag")) +
  theme_minimal()+
  geom_point(aes(x = Days.Elapsed, y = conc, group = Reactor), size = 0.5) +
  geom_line(aes(x = Days.Elapsed, y = conc, group = Reactor, color = biotype)) +
  scale_y_log10() +
  facet_wrap(~PFAS)


########## diPAP.6.2.8.2 ##########

chosen_PFAS <- "PFHxA"

single_PFAS <- filter(summary_plot, PFAS == chosen_PFAS)

levels(single_PFAS$phase) <- c("lag", "peak","decline")

  #####fitting truncated distribution
fitting_values <- tibble(raw_measure = single_PFAS$raw_measure) %>%
  mutate(index1 = seq_len(nrow(.))) %>%
  sample_n(size = nrow(.)) %>%
  mutate(index2 = seq_len(nrow(.)))

lower_bound <- LOQs[chosen_PFAS]/10
sqrtloq <- round(LOQs[chosen_PFAS]/sqrt(2), 5)
cutoff <- LOQs[chosen_PFAS]

fitting_values <-  fitting_values %>%
  arrange(raw_measure) %>%
  mutate(raw_measure = round(raw_measure, 5)) %>%
  mutate(left = case_when(raw_measure == sqrtloq ~ lower_bound,
                          raw_measure == lower_bound ~ as.numeric(NA),
                          TRUE ~ raw_measure)) %>%
  mutate(right = case_when(raw_measure == lower_bound ~ lower_bound,
                           raw_measure == sqrtloq ~ sqrtloq,
                           TRUE ~ raw_measure)) %>%
  arrange(raw_measure)

censdata <- fitting_values[c("left","right")] %>% as.data.frame()

trunc_dist <- fitdistcens(log(censdata), distr = "norm")

trunc_quantiles <- unlist(quantile(trunc_dist, probs = seq(0,100,by=(100/(max(fitting_values$index1)+2)))/100)[1])[-1]

lowestcens <- sum(fitting_values$raw_measure == lower_bound)
lowercens <- sum(fitting_values$raw_measure == sqrtloq)

lowestcens_quatiles <- trunc_quantiles[1:lowestcens]
lowercens_quantiles <- trunc_quantiles[(lowestcens+1):(lowestcens+lowercens)]

inputions <- c(lowestcens_quatiles,lowercens_quantiles)

length(inputions) <- nrow(fitting_values)

inputed_frame <- cbind(fitting_values, inputions) %>%
  arrange(index1) %>%
  mutate(quantiles = ifelse(is.na(inputions), log(raw_measure), inputions)) %>%
  dplyr::select(quantiles) %>%
  cbind(single_PFAS) %>%
  mutate(conc = ifelse(raw_measure < cutoff, exp(quantiles) * dil, raw_measure*dil),
         Days.Elapsed = Days.Elapsed/10)
#####

model <- lmerTest::lmer(data = inputed_frame,
              log(conc) ~ Days.Elapsed + biotype + (Days.Elapsed|Reactor))

model2 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + (phase|Reactor))

model3 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + (1|Reactor))

model4 <- lmer(data = inputed_frame,
               log(conc) ~ Days.Elapsed + biotype + (1|Reactor))

model5 <- lmer(data = inputed_frame,
               log(conc) ~ Days.Elapsed + biotype + (1|Season) + (Days.Elapsed|Reactor))

model6 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + (1|Season) + (phase|Reactor))

model7 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (phase|Reactor))

model8 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (Days.Elapsed|Reactor))

model9 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (1|Reactor))

model10 <- lmer(data = inputed_frame,
                log(conc) ~ phase + biotype + (1|Reactor) + (Days.Elapsed|Reactor))

model11 <- lmer(data = inputed_frame,
                log(conc) ~ Days.Elapsed + I(Days.Elapsed^2) + biotype + (1|Season/Reactor))


anova(model,model2,model3,model4,model5,model6,model7,model8,model9,model10)

preds <- predict(model3, newdata = new.data, type = "response")

inputed_frame <- left_join(inputed_frame, codings)

ggplot(inputed_frame, aes(x = Days.Elapsed*10, group = Reactor)) +
  theme_classic()+
  guides(color = FALSE, fill =FALSE, linetype = FALSE)+
  geom_point(aes(y = conc,
                 fill = phase,
                 color = as.factor(dummy_rep))) +
  #geom_line(aes(x= Days.Elapsed, y = predict(model), color= phase, group = Reactor)) + 
  #geom_line(aes(x= Days.Elapsed, y = predict(model2), color= phase, group = Reactor)) + 
  geom_line(aes(y = exp(predict(model)),
                color = as.factor(dummy_rep))) +
  #geom_line(aes(x= Days.Elapsed, y = exp(predict(model11)), color= phase, group = Reactor)) +
  #geom_line(aes(x= Days.Elapsed, y = predict(model5), color= phase, group = Reactor)) +
  facet_grid(biotype~Season)

summary(model)

(model)

########## PFOA ##########

chosen_PFAS <- "PFOA"

single_PFAS <- filter(summary_plot, PFAS == chosen_PFAS) 

single_PFAS$phase <- ordered(single_PFAS$phase, c("lag","peak","decline"))

#####fitting truncated distribution

fitting_values <- tibble(raw_measure = single_PFAS$raw_measure) %>%
  mutate(index1 = seq_len(nrow(.))) %>%
  sample_n(size = nrow(.)) %>%
  mutate(index2 = seq_len(nrow(.)))

lower_bound <- LOQs[chosen_PFAS]/10
sqrtloq <- round(LOQs[chosen_PFAS]/sqrt(2), 5)
cutoff <- LOQs[chosen_PFAS]

fitting_values <-  fitting_values %>%
  arrange(raw_measure) %>%
  mutate(raw_measure = round(raw_measure, 5)) %>%
  mutate(left = case_when(raw_measure == sqrtloq ~ lower_bound,
                          raw_measure == lower_bound ~ as.numeric(NA),
                          TRUE ~ raw_measure)) %>%
  mutate(right = case_when(raw_measure == lower_bound ~ lower_bound,
                           raw_measure == sqrtloq ~ sqrtloq,
                           TRUE ~ raw_measure)) %>%
  arrange(raw_measure)

censdata <- fitting_values[c("left","right")] %>% as.data.frame()

trunc_dist <- fitdistcens(log(censdata), distr = "norm")

trunc_quantiles <- unlist(quantile(trunc_dist, probs = seq(0,100,by=(100/(max(fitting_values$index1)+2)))/100)[1])[-1]

lowestcens <- sum(fitting_values$raw_measure == lower_bound)
lowercens <- sum(fitting_values$raw_measure == sqrtloq)

lowestcens_quatiles <- trunc_quantiles[1:lowestcens]
lowercens_quantiles <- trunc_quantiles[(lowestcens+1):(lowestcens+lowercens)]

inputions <- c(lowestcens_quatiles,lowercens_quantiles)

length(inputions) <- nrow(fitting_values)

inputed_frame <- cbind(fitting_values, inputions) %>%
  arrange(index1) %>%
  mutate(quantiles = ifelse(is.na(inputions), log(raw_measure), inputions)) %>%
  dplyr::select(quantiles) %>%
  cbind(single_PFAS) %>%
  mutate(conc = ifelse(raw_measure < cutoff, exp(quantiles) * dil, raw_measure*dil))
#####

model <- lmer(data = inputed_frame,
              log(conc) ~ Days.Elapsed + biotype + (Days.Elapsed|Reactor))

model2 <- lmer(data = inputed_frame,
              log(conc) ~ phase + biotype + (phase|Reactor))

model3 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + (1|Reactor))

model4 <- lmer(data = inputed_frame,
               log(conc) ~ Days.Elapsed + biotype + (1|Reactor))

model5 <- lmer(data = inputed_frame,
               log(conc) ~ Days.Elapsed + biotype + (1|Season) + (Days.Elapsed|Reactor))

model6 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + (1|Season) + (phase|Reactor))

model7 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (phase|Reactor))

model8 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (Days.Elapsed|Reactor))

model9 <- lmer(data = inputed_frame,
               log(conc) ~ phase + biotype + Season + (1|Reactor))

model10 <- lmer(data = inputed_frame,
                log(conc) ~ phase + biotype + (1|Reactor) + (Days.Elapsed|Reactor))

model11 <- lmer(data = inputed_frame,
                log(conc) ~ Days.Elapsed + I(Days.Elapsed^2) + biotype + (1|Season/Reactor))
                

anova(model,model2,model3,model4,model5,model6,model7,model8,model9,model10)

preds <- predict(model3, newdata = new.data, type = "response")

inputed_frame <- left_join(inputed_frame, codings)

ggplot(inputed_frame, aes(x = Days.Elapsed*10, group = Reactor)) +
  theme_classic()+
  guides(color = FALSE, fill =FALSE, linetype = FALSE)+
  geom_point(aes(y = conc,
                 fill = phase,
                 color = as.factor(dummy_rep))) +
  #geom_line(aes(x= Days.Elapsed, y = predict(model), color= phase, group = Reactor)) + 
  #geom_line(aes(x= Days.Elapsed, y = predict(model2), color= phase, group = Reactor)) + 
  geom_line(aes(y = exp(predict(model)),
                color = as.factor(dummy_rep))) +
  #geom_line(aes(x= Days.Elapsed, y = exp(predict(model11)), color= phase, group = Reactor)) +
  #geom_line(aes(x= Days.Elapsed, y = predict(model5), color= phase, group = Reactor)) +
  facet_grid(biotype~Season)

summary(model)

