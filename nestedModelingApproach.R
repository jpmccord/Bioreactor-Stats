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

write_csv(mutate(summaries[subset2], Summary_Stat = rownames(summaries)), "Modeled_PFAS.csv")

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
  facet_wrap(~PFAS) +
  xlab("Days Elapsed") +
  ylab("Measured Concentration (ng/L)")

inpute_data <- function(data) {
  
  chosen_PFAS <- unique(data$PFAS)
  
  levels(data$phase) <- c("lag", "peak","decline")
  
  #####fitting truncated distribution
  fitting_values <- tibble(raw_measure = data$raw_measure) %>%
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
    cbind(data) %>%
    mutate(conc = ifelse(raw_measure < cutoff, exp(quantiles) * dil, raw_measure*dil),
           Days.Elapsed.Trans = Days.Elapsed/10)

  return(inputed_frame)
}

model_PFAS <- function(data) {
  model <- lmerTest::lmer(data = data, log(conc) ~ Days.Elapsed.Trans + biotype + Season + (Days.Elapsed.Trans|Reactor))
  
  return(model)
}

plot_results <- function(rownum) {

inputed_frame <- left_join(nested_analysis$inputed_data[[rownum]], codings)

plot <- ggplot(inputed_frame, aes(x = Days.Elapsed, group = Reactor)) +
  theme_classic()+
  guides(color = FALSE, fill =FALSE, linetype = FALSE)+
  geom_point(aes(y = conc,
                 fill = phase,
                 color = as.factor(dummy_rep))) +
  geom_line(aes(y = exp(predict(nested_analysis$model[[rownum]])),
                color = as.factor(dummy_rep))) +
  facet_grid(biotype~Season) +
  ggtitle(paste("Modeled Concentrations for",nested_analysis$PFAS_holder[rownum]))+
  xlab("Days Elapsed") +
  ylab("Concentration (ng/L)") +
  theme(plot.title = element_text(hjust=0.5))

return(plot)

}

nested_analysis <- summary_plot %>%
  mutate(PFAS_holder = PFAS) %>%
  group_by(PFAS_holder) %>%
  nest() %>%
  mutate(inputed_data = map(.x = data, .f = inpute_data),
         model = map(.x = inputed_data, .f = model_PFAS),
         model_sum = map(.x = model, .f = summary))


pdf("Individual_Models.pdf", paper = "USr")
for(rownum in 1:length(nested_analysis$PFAS_holder)) {
  print(plot_results(rownum))
}
dev.off()

grabSummary <- function(model) {
  coefficients <- as.tibble(model$coefficients) %>%
    mutate(vars = rownames(model$coefficients))
}

grabRans <- function(model) {
    effs <- ranef(model) %>%
      mutate(vars = rownames(model$coefficients))
}

coef_summary <- nested_analysis %>%
  mutate(coefs = map(.x= model_sum, .f= grabSummary)) %>%
  unnest(coefs) 


formatted_table <- select(coef_summary, PFAS_holder, Estimate, vars) %>%
  spread(PFAS_holder, Estimate) %>%
  mutate_at(.vars = vars(-vars), signif, digits = 1) %>%
  t() %>%
  as_tibble() %>%
  slice(-1) %>%
  setNames(unique(coef_summary$vars)) %>%
  mutate(Days.Elapsed.Trans = as.numeric(Days.Elapsed.Trans)/10) %>%
  rename(Days.Elapsed = Days.Elapsed.Trans) %>%
  t() %>% t()

formatted_pval <- select(coef_summary, PFAS_holder, `Pr(>|t|)`, vars) %>%
  spread(PFAS_holder, `Pr(>|t|)`) %>%
  mutate_at(.vars = vars(-vars), signif, digits = 1) %>%
  t() %>%
  as_tibble() %>%
  slice(-1) %>%
  setNames(unique(coef_summary$vars)) %>%
  t() %>% t()


for_paper = matrix(paste0(formatted_table,"(",formatted_pval,")"), ncol = 6) %>%
  as.tibble()

colnames(for_paper) = colnames(formatted_table)

for_paper = mutate(for_paper, PFAS = nested_analysis$PFAS_holder)

write_csv(for_paper, "model_estimates.csv")

