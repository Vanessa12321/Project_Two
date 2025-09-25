library(tidyverse)
library(reactable) 
library(glue)
library(lme4) 
library(boot) 

data = read.csv("running_surface_data.csv") 
data$subID = as.factor(data$subID)
data$surfaceID = as.factor(data$surfaceID) 

# recreating figure 3 
axial_acc_mean = data %>% 
  group_by(surfaceID) %>% 
  summarize(
    mean_axial = mean(TibiaAcc_Axial), 
    sd_axial = sd(TibiaAcc_Axial)
  ) %>% 
  mutate(
    lower95 = mean_axial - 2 * sd_axial, 
    upper95 = mean_axial + 2 * sd_axial, 
    surfaceID_num = as.numeric(surfaceID)
  ) %>% 
  rename(TibiaAcc_Axial = mean_axial)

ggplot(data, aes(x=surfaceID, y=TibiaAcc_Axial)) + 
  geom_point(shape = 18, size=2, aes(color=subID), show.legend = F) + 
  geom_line(aes(group = subID, color=subID), show.legend = F) + 
  geom_ribbon(data = axial_acc_mean,
              aes(ymin = lower95, ymax = upper95, x = surfaceID_num),
              fill = "#717171", alpha = 0.2, inherit.aes = F) +
  geom_point(data = axial_acc_mean) + 
  geom_line(data = axial_acc_mean, group=1) + 
  ylim(0, 15) + 
  theme_classic() + 
  labs(
    x = NULL, 
    y = "Acceleration [g]"
  )
 
# recreating figure 4 
res_acc_mean = data %>% 
  group_by(surfaceID) %>% 
  summarize(
    mean_resultant = mean(TibiaAcc_Resultant), 
    sd_resultant = sd(TibiaAcc_Resultant)
  ) %>% 
  mutate(
    lower95 = mean_resultant - 2 * sd_resultant, 
    upper95 = mean_resultant + 2 * sd_resultant, 
    surfaceID_num = as.numeric(surfaceID)
  ) %>% 
  rename(TibiaAcc_Resultant = mean_resultant)

ggplot(data, aes(x=surfaceID, y=TibiaAcc_Resultant)) + 
  geom_point(shape = 18, size=2, aes(color=subID), show.legend = F) + 
  geom_line(aes(group = subID, color=subID), show.legend = F) + 
  geom_ribbon(data = res_acc_mean,
              aes(ymin = lower95, ymax = upper95, x = surfaceID_num),
              fill = "#717171", alpha = 0.2, inherit.aes = F) +
  geom_point(data = res_acc_mean) + 
  geom_line(data = res_acc_mean, group=1) + 
  ylim(0,25) + 
  theme_classic() + 
  labs(
    x = NULL, 
    y = "Acceleration [g]"
  ) 

# recreating table 1 
create_table_1 = function(data) {
  cols = colnames(
    data %>% 
      select(StrideRate, TibiaAcc_Break, TibiaAcc_Propulsion, TibiaAcc_Axial, 
             TibiaAcc_Medial, TibiaAcc_Lateral, TibiaAcc_Resultant, ShockAtten, ShockAtten_Resultant)
  )
  result = tibble() 
  for (i in 1:length(cols)) { 
    # things needed for first four columns 
    response = cols[i]
    mean_sd = data %>% 
      group_by(surfaceID) %>% 
      summarize(
        mean = mean(!!sym(response)), 
        sd = sd(!!sym(response)) 
      ) 
    dirt_mean = mean_sd[1,2]
    dirt_sd = mean_sd[1,3] 
    gravel_mean = mean_sd[2,2]
    gravel_sd = mean_sd[2,3]
    paved_mean = mean_sd[3,2]
    paved_sd = mean_sd[3,3]

    full_formula = as.formula(glue("{response} ~ surfaceID")) 
    reduced_formula = as.formula(glue("{response} ~ 1"))

    # things need for last two columns 
    full = lm(full_formula, data = data) 
    reduced = lm(reduced_formula, data = data) 
    anova = anova(full, reduced) 
    es = round(1 - (anova$RSS[1] / anova$RSS[2]), 2) 
    p = round(anova$`Pr(>F)`[2], 2) 

    entry = tibble(
      Variable = response, 
      Dirt = sprintf("%.2f \u00B1 %.2f", dirt_mean, dirt_sd), 
      Gravel = sprintf("%.2f \u00B1 %.2f", gravel_mean, gravel_sd), 
      Paved = sprintf("%.2f \u00B1 %.2f", paved_mean, paved_sd), 
      ES = es, 
      p_value = p 
    )

    result = bind_rows(result, entry)
  }
  return(result)
}

table_1 = create_table_1(data)

# modify ES values under 0.01 
table_1 = table_1 %>% 
  mutate(
    ES = as.character(ES), 
    ES = case_when(
      ES == "0" ~ "<0.01", 
      .default = "0.01"
    )
  )

reactable(table_1)

no_runner = lm(TibiaAcc_Axial ~ surfaceID, data = data) 
no_runner_tidy = broom::tidy(summary(no_runner)) %>% mutate_if(is.numeric, round, digits=2)
reactable(no_runner_tidy)
with_runner = lm(TibiaAcc_Axial ~ surfaceID + subID, data = data) 
with_runner_tidy = broom::tidy(summary(with_runner)) %>% mutate_if(is.numeric, round, digits=2)
reactable(with_runner_tidy) 


# bootstrapping confidence intervals 
boot_function = function(data, indices, y_var) {
  formula = as.formula(glue("{y_var} ~ surfaceID + subID"))
  d = data[indices, ] 
  m = lm(formula, data = d) 
  return(coef(m))
} 

# first look at axial acceleration 
axial_acc_results = boot(data = data, statistic = boot_function, 
  strata = data[,"subID"], R = 1000, y_var="TibiaAcc_Axial") 
axial_gravel_sims = data.frame(sims = axial_acc_results$t[,2]) 
axial_paved_sims = data.frame(sims = axial_acc_results$t[,3]) 
axial_gravel_ci = boot.ci(axial_acc_results, type = "bca", index = 2)
axial_paved_ci = boot.ci(axial_acc_results, type = "bca", index = 3) 

# gravel surface 
g1 = ggplot(axial_gravel_sims, aes(x=sims)) + 
  geom_histogram(bins=20, fill="skyblue") + 
  geom_vline(xintercept = axial_gravel_ci$bca[,4], linetype="dashed", color="red") + 
  geom_vline(xintercept = axial_gravel_ci$bca[,5], linetype="dashed", color="red") + 
  theme_bw() + 
  labs(
    x = "Regression Coefficient", 
    y = "Count", 
    title = "1000 Bootstrap Samples of\nGravel Surface Coefficient", 
    subtitle = "95% CI = (0.129, 0.791)"
  )

# paved surface 
g2 = ggplot(axial_paved_sims, aes(x=sims)) + 
  geom_histogram(bins=20, fill="orange") + 
  geom_vline(xintercept = axial_paved_ci$bca[,4], linetype="dashed", color="red") + 
  geom_vline(xintercept = axial_paved_ci$bca[,5], linetype="dashed", color="red") + 
  theme_bw() + 
  labs(
    x = "Regression Coefficient", 
    y = "Count", 
    title = "1000 Bootstrap Samples of\nPaved Surface Coefficient", 
    subtitle = "95% CI = (0.134, 0.842)"
  )

gridExtra::grid.arrange(g1, g2, ncol=2) 

# then look at resultant acceleration
res_acc_results = boot(data = data, statistic = boot_function, 
  strata = data[,"subID"], R = 1000, y_var="TibiaAcc_Resultant") 
res_gravel_sims = data.frame(sims = res_acc_results$t[,2]) 
res_paved_sims = data.frame(sims = res_acc_results$t[,3]) 
res_gravel_ci = boot.ci(res_acc_results, type = "bca", index = 2)
res_paved_ci = boot.ci(res_acc_results, type = "bca", index = 3) 

# same thing as before - gravel then paved 
g3 = ggplot(res_gravel_sims, aes(x=sims)) + 
  geom_histogram(bins=20, fill="skyblue") + 
  geom_vline(xintercept = res_gravel_ci$bca[,4], linetype="dashed", color="red") + 
  geom_vline(xintercept = res_gravel_ci$bca[,5], linetype="dashed", color="red") + 
  theme_bw() + 
  labs(
    x = "Regression Coefficient", 
    y = "Count", 
    title = "1000 Bootstrap Samples of\nGravel Surface Coefficient", 
    subtitle = "95% CI = (-0.398, 0.894)"
  )

# paved surface 
g4 = ggplot(res_paved_sims, aes(x=sims)) + 
  geom_histogram(bins=20, fill="orange") + 
  geom_vline(xintercept = res_paved_ci$bca[,4], linetype="dashed", color="red") + 
  geom_vline(xintercept = res_paved_ci$bca[,5], linetype="dashed", color="red") + 
  theme_bw() + 
  labs(
    x = "Regression Coefficient", 
    y = "Count", 
    title = "1000 Bootstrap Samples of\nPaved Surface Coefficient", 
    subtitle = "95% CI = (0.107, 1.421)"
  )

gridExtra::grid.arrange(g3, g4, ncol=2)
