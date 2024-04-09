LATITUDE # full analysis set

# IPW-adjustment
# calculate IPW
# all --- PPI use within 30 days before and after randomization
fit_wgt <- WeightIt::weightit(
  PPI_WITHIN_BIN ~ AGE_CAT + PS + CAD + GI + VD + PSA_CUT + BPI_CAT + GG_CAT + EOD_CAT + VISCERAL + HGB_CUT + LDH_CUT + ARM + BMA_WITHIN + GC_WITHIN,
  data = LATITUDE, method = "ps", estimand = "ATE")
wgt <- fit_wgt[["weights"]]
ps <- fit_wgt[["ps"]]
LATITUDE <- tibble(LATITUDE, wgt)
LATITUDE <- tibble(LATITUDE, ps)

# propensity score and weight distributions
library(ggplot2)
library(ggsci)
theme_set(theme_minimal() + theme(legend.position = "none"))

# propensity score distributions
psdist <- ggplot() + 
  geom_density(data = filter(LATITUDE, PPI_WITHIN == "yes"),
               aes(x = ps, weight = wgt, fill = PPI_WITHIN), alpha = 0.2) + 
  geom_density(data = filter(LATITUDE, PPI_WITHIN == "no"),
               aes(x = ps, weight = wgt, fill = PPI_WITHIN, y = -..density..), alpha = 0.2) + 
  geom_density(data = filter(LATITUDE, PPI_WITHIN == "yes"),
               aes(x = ps, fill = PPI_WITHIN), alpha = 0.6) + 
  geom_density(data = filter(LATITUDE, PPI_WITHIN == "no"),
               aes(x = ps, fill = PPI_WITHIN, y = -..density..), alpha = 0.6) +
  geom_hline(yintercept = 0, color = "white", size = 0.25) + 
  annotate(geom = "label", x = 0.2, y = 5, label = "PPI users\n(Pseudo-population)",
           fill = "#EE0000FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.5, y = 2.5, label = "PPI users\n(actual)",
           fill = "#EE0000FF", alpha = 0.6, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.15, y = -5, label = "Non-users\n(Pseudo-population)",
           fill = "#3B4992FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = -2.5, label = "Non-users\n(actual)",
           fill = "#3B4992FF", alpha = 0.6, color = "white", hjust = 0) +
  scale_fill_aaas() + 
  scale_y_continuous(label = abs) + 
  labs(x = "Propensity score", y = "Density")

# weight distributions
wgtdist <- ggplot() + 
  geom_histogram(data = filter(LATITUDE, PPI_WITHIN == "yes"),
                 bins = 30, aes(x = wgt, fill = PPI_WITHIN), alpha = 0.3) + 
  geom_histogram(data = filter(LATITUDE, PPI_WITHIN == "no"),
                 bins = 30, aes(x =  wgt, fill = PPI_WITHIN, y = -..count..), alpha = 0.3) + 
  geom_hline(yintercept = 0, color = "white", size = 0.25) + 
  annotate(geom = "label", x = 0.5, y = 0.1, label = "PPI users",
           fill = "#EE0000FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.15, y = -1, label = "Non-users",
           fill = "#3B4992FF", alpha = 0.3, color = "white", hjust = 0) +
  scale_fill_aaas() + 
  scale_y_continuous(label = abs) + 
  labs(x = "Count", y = "Inverse probability weight")

wgtdist <- LATITUDE %>%
  ggplot() +
  geom_histogram(aes(x = wgt, y = ..count..), fill = "#3B4992FF", stat = "bin", bins = 20) + 
  labs(
    x = "Weight", y = "Count",
    title = "Distribution of weights (LATITUDE dataset)"
  )

# Patient characteristics --- PPI use within 30 days before and after randomization
library(tableone)
library(survey)
fct_var <- c("AGE_CAT", "PS", "CAD", "GI", "VD", "PSA_CUT", "BPI_CAT", "GG_CAT", "EOD_CAT", "VISCERAL", "HGB_CUT", "LDH_CUT", "BMA_WITHIN", "GC_WITHIN", "ARM")
all_var <- c("AGE_CAT", "PS", "CAD", "GI", "VD", "PSA_CUT", "BPI_CAT", "GG_CAT", "EOD_CAT", "VISCERAL", "HGB_CUT", "LDH_CUT", "BMA_WITHIN", "GC_WITHIN", "ARM")

# table one
tbl_base <- CreateTableOne(vars = all_var, factorVars = fct_var, strata = "PPI_WITHIN", data = LATITUDE, test = FALSE)
tbl <- print(tbl_base, catDigits = 1, contDigits = 1, explain = TRUE, varLabels = TRUE, noSpaces = TRUE, 
             showAllLevels = TRUE, smd = TRUE)

svyLATITUDE <- svydesign(ids = ~ USUBJID, weights = ~ wgt, data = LATITUDE)
wtbl_base <- svyCreateTableOne(vars = all_var, factorVars = fct_var, strata = "PPI_WITHIN", data = svyLATITUDE, test = FALSE)
wtbl <- print(wtbl_base, catDigits = 1, contDigits = 1, explain = TRUE, varLabels = TRUE, noSpaces = TRUE, 
              showAllLevels = TRUE, smd = TRUE)

# kaplan-meier estimations
pfs_fit <- survfit(Surv(TTPROG, PROG_BIN) ~ PPI_WITHIN, data = LATITUDE, weights = wgt)
pfs_plot_mcspc <- ggsurvplot(pfs_fit, 
                             risk.table = FALSE, font.tickslab = 12, font.x = 12, font.y = 12,
                             title = "Radiographic progression-free survival",
                             legend = "top", legend.labs = c("no", "yes"),
                             legend.title = "PPI use", 
                             censor = FALSE, 
                             xlab = "Time from randomization, days",
                             ylab = "Patients without radiographic progression or death",
                             palette = "aaas", size = 0.6,
                             ggtheme = theme_classic())

os_fit <- survfit(Surv(FOLLOWUP, DEATH_BIN) ~ PPI_WITHIN, data = LATITUDE, weights = wgt)
os_plot_mcspc <- ggsurvplot(os_fit, 
                            risk.table = FALSE, font.tickslab = 12, font.x = 12, font.y = 12,
                            title = "Overall survival",
                            legend = "top", legend.labs = c("no", "yes"),
                            legend.title = "PPI use",
                            censor = FALSE, 
                            xlab = "Time from randomization, days",
                            ylab = "Patients who survived",
                            palette = "aaas", size = 0.6,
                            ggtheme = theme_classic())

# weighted log-rank test
# progression-free survival
pfs_cox <- coxph(Surv(TTPROG, PROG_BIN) ~ PPI_WITHIN, data = LATITUDE, weights = wgt)
pfs_pval <- round(summary(pfs_cox)$sctest[3], digits = 6)

# Overall survival
os_cox <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ PPI_WITHIN, data = LATITUDE, weights = wgt)
os_pval <- round(summary(os_cox)$sctest[3], digits = 10)

# interaction analysis within Cox regression model
# OS
library(finalfit)

# interaction analysis for PFS
int_pfs_age <- coxph(Surv(TTPROG, PROG_BIN) ~ AGE_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_ps <- coxph(Surv(TTPROG, PROG_BIN) ~ PS * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_cad <- coxph(Surv(TTPROG, PROG_BIN) ~ CAD * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_gi <- coxph(Surv(TTPROG, PROG_BIN) ~ GI * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_vd <- coxph(Surv(TTPROG, PROG_BIN) ~ VD * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_psa <- coxph(Surv(TTPROG, PROG_BIN) ~ PSA_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_bpi <- coxph(Surv(TTPROG, PROG_BIN) ~ BPI_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_gg <- coxph(Surv(TTPROG, PROG_BIN) ~ GG_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_eod <- coxph(Surv(TTPROG, PROG_BIN) ~ EOD_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_vis <- coxph(Surv(TTPROG, PROG_BIN) ~ VISCERAL * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_hgb <- coxph(Surv(TTPROG, PROG_BIN) ~ HGB_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_ldh <- coxph(Surv(TTPROG, PROG_BIN) ~ LDH_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_bma <- coxph(Surv(TTPROG, PROG_BIN) ~ BMA_WITHIN * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_gc <- coxph(Surv(TTPROG, PROG_BIN) ~ GC_WITHIN * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_pfs_arm <- coxph(Surv(TTPROG, PROG_BIN) ~ ARM * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()

int_pfs_sum <- bind_rows(int_pfs_age, int_pfs_ps, int_pfs_cad, int_pfs_gi, int_pfs_vd,
                         int_pfs_psa, int_pfs_bpi, int_pfs_gg, int_pfs_eod, int_pfs_vis,
                         int_pfs_hgb, int_pfs_ldh, int_pfs_bma, int_pfs_gc, int_pfs_arm)

# interaction analysis for OS
int_os_age <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ AGE_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_ps <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ PS * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_cad <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ CAD * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_gi <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ GI * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_vd <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ VD * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_psa <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ PSA_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_bpi <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ BPI_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_gg <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ GG_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_eod <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ EOD_CAT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_vis <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ VISCERAL * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_hgb <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ HGB_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_ldh <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ LDH_CUT * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_bma <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ BMA_WITHIN * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_gc <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ GC_WITHIN * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()
int_os_arm <- coxph(Surv(FOLLOWUP, DEATH_BIN) ~ ARM * PPI_WITHIN, data = LATITUDE, weights = wgt) %>% fit2df()

int_os_sum <- bind_rows(int_os_age, int_os_ps, int_os_cad, int_os_gi, int_os_vd,
                        int_os_psa, int_os_bpi, int_os_gg, int_os_eod, int_os_vis,
                        int_os_hgb, int_os_ldh, int_os_bma, int_os_gc, int_os_arm)

