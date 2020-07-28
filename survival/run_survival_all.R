library(survival)
#library(ggplot2)
#library(GGally)
#library(parallel)
library(powerSurvEpi)


args = commandArgs(trailingOnly=TRUE)
input_survival_file = args[1]
input_bin_phe_file = args[2]
output_file = args[3]
adjust_BMI = as.logical(args[4])
adjust_SEX = as.logical(args[5])

# input_survival_file = "~/Kimlab/pan_cancer/revision/TCGA_data/phenotype/TCGA_survival.csv"
# input_bin_phe_file = "60k_biobin-phe-bins.csv"
# adjust_BMI = F
# adjust_SEX = F

options(warn=-1)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dem = read.csv(file = input_survival_file, stringsAsFactors = F)
dem = dem[!is.na(dem$VITAL_STATUS),]
dem$SEQN_ID = paste(dem$SEQN_ID, dem$SEQN_ID, sep = '_')

bbin = read.csv(input_bin_phe_file, stringsAsFactors = F)
cox_regression_min_sample_var = 10

bbin = bbin[-c(1:9),]
bbin = bbin[bbin[,2] == 1,]

cov_bbin = merge(x = dem, y = bbin, by.y = "ID", by.x = "SEQN_ID")
input = cov_bbin[,8:(ncol(cov_bbin))]
input = data.frame(apply(input, 2, as.numeric))
notzeros = apply(input, 2, sum)
notzeros = notzeros > 0
input = input[,notzeros]

run_kp = function(name, cox_pval) {
  x = input[,which(colnames(input) == name)]
  x[which(x > 0)] = 1
  x[which(x == 0)] = 0
  surv = Surv(cov_bbin$SURVIVAL_DAYS * 0.0328767, cov_bbin$VITAL_STATUS == 1)
  sdiff = survdiff(surv ~ x, data = cov_bbin)
  fit = survfit(surv ~ x, data = cov_bbin)
  p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  s = ggsurv(fit, cens.col = gg_color_hue(2)) + guides(linetype = FALSE)+ scale_colour_discrete(name = NULL, breaks = c(0,1),labels = c(paste('No Rare variant: ', fit$n[1], sep =''), paste('Rare variant: ', fit$n[2], sep ='')))
  s = s + theme(legend.justification = c(0,0), legend.position=c(0, 0))
  s = s + annotate("text", x=max(fit$time)-(1/6)*max(fit$time), y=0.95, label= paste(paste('cox p-value: ', format(cox_pval, scientific=T), sep='')))
  s = s + xlab('Time(months)') + ggtitle(name)
  png(paste('out/', name, '.png', sep=''), width=7, height=5, units="in", res=200, type = "cairo")
  plot(s)
  dev.off()
}

run_coxreg = function(x) {
  if (sum(x > 0) >= cox_regression_min_sample_var) {
    formula = "Surv(SURVIVAL_DAYS, VITAL_STATUS == 1) ~ x + AGE"
    if (adjust_SEX) {
      formula = paste(formula, " + SEX", sep = "")
    }
    if (adjust_BMI) {
      formula = paste(formula, " + BMI", sep = "")
    }
    cph = coxph(as.formula(formula), data = cov_bbin)
    cph = summary(cph)
    return(cph$coefficients[1, 5])
  } 
  return(NA)
}

coxreg_power = function(x) {
  if (sum(x > 0) >= cox_regression_min_sample_var) {
    formula = "Surv(SURVIVAL_DAYS, VITAL_STATUS == 1) ~ x + AGE"
    s_formula = "x ~ AGE"
    if (adjust_SEX) {
      formula = paste(formula, " + SEX", sep = "")
      s_formula = paste(s_formula, " + SEX", sep = "")
    }
    if (adjust_BMI) {
      formula = paste(formula, " + BMI", sep = "")
      s_formula = paste(s_formula, " + BMI", sep = "")
    }
    cph = coxph(as.formula(formula), data = cov_bbin)
    s = summary(lm(as.formula(s_formula), data = cov_bbin))
    cph = summary(cph)
    p = powerEpiCont.default(length(x), cph$coefficients[1, 2], var(x), sum(cov_bbin$VITAL_STATUS),  s$r.squared, 0.05)
    return(p)
  }
  return(NA)
}

cat("Start cox regression\n")
obj2 = apply(input, 2, run_coxreg)
power = apply(input, 2 , coxreg_power)
pval_matrix = data.frame(colnames(input), obj2, power)
colnames(pval_matrix) = c("gene", "p_value", "power")
pval_matrix$fdr = p.adjust(pval_matrix$p_value, method = "fdr", n = sum(!is.na(pval_matrix$p_value)))
pval_matrix$bonferroni = p.adjust(pval_matrix$p_value, method = "bonferroni", n = sum(!is.na(pval_matrix$p_value)))
pval_matrix = pval_matrix[order(pval_matrix$p_value),]
write.csv(pval_matrix, file=output_file, row.names = F, quote = F)
cat("End cox regression\n")
