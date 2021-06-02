library(pwr)

# For a one-way ANOVA comparing 3 groups of equal n, calculate the
# sample size needed in each group to obtain a power of
# 0.80, when the effect size is "f" and a
# significance level of 0.05 is employed.

##Calculate effect size (f) from gene expression data

#Log transform Fold Change values following BCG stim
mFC_HU <- mean(log(230), log(40.5), log(39.6), log(32.3), log(20.6), log(10.7))
mFC_HE <- mean(log(140), log(25.7), log(58), log(15.8), log(11.4), log(6.53))
stdFC_HU <- sd(c(log(230), log(40.5), log(39.6), log(32.3), log(20.6), log(10.7)))

effSize <- (mFC_HU - mFC_HE)/stdFC_HU


##Calculate Sample Size using effect size based on gene expression data

#
pwr.anova.test(k=3, f=effSize, sig.level=.05, power = 0.8)
pwr.t.test(n = 30 , d = effSize , sig.level = , power = , type = c("two.sample", "one.sample", "paired"))