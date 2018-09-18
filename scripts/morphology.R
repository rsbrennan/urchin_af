### OASV2 body length (morphometrics) analysis
### April D. Garrett

### MORPHOMETRICS: BODY + ARM ###
pluteus <- read.csv("~/Desktop/Brennan_Garrett_Manuscript/OASV2_ReidGarrettPaper_ForR.csv")
attach(pluteus)

### RUNNING LINEAR MODEL WITH VESSEL_ID AS RANDOM EFFECT ###
library(lmerTest)
p_value <- lmer(Size~pH+(1|Vessel_ID), data=pluteus, family=("gaussian"))
summary(p_value)
anova(p_value)
