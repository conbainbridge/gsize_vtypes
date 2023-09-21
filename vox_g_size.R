###########################################
### Group vocalization size estimations ###
########### Constance Bainbridge ##########
###########################################

setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), sep=""))

#####################
### Initial setup ###
#####################

library(stringr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library('lme4')
library('pbkrtest')
library(emmeans)
library(ggformula)
library(caret) # Confusion matrix
library(psycho) # Signal detection theory indices

# Acquire sensR from archives
library(remotes)
install_version("sensR")
4 # Select correct version
library(sensR) # ROC curves


###
############################
### Study 1 - Numerosity ###
############################
###
numerosity_data <- read.csv(file = "g_numerosity_data.csv")

####################
### Demographics ###
####################
### Sample size
# Drop incomplete
numerosity_data = numerosity_data[numerosity_data$subjGend != "None",]

n = length(unique(numerosity_data$id))
n

# Female
n_fem = length(unique(numerosity_data$id[numerosity_data$subjGend == "Female"]))
n_fem

# Male
n_male = length(unique(numerosity_data$id[numerosity_data$subjGend == "Male"]))
n_male

# Other
n_other = length(unique(numerosity_data$id[numerosity_data$subjGend == "Other"]))
n_other

### Age
ages = numerosity_data$age[match(unique(numerosity_data$id),numerosity_data$id)]
ages = strtoi(ages)
mean(ages)
sd(ages)

### Lm stimulus gender
# Create linear model
gendList <- list(g_size=seq(1,6,by=1), gend=c("female", "male"))
modalityList <- list(g_size=seq(1,6,by=1), modality=c("laughing", "sentence", "talking", "yelling"))

stimGendModel <- lm(user_g_size~g_size*stimGend, data = numerosity_data)
summary(stimGendModel)

# Gender over modality
stimGendModeModel <- lm(user_g_size~modality*stimGend, data = numerosity_data)
summary(stimGendModeModel)

# Linear plot stim gender and group size
stimGendPlot <- emmip(stimGendModel, stimGend ~ g_size, at=gendList,CIs=TRUE) +
  labs(x = 'Vocal modality', y = 'User estimates', title = 'Estimates')
stimGendPlot

# Linear plot stim gender and modality - difference alluded to in discussion section of paper
stimGendModePlot <- emmip(stimGendModeModel, stimGend ~ modality, at=modalityList,CIs=TRUE)
stimGendModePlot

#############
### GLMMs ###
#############
lmmDF <- data.frame(id = numerosity_data$id, modality = numerosity_data$modality, g_size = numerosity_data$g_size, user_g_size = numerosity_data$user_g_size, stimGend = numerosity_data$stimGend, subjGend = numerosity_data$subjGend, diffScore =  numerosity_data$diffScore, trialIndex = numerosity_data$trialIndex)

### Build models
m1 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size, data=lmmDF, REML=FALSE)
m2 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality, data=lmmDF, REML=FALSE)
m3 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + g_size:modality, data=lmmDF, REML=FALSE)
m4 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + (1|stimGend), data=lmmDF, REML=FALSE)
m5 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + g_size:modality + (1|stimGend), data=lmmDF, REML=FALSE)
m6 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + (1|stimGend) + (1|subjGend), data=lmmDF, REML=FALSE)
m7 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + (1|stimGend) + (1|subjGend) + (1|stimGend:subjGend), data=lmmDF, REML=FALSE)  
m8 <- lmer(diffScore ~ (1|id) + (1|trialIndex) + g_size + modality + (1|stimGend) + (1|subjGend) + (1|stimGend:g_size), data=lmmDF, REML=FALSE)  

# Look to compare the models
summary(m1)
summary(m2)
summary(m3) # 2nd lowest BIC scores
summary(m4)
summary(m5) # Lowest AIC and BIC scores
summary(m6)
summary(m7)
summary(m8)

# Calculate difference between second lowest and lowest BIC, m3 - m5)
11884.0-11877.5

# get the KR-approximated degrees of freedom
dfm5.KR <- get_Lb_ddf(m5, fixef(m5))

# extract coefficients
coefsm5 <- data.frame(coef(summary(m5)))

# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefsm5$p.KR <- 2 * (1 - pt(abs(coefsm5$t.value), dfm5.KR))
coefsm5

# use normal distribution to approximate p-value
coefsm5$p.z <- 2 * (1 - pnorm(abs(coefsm5$t.value)))
coefsm5

### Plots
# Model
mEstimates <- lmer(user_g_size ~ (1|id) + (1|trialIndex) + g_size + modality + g_size:modality + (1|stimGend), data=lmmDF, REML=FALSE)
summary(mEstimates)

### Fig. 1
estimatesData <- emmip(mEstimates, modality ~ g_size, at=modalityList,CIs=TRUE, plotit=FALSE)

estimatesGGPlot <- ggplot(data=estimatesData, aes(x=g_size, y=yvar, color=modality)) +
  geom_line() +
  geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=modality), alpha=0.4) +
  geom_line(linetype = "dashed", aes(x = estimatesData$g_size, y = estimatesData$g_size), color = "black") +
  theme(text = element_text(size = 15)) +
  labs(x = 'Actual group size', y = 'Participant estimates', color = "Vocalization type", fill = "Vocalization type")
estimatesGGPlot

# Figure 1
tiff("fig1.tiff", units="in", width=6.3, height=5, res=300)
estimatesGGPlot
dev.off()


###
###############################
### Study 2 - Forced-choice ###
###############################
###
fc_trials <- read.csv(file = "g_fc_data.csv")

############################
### FC accuracy measures ###
############################
### Age data
fc_age_data = strtoi(fc_trials$age[!is.na(fc_trials$age)])
mean(fc_age_data)
sd(fc_age_data)

### Sample size
n = length(unique(fc_trials$id))
n

# Female
n_fem = length(na.omit(unique(fc_trials$id[fc_trials$subjGend == "Female"])))
n_fem

# Male
n_male = length(na.omit(unique(fc_trials$id[fc_trials$subjGend == "Male"])))
n_male

# Other
n_other = length(na.omit(unique(fc_trials$id[fc_trials$subjGend == "Other"])))
n_other

# Get proportions of correct / correct + incorrect
SUMmodecombo <- fc_trials %>%
  group_by(modecombo, bigger, smaller, correctDUMMY, g_size_bigger, g_size_smaller) %>%
  summarize(n())

SUMmodecombo <- spread(SUMmodecombo,
                        key = "correctDUMMY",
                        value = "n()")
SUMmodecombo$prop <- (SUMmodecombo$`1` / (SUMmodecombo$`1`+SUMmodecombo$`0`))
SUMmodecombo

# Re-code which mode is bigger/smaller
for (i in 1:length(SUMmodecombo$modecombo)) {
  if (SUMmodecombo$modecombo[i] == "laugh>sentence"){
    SUMmodecombo$bigger[i] <- "laughing"
    SUMmodecombo$smaller[i] <- "sentence"
  } else if (SUMmodecombo$modecombo[i] == "laugh>talking"){
    SUMmodecombo$bigger[i] <- "laughing"
    SUMmodecombo$smaller[i] <- "talking"
  } else if (SUMmodecombo$modecombo[i] == "laugh>yell"){
    SUMmodecombo$bigger[i] <- "laughing"
    SUMmodecombo$smaller[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "laugh<sentence"){
    SUMmodecombo$smaller[i] <- "laughing"
    SUMmodecombo$bigger[i] <- "sentence"
  } else if (SUMmodecombo$modecombo[i] == "laugh<talking"){
    SUMmodecombo$smaller[i] <- "laughing"
    SUMmodecombo$bigger[i] <- "talking"
  } else if (SUMmodecombo$modecombo[i] == "laugh<yell"){
    SUMmodecombo$smaller[i] <- "laughing"
    SUMmodecombo$bigger[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "sentence>talking"){
    SUMmodecombo$bigger[i] <- "sentence"
    SUMmodecombo$smaller[i] <- "talking"
  } else if (SUMmodecombo$modecombo[i] == "sentence>yell"){
    SUMmodecombo$bigger[i] <- "sentence"
    SUMmodecombo$smaller[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "sentence<talking"){
    SUMmodecombo$smaller[i] <- "sentence"
    SUMmodecombo$bigger[i] <- "talking"
  } else if (SUMmodecombo$modecombo[i] == "sentence<yell"){
    SUMmodecombo$smaller[i] <- "sentence"
    SUMmodecombo$bigger[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "talking>yell"){
    SUMmodecombo$bigger[i] <- "talking"
    SUMmodecombo$smaller[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "talking<yell"){
    SUMmodecombo$smaller[i] <- "talking"
    SUMmodecombo$bigger[i] <- "yelling"
  } else if (SUMmodecombo$modecombo[i] == "laughlaugh"){
    SUMmodecombo$smaller[i] <- "laughing"
    SUMmodecombo$bigger[i] <- "laughing"
  } else if (SUMmodecombo$modecombo[i] == "sentencesentence"){
    SUMmodecombo$smaller[i] <- "sentence"
    SUMmodecombo$bigger[i] <- "sentence"
  } else if (SUMmodecombo$modecombo[i] == "talkingtalking"){
    SUMmodecombo$smaller[i] <- "talking"
    SUMmodecombo$bigger[i] <- "talking"
  } else if (SUMmodecombo$modecombo[i] == "yellyell"){
    SUMmodecombo$smaller[i] <- "yelling"
    SUMmodecombo$bigger[i] <- "yelling"
  } else {
    SUMmodecombo$bigger[i] <- NA
    SUMmodecombo$smaller[i] <- NA
  }
}

SUMCmeans <- SUMmodecombo %>%
  group_by(bigger, smaller) %>%
  summarize(prop = mean(prop))

### Heatmaps
FC_heat = ggplot(data =  SUMCmeans, mapping = aes(x = smaller, y = bigger, fill = prop)) +
  geom_tile(aes(fill = prop), colour = "white") +
  scale_fill_gradient(high = "red", low = "yellow") +
  geom_raster(aes(fill = prop)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size = 15), axis.text.y = element_text(size = 15)) +
  geom_text(aes(label = round(prop, digits = 4)), size = 5) +
  labs(fill = "Accuracy") +
  theme(text = element_text(size = 16), axis.title = element_text(size = 16))
FC_heat

tiff("fig2.tiff", units="in", width=6.3, height=5, res=300)
FC_heat
dev.off()

### Exploritory - Modality + numerosity
SUMmodecomboNum <- fc_trials %>%
  group_by(modecombo, g_size1, g_size2, correctDUMMY) %>%
  summarize(n())

# Initialize
SUMmodecomboNum <- spread(SUMmodecomboNum,
                       key = "correctDUMMY",
                       value = "n()")
SUMmodecomboNum$prop <- NA
SUMmodecomboNum$stimN <- (SUMmodecomboNum$`0` + SUMmodecomboNum$`1`)

# Fix NAs from none wrong
SUMmodecomboNum$`0`[is.na(SUMmodecomboNum$`0`)] <- 0

# Get proportions of correct / correct + incorrect
SUMmodecomboNum$prop <- (SUMmodecomboNum$`1` / (SUMmodecomboNum$`1`+SUMmodecomboNum$`0`))
SUMmodecomboNum

# Heatmap - modality + numerosity
ggplot(data =  SUMmodecomboNum, mapping = aes(x = g_size1, y = g_size2, fill = prop)) +
  geom_tile(aes(fill = prop), colour = "white") +
  scale_fill_gradient(high = "red", low = "yellow") +
  geom_raster(aes(fill = prop)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("First group") +
  ylab("Second group") +
  labs(title = "Forced-choice accuracy across modalities and sizes", fill = "Accuracy")

##################################
### FC Signal Detection Theory ###
##################################
### Indicate hits and false alarms
for (i in 1:length(fc_trials[,1])) {
  if (fc_trials$bigger[i] == "laugh" & fc_trials$correct[i] == "true") {
    fc_trials$sdt[i] <- "hit"
  } else if (fc_trials$smaller[i] == "laugh" & fc_trials$correct[i] == "false") {
    fc_trials$sdt[i] <- "fa"
  } else if (fc_trials$bigger[i] == "sentence" & fc_trials$correct[i] == "true") {
    fc_trials$sdt[i] <- "hit"
  } else if (fc_trials$smaller[i] == "sentence" & fc_trials$correct[i] == "false") {
    fc_trials$sdt[i] <- "fa"
  } else if (fc_trials$bigger[i] == "talking" & fc_trials$correct[i] == "true") {
    fc_trials$sdt[i] <- "hit"
  } else if (fc_trials$smaller[i] == "talking" & fc_trials$correct[i] == "false") {
    fc_trials$sdt[i] <- "fa"
  } else if (fc_trials$bigger[i] == "yell" & fc_trials$correct[i] == "true") {
    fc_trials$sdt[i] <- "hit"
  } else if (fc_trials$smaller[i] == "yell" & fc_trials$correct[i] == "false") {
    fc_trials$sdt[i] <- "fa"
  } else {
    fc_trials$sdt[i] <- NA
  }
}

### Aggregate signal detection
laugh_SDT_df <- NULL
laugh_SDT_df$n_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true")
laugh_SDT_df$n_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false")
laugh_SDT_df$n_miss <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "false")
laugh_SDT_df$n_cr <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "true")
# Hits
laugh_SDT_df$six_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 6)
laugh_SDT_df$five_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 5)
laugh_SDT_df$four_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 4)
laugh_SDT_df$three_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 3)
laugh_SDT_df$two_hit <- sum(fc_trials$bigger == "laugh" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 2)
# FAs
laugh_SDT_df$five_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 5)
laugh_SDT_df$four_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 4)
laugh_SDT_df$three_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 3)
laugh_SDT_df$two_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 2)
laugh_SDT_df$one_fa <- sum(fc_trials$smaller == "laugh" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 1)

sentence_SDT_df <- NULL
sentence_SDT_df$n_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true")
sentence_SDT_df$n_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false")
sentence_SDT_df$n_miss <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "false")
sentence_SDT_df$n_cr <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "true")
# Hits
sentence_SDT_df$six_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 6)
sentence_SDT_df$five_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 5)
sentence_SDT_df$four_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 4)
sentence_SDT_df$three_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 3)
sentence_SDT_df$two_hit <- sum(fc_trials$bigger == "sentence" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 2)
# FAs
sentence_SDT_df$five_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 5)
sentence_SDT_df$four_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 4)
sentence_SDT_df$three_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 3)
sentence_SDT_df$two_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 2)
sentence_SDT_df$one_fa <- sum(fc_trials$smaller == "sentence" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 1)

talking_SDT_df <- NULL
talking_SDT_df$n_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true")
talking_SDT_df$n_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false")
talking_SDT_df$n_miss <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "false")
talking_SDT_df$n_cr <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "true")
# Hits
talking_SDT_df$six_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 6)
talking_SDT_df$five_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 5)
talking_SDT_df$four_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 4)
talking_SDT_df$three_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 3)
talking_SDT_df$two_hit <- sum(fc_trials$bigger == "talking" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 2)
# FAs
talking_SDT_df$five_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 5)
talking_SDT_df$four_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 4)
talking_SDT_df$three_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 3)
talking_SDT_df$two_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 2)
talking_SDT_df$one_fa <- sum(fc_trials$smaller == "talking" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 1)

yell_SDT_df <- NULL
yell_SDT_df$n_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true")
yell_SDT_df$n_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false")
yell_SDT_df$n_miss <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "false")
yell_SDT_df$n_cr <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "true")
# Hits
yell_SDT_df$six_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 6)
yell_SDT_df$five_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 5)
yell_SDT_df$four_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 4)
yell_SDT_df$three_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 3)
yell_SDT_df$two_hit <- sum(fc_trials$bigger == "yell" & fc_trials$correct == "true" & fc_trials$g_size_bigger == 2)
# FAs
yell_SDT_df$five_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 5)
yell_SDT_df$four_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 4)
yell_SDT_df$three_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 3)
yell_SDT_df$two_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 2)
yell_SDT_df$one_fa <- sum(fc_trials$smaller == "yell" & fc_trials$correct == "false" & fc_trials$g_size_smaller == 1)

laugh_SDT_df$modality <- "laugh"
sentence_SDT_df$modality <- "sentence"
talking_SDT_df$modality <- "talking"
yell_SDT_df$modality <- "yell"

SDT_df <- data.frame(modality = c(laugh_SDT_df$modality, sentence_SDT_df$modality, talking_SDT_df$modality, yell_SDT_df$modality),
                     n_hit = c(laugh_SDT_df$n_hit, sentence_SDT_df$n_hit, talking_SDT_df$n_hit, yell_SDT_df$n_hit),
                     n_fa = c(laugh_SDT_df$n_fa, sentence_SDT_df$n_fa, talking_SDT_df$n_fa, yell_SDT_df$n_fa),
                     n_miss = c(laugh_SDT_df$n_miss, sentence_SDT_df$n_miss, talking_SDT_df$n_miss, yell_SDT_df$n_miss),
                     n_cr = c(laugh_SDT_df$n_cr, sentence_SDT_df$n_cr, talking_SDT_df$n_cr, yell_SDT_df$n_cr))

# Get indices
FC_SDT_indices <- psycho::dprime(n_hit = SDT_df$n_hit, n_fa = SDT_df$n_fa, n_miss = SDT_df$n_miss, n_cr = SDT_df$n_cr)

### Table 1 data (minus AUC and plus additional stats)
full_FC_SDT <- cbind(SDT_df, FC_SDT_indices)
full_FC_SDT

# Calculate area under curve for ROCs
AUC(full_FC_SDT$dprime[1]) # Laughing
AUC(full_FC_SDT$dprime[2]) # Sentence
AUC(full_FC_SDT$dprime[3]) # Talking
AUC(full_FC_SDT$dprime[4]) # Yelling

### Exploratory - Look at ROC curves using sensR package
# Laughing
sensR::ROC(full_FC_SDT$dprime[1], col = "red", lwd = 2)
# Sentence
sensR::ROC(full_FC_SDT$dprime[2], col = "green", lwd = 2)
# Talking
sensR::ROC(full_FC_SDT$dprime[3], col = "steelblue2", lwd = 2)
# Yelling
sensR::ROC(full_FC_SDT$dprime[4], col = "purple", lty=2, lwd = 2)
title("ROC curves for detecting bigger group")
legend(x = "topright", legend = c("Laughing", "Sentence", "Talking", "Yelling"), col = c("red", "green", "steelblue2", "purple"), lwd = 2)
