#title: "Latent Factors of Language Disturbance and Relationships to Quantitative Speech Features: Correlations"
#author: "Tang, Hänsel, Cong, Nikzad, Mehta, Cho, Berretta, Behbehani, Pradhan, John, & Liberman"
#date: "8/6/2022"

### Related Scripts ###
# 1) Script Compile Factors
# 2) RMD Factors Expanded
# 3) Script Compile Networks
# 4) RMD Network
# 5) Script Correlations ** This one

### Load Packages
library(effsize) #cohen d, etc.
library(lm.beta) #standardized beta
library(psych) #psych stats package, inlcudes factor analysis
library(gridExtra) #arranges ggplot grobs
library(Hmisc) #works with ggplot and others
library(ggpubr) #pairwise and statistics
library(bestNormalize) #yeojohnson and other transformations
library(mice) #multiple imputation
library(plyr) #data manipulations
library(MatchIt) #simpler matching
library(naniar) #na functions
library(arsenal) #tableby
library(hexbin) #Binning and plotting functions for hexagonal bins - network
library(networkD3) #Creates D3 JavaScript network, tree, dendrogram graphs from R.
library(tidytext) #Text mining
library(plotly) #Interactive graphs
library(repr)#String and binary representations - Network
library(RColorBrewer)
library(ggthemes) #Preset themes for ggplot
library(corrplot) #correlation plots
library(tidyverse) #dplyr, ggplot, stringr, et al

set.seed(123)

### Load Clinical Data
load("crossdx_all.R")
load("crossdx_clin.R")
load("vars_network.R")

### Data Dictionary
#crossdx_all - merged clinical (subset of clin vars) and NLP features
#dat.network.* - datasets for correlations including all network features + subsetting for different segments of the data
#cor.network.* - correlation matrices...


################ Demographics ################
tableby(group ~ assessment_type + demo_age + demo_gender + demo_race + demo_ethnicity + tlc_global,
        data=crossdx.all, numeric.stats=c("mean", "range", "sd")) %>% summary(digits = 2, digits.p=3)


### Different from TLC sample? Yes p=0.04
table(crossdx.clin$group)
table(crossdx.all$group)

mat.groups <- matrix(c(76, 130, 47, 90, 37, 87, 14, 64), nrow = 2, ncol = 4, byrow = T)
rownames(mat.groups) <- c('TLC', 'NLP')
colnames(mat.groups) <- c('HV', 'OD', 'PSY', 'SSD')
mat.groups

chisq.test(mat.groups)

################ Setting up function ################
fun.cor.nlp <- function(dat, grps, astyp, extension) {
      # Correlation matrix with all items
      dat.network <- filter(dat, group %in% grps & assessment_type %in% astyp) %>% select(all_of(vars_network))
      cor.network <- cor(dat.network, method = "spearman", use = "complete.obs")
      
      # Subset unique comparisons
      n <- dim(cor.network)[2]
      cor.network <- cor.network[4:n, 1:3]
      
      # Names
      colnames(cor.network) <- c("Inefficient Sp.", "Incoherent Sp.", "Imp. Expressivity")
      cor.network %>% rownames()
      rownames(cor.network) <- c("QUAN:TypeTokenRat",
                                         "QUAN:UtterLength",
                                         "GPH:SeqClustCoeff",
                                         "GPH:SeqLargestClique",
                                         "GPH:ActionPredClustCoeff",
                                         "POS:Adjectives",
                                         "POS:Adpositions",
                                         "POS:CoorConjunc",
                                         "POS:SubConjunc",
                                         "POS:Determiners",
                                         "POS:Particles",
                                         "ERR:PartialWords",
                                         "ERR:RepeatedWords",
                                         "ERR:FilledPauses",
                                         "COS:LSAcosineDist",
                                         "LEX:AgeofAquisition",
                                         "SENT:PositiveValence",
                                         "SENT:NegativeValence",
                                         "SENT:Arousal",
                                         "VQ:MeanPitch",
                                         "VQ:PitchVari",
                                         "VQ:ShimmerVari",
                                         "TEM:MeanPauseLength",
                                         "TEM:PauseLengthVari",
                                         "TEM:MeanSpeakRate",
                                         "TEM:MinSpeakRate",
                                         "TEM:MeanTurnLatency")
      
      #Correlation plot
      png(height=1200, width=400, file=paste("fig.cor.", extension, ".png", sep=""), type = "cairo")
      corrplot(cor.network, method = 'shade', addCoef.col = 'black', tl.col = 'black', cl.pos = 'n')
      dev.off()
      
}



################ Whole Sample ################
## Cross-dx
dat=crossdx.all
grps=c("HC", "HSC", "OPD", "SSD") 
astyp=c("in-person", "virtual")
extension="crossdx"

fun.cor.nlp(dat, grps, astyp, extension)

## HC
grps=c("HC") 
extension="hv"

fun.cor.nlp(dat, grps, astyp, extension)

## HSC
grps=c("HSC") 
extension="od"

fun.cor.nlp(dat, grps, astyp, extension)

## OPD
grps=c("OPD") 
extension="psy"

fun.cor.nlp(dat, grps, astyp, extension)

## SSD
grps=c("SSD") 
extension="ssd"

fun.cor.nlp(dat, grps, astyp, extension)



################ Assessment type-specific Sample ################
## Assessment location x group
table(crossdx.all$group, crossdx.all$assessment_type)

## HC
dat=crossdx.all
grps=c("HC") 
astyp=c("virtual")
extension="hv-vir"

fun.cor.nlp(dat, grps, astyp, extension)

## HSC
grps=c("HSC") 
astyp=c("in-person")
extension="od-inp"

fun.cor.nlp(dat, grps, astyp, extension)

## OPD
grps=c("OPD") 
astyp=c("in-person")
extension="psy-inp"

fun.cor.nlp(dat, grps, astyp, extension)

## SSD - in-person
grps=c("SSD") 
astyp=c("in-person")
extension="ssd-inp"

fun.cor.nlp(dat, grps, astyp, extension)

## SSD - virtual
grps=c("SSD") 
astyp=c("virtual")
extension="ssd-vir"

fun.cor.nlp(dat, grps, astyp, extension)




################ Demo-Matched Sample: Compile ################
# Matching gender, race, age
# Index group is HV. Not including PSY due to two few people

### setting up
# White/Black Race 
crossdx.all$demo_white <- rep(0, dim(crossdx.all)[1])
crossdx.all$demo_black <- rep(0, dim(crossdx.all)[1])

crossdx.all$demo_white[which(crossdx.all$demo_race == "White")] <- 1
crossdx.all$demo_black[which(crossdx.all$demo_race == "Black")] <- 1

table(crossdx.all$demo_white)
table(crossdx.all$demo_black)

### Matching

## HC + HSC
temp <- filter(crossdx.all, group == "HSC" | group == "HC")
chisq.test(temp$group, temp$demo_white)
chisq.test(temp$group, temp$demo_black)
t.test(temp$demo_age ~ temp$group)
chisq.test(temp$group, temp$demo_gender)

# Checking std. Mean differences prior to Match
matchit(group ~ demo_gender + demo_age + demo_white + demo_black, data = temp, method = NULL, distance = "glm") %>% summary()

# Match: 
mch.age_hsc <- matchit(group ~ demo_age + demo_white + demo_black, data = temp, method = "optimal", distance = "glm", link="probit")
summary(mch.age_hsc)
temp1 <- match.data(mch.age_hsc)
chisq.test(temp1$group, temp1$demo_race)
chisq.test(temp1$group, temp1$demo_gender)
t.test(temp1$demo_age ~ temp1$group)

## HC + SSD 
temp <- filter(crossdx.all, group == "SSD" | group == "HC")
chisq.test(temp$group, temp$demo_white)
chisq.test(temp$group, temp$demo_black)
chisq.test(temp$group, temp$demo_gender)
t.test(temp$demo_age ~ temp$group)

# Checking std. Mean differences prior to Match
matchit(group ~ demo_age + demo_white + demo_black, data = temp, method = NULL, distance = "glm") %>% summary()

# Match: 
mch.age_ssd <- matchit(group ~ demo_age + demo_white + demo_black, data = temp, method = "optimal", distance = "glm", link="probit")
summary(mch.age_ssd)
temp2 <- match.data(mch.age_ssd)
chisq.test(temp2$group, temp2$demo_race)
chisq.test(temp2$group, temp2$demo_gender)
t.test(temp2$demo_age ~ temp2$group)


### Compile
mch.nlp <- temp1
temp2 <- filter(temp2, group == "SSD")
mch.nlp <- rbind(mch.nlp, temp2)

### Demographics
tableby(group ~ assessment_type + demo_age + demo_gender + demo_race + demo_ethnicity + tlc_global,
        data=mch.nlp, numeric.stats=c("mean", "range", "sd")) %>% summary(digits = 2, digits.p=3)

################ Demo-matched sample: correlations ################
## Cross-dx
dat=mch.nlp
grps=c("HC", "HSC", "SSD") 
astyp=c("in-person", "virtual")
extension="mch-crossdx"

fun.cor.nlp(dat, grps, astyp, extension)

## HC
grps=c("HC") 
extension="mch-hv"

fun.cor.nlp(dat, grps, astyp, extension)

## HSC
grps=c("HSC") 
extension="mch-od"

fun.cor.nlp(dat, grps, astyp, extension)

## SSD
grps=c("SSD") 
extension="mch-ssd"

fun.cor.nlp(dat, grps, astyp, extension)


################ Word count ################
# Overall
tableby(group ~ word_count,
        data=crossdx.all, numeric.stats=c("mean", "range", "sd")) %>% summary(digits = 2, digits.p=3)

# Just in-person
temp <- filter(crossdx.all, assessment_type == "in-person")
tableby(group ~ word_count,
        data=temp, numeric.stats=c("mean", "range", "sd")) %>% summary(digits = 2, digits.p=3)

# Just Virtual
temp <- filter(crossdx.all, assessment_type == "virtual")
tableby(group ~ word_count,
        data=temp, numeric.stats=c("mean", "range", "sd")) %>% summary(digits = 2, digits.p=3)
