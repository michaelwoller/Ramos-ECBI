library(ggplot2)
library(lme4)
library(lmerTest)
library(haven)



# Loading in data
dat <- read.csv("PCIT_Long.csv") # long data

dat$RLP <- dat$Race_Language_Predictor 


# Creating variables
dat$Ind <- rep(0, nrow(dat))
dat$Ind[is.na(dat$PDI_ID) == F] <- 1

# Centering time

dat$SID_Cent <- rep(0, length(dat$SessionID))

parts <- unique(dat$ParticipantID)
lp <- length(parts)

for (i in 1:lp) {
  dt <- dat[dat$ParticipantID == parts[i],]
  cdi <- dt$CDI_ID[is.na(dt$CDI_ID) == F]
  
  # Cetnering CDI half
  SID_CDI <- dt$SessionID[dt$Ind == 0]
  
  lcdi <- length( SID_CDI )
  
  CDIStart <- min(SID_CDI)
  CDIEnd <- max(SID_CDI)
  Cmiss <- setdiff(CDIStart:CDIEnd, SID_CDI)
  lCmiss <- length(Cmiss)
  
  
  
  CDICent <- SID_CDI -  CDIStart - lcdi - lCmiss
  
  # Centering PDI half
  SID_PDI <- dt$SessionID[dt$Ind == 1]
  lpdi <- length( SID_PDI )
  
  if (lpdi != 0) {
    PDIStart <- min(SID_PDI)
    PDIEnd <- max(SID_PDI)
    Pmiss <- setdiff(PDIStart:PDIEnd, dt$SessionID)
    lPmiss <- length(Pmiss)
    
    PDICent <- SID_PDI -  PDIStart + 1
    
    # Combining and centering full variable
    SID_Cent <- c(CDICent, PDICent)
    dt$SID_Cent <- SID_Cent
    
  } else {
    # Combining and centering full variable
    SID_Cent <- CDICent
    dt$SID_Cent <- SID_Cent
  }
  dat[dat$ParticipantID == parts[i],] <-  dt
}

# Factoring ParticipantID
dat$ParticipantID <- factor(dat$ParticipantID)

# Factoring Indicator
dat$Ind <- factor(dat$Ind)
# Redoing baseline with reference as PDI slope
dat$Ind2 <- relevel(dat$Ind, ref = "1")


# Factoring RLP with white/eng as comparison
dat$RLP <- factor(dat$RLP, levels = c(0,1,2), labels = c("White/Eng", "Latin/Eng", "Latin/Span"))
# Redoing baseline with reference being Latin/Eng
dat$RLP2 <- relevel(dat$RLP, ref = "Latin/Eng")
# Redoing baseline with reference being Latin/Span
dat$RLP3 <- relevel(dat$RLP, ref = "Latin/Span")

###################################
#### Initial DPICS_P Model #######
###################################

# Setting up the BOBYQA optimzer
ctrl <- lme4::lmerControl(optimizer="bobyqa",
                          optCtrl=list(maxfun=2e5))

# Setting up the model
# Note: use lmer() function from "lmerTest" package to get p-values DONT use lme() from "lme4"

# CDI as reference (Indicator = 0 as baseline)
modelP <- lmerTest::lmer(DPICS_P ~ 
                           # The CDI and PDI slopes
                           SID_Cent + SID_Cent:Ind + 
                           # Race-Language Covariate and Interactions
                           RLP + SID_Cent*RLP + (SID_Cent:Ind)*RLP + 
                           # The random slopes for CDI and PDI for each individual
                           (1 + SID_Cent | ParticipantID),
                         data = dat, REML = F, control = ctrl)

# Note, 'lmerTest' package needs to be loaded to use the ddf argument
# S stands for Satterthwaite df correction
summary(modelP, ddf = "S")

# PDI as reference (Indicator = 1 as baseline) (have to reverse the random effect too)
modelPp <- lmerTest::lmer(DPICS_P ~  SID_Cent + SID_Cent:Ind2 + 
                            RLP + SID_Cent*RLP + (SID_Cent:Ind2)*RLP + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelPp, ddf = "S")



#### Redoing Baseline RLP as LatinX/Eng #######

# CDI slopes
modelP2 <- lmerTest::lmer(DPICS_P ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind)*RLP2 + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelP2, ddf = "S")

# PDI slopes
modelPp2 <- lmerTest::lmer(DPICS_P ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind2)*RLP2 + 
                             (1 + SID_Cent | ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelPp2, ddf = "S")

#### Redoing Baseline RLP as LatinX/Span #######

# CDI slopes
modelP3 <- lmerTest::lmer(DPICS_P ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind)*RLP3 + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelP3, ddf = "S")
# PDI slopes
modelPp3 <- lmerTest::lmer(DPICS_P ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind2)*RLP3 + 
                             (1 + SID_Cent | ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelPp3, ddf = "S")


# White/Eng
dd <- summary(modelP, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[1]

Plines <- data.frame(xcoor,ycoor, group)

# Latin/Eng
dd <- summary(modelP2, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

group <- levels(dat$RLP)[2]

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

lines2 <- data.frame(xcoor,ycoor, group)

Plines <- rbind(Plines, lines2)


# Latin/Span
dd <- summary(modelP3, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[3]
lines3 <- data.frame(xcoor,ycoor, group)

Plines <- rbind(Plines, lines3)

plot <- ggplot(dat, aes(x = SID_Cent, y = DPICS_P, group = ParticipantID)) +
  geom_point() +
  
  # Line segments per RLP group
  geom_segment(data = subset(dat, RLP == "White/Eng"), 
               aes(x = Plines[1,1], y = Plines[1,2], xend = Plines[2,1], yend = Plines[2,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "White/Eng"),
               aes(x = Plines[3,1], y = Plines[3,2], xend = Plines[4,1], yend = Plines[4,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"), 
               aes(x = Plines[5,1], y = Plines[5,2], xend = Plines[6,1], yend = Plines[6,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"),
               aes(x = Plines[7,1], y = Plines[7,2], xend = Plines[8,1], yend = Plines[8,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"), 
               aes(x = Plines[9,1], y = Plines[9,2], xend = Plines[10,1], yend = Plines[10,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"),
               aes(x = Plines[11,1], y = Plines[11,2], xend = Plines[12,1], yend = Plines[12,2], color = "black")) +
  
  labs(x = "Session", y = "DPICS Positive", title = "Figure 7. Average Piecewise for DPICS Positive per Race Language Group") +
  guides(col = FALSE) +
  facet_grid(. ~ RLP) +
  theme_minimal()
plot


###################################
#### Initial DPICS_N Model #######
###################################

# Setting up the BOBYQA optimzer
ctrl <- lme4::lmerControl(optimizer="bobyqa",
                          optCtrl=list(maxfun=2e5))

# Setting up the model
# Note: use lmer() function from "lmerTest" package to get p-values

# CDI slopes
modelN <- lmerTest::lmer(DPICS_N ~ 
                           # The CDI and PDI slopes
                           SID_Cent + SID_Cent:Ind + 
                           # Race-Language Covariate and Interactions
                           RLP + SID_Cent*RLP + (SID_Cent:Ind)*RLP + 
                           # The random slopes for CDI and PDI for each individual
                           (1 + SID_Cent | ParticipantID),
                         data = dat, REML = F, control = ctrl)

# Note, 'lmerTest' package needs to be loaded to use the ddf argument
# S stands for Satterthwaite df correction
summary(modelN, ddf = "S")

# PDI as reference (Indicator = 1 as baseline)
modelNp <- lmerTest::lmer(DPICS_N ~  SID_Cent + SID_Cent:Ind2 + 
                            RLP + SID_Cent*RLP + (SID_Cent:Ind2)*RLP + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelNp, ddf = "S")


#### Redoing Baseline RLP as LatinX/Eng #######

# CDI slopes
modelN2 <- lmerTest::lmer(DPICS_N ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind)*RLP2 + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelN2, ddf = "S")
# PDI as baseline
modelNp2 <- lmerTest::lmer(DPICS_N ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind2)*RLP2 + 
                             (1 + SID_Cent | ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelNp2, ddf = "S")

#### Redoing Baseline RLP as LatinX/Span #######
# CDI slopes
modelN3 <- lmerTest::lmer(DPICS_N ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind)*RLP3 + 
                            (1 + SID_Cent | ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelN3, ddf = "S")
# PDI slopes
modelNp3 <- lmerTest::lmer(DPICS_N ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind2)*RLP3 + 
                             (1 + SID_Cent | ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelNp3, ddf = "S")


# White/Eng
dd <- summary(modelN, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[1]

Nlines <- data.frame(xcoor,ycoor, group)

# Latin/Eng
dd <- summary(modelN2, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[2]
lines2 <- data.frame(xcoor,ycoor, group)

Nlines <- rbind(Nlines, lines2)


# Latin/Span
dd <- summary(modelN3, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[3]
lines3 <- data.frame(xcoor,ycoor, group)

Nlines <- rbind(Nlines, lines3)

####### Facet Plot ######

plot <- ggplot(dat, aes(x = SID_Cent, y = DPICS_N, group = ParticipantID)) +
  geom_point() +
  
  # Line segments per RLP group
  geom_segment(data = subset(dat, RLP == "White/Eng"), 
               aes(x = Nlines[1,1], y = Nlines[1,2], xend = Nlines[2,1], yend = Nlines[2,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "White/Eng"),
               aes(x = Nlines[3,1], y = Nlines[3,2], xend = Nlines[4,1], yend = Nlines[4,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"), 
               aes(x = Nlines[5,1], y = Nlines[5,2], xend = Nlines[6,1], yend = Nlines[6,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"),
               aes(x = Nlines[7,1], y = Nlines[7,2], xend = Nlines[8,1], yend = Nlines[8,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"), 
               aes(x = Nlines[9,1], y = Nlines[9,2], xend = Nlines[10,1], yend = Nlines[10,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"),
               aes(x = Nlines[11,1], y = Nlines[11,2], xend = Nlines[12,1], yend = Nlines[12,2], color = "black")) +
  
  labs(x = "Session", y = "DPICS Negative", title = "Figure 8. Average Piecewise for DPICS Negative per Race Language Group") +
  guides(col = FALSE) +
  facet_grid(. ~ RLP) +
  theme_minimal()
plot


###################################
#### Initial ECBI Model #######
###################################

# Setting up the BOBYQA optimzer
ctrl <- lme4::lmerControl(optimizer="bobyqa",
                          optCtrl=list(maxfun=2e5))

# Setting up the model
# Note: use lmer() function from "lmerTest" package to get p-values
# CDI slopes
modelE <- lmerTest::lmer(ECBI ~ 
                           # The CDI and PDI slopes
                           SID_Cent + SID_Cent:Ind + 
                           # Race-Language Covariate and Interactions
                           RLP + SID_Cent*RLP + (SID_Cent:Ind)*RLP + 
                           # The random slopes for CDI and PDI for each individual
                           (1 + SID_Cent + SID_Cent:Ind| ParticipantID),
                         data = dat, REML = F, control = ctrl)

# Note, 'lmerTest' package needs to be loaded to use the ddf argument
# S stands for Satterthwaite df correction
summary(modelE, ddf = "S")

# PDI as reference
modelEp <- lmerTest::lmer(ECBI ~ 
                            SID_Cent + SID_Cent:Ind2 + 
                            RLP + SID_Cent*RLP + (SID_Cent:Ind2)*RLP + 
                            (1 + SID_Cent + SID_Cent:Ind2| ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelEp, ddf = "S")

#### Redoing Baseline RLP as LatinX/Eng #######
# CDI slopes
modelE2 <- lmerTest::lmer(ECBI ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind)*RLP2 + 
                            (1 + SID_Cent + SID_Cent:Ind| ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelE2, ddf = "S")
# PDI as reference
modelEp2 <- lmerTest::lmer(ECBI ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP2 + SID_Cent*RLP2 + (SID_Cent:Ind2)*RLP2 + 
                             (1 + SID_Cent + SID_Cent:Ind2| ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelEp2, ddf = "S")

#### Redoing Baseline RLP as LatinX/Span #######
# CDI slopes
modelE3 <- lmerTest::lmer(ECBI ~ 
                            SID_Cent + SID_Cent:Ind + 
                            RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind)*RLP3 + 
                            (1 + SID_Cent + SID_Cent:Ind| ParticipantID),
                          data = dat, REML = F, control = ctrl)
summary(modelE3, ddf = "S")
# PDI as reference
modelEp3 <- lmerTest::lmer(ECBI ~ 
                             SID_Cent + SID_Cent:Ind2 + 
                             RLP3 + SID_Cent*RLP3 + (SID_Cent:Ind2)*RLP3 + 
                             (1 + SID_Cent + SID_Cent:Ind2| ParticipantID),
                           data = dat, REML = F, control = ctrl)
summary(modelEp3, ddf = "S")


# White/Eng
dd <- summary(modelN, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[1]

Nlines <- data.frame(xcoor,ycoor, group)

# Latin/Eng
dd <- summary(modelN2, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[2]
lines2 <- data.frame(xcoor,ycoor, group)

Nlines <- rbind(Nlines, lines2)


# Latin/Span
dd <- summary(modelN3, ddf = "S")
coefsMean <- coef(dd)[,1]
int <- coefsMean[1]
slope1 <- coefsMean[2]
slope2 <- coefsMean[2] + coefsMean[5]
slopes <- c(slope1, slope2)

xstart <- c(min(dat$SID_Cent), 0)
xend <- c(0, max(dat$SID_Cent))
xcoor <- c(xstart, xend)

ystart <- slopes*xstart + int
yend <- slopes*xend + int
ycoor <- c(ystart,yend)

group <- levels(dat$RLP)[3]
lines3 <- data.frame(xcoor,ycoor, group)

Nlines <- rbind(Nlines, lines3)

####### Facet Plot ######

plot <- ggplot(dat, aes(x = SID_Cent, y = DPICS_N, group = ParticipantID)) +
  geom_point() +
  
  # Line segments per RLP group
  geom_segment(data = subset(dat, RLP == "White/Eng"), 
               aes(x = Nlines[1,1], y = Nlines[1,2], xend = Nlines[2,1], yend = Nlines[2,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "White/Eng"),
               aes(x = Nlines[3,1], y = Nlines[3,2], xend = Nlines[4,1], yend = Nlines[4,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"), 
               aes(x = Nlines[5,1], y = Nlines[5,2], xend = Nlines[6,1], yend = Nlines[6,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Eng"),
               aes(x = Nlines[7,1], y = Nlines[7,2], xend = Nlines[8,1], yend = Nlines[8,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"), 
               aes(x = Nlines[9,1], y = Nlines[9,2], xend = Nlines[10,1], yend = Nlines[10,2], color = "black")) +
  geom_segment(data = subset(dat, RLP == "Latin/Span"),
               aes(x = Nlines[11,1], y = Nlines[11,2], xend = Nlines[12,1], yend = Nlines[12,2], color = "black")) +
  
  labs(x = "Session", y = "DPICS Negative", title = "Figure 8. Average Piecewise for DPICS Negative per Race Language Group") +
  guides(col = FALSE) +
  facet_grid(. ~ RLP) +
  theme_minimal()
plot


###################################### 
##### Followup Mean Comparisons ##### 
###################################### 


df <- read_sav("PCIT_Comparisons.sav")

df <- as.data.frame(df)
df$RLP <- factor(df$Race_Language_Predictor, levels = c(0,1,2), labels = c("White/Eng", "Latin/Eng", "Latin/Span"))
df2 <- subset(df, !ParticipantID %in% c(1, 8,12, 14, 15, 17, 18, 19, 21, 22, 31, 32, 33, 35, 36, 43, 45, 46, 48, 49, 50, 51, 53, 54, 60, 65, 67, 68, 72, 74, 75, 77, 80, 82, 83, 85, 88, 92, 93, 102, 103))


########## number of weeks in the Clinic (Total_Weeks_Clinic)
##### all participants
## Kruskal-Wallace test 

kruskal.test(Total_Weeks_Clinic~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(Total_Weeks_Clinic~ RLP, data = df2)

########## of no-shows or cancellations (Total_NoShows_Cancellations)

##### all participants
## Kruskal-Wallace test 

kruskal.test(Total_NoShows_Cancellations~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(Total_NoShows_Cancellations~ RLP, data = df2)

##### complete participants
## Kruskal-Wallace test 

########## of weeks in CDI (CDI_Weeks)

##### all participants
## Kruskal-Wallace test 

kruskal.test(CDI_Weeks~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(CDI_Weeks~ RLP, data = df2)


########## of actual CDI sessions attended (CDI_Sessions_Attended)

##### all participants
## Kruskal-Wallace test 

kruskal.test(CDI_Sessions_Attended~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(CDI_Sessions_Attended~ RLP, data = df2)

## Further Man-Whitney U test for pairwise comparisons
# adding small noise to prevent exact matches
mean(df2$CDI_Sessions_Attended[df2$RLP == "White/Eng"])
mean(df2$CDI_Sessions_Attended[df2$RLP == "Latin/Eng"])
mean(df2$CDI_Sessions_Attended[df2$RLP == "Latin/Span"])


df2$CDI_Sessions_Attended[1] <- df2$CDI_Sessions_Attended[1] + .001
df2$CDI_Sessions_Attended[2] <- df2$CDI_Sessions_Attended[2] + .001
df2$CDI_Sessions_Attended[3] <- df2$CDI_Sessions_Attended[3] + .001
df2$CDI_Sessions_Attended[4] <- df2$CDI_Sessions_Attended[4] + .001
df2$CDI_Sessions_Attended[5] <- df2$CDI_Sessions_Attended[5] + .001
df2$CDI_Sessions_Attended[6] <- df2$CDI_Sessions_Attended[6] + .001
df2$CDI_Sessions_Attended[7] <- df2$CDI_Sessions_Attended[7] + .001
df2$CDI_Sessions_Attended[8] <- df2$CDI_Sessions_Attended[8] - .001
df2$CDI_Sessions_Attended[9] <- df2$CDI_Sessions_Attended[9] - .002
df2$CDI_Sessions_Attended[10] <- df2$CDI_Sessions_Attended[10] + .002
df2$CDI_Sessions_Attended[12] <- df2$CDI_Sessions_Attended[12] - .002
df2$CDI_Sessions_Attended[13] <- df2$CDI_Sessions_Attended[13] + .002
df2$CDI_Sessions_Attended[18] <- df2$CDI_Sessions_Attended[18] - .001



pairwise_result <- pairwise.wilcox.test(df2$CDI_Sessions_Attended, df2$RLP, 
                                        alternative = "two.sided",  
                                        p.adjust.method = "holm", 
                                        paired = FALSE)
WE <- df2[df2$RLP == "White/Eng",]
LE <- df2[df2$RLP == "Latin/Eng",]
LS <- df2[df2$RLP == "Latin/Span",]

w <- wilcox.test(WE$CDI_Sessions_Attended, LE$CDI_Sessions_Attended,  
                 alternative = "two.sided",
                 p.adjust.method = "holm",
                 paired = FALSE)

w2 <- wilcox.test(WE$CDI_Sessions_Attended, LS$CDI_Sessions_Attended,  
                  alternative = "two.sided",
                  p.adjust.method = "holm",
                  paired = FALSE)
w3 <- wilcox.test(LE$CDI_Sessions_Attended, LS$CDI_Sessions_Attended,  
                  alternative = "two.sided",
                  p.adjust.method = "holm",
                  paired = FALSE)

w$statistic
w2$statistic
w3$statistic

# So there is a significant difference between LatinX/Spanish and BOTH White/English and LatinX/English

########## weeks in PDI (PDI_Weeks)

##### all participants
## Kruskal-Wallace test 

kruskal.test(PDI_Weeks~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(PDI_Weeks~ RLP, data = df2)

########### actual PDI sessions attended (PDI_Sessions_Attended)
##### all participants
## Kruskal-Wallace test 

kruskal.test(PDI_Sessions_Attended~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(PDI_Sessions_Attended~ RLP, data = df2)

########## ECBI score by the end of CDI (ECBI_CDI)
##### all participants
## Kruskal-Wallace test 

kruskal.test(ECBI_CDI~ RLP, data = df)

##### complete participants
## Kruskal-Wallace test 

kruskal.test(ECBI_CDI~ RLP, data = df2)

################################################################################
# for the dichotomous variables, I ended up doing both a fisher exact test and a chi square
# though I think chi-square would be fine since the groups have at LEAST N = 5 
# which is the blanket recommended

########## % of participants who completed CDI (CDI_Mastery) 

##### all participants

contingency_table <- table(df$CDI_Mastery, df$RLP)

chisq.test(contingency_table, simulate.p.value = TRUE)
fisher.test(contingency_table)

##### complete participants
# cannot be done because ALL participants who went to completion reached mastery


########## % of participants who completed PDI (PDI_Mastery) 
##### all participants
contingency_table <- table(df$PDI_Mastery, df$RLP)

chisq.test(contingency_table)
fisher.test(contingency_table)

##### complete participants
# cannot be done because ALL participants who went to completion reached mastery


########## % of participants who dropped out (Treatment_DropOut) 
##### all participants
contingency_table <- table(df$Treatment_DropOut, df$RLP)

chisq.test(contingency_table)
fisher.test(contingency_table)

##### complete participants
# cannot be done because obviously this group did not drop out


# Creating variables
dat$Ind <- rep(0, nrow(dat))
dat$Ind[is.na(dat$PDI_ID) == F] <- 1

# Centering time

dat$SID_Cent <- rep(0, length(dat$SessionID))

dat <- dat[, c(1,2, 10, 3, 4, 5, 6, 7, 8, 9)]

parts <- unique(dat$ParticipantID)
lp <- length(parts)

for (i in 1:lp) {
  dt <- dat[dat$ParticipantID == parts[i],]
  cdi <- dt$CDI_ID[is.na(dt$CDI_ID) == F]
  
  # Cetnering CDI half
  SID_CDI <- dt$SessionID[dt$Ind == 0]
  
  lcdi <- length( SID_CDI )
  
  CDIStart <- min(SID_CDI)
  CDIEnd <- max(SID_CDI)
  Cmiss <- setdiff(CDIStart:CDIEnd, SID_CDI)
  lCmiss <- length(Cmiss)
  
  
  
  CDICent <- SID_CDI -  CDIStart - lcdi - lCmiss
  
  # Centering PDI half
  SID_PDI <- dt$SessionID[dt$Ind == 1]
  lpdi <- length( SID_PDI )
  
  if (lpdi != 0) {
    PDIStart <- min(SID_PDI)
    PDIEnd <- max(SID_PDI)
    Pmiss <- setdiff(PDIStart:PDIEnd, dt$SessionID)
    lPmiss <- length(Pmiss)
    
    PDICent <- SID_PDI -  PDIStart + 1
    
    # Combining and centering full variable
    SID_Cent <- c(CDICent, PDICent)
    dt$SID_Cent <- SID_Cent
    
  } else {
    # Combining and centering full variable
    SID_Cent <- CDICent
    dt$SID_Cent <- SID_Cent
  }
  dat[dat$ParticipantID == parts[i],] <-  dt
}
parts <- unique(dat$ParticipantID)

# beginnings
lp <- length(parts)

starts <- data.frame("ParticipantID" = 0, "SessionID" = 0,
                     "DPICS_P" = 0, "DPICS_N" = 0, "ECBI" = 0,
                     "RLP" = 0)
# Extracts relevent data groups
for (i in 1:lp) {
  dt <- dat[dat$ParticipantID == parts[i],]
  starts <- rbind(starts, c(dt$ParticipantID[1], dt$SessionID[1], dt$DPICS_P[1],
                            dt$DPICS_N[1], dt$ECBI[1], dt$RLP[1]))
}
starts <- starts[-1,]
starts$RLP <- factor(starts$RLP, levels = c(1,2,3), labels = c("White/Eng", "Latin/Eng", "Latin/Span"))



#####################################################################
# Means and SD's                
#####################################################################


################### Whole Sample
# ECBI
mean(dat$ECBI, na.rm = T)
sd(dat$ECBI, na.rm = T)
# 129.322 (40.27)

# DPICS_P
mean(dat$DPICS_P, na.rm = T)
sd(dat$DPICS_P, na.rm = T)
# 26.388 (12.02)

# DPICS_N
mean(dat$DPICS_N, na.rm = T)
sd(dat$DPICS_N, na.rm = T)
# 7.671 (9.126)



#################### White/English
WE <- dat[dat$RLP == 0,]

# ECBI
mean(WE$ECBI, na.rm = T)
sd(WE$ECBI, na.rm = T)
# 137.892 (38.413)

# DPICS_P
mean(WE$DPICS_P, na.rm = T)
sd(WE$DPICS_P, na.rm = T)
# 26.117 (12.007)

# DPICS_N
mean(WE$DPICS_N, na.rm = T)
sd(WE$DPICS_N, na.rm = T)
# 5.264 (5.587)



#################### Latin/English
LE <- dat[dat$RLP == 1,]

# ECBI
mean(LE$ECBI, na.rm = T)
sd(LE$ECBI, na.rm = T)
# 129.383 (41.521)

# DPICS_P
mean(LE$DPICS_P, na.rm = T)
sd(LE$DPICS_P, na.rm = T)
# 27.441 (12.159)

# DPICS_N
mean(LE$DPICS_N, na.rm = T)
sd(LE$DPICS_N, na.rm = T)
# 8.105 (9.966)



#################### Latin/Spanish
LS <- dat[dat$RLP == 2,]

# ECBI
mean(LS$ECBI, na.rm = T)
sd(LS$ECBI, na.rm = T)
# 117.815 (37.250)

# DPICS_P
mean(LS$DPICS_P, na.rm = T)
sd(LS$DPICS_P, na.rm = T)
# 24.524 (11.600)

# DPICS_N
mean(LS$DPICS_N, na.rm = T)
sd(LS$DPICS_N, na.rm = T)
# 9.790 (10.196)


############### Last sessions for those who complete CDI and PDI 

### CDI
CDImast <- df[df$CDI_Mastery == 1,]

datCDImast <- dat[dat$ParticipantID %in% unique(CDImast$ParticipantID) == T,]
CDImastCDIend <- datCDImast[datCDImast$SID_Cent == -1,]
CDImastCDIendWE <- CDImastCDIend[CDImastCDIend$RLP == 0,]
CDImastCDIendLE <- CDImastCDIend[CDImastCDIend$RLP == 1,]
CDImastCDIendLS <- CDImastCDIend[CDImastCDIend$RLP == 2,]

### PDI
PDImast <- df[df$PDI_Mastery == 1,]

datPDImast <- dat[dat$ParticipantID %in% unique(PDImast$ParticipantID) == T,]

parts <- unique(datPDImast$ParticipantID)
lp <- length(parts)

PDImastPDIend <- matrix(rep(0,10), ncol = 10)
PDImastPDIend <- as.data.frame(PDImastPDIend)
colnames(PDImastPDIend) <- colnames(datPDImast)

for (i in 1:lp) {
  dt <- datPDImast[datPDImast$ParticipantID == parts[i],]
  maxline <- dt[dt$PDI_ID == tail(dt$PDI_ID,1),]
  end <- tail(maxline, 1)
  PDImastPDIend <- rbind(PDImastPDIend, end)
}
PDImastPDIend <-  PDImastPDIend[-1,]
PDImastPDIendWE <- PDImastPDIend[PDImastPDIend$RLP == 0,]
PDImastPDIendLE <- PDImastPDIend[PDImastPDIend$RLP == 1,]
PDImastPDIendLS <- PDImastPDIend[PDImastPDIend$RLP == 2,]

########## Last session of CDI for those that 
########## achieved CDI criteria

### Whole Sample
# ECBI
mean(CDImastCDIend$ECBI, na.rm = T)
sd(CDImastCDIend$ECBI, na.rm = T)

# DPICS_P
mean(CDImastCDIend$DPICS_P, na.rm = T)
sd(CDImastCDIend$DPICS_P, na.rm = T)

# DPICS_N
mean(CDImastCDIend$DPICS_N, na.rm = T)
sd(CDImastCDIend$DPICS_N, na.rm = T)

### White English
# ECBI
mean(CDImastCDIendWE$ECBI, na.rm = T)
sd(CDImastCDIendWE$ECBI, na.rm = T)

# DPICS_P
mean(CDImastCDIendWE$DPICS_P, na.rm = T)
sd(CDImastCDIendWE$DPICS_P, na.rm = T)

# DPICS_N
mean(CDImastCDIendWE$DPICS_N, na.rm = T)
sd(CDImastCDIendWE$DPICS_N, na.rm = T)


### Latin English
# ECBI
mean(CDImastCDIendLE$ECBI, na.rm = T)
sd(CDImastCDIendLE$ECBI, na.rm = T)

# DPICS_P
mean(CDImastCDIendLE$DPICS_P, na.rm = T)
sd(CDImastCDIendLE$DPICS_P, na.rm = T)

# DPICS_N
mean(CDImastCDIendLE$DPICS_N, na.rm = T)
sd(CDImastCDIendLE$DPICS_N, na.rm = T)

### Latin Spanish
# ECBI
mean(CDImastCDIendLS$ECBI, na.rm = T)
sd(CDImastCDIendLS$ECBI, na.rm = T)

# DPICS_P
mean(CDImastCDIendLS$DPICS_P, na.rm = T)
sd(CDImastCDIendLS$DPICS_P, na.rm = T)

# DPICS_N
mean(CDImastCDIendLS$DPICS_N, na.rm = T)
sd(CDImastCDIendLS$DPICS_N, na.rm = T)


##### Kruskal Wallace Tests

kruskal.test(ECBI ~ RLP, data = CDImastCDIend)
kruskal.test(DPICS_P ~ RLP, data = CDImastCDIend)
kruskal.test(DPICS_N ~ RLP, data = CDImastCDIend)


########## Last session of PDI for those that 
########## achieved PDI criteria


### Whole Sample
# ECBI
mean(PDImastPDIend$ECBI, na.rm = T)
sd(PDImastPDIend$ECBI, na.rm = T)

# DPICS_P
mean(PDImastPDIend$DPICS_P, na.rm = T)
sd(PDImastPDIend$DPICS_P, na.rm = T)

# DPICS_N
mean(PDImastPDIend$DPICS_N, na.rm = T)
sd(PDImastPDIend$DPICS_N, na.rm = T)

### White English
# ECBI
mean(PDImastPDIendWE$ECBI, na.rm = T)
sd(PDImastPDIendWE$ECBI, na.rm = T)

# DPICS_P
mean(PDImastPDIendWE$DPICS_P, na.rm = T)
sd(PDImastPDIendWE$DPICS_P, na.rm = T)

# DPICS_N
mean(PDImastPDIendWE$DPICS_N, na.rm = T)
sd(PDImastPDIendWE$DPICS_N, na.rm = T)


### Latin English
# ECBI
mean(PDImastPDIendLE$ECBI, na.rm = T)
sd(PDImastPDIendLE$ECBI, na.rm = T)

# DPICS_P
mean(PDImastPDIendLE$DPICS_P, na.rm = T)
sd(PDImastPDIendLE$DPICS_P, na.rm = T)

# DPICS_N
mean(PDImastPDIendLE$DPICS_N, na.rm = T)
sd(PDImastPDIendLE$DPICS_N, na.rm = T)

### Latin Spanish
# ECBI
mean(PDImastPDIendLS$ECBI, na.rm = T)
sd(PDImastPDIendLS$ECBI, na.rm = T)

# DPICS_P
mean(PDImastPDIendLS$DPICS_P, na.rm = T)
sd(PDImastPDIendLS$DPICS_P, na.rm = T)

# DPICS_N
mean(PDImastPDIendLS$DPICS_N, na.rm = T)
sd(PDImastPDIendLS$DPICS_N, na.rm = T)


##### Kruskal Wallace Tests

kruskal.test(ECBI ~ RLP, data = PDImastPDIend)
kruskal.test(DPICS_P ~ RLP, data = PDImastPDIend)
kruskal.test(DPICS_N ~ RLP, data = PDImastPDIend)


########## Number of sessions attended 

dd <- df

dd$SA_whole <- df$CDI_Sessions_Attended + df$PDI_Sessions_Attended


### whole sample
mean(SA_whole)
sd(SA_whole)


### Race Language Groups
WE <- dd[dd$Race_Language_Predictor == 0,]
mean(WE$SA_whole)
sd(WE$SA_whole)

LE <- dd[dd$Race_Language_Predictor == 1,]
mean(LE$SA_whole)
sd(LE$SA_whole)


LS <- dd[dd$Race_Language_Predictor == 2,]
mean(LS$SA_whole)
sd(LS$SA_whole)

# Kruskal Wallace Test
kruskal.test(SA_whole ~ Race_Language_Predictor, data = dd)


##### Mastery Achieved

### CDI Mastery
CDImast2 <- dd[dd$CDI_Mastery == 1,]

# Whole Sample
mean(CDImast2$SA_whole)
sd(CDImast2$SA_whole)

# Race Language Groups
WE <- CDImast2[CDImast2$Race_Language_Predictor == 0,]
mean(WE$SA_whole)
sd(WE$SA_whole)

LE <- CDImast2[CDImast2$Race_Language_Predictor == 1,]
mean(LE$SA_whole)
sd(LE$SA_whole)


LS <- CDImast2[CDImast2$Race_Language_Predictor == 2,]
mean(LS$SA_whole)
sd(LS$SA_whole)


# Kruskal Wallace Test
kruskal.test(SA_whole ~ Race_Language_Predictor, data = CDImast2)


### PDI Mastery
PDImast2 <- dd[dd$PDI_Mastery == 1,]

# Whole Sample
mean(PDImast2$SA_whole)
sd(PDImast2$SA_whole)

# Race Language Groups
WE <- PDImast2[PDImast2$Race_Language_Predictor == 0,]
mean(WE$SA_whole)
sd(WE$SA_whole)

LE <- PDImast2[PDImast2$Race_Language_Predictor == 1,]
mean(LE$SA_whole)
sd(LE$SA_whole)


LS <- PDImast2[PDImast2$Race_Language_Predictor == 2,]
mean(LS$SA_whole)
sd(LS$SA_whole)


# Kruskal Wallace Test
kruskal.test(SA_whole ~ Race_Language_Predictor, data = PDImast2)


