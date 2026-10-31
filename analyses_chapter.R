library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)
library(gridExtra)
library(RColorBrewer)
library(landmarkregUtils)
library(lmerTest)
library(ggpubr)
library(MuMIn)

#####
##### 0. INTRO
#####

# this is the directory of my df, you can replace with your directory
setwd("//phonetik.uni-muenchen.de/vdata/users_share_tmp/Quebec&FrenchNasals/Merged_Data")

D = readRDS("merged_data_dualchannel.RDS")

# This script represents a special version of "exLand.R" by Michele Gubian which you can find here: 
# https://github.com/uasolo/FPCA-phonetics-workshop/blob/master/scripts/

# for this analysis we had to exclude two speakers, nmf23 and nmf03:

levels(factor(D$speaker)) %>% length()

D = D %>% filter(!speaker %in%c("NMF23", "NMF03"))
#check, they should be now 55
levels(factor(D$speaker)) %>% length() # all right! 

# some graphical settings, which you can modify at any time

mytheme <- theme_light() +
  theme(text = element_text(size = 12)) # in Michele's script it's size = 16

Category.colors <- c("gold3", "darkgreen") # blue for NMF, green for QF

#####
##### 1. EXPLORE the data and CREATE the DFs
#####

names(D) %>% as.tibble() # how many and which variables?
head(D, n=3)

# In this analysis, we have 2 curves and duration, so 3 dims!

# citing Michele in his e-mail: (if you do not want to do the landmarking) "The curves are encoded using a 
# linearly time-normalized axis, the landmarks determine start and end in non-normalized units."

# for the raw curves we need $token (curveID), $timepoint_percent (although in Michele's script there are real timepoints) $nasalance (our y-value)
# change names of columns we want so that it fits to Michele's script, BUT I'll leave "nasalance"

# NB we have 2 curves per token! so we keep the column "token" and create curveID separating the oral from the nasal curve

raw_curves = D
raw_curves$curveId <- paste(raw_curves$token, raw_curves$channel, sep = "_")

#re-convert the normalized durations into real seconds (otherwise Michele's script does not work...)

raw_curves$time <- raw_curves$timepoint_percent * raw_curves$duration / 100

# NB Michele (23.10.2026) suggested to normalise the y-axis ("pascal") from 0 to 1000, to avoid huge stdev later
# at the end we discovered it does not change anything, so we won't use it for the chapter
raw_curves <- raw_curves %>%
  mutate(pascal_norm = (pascal - min(pascal, na.rm = TRUE)) /
           (max(pascal, na.rm = TRUE) - min(pascal, na.rm = TRUE)) * 10)

# for the duration, we need the df "land" with $token and the $landmarks: onset + $duration (in seconds)
# variable names "l1" (onset) and "l2" (offset)
# we need just for each token, so we select only "oral"

land = raw_curves %>% filter (time ==0.00) %>% rename("l1" = time, "l2" = duration)

# column l1 before l2 (otherwise landmark reg will not work)
land <- land %>% relocate(l1, .before = l2)

# "land" has now 1 row per per curveID (2 rows per token)

#####
##### 2. HAVE A LOOK AT 'raw_curves' and 'land' (it's just plots)
#####

# plot a few curves
set.seed(123)
subset_curveId <- raw_curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(raw_curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = pascal, group = curveId, color = dialect) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  facet_grid(~channel)+
  theme(legend.position = "bottom")

# in pascal_norm:
ggplot(raw_curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = pascal_norm, group = curveId, color = dialect) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  facet_grid(~channel)+
  theme(legend.position = "bottom")
# they look the same!

# with landmark position
land_y <- land %>%
  pivot_longer(l1:l2, names_to = "landmark", values_to = "time") %>% 
  inner_join(subset_curveId, by = "curveId") %>% 
  group_by(curveId) %>% 
  mutate(y = {
    curve <- raw_curves %>%
      inner_join(cur_group(), by = 'curveId') 
    approx(curve$time, curve$pascal_norm, time, rule = 2)$y
  }) 


ggplot(raw_curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = pascal_norm, group = curveId) +
  geom_line(linewidth = 0.6, color = 'slategray4') +
  geom_point(data = land_y,
             mapping = aes(time, y, color = landmark, group = curveId),
             inherit.aes = FALSE,
             size = 2) +
  scale_color_brewer(palette = "Dark2") + 
  mytheme  +
  facet_grid(~channel)+
  theme(legend.position = "bottom")


#####
##### 3. LANDMARK REGISTRATION (i.e. TIME-WARPING)
#####

# 1. compute time warping curves h(t), based on landmarks (not on curves!) 
reg <- landmarkreg_nocurves(inputMarks = land %>% select(starts_with("l")), 
                            targetMarks = c(l1 = 0, l2 = 1),  # duration from 0 to 1
                            njobs = 2)
# 2. apply time warping to curves
curvesReg <- applyReg(dat = raw_curves %>% filter(!is.na(pascal)), # remove NAs otherwise applyReg fails
                      reg = reg,
                      grid = seq(0, last(reg$landmarks), length.out = 100), # was also in Michele's script
                      id = "curveId", time = "time", value = "pascal") %>% 
  # to get back the other useful categorical variables (the ones you prefer)
  left_join(raw_curves %>% distinct(curveId, token, channel, dialect, speaker, vowel, word, 
                                    prec, foll, foll_voicing, syllable), by = "curveId")

# Explore reg
reg$h %>% plot()
reg$lograte %>% plot() #!!! these are the time-warped curves we are going to use to reconstruct durations!!!! #
reg$landmarks # position of aligned landmarks after time-warping, they are mean landmark positions
              # in this case 0 and 1 because I set them like this

# plot a few reg curves
ggplot(curvesReg %>% inner_join(subset_curveId, by = "curveId")) +
  aes(time, pascal, group = curveId, color = dialect) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  mytheme +
  facet_grid(~channel)+
  theme(legend.position = "bottom")


#####
##### 4. mFPCA - 3D: nasalance curves *and* segment DURATION
#####

# curve dimension - # build a funData object
# in this case, we build two separate objects - 'cause we have two curves

curvesReg_or = curvesReg %>% filter(channel == "oral")
curvesFun_or <- long2irregFunData(curvesReg_or, id = "curveId", time = "time", value = "pascal") %>% 
  as.funData()

curvesReg_nas = curvesReg %>% filter(channel == "nasal")
curvesFun_nas <- long2irregFunData(curvesReg_nas, id = "curveId", time = "time", value = "pascal") %>% 
  as.funData()

# check that order of raws in "raw_curves"2" is the same in "land", otherwise 
#the results are not good
ids_input <- raw_curves %>%
  filter(!is.na(pascal)) %>%
  distinct(curveId) %>%
  pull(curveId) %>%
  as.character()

identical(ids_input, as.character(land$curveId))
#TRUE -> good!

# Get token order in the two channels
ids_or <- curvesReg_or %>%
  ungroup() %>%
  distinct(curveId, token) %>%
  pull(token) %>%
  as.character()

ids_nas <- curvesReg_nas %>%
  ungroup() %>%
  distinct(curveId, token) %>%
  pull(token) %>%
  as.character()

# Check that the tokens are in the same order
identical(ids_or, ids_nas)
# true

# Check that each token occurs only once per channel
anyDuplicated(ids_or)
anyDuplicated(ids_nas)

# 0

# Check correspondence with funData observations
length(ids_or) == nObs(curvesFun_or)
length(ids_nas) == nObs(curvesFun_nas)

# all tests fine!

# now, create a 3D MultiFunData object in 6 further steps,
# being sure that order of the tokens is correct and
# that the durations are there only once per token (so no double rows)

# 1. Create funData with ALL lograte functions (similarly to Michele's orginal script)
lograteFun_all <- fd2funData(
  reg$lograte,
  argvals = curvesFun_or@argvals[[1]] # we took one of the two channels
)

# 2. Find the first occurrence of each token in land
idx <- match(ids_or, as.character(land$token))

# 3. another test ( you could also skipthis)

stopifnot(
  !anyNA(idx),
  nObs(lograteFun_all) == nrow(land)
)

# 4. Select only one lograte function per token
lograteFun <- funData::extractObs(
  lograteFun_all,
  obs = idx
)

# 5. Check the number of observations, should be all the same
c(
  oral = nObs(curvesFun_or),
  nasal = nObs(curvesFun_nas),
  duration = nObs(lograteFun)
)

# Fine!

# 6. put together a 3D fd object: 2 curves and lograte (=duration)

yRegMult <- multiFunData(list(curvesFun_or, curvesFun_nas, lograteFun))

# have a look at one curve/lograte curve
yRegMult[140] %>% plot()  

# compute integral of uni-dim variances, will be used to weight dimensions. # NB: NOT IN THE TUTORIAL
uniDvars <- sapply(c(curvesFun_or, curvesFun_nas, lograteFun), \(f) {
  f %>% PACE() %>% `$`("estVar") %>% integrate()
}, simplify = TRUE) # still great variances though, in spite of normalising Pascals

# multidim FPCA --- video-tutorial from 2:39:40, some differences between updated code and tutorial version
mfpca <- MFPCA(yRegMult,
               M = 5, # n. of PCs to calculate: set M > nPC only to get a good approx of explained var below
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"), list(type = "uFPCA")),
               weights = 1/uniDvars # each dim divided by its variance,  # added later, NOT IN THE TUTORIAL
               # see Happ & Greven 2018, section 5
)
# Prop of explained var
mfpca$values  / sum( mfpca$values)
#[1] 0.41750445 0.25402421 0.22506876 0.05466110 0.04874147  # same as before normalising Pascal scale

#so the first 2 explain around 76% of variance in the data

nPC <- 3 # we will look at the first 3 nPC components
# scores st. dev.
sdFun <- mfpca$values %>% sqrt()
fractionOfStDev <-  seq(-1, 1, by=.25)
# PC curves to be plotted
DimCurves <- c(1,2)
DimLograte <- 3 # lograte is the third dimension in mfpca

#here you build DF to plot PCs. Michele left it "open" because you could adjust it to your taste
#code lines written in tidyverse
PCcurves <- expand_grid(PC = 1:nPC,
                        Dim = DimCurves , #c(DimCurves), # DimLograte),
                        fractionOfStDev = fractionOfStDev) %>%
  group_by(PC, Dim, fractionOfStDev) %>%
  reframe(time = mfpca$meanFunction[[Dim]]@argvals[[1]],
          y = (mfpca$meanFunction[[Dim]] + # it's the mean curve, plus coefficient multiplied to PCs (line below)
                 fractionOfStDev * sdFun[PC] * mfpca$functions[[Dim]][PC])@X %>% as.numeric()
  )

# Plot
DimLabels <- c(`1` = "oral", `2` = "nasal")

PCcurves_plot <- ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line(linewidth = 0.8) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0),
            color = "black", linewidth = 1.2) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1, suffix = "")) +
  facet_grid(Dim ~ PC,
             scales = "free_y",
             labeller = labeller(PC  = ~ str_glue("PC{.x}"),
                                 Dim = DimLabels)) +
  labs(color = expression(paste(s[k], ", z-scores"))) +
  xlab("Normalized time (%)") +
  ylab("Channel (Pa)") +
  mytheme +
  theme(legend.position = "top",
        strip.text.x = element_text(color = "black", size = 12),
        strip.text.y = element_text(color = "black", size = 12)) +
  ggtitle ("PC1: 42% variance, PC2: 25%, PC3: 22%")

PCcurves_plot

#################################################################################################
# OTHER TYPE of VISUALIZATION - curves + duration in one plot - used in the paper

# Plot PC curves against original time axis + landmark positions
# (this is an example for 1-dim curves)

PCdur <- expand_grid(PC = 1:nPC,
                     fractionOfStDev = fractionOfStDev,
                     landmarks2long(reg$landmarks)
) %>%
  group_by(PC, fractionOfStDev, leftBoundary) %>%
  mutate(duration = lograte2duration(  #function by package landmarkregUtils, converting curve to duration again
    lograte = mfpca$meanFunction[[DimLograte]] +
      fractionOfStDev * sdFun[PC] * mfpca$functions[[DimLograte]][PC],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))

# add landmark position (except for the first one at time = 0)
PCdur <- PCdur %>% 
  group_by(PC, fractionOfStDev) %>% 
  mutate(land = cumsum(duration))

# compute original time samples from registered time samples
PCcurves <- PCcurves %>% 
  # inner_join(PCdur %>% distinct(PC, fractionOfStDev), by = c("PC", "fractionOfStDev")) %>% 
  group_by(PC, Dim, fractionOfStDev) %>% 
  mutate(origTime = {
    landmarks <- c(0, # first landmark at time = 0
                   PCdur %>% 
                     inner_join(cur_group(), by = c("PC", "fractionOfStDev")) %>% 
                     pull(land))
    landmarkreg_timeSamples(timeSamples = time,
                            inputMarks = reg$landmarks,
                            targetMarks = landmarks)
  })

# compute (time, y) coordinates of landmark points on predicted curves
PCland <- PCdur %>% 
  group_by(PC, fractionOfStDev) %>% 
  select(PC, fractionOfStDev, land) %>% 
  reframe(land = list(c(0, land))) %>% 
  inner_join(PCcurves %>% select(PC, fractionOfStDev, origTime, y),
             by =  c("PC", "fractionOfStDev"),
             multiple = "all") %>% 
  nest(.by = c(PC, fractionOfStDev, land), .key = "curves") %>% 
  group_by(PC, fractionOfStDev) %>%
  mutate(y = list(approx(curves[[1]]$origTime, curves[[1]]$y, land[[1]], rule = 2)$y)) %>% 
  select(!curves) %>% 
  unnest(cols = c(land, y))

PCcurves_origTime_plot <- ggplot(PCcurves) +
  aes(x = origTime, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line(size=0.8) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid( Dim ~ PC,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"),
                                 Dim = DimLabels)) +
  #geom_point(data = PCland,
     #        mapping = aes(land, y, color = fractionOfStDev),
      #       inherit.aes = FALSE,
      #       size = 3
  #) +
  labs(color = expression(paste(s[k], ", z-scores"))) +
  xlab("Time (seconds)") + 
  ylab("Channel (Pa)") +
  mytheme +
  theme(legend.position = "top",
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12))

PCcurves_origTime_plot
######################################################################


# PC durations to be plotted


PCdur_plot <- ggplot(PCdur %>% filter(fractionOfStDev %in% seq(-1, 1, by = .5))) +
  aes(x = fractionOfStDev %>% factor(labels = ""),
      y = duration, fill = fractionOfStDev) + 
  geom_bar(stat="identity", color = "white", alpha = 0.7, width = 0.98) +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 4, color="white",
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_grid(~ PC, labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "orangered") +
  labs(color = expression(frac(s[k], sigma[s[k]]))) +
  ylab("Time (seconds)") +
  xlab("Duration")+
  coord_flip() + 
  mytheme +
  theme(legend.position = "none",
        strip.text.x = element_text(color = "black", size=12),
        strip.text.y = element_text(color = "black", size=12),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank())

PCdur_plot

#grid.arrange(PCcurves_plot, PCdur_plot, nrow=2) ---- alternative FOR THE ICPHS PAPER:

PCcurves_plot <- PCcurves_plot +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

PCdur_plot <- PCdur_plot +
  theme(plot.margin = margin(t = 8, r = 5, b = 5, l = 5))

ggarrange(PCcurves_plot, PCdur_plot,
          nrow = 2,
          heights = c(1.3, 0.8),
          align = "v")

#####
##### 5. Analysis of PC-scores
#####

# collect PC scores
PCscores <- mfpca$scores[, 1:nPC, drop=FALSE] %>%
  `colnames<-`( paste0("s", 1:nPC)) %>%
  as_tibble() %>%
  bind_cols(raw_curves %>% distinct(token, dialect, speaker, vowel, word,  #NB this time *token*, not curveID
                                    prec, foll, foll_voicing, syllable), .)
# variables you can explore:
names(PCscores) %>% as.tibble()

## SYLLABLE

# scatterplot PC scores s1 and s2 by dialect and syllable (close vs open)
ggplot(PCscores) +
  aes(x = s1, y = s2, color = dialect) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  facet_grid(~syllable)+
  ylim(-2.5, 2.5)+
  xlim(-3, 3)+
  mytheme +
  theme(legend.position = "right")

# scatterplot PC scores s1 and s3 by dialect and syllable (close vs open)
ggplot(PCscores) +
  aes(x = s1, y = s3, color = dialect) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  facet_grid(~syllable)+
  ylim(-2.5, 2.5)+
  xlim(-3, 3)+
  mytheme +
  theme(legend.position = "right")

# boxplots PC scores by by dialect and syllable
PCscores %>% 
  pivot_longer(cols = s1:all_of(str_glue("s{nPC}")), 
               names_to = "PC", values_to = "score") %>% 
  ggplot() +
  aes(x = dialect, y = score, color = dialect) +
  geom_violin() +
  geom_boxplot(width = 0.4) + #set "outlier.shape = NA" to delete outlier dots
  facet_grid(syllable~PC) +
  ylim(-3, 3)+
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")

### FOLLOWING VOICING (only closed syllables) #- NO DIFFS

ggplot(PCscores %>% filter(syllable=="closed")) +
  aes(x = s1, y = s2, color = dialect) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  facet_grid(syllable~foll_voicing)+
  ylim(2.5, -2.5)+
  mytheme +
  theme(legend.position = "right")

PCscores %>% filter(syllable=="closed") %>%
  pivot_longer(cols = s1:all_of(str_glue("s{nPC}")), 
               names_to = "PC", values_to = "score") %>% 
  ggplot() +
  aes(x = dialect, y = score, color = dialect) +
  geom_violin() +
  geom_boxplot(width = 0.4) + #set "outlier.shape = NA" to delete outlier dots
  facet_grid(PC ~ foll_voicing) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")


#####
##### 6. Statistics of PC-scores
#####

# NB: we run here an anova on output from lmerTest, BUT to reconstruct the curves we use lme4!
# HOWEVER I tested this and it looks like the output from emmeans does not vary

##################
### syllable type

### s1

s1_syll.lmer= lmer(s1 ~ dialect * syllable +
                 (1 | word) + (syllable|speaker), 
                 data = PCscores)

anova(s1_syll.lmer) # all highly sig.

# Type III Analysis of Variance Table with Satterthwaite's method
#                  Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# dialect          19.667  19.667     1  81.449  41.535 7.638e-09 ***
# syllable         22.209  22.209     1  26.040  46.904 2.839e-07 ***
# dialect:syllable 13.094  13.094     1 158.611  27.654 4.631e-07 ***

s1_syll.lmer %>% emmeans(.,  ~ syllable|dialect) %>%
  pairs(., simple = "each", combine=T) %>%
  as.data.frame() %>%
  filter(p.value <= 0.05) %>%
  as_tibble() %>%
  mutate(p.value = round(p.value, 4)) %>%
  print(., n = nrow(.), max_extra_cols = ncol(.))

# A tibble: 4 × 8
# dialect syllable contrast      estimate    SE    df t.ratio p.value
# <chr>   <chr>    <chr>            <dbl> <dbl> <dbl>   <dbl>   <dbl>
#   1 NMF     .        closed - open    0.432 0.167  60.8    2.58  0.0487
# 2 QF      .        closed - open    1.34  0.145  45.2    9.27  0     
# 3 .       closed   NMF - QF        -1.35  0.161 103.    -8.35  0     
# 4 .       open     NMF - QF        -0.440 0.168 100.0   -2.62  0.0401


s1_syll.lmer %>% 
  emmeans(~ syllable*dialect) %>% 
  contrast(interaction = "pairwise")

# ALSO TRY THIS (suggested by Josie) for more complete results
s1_syll.lmer %>% emmeans(., pairwise~ syllable*dialect)

#### s2# to be CONTINUED

s2_syll.lmer = lmer(s2 ~ dialect * syllable + #or lme4::lmer() etc
                      (1 | word) + (syllable|speaker), 
                    data = PCscores)

anova(s2_syll.lmer) 

s2_syll.lmer %>% emmeans(.,  ~ syllable|dialect) %>%
  pairs(., simple = "each", combine=T) %>%
  as.data.frame() #%>%
 # filter(p.value <= 0.05) %>%
#  as_tibble() %>%
 # mutate(p.value = round(p.value, 4)) %>%
 # print(., n = nrow(.), max_extra_cols = ncol(.))



s2_syll.lmer %>% 
  emmeans(~ syllable*dialect) %>% 
  contrast(interaction = "pairwise")

# ALSO TRY THIS (suggested by Josie) for more complete results
s2_syll.lmer %>% emmeans(., pairwise~ syllable*dialect)

#### s3

s3_syll.lmer = lmer(s3 ~ dialect * syllable + #or lme4::lmer() etc
                      (1 | word) + (syllable|speaker), 
                    data = PCscores)

anova(s3_syll.lmer) 


s3_syll.lmer %>% emmeans(., pairwise~ syllable*dialect)




# Pseudo-R2:
r.squaredGLMM(s1_syll.lmer)
#            R2m      R2c
#[1,] 0.4110373 0.643133

r.squaredGLMM(s2_syll.lmer)
#           R2m       R2c
# [1,] 0.2398415 0.5046084

r.squaredGLMM(s3_syll.lmer)
#            R2m       R2c
# [1,] 0.01912389 0.7390992

### to be continued
############################### 

#####
##### 7. Reconstruct predicted curves on the original time axis
#####

# based on dialect and syllable type, using s1 scores only + duration

## first predicted curves 

predCurves <- emmeans(s1_syll.lmer, pairwise ~ syllable|dialect)$emmeans %>%
  as_tibble() %>%
  expand_grid(Dim = DimCurves) %>% 
  group_by(syllable, dialect, Dim) %>% 
  reframe(bind_cols(
    funData2long1(mfpca$meanFunction[[Dim]] +  #again, mean curve plus scores predicted by LMER multipl. by the PCs
                    emmean * mfpca$functions[[Dim]][1], value = "y"),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    lower.CL * mfpca$functions[[Dim]][1], value = "yl") %>% #l=lower CL
      select(yl),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    upper.CL * mfpca$functions[[Dim]][1], value = "yu") %>% #u=upper CL
      select(yu)
  ))



## then predicted durations

predDur <- expand_grid(emmeans(s1_syll.lmer, pairwise ~ syllable|dialect)$emmeans %>% 
                         as_tibble()
                       , landmarks2long(reg$landmarks)) %>% 
  group_by(syllable, dialect, leftBoundary) %>% 
  mutate(duration = lograte2duration(
    lograte = mfpca$meanFunction[[DimLograte]] +
      emmean * mfpca$functions[[DimLograte]][1],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))


# add landmark position (except for the first one at time = 0)
predDur <- predDur %>% 
  group_by(syllable, dialect) %>% 
  mutate(land = cumsum(duration)) # cumulative sum # in my case, the same as duration

# compute original time samples from registered time samples
predCurves <- predCurves %>% 
  group_by(syllable, dialect) %>% 
  mutate(origTime = {
    landmarks <- c(0, # first landmark at time = 0
                   predDur %>% 
                     inner_join(cur_group(), by = c("dialect", "syllable")) %>% 
                     pull(land))
    landmarkreg_timeSamples(timeSamples = time,
                            inputMarks = reg$landmarks,
                            targetMarks = landmarks)
  })

# compute (time, y) coordinates of landmark points on predicted curves
predLand <- predDur %>% 
  group_by(dialect, syllable) %>% 
  select(dialect, syllable, land) %>% 
  reframe(land = list(c(0, land))) %>% 
  inner_join(predCurves %>% select(dialect, syllable, origTime, y),
             by = c("dialect", "syllable"),
             multiple = "all") %>% 
  nest(.by = c(dialect, syllable, land), .key = "curves") %>% 
  group_by(dialect, syllable) %>%
  mutate(y = list(approx(curves[[1]]$origTime, curves[[1]]$y, land[[1]], rule = 2)$y)) %>% 
  select(!curves) %>% 
  unnest(cols = c(land, y))

####

syll.colors=c("purple", "slategrey")

ggplot(predCurves) +
  aes(origTime, y, color = syllable) +
  geom_line(linewidth = 0.8) +
  geom_point(data = predLand,
             mapping = aes(land, y, color = syllable),
             inherit.aes = FALSE,
             size = 3
  ) +
  facet_grid(~dialect)+
  scale_color_manual(values=syll.colors) +
  xlab("Time (s)") +
  ylab("Nasalance")+
  # ggtitle(str_glue("Reconstructed curves according to regr model: s{s} ~ Category")) +
  mytheme +
  theme(legend.position = "bottom")

#####
##### 8. Reconstruct predicted curves on the original time axis - 2 PCs
#####

## see CLAUDE-AI script "ricostruzione_curve_PC1_PC2.R"

### NB: IT SHOULD BE VERIFIED BY MICHELE, but to me it looks fine

############################### 

#####
##### 8. Numbers for the method
#####


levels(factor(D$word)) %>%length() #27

table(D$word, D$syllable)

levels(factor(D$token)) %>% length()
#[1] 1913

EF = D %>% filter(dialect == "EF"
levels(factor(EF$token)) %>% length()
#[1] 776
1913-776
#1137
levels(factor(D$speaker)) %>% length()
#57
levels(factor(EF$speaker)) %>% length()
#27

# how many open/closed syllables by dialect?
D %>% 
  filter(dialect == "EF" & syllable == "open") %>% 
  distinct(token) %>% 
  nrow()
#[1] 381

D %>% 
  filter(dialect == "EF" & syllable == "closed") %>% 
  distinct(token) %>% 
  nrow()
#[1] 395

D %>% 
  filter(dialect == "QF" & syllable == "open") %>% 
  distinct(token) %>% 
  nrow()
#[1] 418

D %>% 
  filter(dialect == "QF" & syllable == "closed") %>% 
  distinct(token) %>% 
  nrow()

#[1] 719
