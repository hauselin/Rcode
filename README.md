# R helper functions

Helper functions to make it easier to analyse and summarise data and results in R. Use at your own risk!!! If anything doesn't work or is wrong, let me know and I'll try to fix it as soon as possible (create an issue or email me at hauselin@gmail.com).

- [Summarise statistical models plus effect sizes](#summarise-statistical-models-plus-effect-sizes)
- [Convert between effect sizes](#convert-between-effect-sizes)
- [Detect and remove outliers](#detect-and-remove-outliers)
- [Compute between- and within-subjects standard errors and confidence intervals](#compute-between--and-within-subjects-standard-errors-and-confidence-intervals)
- [Fit Wagenmaker's EZ-diffusion model for two-choice response-time tasks](#fit-ez-diffusion-model-for-two-choice-response-time-tasks)
- [Fit drift-diffusion model for two-choice response time tasks using maximum likelihood (estimates a, v, t0, z)](https://github.com/hauselin/Rcode#fit-drift-diffusion-model-for-two-choice-response-time-tasks-using-maximum-likelihood-estimation-parameters-a-v-t0-z)

## Summarise statistical models plus effect sizes

When fitting models, use ```summaryh()``` in place of ```summary()``` to get APA (American Psychological Association) formatted output that also includes **effect size estimates for each effect** (*r* effect size).

Currently accepts models fitted with these functions: ```lm```, ```anova```, ```aov```, ```chisq.test```, ```cor.test```, ```glm```, ```lmer```, `glmer`, ```lme```, ```t.test```. For ```lmer``` models, p-values must have been computed with ```lmerTest```.

To use ```summaryh```, run this line of code each time you start a new R session: ```source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R")```. The first time you run this line of code, it will take some time because it's going to install a few R packages. Subsequently, it should load the functions much faster.

Arguments in ```summaryh(model, decimal = 2, showTable = F, showEffectSizesTable = F)```

* **model** (required): fitted model
* **decimal** (default = 2): decimal places of output
* **showTable** (default = F): show the results in table format
* **showEffectSizesTable** (default = F): show other effect sizes computed using ```es``` function (see sections below) (d, r, R<sup>2</sup>, f, odds ratio, log odds ratio, area under curve)

Example outputs (output is data.table class)

- ```summaryh(lm(mpg ~ qsec, mtcars))```: b = 1.41, SE = 0.56, t(30) = 2.53, p = .017, r = 0.42
- ```summaryh(aov(mpg ~ gear, mtcars))```: F(1, 30) = 9.00, p = .005, r = 0.48
- ```summaryh(cor.test(mtcars$mpg, mtcars$gear))```: r(30) = 0.48, p = .005
- ```summaryh(t.test(mpg ~ vs, mtcars))```: t(23) = −4.67, p < .001, r = 0.70
- ```summaryh(glm(vs ~ 1, mtcars, family = "binomial"))```: b = −0.25, SE = 0.36, z(31) = −0.71, p = .481, r = −0.07

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R")

# linear regression
model_lm <- lm(mpg ~ cyl, mtcars) 
summary(model_lm) # base R summary()
summaryh(model_lm) # returns APA-formatted output in a data.table (output below)
##           term                                                 results
## 1: (Intercept) b = 37.88, SE = 2.07, t(30) = 18.27, p < .001, r = 0.96
## 2:         cyl b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = 0.85 
summaryh(model_lm, decimal = 5, showTable = T, showEffectSizesTable = T) # optional arguments

# linear mixed effects regression
library(lme4); library(lmerTest) # load packages to fit mixed effects models
model <- lmer(weight ~ Time * Diet  + (1 + Time | Chick), data = ChickWeight)
summary(model)
summaryh(model, decimal = 4, showTable = T, showEffectSizesTable = T) # optional arguments

# ANOVA
summaryh(aov(mpg ~ gear, mtcars))

# correlation
cor.test(mtcars$mpg, mtcars$cyl)
summaryh(cor.test(mtcars$mpg, mtcars$cyl))
summaryh(cor.test(mtcars$mpg, mtcars$cyl), 3, T, T)
```

## Convert between effect sizes

The ```es``` function converts one effect size into other effect sizes (e.g., d, r, R<sup>2</sup>, f, odds ratio, log odds ratio, area-under-curve [AUC]). Note that AUC calculations may not be correct!

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/es.R")

es(d = 0.3) # convert Cohen's d to other effect sizes
es(r = c(0.1, 0.3, 0.5)) # convert multiple effect sizes (r) to other effect sizes
```

## Detect and remove outliers
We can identify and remove outliers in our data by identifying data points that are too extreme—either too many standard deviations (SD) away from the mean or too many median absolute deviations (MAD) away from the median. The SD approach might not be ideal with extreme outliers, whereas the MAD approach is much more robust (for comparison of both approaches, see [Leys et al., 2013, Journal of Experimental Social Psychology](https://s3.amazonaws.com/academia.edu.documents/32918779/Leys_MAD_final_copy.pdf?AWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1522001231&Signature=RL2IhDQKFM8W9z32xELWJE%2BGWBM%3D&response-content-disposition=inline%3B%20filename%3DDetecting_outliers_Do_not_use_standard_d.pdf)).

```detectOutliers``` script contains two functions: ```outliersZ``` and ```outliersMAD```

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/detectOutliers.R")

example <- c(1, 3, 3, 6, 8, 10, 10, 1000) # 1000 is an outlier
outliersZ(example) # SD approach
outliersMAD(example) # MAD approach

# changing a few default values (note that if zCutOff is 3, 1000 ISN'T considered an outlier!)
outliersZ(x = example, zCutOff = 3, replaceOutliersWith = -999) # common to use 1.96, 2.5, 3 for Z cutoff
outliersMAD(x = example, MADCutOff = 3, replaceOutliersWith = -888) # Leys et al. (2003) recommends 2.5 to 3 for MAD cutoff
```

## Compute between- and within-subjects standard errors and confidence intervals

Code adapted from <a href="http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/">Cookbook for R</a>

When using the functions ```se``` (between-subjects) or ```seWithin``` (within-subjects), you can specify more than one outcome variable via the ```measurevar``` argument. If you specify more than one outcome variable, the output will be a list that has length of the number of outcome variables provided.

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/se.R")

# between-subjects standard error and confidence intervals (95% default)
se(data = mtcars, measurevar = "disp", groupvars = "cyl")
se(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "am"))
se(data = ChickWeight, measurevar = "weight", groupvars = "Diet")

# within-subjects standard error and confidence intervals (95% default)
seWithin(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
```

## Fit EZ-diffusion model for two-choice response time tasks
```fit_ezddm``` function fits [Wagenmaker et al.'s (2007)](https://link.springer.com/article/10.3758/BF03194023) EZ-diffusion model for two-choice response time tasks. To use the function, ensure your dataframe is in long form, has single-trial reaction time (in seconds) and accuracy (coded as 0 or 1) on each row. You can use the function to fit the EZ-diffusion model to just a single subject or multiple subjects, and separately for each experimental condition (see below for examples).

Assumptions of EZ-diffusion model

- error and correct reaction-time distributions are identical (often violated!)
- z = .5: starting point is equidistant from the response boundaries
- sv = 0: inter-trial variability in drift rate is negligible
- sz = 0: inter-trial variability in starting point is negligible
- st = 0: inter-trial range in nondecision time is negligible

To use/download ```fit_ezddm```, run this line of code: ```source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ezddm.R")```. The first time you run this line of code, it will take some time because; subsequently, it should load the functions much faster.

Arguments in ```fit_ezddm(data, rts, responses, id = NULL, group = NULL)```

* **data** (required): data object with reaction time and accuracy variables (long form data expected)
* **rts** (required; in seconds): specify in characters the name of the reactiontime column
* **responses** (required; coded as 0/1): specify in characters the name of the accuracy column
* **id** (default = NULL): specify in characters the name of your subject/id column (if not specified, assumes data [all rows] belong to a single subject)
* **group** (default = NULL): specify in characters the name of your column(s) indicating various conditions

Output (data.table class)

* subject: returns this variable only there's more than one subject
* group/condition names: returns these variables only if you specify grouping variables
* n: number of trials
* a_threshold: boundary/threshold
* v_drift: drift rate/evidence accumulation rate
* ndt_Ter: non-decision time
* rt_correct: mean reaction time for correct trials (used by ezddm to compute parameters)
* rtVar_correct: reaction time variance for correct trials (used by ezddm to compute parameters)
* acc: mean accuracy or proportion of upper bound (1) responses (used by ezddm to compute parameters)
* acc_adjust: indicates if mean accuracies (acc) have been adjusted; ezddm can't estimate parameters if mean accuracy is exactly 0.5 or 1.0; if acc_adjust is 0, no adjustments have been made; if acc_adjust is 1, minor adjustments have been made

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ezddm.R")

# simulate some data
library(rtdists)
data1 <- rdiffusion(n = 100, a = 2, v = 0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
data2 <- rdiffusion(n = 100, a = 2, v = -0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
dataAll <- rbind(data1, data2) # join data
dataAll$response <- ifelse(dataAll$response == "upper", 1, 0) # convert responses to 1 and 0
dataAll$subject <- rep(c(1, 2), each = 100) # assign subject id
dataAll$cond1 <- sample(c("a", "b"), 200, replace = T) # randomly assign conditions a/b
dataAll$cond2 <- sample(c("y", "z"), 200, replace = T) # randomly assign conditions y/z

# fit model to just entire data set (assumes all data came from 1 subject)
fit_ezddm(data = dataAll, rts = "rt", responses = "response")
# fit model to each subject (no conditions)
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject") 
# fit model to each subject by cond1
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1") 
# fit model to each subject by cond1,cond2
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2"))
```

## Fit drift-diffusion model for two-choice response time tasks using maximum likelihood estimation (parameters: a, v, t0, z)

`fit_ddm` function fits four-parameter (a, v, t0, z) drift diffusion model (also known as Wiener diffusio model) to two-choice response time tasks using maximum likelihood estimation (R `ucminf` optimization). Assumes no or negligible inter-trial variability in drift rate (sv), starting point (sz), and non-decision time (st)—these parameters aren't not estimated by `fit_ddm`.

To use/download ```fit_ddm```, run this line of code: ```source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ddm.R")```. The first time you run this line of code, it will take some time because; subsequently, it should load the functions much faster.

Arguments in ```fit_ddm(data, rts, responses, id = NULL, group = NULL)```

- **data** (required): data object with reaction time and accuracy variables (long form data expected)
- **rts** (required; in seconds): specify in characters the name of the reactiontime column
- **responses** (required; coded as 0/1 or "lower"/"upper"): specify in characters the name of the accuracy column
- **id** (default = NULL): specify in characters the name of your subject/id column (if not specified, assumes data [all rows] belong to a single subject)
- **group** (default = NULL): specify in characters the name of your column(s) indicating various conditions
- **startParams** (default = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)): starting parameters for likelihood estimation with `ucminf`

Output (data.table class)

- subject: returns this variable only there's more than one subject
- group/condition names: returns these variables only if you specify grouping variables
- n: number of trials
- a: boundary/threshold
- v: drift rate/evidence accumulation rate
- t0: non-decision time
- z: starting-point bias (0.5 is no bias)
- convergence: reason for optimization termination (see `?ucminf`)
- value: objective function value at computed miminizer (see `?ucminf`)

```R
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ddm.R")

# simulate some data
library(rtdists)
data1 <- rdiffusion(n = 100, a = 2, v = 0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
data2 <- rdiffusion(n = 100, a = 2, v = -0.3, t0 = 0.5, z = 0.5 * 2) # simulate data
dataAll <- rbind(data1, data2) # join data
dataAll$subject <- rep(c(1, 2), each = 100) # assign subject id
dataAll$cond1 <- sample(c("a", "b"), 200, replace = T) # randomly assign conditions a/b
dataAll$cond2 <- sample(c("y", "z"), 200, replace = T) # randomly assign conditions y/z

# fit model to just entire data set (assumes all data came from 1 subject)
fit_ddm(data = dataAll, rts = "rt", responses = "response")
# fit model to each subject (no conditions)
fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject") 
# fit model to each subject by cond1
fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1") 
# fit model to each subject by cond1,cond2
fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2"))
```