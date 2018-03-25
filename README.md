# R helper functions

Helper functions to make it easier to analyse and summarise data and results in R. Use at your own risk!!! If anything doesn't work or is wrong, let me know and I'll try to fix it as soon as possible (create an issue or email me at hauselin@gmail.com).

- [Summarise statistical models plus effect sizes](#summarise-statistical-models-plus-effect-sizes)
- [Convert between effect sizes](#convert-between-effect-sizes)
- [Detect and remove outliers](#detect-and-remove-outliers)
- [Compute between- and within-subjects standard errors and confidence intervals](#compute-between--and-within-subjects-standard-errors-and-confidence-intervals)

## Summarise statistical models plus effect sizes

When fitting models, use ```summaryh()``` instead of ```summary()``` to get APA (American Psychological Association) formatted output that also includes **effect size estimates for each effect** (*r* effect size).

Currently accepts models fitted with these functions: ```lm```, ```anova```, ```aov```, ```chisq.test```, ```cor.test```, ```glm```, ```lmer```, ```lme```, ```t.test```.

To use/download ```summaryh```, run this line of code: ```source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R")```. The first time you run this line of code, it will take some time because it's going to install a few useful R packages. Subsequently, it should load the functions much faster.

Example outputs

* ```summaryh(lm(mpg ~ qsec, mtcars))```: b = 1.41, SE = 0.56, t(30) = 2.53, p = .017, r = 0.42
* ```summaryh(aov(mpg ~ gear, mtcars))```: F(1, 30) = 9.00, p = .005, r = 0.48)
* ```summaryh(cor.test(mtcars$mpg, mtcars$gear))```: r(30) = 0.48, p = .005
* ```summaryh(t.test(mpg ~ vs, mtcars))```: t(23) = −4.67, p < .001, r = 0.70

Arguments in ```summaryh(model, decimal = 2, showTable = F, showEffectSizesTable = F)```

* model: fitted model (required)
* decimal (default = 2): decimal places of output
* showTable (default = F): show the results in table format
* showEffectSizesTable (default = F): show other effect sizes computed using ```es``` function (see sections below) (d, r, R<sup>2</sup>, f, odds ratio, log odds ratio, area under curve)

```
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R")

# linear regression
summary(lm(mpg ~ cyl, mtcars)) # base R summary()
summaryh(lm(mpg ~ cyl, mtcars)) # returns APA-formatted output
summaryh(lm(mpg ~ cyl, mtcars), decimal = 5, showTable = T, showEffectSizesTable = T) # optional arguments

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

```
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/es.R")

es(d = 0.3) # convert Cohen's d to other effect sizes
es(r = c(0.1, 0.3, 0.5)) # convert multiple effect sizes (r) to other effect sizes
```

## Detect and remove outliers
We can identify and remove outliers in our data by identifying data points that are too extreme—either too many standard deviations (SD) away from the mean or too many median absolute deviations (MAD) away from the median. The SD approach might not be ideal with extreme outliers, whereas the MAD approach is much more robust (for comparison of both approaches, see [Leys et al., 2013, Journal of Experimental Social Psychology](https://s3.amazonaws.com/academia.edu.documents/32918779/Leys_MAD_final_copy.pdf?AWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1522001231&Signature=RL2IhDQKFM8W9z32xELWJE%2BGWBM%3D&response-content-disposition=inline%3B%20filename%3DDetecting_outliers_Do_not_use_standard_d.pdf)).

```detectOutliers``` script contains two functions: ```outliersZ``` and ```outliersMAD```

```
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

```
# load functions from my github site
source("https://raw.githubusercontent.com/hauselin/Rcode/master/se.R")

# between-subjects standard error and confidence intervals (95% default)
se(data = mtcars, measurevar = "disp", groupvars = "cyl")
se(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "am"))
se(data = ChickWeight, measurevar = "weight", groupvars = "Diet")

# within-subjects standard error and confidence intervals (95% default)
seWithin(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
```
