# R helper functions

- [Summarise statistical models plus effect sizes](#summarise-statistical-models-plus-effect-sizes)
- [Convert between effect sizes](#convert-between-effect-sizes)
- [Detect and remove outliers](#detect-and-remove-outliers)
- [Compute between- and within-subjects standard errors and confidence intervals](#compute-between--and-within-subjects-standard-errors-and-confidence-intervals)

## Summarise statistical models plus effect sizes

When fitting models, use ```summaryh()``` instead of ```summary()``` to get APA (American Psychological Association) formatted output that includes effect size estimates for each effect. The first time you ```source()``` the the functions, it might take some time because it will install a few R packages (e.g., dplyr). Subsequently, ```source()``` should load the functions much faster.

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R") # load functions from my github site

# linear regression
summary(lm(mpg ~ cyl, mtcars)) # base R summary()
summaryh(lm(mpg ~ cyl, mtcars)) # returns APA-formatted output
summaryh(lm(mpg ~ cyl, mtcars), decimal = 5, showTable = T, showEffectSizesTable = T) # optional arguments

# linear mixed effects regression
library(lme4); library(lmerTest) # load packages to fit mixed effects models
model <- lmer(weight ~ Time * Diet  + (1 + Time | Chick), data = ChickWeight)
summary(model)
summaryh(model, decimal = 4, showTable = T, showEffectSizesTable = T) # optional arguments

# correlation
cor.test(mtcars$mpg, mtcars$cyl)
summaryh(cor.test(mtcars$mpg, mtcars$cyl))
summaryh(cor.test(mtcars$mpg, mtcars$cyl), 3, T, T)
```

## Convert between effect sizes

The ```es``` function converts one effect size into other effect sizes (e.g., d, r, R<sup>2</sup>, f, odds ratio, log odds ratio, area-under-curve [AUC]). Note that AUC calculations may not be correct!

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/es.R") # load functions from my github site

es(d = 0.3) # Cohen's d
es(r = c(0.1, 0.3, 0.5)) # r
```

## Detect and remove outliers
We can identify and remove outliers in our data by identifying data points that are too extreme—either too many standard deviations (SD) away from the mean or too many median absolute deviations (MAD) away from the median. The SD approach might not be ideal with extreme outliers, whereas the MAD approach is much more robust (for comparison of both approaches, see Leys et al., 2013, Journal of Experimental Social Psychology).

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/detectOutliers.R") # load functions from my github site

example <- c(1, 3, 3, 6, 8, 10, 10, 1000) # 1000 is an outlier
outliersZ(example) # SD approach
outliersMAD(example) # MAD approach

# changing a few default values (note that if zCutOff is 3, 1000 ISN'T considered an outlier!!!)
outliersZ(x = example, zCutOff = 3, replaceOutliersWith = -999) # common to use 1.96, 2.5, 3 for Z cutoff
outliersMAD(x = example, MADCutOff = 3, replaceOutliersWith = -888) # Leys et al. (2003) recommends 2.5 to 3 for MAD cutoff
```

## Compute between- and within-subjects standard errors and confidence intervals

Code adapted from [Cookbook for R](http://www.cookbook-r.com/Graphs/)

When using the functions ```se``` (between-subjects) or ```seWithin``` (within-subjects), you can specify more than one outcome variable via the ```measurevar``` argument. If you specify more than one outcome variable, the output will be a list that has length of the number of outcome variables provided.

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/se.R") # load functions from my github site

# between-subjects
se(data = mtcars, measurevar = "disp", groupvars = "cyl")
se(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "am"))
se(data = ChickWeight, measurevar = "weight", groupvars = "Diet")

# within-subjects
seWithin(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
```
