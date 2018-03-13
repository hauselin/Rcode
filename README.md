# R helper functions

Helper functions for R.

## Fitting statistical models in R

When fitting models, use summaryh() instead of summary() to get APA (American Psychological Association) formatted output. The first time you source() the the functions, it might take some time because it will install a few R packages (e.g., dplyr). Subsequently, source() should load the functions much faster.

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
