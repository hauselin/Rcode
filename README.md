# R helper functions

Helper functions for R.

## Fitting statistical models in R

When fitting models, using summaryh() instead of summary() to get APA (American Psychological Association) formatted output. The first time you source() the the functions, it might take some time because it will install a few R packages (e.g., dplyr). Subsequently, source() should load the functions much faster.

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R") # load functions from my github site
summary(lm(mpg ~ cyl, mtcars)) # base R summary()
summaryh(lm(mpg ~ cyl, mtcars)) # returns APA-formatted output
```
