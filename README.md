# R helper functions

Helper functions for R.

## Fitting statistical models in R

When fitting models, using summaryh() instead of summary() to get APA (American Psychological Association) formatted output.

```
source("https://raw.githubusercontent.com/hauselin/Rcode/master/summaryh.R") # load functions from my github site
summary(lm(mpg ~ cyl, mtcars)) # base R summary()
summaryh(lm(mpg ~ cyl, mtcars)) # returns APA-formatted output
```
