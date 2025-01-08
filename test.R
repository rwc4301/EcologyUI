library(ggplot2)

p <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point() + theme_minimal()