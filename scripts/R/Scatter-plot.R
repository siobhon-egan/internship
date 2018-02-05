### -- Scatter plot -- ##
# Remove current environment variables
rm(list=ls())
### load libraries ###
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(forcats)
## Set working directory
setwd("~/Documents/Uni/Cryptick lab")
# You can check your working directory
getwd()
## Load data, you can do this using the GUI features in R studio and easily change factors such as 'first row equals variable names"
sample_data <- read_csv("sample_data.csv")
## check data has loaded correctly
View(sample_data)
str(data)
## begin constructig chart
p <- ggplot(sample_data, aes(measurement_one, measurement_two))
p1 <- p + geom_point()
p1
# Code sample points by a factor using colour
p2 <- p1 + geom_point(aes(colour = factor(species)))
p2
# Want to code sample points by a factor using shape? Do this
#> p1 + geom_point(aes(shape = factor(species)))
# You can even use both colour on the same factor, or different factors
#> p1 + geom_point(aes(shape = factor(species), colour = factor(species)))
# changing axis
p3 <- p2 + theme(axis.text.x=element_text(size=12),
                 axis.title.x=element_text(size=14),
                 axis.title.y=element_text(size=14),
                 axis.text.y=element_text(size=12),
                 strip.text.x=element_text(size=10,face="italic"))
p3
p4 <- p3 + xlab ("Fluffiness (kg)") + ylab ("Cuteness (cm)")
p4
# need to add scientific notaion to your axis
# p4 + xlab(expression(paste("Fluffiness, ", mu, sigma)))
# see here for more details on symbols https://stats.idre.ucla.edu/r/codefragments/greek_letters/
## changing legend
p5 <- p4 + theme(legend.text=element_text(size=12),
                  legend.title=element_text(size=14))
p5
# make legend appear on top instead
#> p5 + theme(legend.position="top")
# Change the axis limits
p6 <- p5 + ylim(0,145) + xlim(0,145)
p6
# Want to emove the plot legend?
#> p6 + theme(legend.position='none')
# save plot
ggsave("p6.png")