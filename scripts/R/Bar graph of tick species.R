### -- Bar chart -- ##
### load libraries ###
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(forcats)
rm(list=ls())
## Set working directory
setwd("~/Documents/Uni/Honours/Results/R data/R spreadsheets input")
## Load data
data <- read.csv("Summary_table1.csv")
## check data has loaded correctly
str(data)
## begin constructig chart
g <- ggplot(data = data, aes(Host))
g + geom_bar(aes(weight = Sum, fill = GenusSpecies))
## reverse orientation of chart
g1 <- g +
  geom_bar(aes(weight = Sum, fill = GenusSpecies), position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme(legend.position = "top")
g1
## order the x axis from largest to smallest
g2 <- g1 + aes(fct_infreq(factor(Host)))
g2
## changing axis
g3 <- g2 + theme(axis.text.x=element_text(size=12),
                 axis.title.x=element_text(size=14),
                 axis.title.y=element_text(size=14),
                 axis.text.y=element_text(size=12),
                 strip.text.x=element_text(size=10,face="italic"))
g3
g4 <- g3 + xlab ("Bandicoot host") + ylab ("No. of ticks")
g4
## changing legend
g5 <- g4 + theme(legend.text=element_text(size=12),
                  legend.title=element_text(size=14))
g5
g6 <- g5 + guides(fill=guide_legend(title="Tick species"))
g6
g7 <- g6 + theme(legend.text = element_text(face = "italic"))
g7
###this is not working - trying to alter order of bars and add black line around the different colours
# g7 <- g6 + ggplot(mydata, aes(x=reorder(Host, - TickSum), y=TickSum)) + geom_bar()
# g7
# g8 <- g7 + ggplot(fill = "identity", colour="black") + geom_bar()
# g8
g8 <- g7 + ylim(0,145)
g8
print(g7)
ggsave("g8.png", width = 12, height = 7)
ggsave("g7.png", plot=g7, width=15, height=5, units="in")
pdf("g7.pdf")
print(g7.pdf(...))
g7 + coord_cartesian(
  xlim = c(0, 100), ylim = c(10, 20))
