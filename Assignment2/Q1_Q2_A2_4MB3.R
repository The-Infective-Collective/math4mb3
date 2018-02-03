## QUESTION 1
##Import data. Note: before running this, set working directory to location under 'session'
datafile <- "pim_us_phila_city_1918_dy.csv"
philadata <- read.csv(datafile)
philadata$date <- as.Date(philadata$date) #encodes date character strings as dates

## first make the box with no annotation or curves
plot(philadata$date, philadata$pim, type="n", # don't actually plot anything
     bty="L", # no upper or right box lines
     ylim=c(0,800), # axis limits
     yaxp=c(0,800,4), #first two numbers is coordinates of extreme tick marks, third number is the number of marks
     xaxt='n', #supress x ticks and labels
     xlab="",
     ann=FALSE, # no axis annotation (i.e., no title or axis labels)
     xaxs="i", #first tick on x axis is the y axis
     las=1 # axis label style: always horizontal
)
month <- c(9,10,11,12) #want sept, oct, nov, dec labels
ticks <- seq(philadata$date[1], philadata$date[length(philadata$date)], by="months") 
axis(1, at = ticks, labels = month.abb[month], tcl = -0.3) #put ticks where we want them and labels

## creat label for axes
mtext("Date", side=1, adj=1, line=1.5, font=1, cex=1.75, col="blue") #putting x label 'Date' on
mtext("P&I Deaths", side=2, at=900, line=-4, font=1, las=1, cex=1.75, col="blue") #putting y label on
#at =900 is because it is at around 900 on the y axis plot (just above the top which is 800)
#font=1 means normal font (not italics or bold)

## plot data
lines(philadata, col="grey", lwd=1.75) #putting the grey line, lwd gives line thickness relativ to default
points(philadata$date, philadata$pim, pch=21, bg="red") #putting the red points



## QUESTION 2
library(tidyverse)
library(ggplot2)

x = philadata %>%
  filter(date <= "1918-10-7" & date >= "1918-09-15") #filter data

ggplot(x, aes(x=date,y=pim))+ #plot log10 plot of filtered data
  geom_point()+
  scale_y_log10() +
  geom_smooth(method = "lm", se=F)

lm(formula ="date~log(pim)", data= x) #generates lm and gives slope in base e


       


