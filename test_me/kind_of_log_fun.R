x<-1:1000
df <- data.frame(x)
library(ggplot2)

myfun <- function(x) 3+x^(1/3.1)

ggplot(df,aes(x))+
 geom_function(fun=myfun)+
  geom_hline(yintercept=12.5, linetype="dashed", color = "red")+
  annotate("text", x = 4, y = 13, label = "sill", color = "red")+
  #geom_hline(yintercept=4, linetype="dashed", color = "red")+
  ylab("gamma(h)")+
  xlab("h")+
  scale_x_continuous(breaks = seq(0,1000,100))+
  scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        #axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        #axis.ticks.y=element_blank(), #remove y axis ticks
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white")
        )

