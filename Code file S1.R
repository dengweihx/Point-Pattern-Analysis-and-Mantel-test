# Code file S1. R code for analysis


######Point Pattern Analysis of Arthrobotrys oligospora population

###Setting up the work path
my.path= "E:/data/"

my.fn="site.csv"

my.fn=paste0(my.path, my.fn)

#### Read data
mydata=read.csv(my.fn)#Contains latitude and longitude information

####Loading R packages
library(SoDA)
library(spatstat)

###
car <- geoXY(mydata$lat, mydata$long, unit = 1000)
write.csv(car, "car1.csv")
sppdata = read.csv("car1.csv", header= T)

spp=ppp(sppdata$X, sppdata$Y,
        window=owin(xrange = c(0,800),yrange=c(0,800)))
cycl.envelope<-envelope(spp,fun=Kest,nsim=199)##K-function
write.csv(cycl.envelope, "cyclenvelope.csv")

####chi-square test
chisq.test(cycl.envelope$obs,cycl.envelope$theo)


##### Mantel test 
install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
packageVersion("linkET")
install.packages("FD")

library(dplyr)
library(vegan)
library(linkET)
library(ggplot2)
library(FD)

spe <- read.csv("spe.csv",row.names = 1)
env <- read.csv("mat.csv",row.names = 1)

mantel <- mantel_test(spe,env,
                      spec_select = list(pp = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#> `mantel_test()` using 'bray' dist method for 'spec'.
#> `mantel_test()` using 'euclidean' dist method for 'env'.

qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

