##Rscript for interpolating values from standard curves, specifically ELISAs

library(drc)
#Make dataframe with known concentrations and Absorbances
CONC <- c(0.0000640, 0.0003200, 0.0016, 0.008, 0.04, 0.2, 1)
BD <- c(0, 0.031063, 0.072063, 0.172563, 0.547063, 1.387063, 2.048063)
stdCurve <- data.frame(CONC, BD)


model1 <- drm(BD~CONC, 
              fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")),data=stdCurve)
plot(model1)

#Here are some observed abosorbances from our unknown samples
titer1in50 <- c(1.66, 1.15, 0.98, 1.68, 2.33, 2.54, 1.68, 0.25, 2.12, 1.12, 0.52, 0.48)   
response <- titer1in50

#The ED function is used to estimate concentration using our model
CONCx <- ED(model1, response)


# type="absolute" gives you the ability to use absolute values for the response, to 
# estimate the DOSE
response<-0.5   #lets use 0.5 for the response
CONCx<-ED(model1,response,type="absolute",display=F)[1:5] # the estimated DOSE
points(y=rep(response,5),x=CONCx,col="blue",pch=1:5)
