

#install.packages("")
library(ODEsensitivity)
library(ODEnetwork)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(readr)
library(splines)
library(reshape2)


theq<-0.341196314
i1<- 3.45125E-05



mymodel1 <- function(t, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dth1<- th2
    dth2<- (-K*(th1-theq)-B*th2)/I
    
    
    
    return(list(c(dth1,dth2)))
  })
}



# The parameters to be included in the sensitivity analysis and their lower
# and upper boundaries:
modelpars  <- c("K", "B","I")
modellwM <- c(0.1,0.1,1.7E-7)
modelupM <- c(10,2,i1)
# The initial values of the state variables:
modelinit  <- c(th1=0.753809808,
                th2=0)
# The timepoints of interest:
modeltimes <- seq(0.008,0.3933,by=0.01)


set.seed(3)
# Sobol's sensitivity analysis (here only with n = 500, but n = 1000 is
# recommended):

# Warning: The following code might take very long!
mymodel3_sobol <- ODEsobol(mod = mymodel1,
                           pars = modelpars,
                           state_init = modelinit,
                           t = modeltimes,
                           n = 5000,
                           rfuncs = "runif",
                           rargs = paste0("min = ", modellwM,
                                          ", max = ", modelupM),
                           sobol_method = "Martinez",
                           ode_method = "lsoda",atol = 1e-5,rtol = 1e-5)


Th1<-mymodel3_sobol$th1$T
Th1<-as.data.frame(t(Th1))


colours<- c("K"="yellow","B"="red","I"="blue")
#Plot
plot1<-ggplot(data = Th1)+geom_line(aes(x=time,y=K,colour="K"),size=2)+
  geom_line(aes(x=time,y=B,colour="B"),size=2)+geom_line(aes(x=time,y=I,colour="I"),size=2)+ggtitle("Sensitivity Analysis")+labs(x="Time (s)", y="Total Order Sobol Indices", colour="Legend")+
  theme_pubr()+theme(plot.title = element_text(hjust = 0.5))

ggsave("Sensitivity Chapter7.png",dpi=600,path="C:/Users/u1857308/OneDrive - University of Warwick/PhD/Sensitivity Analysis/Linda R/")

th1_f<-as.data.frame(t(mymodel3_sobol$th1$S))
F_Sobol<-colMeans(th1_f)[-1]
tTh1<-colMeans(Th1)[-1]
if (tTh1[3]<0.01)
  {
  tTh1[3]<-0.005
}

if (F_Sobol[3]<0.01)
{
  F_Sobol[3]<-0.001
}
Parameters<-colnames(Th1)[-1]
soboldata <- data.frame(Parameters, tTh1, F_Sobol)
colnames(soboldata)[2] <- "Total Order Sobol Indices"
colnames(soboldata)[3] <- "First Order Sobol Indices"

# Melt data for ggplot
soboldata_long <- melt(soboldata, id.vars = "Parameters", variable.name = "Type", value.name = "Value")

# Plotting
plot2 <- ggplot(data = soboldata_long, aes(x = Parameters, y = Value, fill = Type)) +
  geom_col(position = position_dodge(width = 0.9), size = 1) +
  ylab(expression("Sobol Indices")) +
  xlab("Parameters") +
  theme_pubr() +
  ggtitle("Sobol Indices") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot
ggsave("Sensitivity_Chapter7_Bars.png", plot = plot2, dpi = 600, path = "C:/Users/u1857308/OneDrive - University of Warwick/PhD/Sensitivity Analysis/Linda R/")