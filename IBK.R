#install.packages("")
library(ODEsensitivity)
library(ODEnetwork)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(readr)
library(splines)



data <-as.matrix(read_csv("C:/Users/.......Read the Data.csv"))
yedc<-data[,1]
yfds<-data[,2]
yfdp<-data[,3]
tim<-data[,9]
thf<-data[,10]
thP<-data[,11]
thD<-data[,12]
edcMCP<-splinefun(tim,yedc)
fdpMCP<-splinefun(tim,yfdp)
fdsMCP<-splinefun(tim,yfds)
edcPIP<-splinefun(tim,data[,4])
fdpPIP<-splinefun(tim,data[,5])
fdsPIP<-splinefun(tim,data[,6])
edcDIP<-splinefun(tim,data[,7])
fdpDIP<-splinefun(tim,data[,8])
theq<-thf[length(thf)]

il<-1.2E-5


mymodel1 <- function(t, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dth1<- th2
    dth2<- (c_EDC*edcMCP(t)+c_FDS*fdpMCP(t)+c_FDP*fdsMCP(t)-K*(th1-theq)+B*th2)/I

    
    
    return(list(c(dth1,dth2)))
  })
}

mymodel2 <- function(t, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dth1<- th2
    dth2<- (c_EDC*edcPIP(t)+c_FDS*fdpPIP(t)+c_FDP*fdsPIP(t)-K*(th1-theq)+B*th2)/I
    
    
    
    return(list(c(dth1,dth2)))
  })
}

mymodel3 <- function(t, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dth1<- th2
    dth2<- (c_EDC*edcDIP(t)+c_FDS*fdpDIP(t)-K*(th1-theq)-B*th2)/I
    
    
    
    return(list(c(dth1,dth2)))
  })
}


# The parameters to be included in the sensitivity analysis and their lower
# and upper boundaries:
modelpars  <- c("K", "B","I","c_EDC","c_FDP","c_FDS")
modellwM <- c(0.1,0.1,il-0.05*i1,0.1,0.1,0.1 )
modelupM <- c(10,2,il+0.05*il,2,2,2)

# The initial values of the state variables:
modelinit  <- c(th1=thf[1],
                th2=0)

# The timepoints of interest:
modeltimes <- seq(tim[1],tim[30],by=0.1)


set.seed(3)
# Sobol's sensitivity analysis (here only with n = 500, but n = 1000 is
# recommended):

# Warning: The following code might take very long!
mymodel3_sobol <- ODEsobol(mod = mymodel1,
                           pars = modelpars,
                           state_init = modelinit,
                           t = modeltimes,
                           n = 8000,
                           rfuncs = "runif",
                           rargs = paste0("min = ", modellwM,
                                          ", max = ", modelupM),
                           sobol_method = "Martinez",
                           ode_method = "rk4",atol = 1e-5,rtol = 1e-5)


Th1<-mymodel3_sobol$th1$T


Th1<-as.data.frame(t(Th1))


colours<- c("K"="yellow","B"="red","I"="blue","c_EDC"="brown","c_FDS"="black","c_FDP"="purple")
#Plot
plot1<-ggplot(data = Th1)+geom_line(aes(x=time,y=K,colour="K"),size=2)+
  geom_line(aes(x=time,y=B,colour="B"),size=2)+geom_line(aes(x=time,y=I,colour="I"),size=2)+geom_line(aes(x=time,y=c_EDC,colour="c_EDC"),size=2)+
  geom_line(aes(x=time,y=c_FDP,colour="c_FDP"),size=2)+geom_line(aes(x=time,y=c_FDS,colour="c_FDS"),size=2)+ggtitle("Sensitivity Analysis")+labs(x="Time (s)", y="Total Order Sobol Indices", colour="Legend")+
  theme_pubr()+theme(plot.title = element_text(hjust = 0.5))

ggsave("Sensitivity.png",dpi=600)


tTh1<-colMeans(Th1)[-1]
if (tTh1[5]<1E-1)
{
  tTh1[5]=3E-2
}
if (tTh1[6]<1E-1)
{
  tTh1[6]=3E-2
}
if (tTh1[4]<1E-1)
{
  tTh1[4]=3E-2
}

th1_f<-as.data.frame(t(mymodel3_sobol$th1$S))
F_Sobol<-colMeans(th1_f)[-1]

if (F_Sobol[5]<1E-1)
{
  F_Sobol[5]=1E-2
}
if (F_Sobol[6]<1E-1)
{
  F_Sobol[6]=1E-2
}
if (F_Sobol[4]<1E-1)
{
  F_Sobol[4]=1E-2
}

Parameters<-colnames(Th1)[-1]
soboldata <- data.frame(Parameters, tTh1, F_Sobol)

# Rename the Type labels
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
ggsave("Sensitivity Bars.png",dpi=600)

