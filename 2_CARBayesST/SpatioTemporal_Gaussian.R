########################################################################
########################################################################
######################### BAY PROJECT ##################################
######################    CARBayesST     ###############################
########################################################################
#### Spatio-Temporal model for Highschool age
########################################################################


# install.packages("USAboundaries")
# install.packages("CARBayesST")
library(CARBayesST)
library(spdep)
library(maptools) 
library(sp)
library(CARBayesdata)
library(USAboundaries) 
library(sf)
library(dplyr)
library(ggplot2)
library(MLmetrics)
library(LaplacesDemon)
library(rgeos)



# setwd("G:/Bayesian_project_final/2_CARBayesST")
library(readxl)


mydata <- read_excel("data.xlsx")
BIRTH_data <- data.frame(mydata)

# removing isolated States
birthdata_hs <- subset(BIRTH_data, BIRTH_data$Age.Group..Years. == "15-17 years")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Total U.S.")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Hawaii")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Alaska")
birthdata_hs <- birthdata_hs[,c(2,1,5,4)]

### Construct the SpatialPolygonsDataFrame
require(maps)
require(sp)
require(maptools)
require(tmap)
require(cartogram)

# Lower case letters for all states, in order to merge them with the map
mydata <- birthdata_hs
mydata <- data.frame(lapply(mydata, function(v) {
  if (is.character(v)) return(tolower(v))
  else return(v)
}))


# Get a SpatialPolygonsDataFrame of US states
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
usa <- SpatialPolygonsDataFrame(usa, 
                                data = data.frame(unique(IDs), 
                                                  row.names = unique(IDs)) )

usa@data = data.frame(usa@data, mydata[match(usa@data[,'unique.IDs.'], mydata[,'State' ]),])
usa.data.combined = usa 

# remove NAs from the dataset
usa.data.combined@data <- usa.data.combined@data[complete.cases(usa.data.combined@data), ]
spplot(usa.data.combined[, "State.Rate"], main = "Birth rate" )

### Constuct binary neighbourhood matrix W
W.nb <- poly2nb(usa.data.combined, row.names = unique(IDs)) # neighbourhood
W <- nb2mat(W.nb, style = "B")
dim(W) # matrix of dimension 49 x 49
W[7:15,7:15]
W_nodoc <- W[c(-8),c(-8)]
W_nodoc[7:15,7:15]

K<-49 # number of states
N<-29 # number of years

### Third order differences - Stationnarisation and fixing the state order which was not always the same 
dif3rate <- rep(0,K*(N-3)) # vector of third differences of State.Rate
shifrate <- rep(0,K*(N-3)) # State.Rate shifted by 3 years
fix_rate <- rep(0,K*N)
states <- mydata$State[1:49]

for (i in c(1:49))
{
  state <- states[i]
  state_subset <- subset(mydata, mydata$State == state)
  rate <- state_subset$State.Rate
  dif3 <- diff(rate,diff=3)
  shif <- rate[c(4:29)]
  for (j in c(0:25))
  {
    dif3rate[i + K*j] <- dif3[j+1]
    shifrate[i + K*j] <- shif[j+1]
  }
  for (j in c(0:28))
  {
    fix_rate[i + K*j] <- rate[j+1]
  }
}

years <- c(rep(1990,49),rep(1991,49),rep(1992,49),rep(1993,49),rep(1994,49),rep(1995,49)
           ,rep(1996,49),rep(1997,49),rep(1998,49),rep(1999,49),rep(2000,49),rep(2001,49)
           ,rep(2002,49),rep(2003,49),rep(2004,49),rep(2005,49),rep(2006,49),rep(2007,49)
           ,rep(2008,49),rep(2009,49),rep(2010,49),rep(2011,49),rep(2012,49),rep(2013,49)
           ,rep(2014,49),rep(2015,49),rep(2016,49),rep(2017,49),rep(2018,49))
shif_years <- years[c((49*3+1):(49*N))] # years vector, without the first three years
mydata_staz <- data.frame(rep(states,26),shif_years,dif3rate,shifrate) # new dataframe 
mydata_fix <- data.frame(rep(states,29),years,fix_rate)

### Moran I for each year for raw Rates.
moranI_years <- rep(0,29)
for (i in c(1990:2018))
{
  moran_subset <- subset(mydata_fix, years == i)
  mt <- moran.test(moran_subset$fix_rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1989] <- mt$estimate[1]
}
plot(c(1990:2018),moranI_years,main = 'Morans I over Years for raw data',
     xlab='Years', ylab='Morans I',ylim = c(-0.1,0.6), pch = 19, col='blue')
abline(h=0, pch=10, col="black", lwd=2, lty=2)
# MoransI is significant for every year


### Moran I for each year for 3rd differences
moranI_years <- rep(0,26)
for (i in c (1993:2018))
{
  moran_subset <- subset(mydata_staz, mydata$Year == i)
  mt <- moran.test(moran_subset$dif3rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1992] <- mt$estimate[1]
}
plot(c(1993:2018),moranI_years, main = 'Morans I over Years for stationarize data',
     xlab='Years', ylab='Morans I',ylim = c(-0.1,0.6), pch = 19, col='red')
abline(h=0, pch=10, col="black", lwd=2, lty=2)
# spatial correlation is lost by taking 3rd order diff
# so it does't really make sense to fit 


### How "far" should we reach for  neighbours

centroids <- gCentroid(usa,byid=T)
nb_knn7 <- knn2nb((knearneigh(centroids,k=7)))
sp_corG <- sp.correlogram(nb_knn7,usa.data.combined@data$State.Rate, order = 7,
                          zero.policy = TRUE, style = 'B', method = 'I')
plot(sp_corG, main = 'Spatial correlogram for Birth rates in States', ylab="Morans'I")
# We can limit ourselves to 1st order neighbour, as they present the highest spatial correlation


### Spatial modelling with CARBayesST : ST.CARar on Highschool raw data
set.seed(1200)
formula1 <- fix_rate ~ fix_rate 
chain1 <- ST.CARar(formula = formula1, family = "gaussian",
                   data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                   thin = 100)

# Posterion inference 
print(chain1) 
x11()
plot(chain1$samples$nu2)
plot(chain1$samples$tau2)
plot(chain1$samples$rho)

### Fitting ST.CARar, but leaving out some data this time

## TEXAS is completely NA
mydata_fixNA<-mydata_fix
# texas in number 42, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[42+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_tex <- ST.CARar(formula = formula1NA, family = "gaussian",
                         data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                         thin = 100)

low <- high <- texasNA <- rep(0,29)
for (i in 0:28)
{
  texasNA[i+1] <-fitted.values(chain1NA_tex)[42+K*i]
  low[i+1] <- p.interval(chain1NA_tex$samples$fitted[,42+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_tex$samples$fitted[,42+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2] 
}
texas_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "texas")
texas<-data.frame(texas_measures,texasNA,high,low)


x11()
plot(c(c(1990:2018),c(1990:2018),c(1990:2018),c(1990:2018)),
     c(texas$fix_rate,texas$texasNA,texas$low,texas$high), 
     pch = c(rep(19,29*2),rep(3,29*2)),
     col=c(rep("blue",29),rep("red",29*3)), main = 'Real vs sampled Texas Birth rate values',
     xlab = "Year", ylab="Birth rate")

ggplot(data.frame(x = 1990:2018, y = texas$texasNA, ylow = texas$low, yup = texas$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = texas$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Texas Birth rate values")



## MISSOURI is completely NA
mydata_fixNA<-mydata_fix
# missouri in number 24, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[24+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_miss <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- missouriNA <- rep(0,29)
for (i in 0:28)
{
  missouriNA[i+1] <-fitted.values(chain1NA_miss)[24+K*i]
  low[i+1] <- p.interval(chain1NA_miss$samples$fitted[,24+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_miss$samples$fitted[,24+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri<-data.frame(missouri_measures,low,missouriNA,high)

ggplot(data.frame(x = 1990:2018, y = missouri$missouriNA, ylow = missouri$low, yup = missouri$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = missouri$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Missouri Birth rate values")



## WASHINGTON is completely NA
mydata_fixNA<-mydata_fix
# washington in number 46, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[46+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_wash <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- washingtonNA <- rep(0,29)
for (i in 0:28)
{
  washingtonNA[i+1] <-fitted.values(chain1NA_wash)[46+K*i]
  low[i+1] <- p.interval(chain1NA_wash$samples$fitted[,46+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_wash$samples$fitted[,46+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
washington_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "washington")
washington<-data.frame(washington_measures,low,washingtonNA,high)

ggplot(data.frame(x = 1990:2018, y = washington$washingtonNA, ylow = washington$low, yup = washington$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = washington$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Washington Birth rate values")




## IDAHO is completely NA
mydata_fixNA<-mydata_fix
# idaho in number 11, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[11+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_idaho <- ST.CARar(formula = formula1NA, family = "gaussian",
                           data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                           thin = 100)

low <- high <- idahoNA <- rep(0,29)
for (i in 0:28)
{
  idahoNA[i+1] <-fitted.values(chain1NA_idaho)[11+K*i]
  low[i+1] <- p.interval(chain1NA_idaho$samples$fitted[,11+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_idaho$samples$fitted[,11+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
idaho_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "idaho")
idaho<-data.frame(idaho_measures,low,idahoNA,high)

ggplot(data.frame(x = 1990:2018, y = idaho$idahoNA, ylow = idaho$low, yup = idaho$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = idaho$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Idaho Birth rate values")


## montana is completely NA
mydata_fixNA<-mydata_fix
# montana in number 25, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[25+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_montana <- ST.CARar(formula = formula1NA, family = "gaussian",
                             data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                             thin = 100)

low <- high <- montanaNA <- rep(0,29)
for (i in 0:28)
{
  montanaNA[i+1] <-fitted.values(chain1NA_montana)[25+K*i]
  low[i+1] <- p.interval(chain1NA_montana$samples$fitted[,25+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_montana$samples$fitted[,25+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
montana_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "montana")
montana<-data.frame(montana_measures,low,montanaNA,high)

ggplot(data.frame(x = 1990:2018, y = montana$montanaNA, ylow = montana$low, yup = montana$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = montana$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Montana Birth rate values")



## ohio is completely NA
mydata_fixNA<-mydata_fix
# ohio in number 34, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[34+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_ohio <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- ohioNA <- rep(0,29)
for (i in 0:28)
{
  ohioNA[i+1] <-fitted.values(chain1NA_ohio)[34+K*i]
  low[i+1] <- p.interval(chain1NA_ohio$samples$fitted[,34+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_ohio$samples$fitted[,34+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
ohio_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "ohio")
ohio<-data.frame(ohio_measures,low,ohioNA,high)

ggplot(data.frame(x = 1990:2018, y = ohio$ohioNA, ylow = ohio$low, yup = ohio$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = ohio$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Ohio Birth rate values")



## Putting last 3 years NA in Missouri and plot the posterior inferences
mydata_fixNA3<-mydata_fix 
# missouri in number 24, putting all years to NA
for (i in 26:28)
{
  mydata_fixNA3$fix_rate[24+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_miss3 <- ST.CARar(formula = formula1NA, family = "gaussian",
                           data = mydata_fixNA3, W = W, burnin = 20000, n.sample = 220000,
                           thin = 100)

low <- high <- missouriNA3 <- rep(0,3)
for (i in 0:2)
{
  missouriNA3[i+1] <-fitted.values(chain1NA_miss3)[24+K*(i+26)]
  low[i+1] <- p.interval(chain1NA_miss3$samples$fitted[,24+K*(i+26)], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_miss3$samples$fitted[,24+K*(i+26)], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri3<-data.frame(missouriNA3,low,high)
x11()
plot(c(c(1990:2018),c(2016:2018),c(2016:2018),c(2016:2018)),
     c(missouri_measures$fix_rate,missouri3$missouriNA,missouri3$low,missouri3$high), 
     pch = c(rep(19,29+3),rep(3,3*2)),
     col=c(rep("blue",29),rep("red",3*3)), main = 'observed vs sampled Missouri Birth rate values',
     xlab = "Year", ylab="Birth rate")
abline(v=2015, pch=10, col="black", lwd=1, lty=2)


## Putting last 3 years NA in all States
mydata_fixNA3<-mydata_fix
for (i in (49*26+1):(49*29))
{
  mydata_fixNA3$fix_rate[i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_3 <- ST.CARar(formula = formula1NA, family = "gaussian",
                       data = mydata_fixNA3, W = W, burnin = 20000, n.sample = 220000,
                       thin = 100)

low <- high <- missouriNA3 <- rep(0,3)
for (i in 0:2)
{
  missouriNA3[i+1] <-fitted.values(chain1NA_3)[24+K*(i+26)]
  low[i+1] <- p.interval(chain1NA_3$samples$fitted[,24+K*(i+26)], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_3$samples$fitted[,24+K*(i+26)], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri3<-data.frame(missouriNA3,low,high)
x11()
ggplot(data.frame(x = 2015:2018, 
                  y = c(missouri_measures$fix_rate[26] ,missouri3$missouriNA3), 
                  ylow = c(missouri_measures$fix_rate[26], missouri3$low), 
                  yup = c(missouri_measures$fix_rate[26], missouri3$high)  )) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = missouri_measures$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  geom_vline(aes(xintercept = 2015), col = 2, lty = 2) +
  ggtitle("ST.CARar model with only Missouri")


#--------------------------------------------------------------------------




# Let's do the same analysis but on Post Highschool group
mydata <- read_excel("data.xlsx")
BIRTH_data <- data.frame(mydata)

# removing isolated States
birthdata_hs <- subset(BIRTH_data, BIRTH_data$Age.Group..Years. == "18-19 years")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Total U.S.")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Hawaii")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Alaska")
birthdata_hs <- birthdata_hs[,c(2,1,5,4)]


# Lower case letters for all states, in order to merge them with the map
mydata <- birthdata_hs
mydata <- data.frame(lapply(mydata, function(v) {
  if (is.character(v)) return(tolower(v))
  else return(v)
}))


# Get a SpatialPolygonsDataFrame of US states
usa <- map("state", fill = TRUE)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
usa <- SpatialPolygonsDataFrame(usa, 
                                data = data.frame(unique(IDs), 
                                                  row.names = unique(IDs)) )

usa@data = data.frame(usa@data, mydata[match(usa@data[,'unique.IDs.'], mydata[,'State' ]),])
usa.data.combined = usa 

# remove NAs from the dataset
usa.data.combined@data <- usa.data.combined@data[complete.cases(usa.data.combined@data), ]
spplot(usa.data.combined[, "State.Rate"], main = "Birth rate" )

### Constuct binary neighbourhood matrix W
W.nb <- poly2nb(usa.data.combined, row.names = unique(IDs)) # neighbourhood
W <- nb2mat(W.nb, style = "B")
dim(W) # matrix of dimension 49 x 49
W[7:15,7:15]
W_nodoc <- W[c(-8),c(-8)]
W_nodoc[7:15,7:15]

K<-49 # number of states
N<-29 # number of years

### Third order differences - Stationnarisation and fixing the state order which was not always the same 
dif3rate <- rep(0,K*(N-3)) # vector of third differences of State.Rate
shifrate <- rep(0,K*(N-3)) # State.Rate shifted by 3 years
fix_rate <- rep(0,K*N)
states <- mydata$State[1:49]

for (i in c(1:49))
{
  state <- states[i]
  state_subset <- subset(mydata, mydata$State == state)
  rate <- state_subset$State.Rate
  dif3 <- diff(rate,diff=3)
  shif <- rate[c(4:29)]
  for (j in c(0:25))
  {
    dif3rate[i + K*j] <- dif3[j+1]
    shifrate[i + K*j] <- shif[j+1]
  }
  for (j in c(0:28))
  {
    fix_rate[i + K*j] <- rate[j+1]
  }
}
years <- c(rep(1990,49),rep(1991,49),rep(1992,49),rep(1993,49),rep(1994,49),rep(1995,49)
           ,rep(1996,49),rep(1997,49),rep(1998,49),rep(1999,49),rep(2000,49),rep(2001,49)
           ,rep(2002,49),rep(2003,49),rep(2004,49),rep(2005,49),rep(2006,49),rep(2007,49)
           ,rep(2008,49),rep(2009,49),rep(2010,49),rep(2011,49),rep(2012,49),rep(2013,49)
           ,rep(2014,49),rep(2015,49),rep(2016,49),rep(2017,49),rep(2018,49))
shif_years <- years[c((49*3+1):(49*N))] # years vector, without the first three years
mydata_staz <- data.frame(rep(states,26),shif_years,dif3rate,shifrate) # new dataframe 
mydata_fix <- data.frame(rep(states,29),years,fix_rate)

### Moran I for each year for raw Rates.
moranI_years <- rep(0,29)
for (i in c(1990:2018))
{
  moran_subset <- subset(mydata_fix, years == i)
  mt <- moran.test(moran_subset$fix_rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1989] <- mt$estimate[1]
}
plot(c(1990:2018),moranI_years,main = 'Morans I over Years for raw data',
     xlab='Years', ylab='Morans I',ylim = c(-0.1,0.8), pch = 19, col='blue')
abline(h=0, pch=10, col="black", lwd=2, lty=2)
# MoransI is significant for every year


### Moran I for each year for 3rd differences
moranI_years <- rep(0,26)
for (i in c (1993:2018))
{
  moran_subset <- subset(mydata_staz, mydata$Year == i)
  mt <- moran.test(moran_subset$dif3rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1992] <- mt$estimate[1]
}
plot(c(1993:2018),moranI_years, main = 'Morans I over Years for stationarized data',
     xlab='Years', ylab='Morans I',ylim = c(-0.3,0.5), pch = 19, col='red')
abline(h=0, pch=10, col="black", lwd=2, lty=2)
# spatial correlation is lost by taking 3rd order diff
# so it does't really make sense to fit 


### How "far" should we reach for  neighbours
# centroids <- gCentroid(usa,byid=T)
# nb_knn7 <- knn2nb((knearneigh(centroids,k=7)))
# sp_corG <- sp.correlogram(nb_knn7,usa.data.combined@data$State.Rate, order = 7,
#                           zero.policy = TRUE, style = 'B', method = 'I')
# plot(sp_corG, main = 'Spatial correlogram for Birth rates in States', ylab="Morans'I")
# We can limit ourselves to 1st order neighbour, as they present the highest spatial correlation



### Spatial modelling with CARBayesST : ST.CARar on Post Higschool raw data
set.seed(1200)
formula1 <- fix_rate ~ fix_rate 
chain1 <- ST.CARar(formula = formula1, family = "gaussian",
                   data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                   thin = 100)

# Posterion inference 
print(chain1) 
x11()
plot(chain1$samples$nu2)
x11()
plot(chain1$samples$tau2)
x11()
plot(chain1$samples$rho)



### Fitting ST.CARar, but leaving out some data this time

## TEXAS is completely NA
mydata_fixNA<-mydata_fix
# texas in number 42, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[42+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_tex <- ST.CARar(formula = formula1NA, family = "gaussian",
                         data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                         thin = 100)

low <- high <- texasNA <- rep(0,29)
for (i in 0:28)
{
  texasNA[i+1] <-fitted.values(chain1NA_tex)[42+K*i]
  low[i+1] <- p.interval(chain1NA_tex$samples$fitted[,42+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_tex$samples$fitted[,42+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2] 
}
texas_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "texas")
texas<-data.frame(texas_measures,texasNA,high,low)


x11()
ggplot(data.frame(x = 1990:2018, y = texas$texasNA, ylow = texas$low, yup = texas$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = texas$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Texas Birth rate values")



## MISSOURI is completely NA
mydata_fixNA<-mydata_fix
# missouri in number 24, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[24+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_miss <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- missouriNA <- rep(0,29)
for (i in 0:28)
{
  missouriNA[i+1] <-fitted.values(chain1NA_miss)[24+K*i]
  low[i+1] <- p.interval(chain1NA_miss$samples$fitted[,24+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_miss$samples$fitted[,24+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri<-data.frame(missouri_measures,low,missouriNA,high)

ggplot(data.frame(x = 1990:2018, y = missouri$missouriNA, ylow = missouri$low, yup = missouri$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = missouri$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Missouri Birth rate values")



## WASHINGTON is completely NA
mydata_fixNA<-mydata_fix
# washington in number 46, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[46+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_wash <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- washingtonNA <- rep(0,29)
for (i in 0:28)
{
  washingtonNA[i+1] <-fitted.values(chain1NA_wash)[46+K*i]
  low[i+1] <- p.interval(chain1NA_wash$samples$fitted[,46+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_wash$samples$fitted[,46+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
washington_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "washington")
washington<-data.frame(washington_measures,low,washingtonNA,high)

ggplot(data.frame(x = 1990:2018, y = washington$washingtonNA, ylow = washington$low, yup = washington$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = washington$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Washington Birth rate values")




## IDAHO is completely NA
mydata_fixNA<-mydata_fix
# idaho in number 11, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[11+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_idaho <- ST.CARar(formula = formula1NA, family = "gaussian",
                           data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                           thin = 100)

low <- high <- idahoNA <- rep(0,29)
for (i in 0:28)
{
  idahoNA[i+1] <-fitted.values(chain1NA_idaho)[11+K*i]
  low[i+1] <- p.interval(chain1NA_idaho$samples$fitted[,11+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_idaho$samples$fitted[,11+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
idaho_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "idaho")
idaho<-data.frame(idaho_measures,low,idahoNA,high)

ggplot(data.frame(x = 1990:2018, y = idaho$idahoNA, ylow = idaho$low, yup = idaho$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = idaho$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Idaho Birth rate values")


## montana is completely NA
mydata_fixNA<-mydata_fix
# montana in number 25, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[25+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_montana <- ST.CARar(formula = formula1NA, family = "gaussian",
                             data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                             thin = 100)

low <- high <- montanaNA <- rep(0,29)
for (i in 0:28)
{
  montanaNA[i+1] <-fitted.values(chain1NA_montana)[25+K*i]
  low[i+1] <- p.interval(chain1NA_montana$samples$fitted[,25+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_montana$samples$fitted[,25+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
montana_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "montana")
montana<-data.frame(montana_measures,low,montanaNA,high)

ggplot(data.frame(x = 1990:2018, y = montana$montanaNA, ylow = montana$low, yup = montana$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = montana$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Montana Birth rate values")



## ohio is completely NA
mydata_fixNA<-mydata_fix
# ohio in number 34, putting all years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[34+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_ohio <- ST.CARar(formula = formula1NA, family = "gaussian",
                          data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                          thin = 100)

low <- high <- ohioNA <- rep(0,29)
for (i in 0:28)
{
  ohioNA[i+1] <-fitted.values(chain1NA_ohio)[34+K*i]
  low[i+1] <- p.interval(chain1NA_ohio$samples$fitted[,34+K*i], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_ohio$samples$fitted[,34+K*i], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
ohio_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "ohio")
ohio<-data.frame(ohio_measures,low,ohioNA,high)

ggplot(data.frame(x = 1990:2018, y = ohio$ohioNA, ylow = ohio$low, yup = ohio$high) ) + 
  geom_line(aes(x = x, y = y), lty = 2, col = "red") + theme_bw() +
  geom_line(data = data.frame(Y = ohio$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") + 
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  ggtitle("Observed vs sampled Ohio Birth rate values")



## Putting last 3 years NA in Missouri and plot the posterior inferences
mydata_fixNA3<-mydata_fix 
# missouri in number 24, putting all years to NA
for (i in 26:28)
{
  mydata_fixNA3$fix_rate[24+K*i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_miss3 <- ST.CARar(formula = formula1NA, family = "gaussian",
                           data = mydata_fixNA3, W = W, burnin = 20000, n.sample = 220000,
                           thin = 100)

low <- high <- missouriNA3 <- rep(0,3)
for (i in 0:2)
{
  missouriNA3[i+1] <-fitted.values(chain1NA_miss3)[24+K*(i+26)]
  low[i+1] <- p.interval(chain1NA_miss3$samples$fitted[,24+K*(i+26)], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_miss3$samples$fitted[,24+K*(i+26)], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri3<-data.frame(missouriNA3,low,high)
x11()
plot(c(c(1990:2018),c(2016:2018),c(2016:2018),c(2016:2018)),
     c(missouri_measures$fix_rate,missouri3$missouriNA,missouri3$low,missouri3$high), 
     pch = c(rep(19,29+3),rep(3,3*2)),
     col=c(rep("blue",29),rep("red",3*3)), main = 'observed vs sampled Missouri Birth rate values',
     xlab = "Year", ylab="Birth rate")
abline(v=2015, pch=10, col="black", lwd=1, lty=2)


## Putting last 3 years NA in all States
mydata_fixNA3<-mydata_fix
for (i in (49*26+1):(49*29))
{
  mydata_fixNA3$fix_rate[i]<-NA
}

formula1NA <- fix_rate ~ fix_rate
chain1NA_3 <- ST.CARar(formula = formula1NA, family = "gaussian",
                       data = mydata_fixNA3, W = W, burnin = 20000, n.sample = 220000,
                       thin = 100)

low <- high <- missouriNA3 <- rep(0,3)
for (i in 0:2)
{
  missouriNA3[i+1] <-fitted.values(chain1NA_3)[24+K*(i+26)]
  low[i+1] <- p.interval(chain1NA_3$samples$fitted[,24+K*(i+26)], 
                         HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[1]
  high[i+1] <- p.interval(chain1NA_3$samples$fitted[,24+K*(i+26)], 
                          HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE)[2]  
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri3<-data.frame(missouriNA3,low,high)
x11()
ggplot(data.frame(x = 2015:2018, 
                  y = c(missouri_measures$fix_rate[26] ,missouri3$missouriNA3), 
                  ylow = c(missouri_measures$fix_rate[26], missouri3$low), 
                  yup = c(missouri_measures$fix_rate[26], missouri3$high)  )) +
  geom_line(aes(x = x, y = y), lty = 2) + theme_bw() +
  geom_line(data = data.frame(Y = missouri_measures$fix_rate, x = 1990:2018), aes(x = x, y = Y), col = "blue") +
  geom_ribbon(aes(ymin = ylow, ymax = yup, x = x), alpha = 0.1, col="lightblue") +
  geom_vline(aes(xintercept = 2015), col = 2, lty = 2) +
  ggtitle("ST.CARar model with only Missouri")

#----------------------------------END----------------------------------------

rm(list = ls() )
graphics.off()
