########################################################################
########################################################################
######################### BAY PROJECT ##################################
######################    CARBayesST     ###############################
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

setwd("C:/Users/Katarina/Desktop/Magistrale polimi/Ing Mat/M2/Bayesiana/project/data")
library(readxl)
mydata <- read_excel("datajay.xlsx")

BIRTH_data <- data.frame(mydata)
head(BIRTH_data)
dim(BIRTH_data)

birthdata_hs <- subset(BIRTH_data, BIRTH_data$Age.Group..Years. == "15-17 years")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Total U.S.")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Hawaii")
#birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "District of Columbia")
birthdata_hs <- subset(birthdata_hs, birthdata_hs$State != "Alaska")
birthdata_hs <- birthdata_hs[,c(2,1,5,4)]

birthdata_posths <- subset(BIRTH_data, BIRTH_data$Age.Group..Years. == "18-19 years")
birthdata_posths <- subset(birthdata_posths, birthdata_posths$State != "Total U.S.")
birthdata_posths <- subset(birthdata_posths, birthdata_posths$State != "Hawaii")
#birthdata_posths <- subset(birthdata_posths, birthdata_posths$State != "District of Columbia")
birthdata_posths <- subset(birthdata_posths, birthdata_posths$State != "Alaska")
birthdata_posths <- birthdata_posths[,c(2,1,5,4)]

####################################################################################
### Construct the SpatialPolygonsDataFrame
require(maps)
require(sp)
require(maptools)
require(tmap)
require(cartogram)


# Lower case letters for all states, in order to merde them with the map
mydata <- birthdata_posths # Change here to see birthdata_hs.

mydata <- data.frame(lapply(mydata, function(v) {
  if (is.character(v)) return(tolower(v))
  else return(v)
}))


# Get a SpatialPolygonsDataFrame of US states
usa <- map("state", fill = TRUE) # for some reason if I put false, it wont make the spatial polygons wtf
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
usa <- SpatialPolygonsDataFrame(usa, 
                                data = data.frame(unique(IDs), 
                                                  row.names = unique(IDs)) )

usa@data = data.frame(usa@data, mydata[match(usa@data[,'unique.IDs.'], mydata[,'State' ]),])
usa.data.combined = usa 
# for some reason there is only the yeas 1990,ok... in any case, we only need this for the W matrix
# remove NAs from the dataset
usa.data.combined@data <- usa.data.combined@data[complete.cases(usa.data.combined@data), ]

spplot(usa.data.combined[, "State.Rate"], main = "Birth rate" )
# if I remove district of columbia, It doesn't draw Wyoming, thus we leave 49 states, for now...

################################################################
### Constuct binary neighbourhood matrix W
W.nb <- poly2nb(usa.data.combined, row.names = unique(IDs)) # neighbourhood
is(W.nb)
W <- nb2mat(W.nb, style = "B")
is(W)
dim(W) # matrix of dimension 49 x 49
W[7:15,7:15]
W_nodoc <- W[c(-8),c(-8)]
W_nodoc[7:15,7:15] # weight matrix withoud district of columbia

# Moran I for each year.
moranI_years <- rep(0,29)
for (i in c (1990:2018))
{
  moran_subset <- subset(mydata, mydata$Year == i)
  mt <- moran.test(moran_subset$State.Rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1989] <- mt$estimate[1]
}
plot(c(1990:2018),moranI_years)
# strangely low dependencies for 1991, 2004 and 2008 (in the case of post highschool)
# in 1991 bad ordering for the states, missouri is first
# in 2004 bad ordering, oregon is first
# in 2008 bad ordering, nevada is first
# this should be fixed because it will give bad results in CARar.

########################################################################
### How "far" should we reach for  neighbours ( takes into account only 1990)
library(rgeos)
centroids <- gCentroid(usa,byid=T)
nb_knn7 <- knn2nb((knearneigh(centroids,k=7)))
sp_corG <- sp.correlogram(nb_knn7,usa.data.combined@data$State.Rate, order = 7,
                          zero.policy = TRUE, style = 'B', method = 'I')
plot(sp_corG, mail = 'Birth Rate correlogram')
# We can limit ourselves to 1st order neighbour, as they present the highest spatial correlation

########################################################################
# Boxplot of Birth rates by years
# a slight decrease in the variance among different states is observed over the years
ggplot(mydata, aes(x = factor(Year), y = State.Rate)) +
  geom_boxplot(fill="red", alpha=0.7) +
  scale_x_discrete(name = "Year") +
  scale_y_continuous(name = "Teenage birth rate") +
  theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold"))

########################################################################
### Spatial modelling with CARBayesST
##########################################
K<-49
N<-29
ones <- rep(1,K*N)


########################################################################
### Third order differences - Stationnarisation
### and fixing the state order which was not always the same (as noticed in MoransI)
################################
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
# plot(subset(mydata_staz, mydata_staz$rep.states..26.=='missouri')$dif3rate)
# checking if its the same as in 2nd prez. ok.

# check if things changed in morans I
# Moran I for each year for raw Rates.
moranI_years <- rep(0,29)
for (i in c(1990:2018))
{
  moran_subset <- subset(mydata_fix, years == i)
  mt <- moran.test(moran_subset$fix_rate,nb2listw(W.nb,style = "B"))
  moranI_years[i-1989] <- mt$estimate[1]
}
plot(c(1990:2018),moranI_years)


# Moran I for each year for 3rd differences
#moranI_years <- rep(0,26)
#for (i in c (1993:2018))
#{
#  moran_subset <- subset(mydata_staz, mydata$Year == i)
#  mt <- moran.test(moran_subset$dif3rate,nb2listw(W.nb,style = "B"))
#  moranI_years[i-1992] <- mt$estimate[1]
#}
#plot(c(1993:2018),moranI_years)
# spatial correlation is lost by taking 3rd order diff
# so it does't really make sense to fit 


########################################################################
### Spatial modelling with CARBayesST
## on stationnarized data
##########################################

#formula2 <- dif3rate ~ dif3rate 
#chain2 <- ST.CARar(formula = formula2, family = "gaussian",
#                   data = mydata_staz, W = W, burnin = 20000, n.sample = 220000,
#                   thin = 100)

# Posterion inference 
#print(chain2) 
#x11()
#plot(chain2$samples$tau2)
#plot(chain2$samples$rho[,1])
#plot(chain2$samples$rho[,2])
# very ugly.

########################################################################
### Spatial modelling with CARBayesST
##  on not stationnarized data
##########################################

?ST.CARar

formula1 <- fix_rate ~ fix_rate 
chain1posths <- ST.CARar(formula = formula1, family = "gaussian",
                         data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                         thin = 100)

# Posterion inference 

print(chain1posths) 
x11()
plot(chain1posths$samples$tau2)
plot(chain1posths$samples$rho)

# chain1posths$fitted.values[1100:1150]
# A vector of fitted values for each area and time period.
# x11()
# plot(chain1posths$samples$phi[,1100])
# thus the phis are not the fitted values, but we sum beta and phi to obtain the fitted values.
# chain1posths$fitted.values[1100:1150]-mydata_fix$fix_rate[1100:1150]
# it is very close to our real data, yaay.

# some results for post high school age group
# very nice mcmc chains for chain3 :)
#Median    2.5%   97.5% n.effective Geweke.diag
#(Intercept) 66.4455 66.3986 66.4924      2000.0        -0.5
#tau2        40.3406 36.2776 45.0259      1487.3        -0.4
#nu2          0.8108  0.3653  1.2820      1240.2        -0.6
#rho.S        0.9547  0.9286  0.9732      2000.0        -1.0
#rho.T        0.9714  0.9568  0.9858      2000.0        -0.8

#DIC =  4854.274       p.d =  1179.973       LMPL =  -2807.509 

# and the results for the high schoolaged group
#Median    2.5%   97.5% n.effective Geweke.diag
#(Intercept) 22.3147 22.3099 22.3197      2000.0         0.4
#tau2        14.9766 13.9675 16.1623      2000.0         0.1
#nu2          0.0066  0.0020  0.0252       361.1        -1.7
#rho.S        0.9546  0.9304  0.9740      2000.0        -0.7
#rho.T        0.9426  0.9244  0.9600      2000.0        -1.1
#DIC =  -1952.891       p.d =  1114.852       LMPL =  396.1347 

########################################################################
### Spatial modelling with CARBayesST
##  on not stationnarized data, without district of columbia
##########################################

mydata_nodoc <- subset(mydata_fix, mydata_fix$rep.states..29. != "district of columbia")

formula2 <- fix_rate ~ fix_rate
chain2 <- ST.CARar(formula = formula2, family = "gaussian",
                   data = mydata_nodoc, W = W_nodoc, burnin = 20000, n.sample = 220000,
                   thin = 100)

# Posterion inference 
print(chain2) 
x11()
plot(chain2$samples$beta)
plot(chain2$samples$tau2)
plot(chain2$samples$nu2)
plot(chain2$samples$rho)

# doesn't seem to get better without doc
#################################################################

# Since our goal is to predict the number of pregnancies in the upcoming years
# We should make some predictions for 2019 for example, shaow that on some nice plots and thnx bye

## Maybe make a linear regression on health data and then predict ?

##################################################################
### Fitting ST.CARadaptive model
##  takes much longer to fit than the ST.CARar

formula3 <- fix_rate ~ fix_rate
chain3 <- ST.CARadaptive(formula = formula3, family = "gaussian",
                         data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                         thin = 100)
# Posterion inference 
print(chain3) 
x11()
plot(chain3$samples)
plot(chain3$samples$tau2)
plot(chain3$samples$nu2)
plot(chain3$samples$rho)
# slightly better indexes, slightly worse n.effective, but I don't understand the model.

###################################################################
### Fitting ST.CARadaptive model
##  takes much longer to fit than the ST.CARar

mydata_fixNA<-mydata_fix
# missouri in number 24, putting last 3 years to NA
for (i in 0:28)
{
  mydata_fixNA$fix_rate[24+K*i]<-NA
}


formula1NA <- fix_rate ~ fix_rate
chain1NA <- ST.CARar(formula = formula1NA, family = "gaussian",
                     data = mydata_fixNA, W = W, burnin = 20000, n.sample = 220000,
                     thin = 100)
# Posterion inference 
print(chain1NA) 
x11()
plot(chain1NA$samples)
plot(chain1NA$samples$tau2)
plot(chain1NA$samples$nu2)
plot(chain1NA$samples$rho)
# fitted.values(chain1NA)[24+K*28]
missouriNA <- rep(0,29)
missouri_measures <- rep(0,29)

for (i in 0:28)
{
  missouriNA[i+1] <-fitted.values(chain1NA)[24+K*i]
}
missouri_measures <- subset(mydata_fix, mydata_fix$rep.states..29. == "missouri")
missouri<-data.frame(missouri_measures,missouriNA)
# when we take out only the last 3 values ,the values that were NA are underestimated, 
# this happens because our time series is not stationnary
# on the other hand by stationarizing it we loose the spatial correlation
# by making NA the whole missouri column, they are overestimated
# this happens thus becaure we take out only an intercept out of our data (phi)
# try using a covariate first...

##################################################################
### Fitting ST.CARar model
##  but this time, take out a linear regression, not only intercept
##  think about it.

formula4 <- fix_rate ~ fix_rate + years
chain4 <- ST.CARadaptive(formula = formula4, family = "gaussian",
                         data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                         thin = 100)
# Posterion inference 
print(chain4) 

years <- c(rep(1990,49),rep(1991,49),rep(1992,49),rep(1993,49),rep(1994,49),rep(1995,49)
           ,rep(1996,49),rep(1997,49),rep(1998,49),rep(1999,49),rep(2000,49),rep(2001,49)
           ,rep(2002,49),rep(2003,49),rep(2004,49),rep(2005,49),rep(2006,49),rep(2007,49)
           ,rep(2008,49),rep(2009,49),rep(2010,49),rep(2011,49),rep(2012,49),rep(2013,49)
           ,rep(2014,49),rep(2015,49),rep(2016,49),rep(2017,49),rep(2018,49))

formula5 <- fix_rate ~ fix_rate + years
chain5 <- ST.CARar(formula = formula5, family = "gaussian",
                   data = mydata_fix, W = W, burnin = 20000, n.sample = 220000,
                   thin = 100)
# Posterion inference 
print(chain5) 
x11()
plot(chain5$samples$tau2)
plot(chain5$samples$nu2)
plot(chain5$samples$rho)
