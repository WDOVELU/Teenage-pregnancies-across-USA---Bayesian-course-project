# Teenage pregnanciesacross the United States (Bayesian_Statistics_Project_2020Fall)
This is a project in Bayesian Statistics at Politecnico di Milano fall of 2020. Here you will find the main code and data that is used for our analysis and a [paper](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/Paper.pdf) summarizing the entire work.

Our dataset is made out of observations across 50 States of America, over 29 years (1990-2018) for 2 separate age groups: 
+ girls aged 15-17, also called high-school group 
+ girls aged 18-19, also called post high-school group.   
This counts for a total of 50x29x2 = 2900 observations of Birth rates.   
<img src="https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/3_Pictures/Markdown_pic/dataset_1.jpg" width="50%" height="50%">


The project is based on two steps, first time series analysis from bayesian view and then Spatio-Temporal analysis based on CARBayesST package. 


## Table of contents
* [Structure](#Structure)
* [R packages requirement and  Installation](#R packages requirement and  Installation)
* [First part: Time series analysis](#First part: Time series analysis)
* [Second part: Spatio-Temporal analysis](#Second part: Spatio-Temporal analysis)
* [Conclusion](#Conclusion)

## Structure
Useful hints about our Github.
In [0_Dataset](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/tree/master/0_Dataset) you will find the dataset used in the study.

In [1_Time_series_analysis](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/tree/master/1_Time_series_analysis) folder you will find the following codes, all regarding time series analysis where contains 4 ARIMA models.

* [full_script.R](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/1_Time_series_analysis/full_script.R) - Main code file. 
* [ARIMA_2_3_3.stan](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/1_Time_series_analysis/ARIMA_2_3_3.stan) - Stan file related to ARIMA(2,3,3) model.
* [StSp_model_2ages.stan](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/1_Time_series_analysis/StSp_model_2ages.stan) - Stan file related to ARIMA(1,3,0) model shared a common baseline. 
* [StSp_model_2ages_second_version.stan](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/1_Time_series_analysis/StSp_model_2ages_second_version.stan) - Stan file related to ARIMA(1,3,0) model shared a common temporal depence parameter. 
* [StSp_model_2ages_third_version.stan](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/1_Time_series_analysis/StSp_model_2ages_third_version.stan) - Stan file related to ARIMA(1,3,0) model shared nothing, which means the baseline parameter mu0_j and temporal dependence parameters phi_j is different for 2 groups.


In [2_CARBayesST](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/tree/master/2_CARBayesST) folder you will find the following codes:
* [SpatioTemporal_Gaussian.R](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/2_CARBayesST/SpatioTemporal_Gaussian.R) - Main code file for ST.carar model using gaussian likehood.


* [SpatioTemporal_Highschool_Possion_model.R](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/2_CARBayesST/SpatioTemporal_Highschool_Possion_model.R) - Code file for ST.carar model using the integer part of brith rate and  possion likehood.
* [SpatioTemporal_Highschool_Possion_model_data10.R](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/blob/master/2_CARBayesST/SpatioTemporal_Highschool_Possion_model_data10.R) - code file for ST.carar model multipling brith rate by 10 and  using possion likehood.

In [3_Pictures](https://github.com/WDOVELU/Teenage-pregnancies-across-USA---Bayesian-course-project/tree/master/3_Pictures), there are pictures using in paper and makedown documents.

## R packages requirement and  Installation
Before running our code, please install all R packages we used using R code below:
```R
needed_packages  <- c("rstan","coda","tseries","ggplot2","tidyr","dplyr","purrr","ggsci","readxl","Metrics","MLmetrics","CARBayesST","spdep","maptools", "sp","CARBayesdata","USAboundaries", "sf","dplyr","LaplacesDemon","rgeos")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])] if (length(new_packages))
install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)
```

## First part: Time series analysis








## Second part: Spatio-Temporal analysis





## Conclusion

With this work we were able to apply some of the techniques offered by the Bayesian statistical approach to a real case study. We finally were able to successfully fit two models, the ARIMA(1,3,0) and ST.CARar model.

We have implemented the ARIMA(1,3,0) in STAN language. Our model takes into account only one state, and introduces random effects on the two age groups. It can be used to forecast future Birth rates in a state for the upcoming years, for both age groups. Before fitting the model, data is stationnarized by means of third order differences.

As a second model, the ST.CARar model from CARBayesST was used, in order to work with areal data. This model works with raw (not-stationnarized) data, and takes into account only one age group. We can imagine to use this model in order to assess whether the data observed in one state seems reasonable with respect to what would be expected regarding its neighbouring states, to obtain Birth rates for those states that hide their Birth rate data, or
maybe even detect fraudulent data. 

This being said, no model is better than the other. The two models are complementary and give answers to different questions.



























Made by:  
Radisic Katarina   
Santamaria Andrea   
Lu Jiajie
