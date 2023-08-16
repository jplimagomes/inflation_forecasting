# Inflation forecasting exercises for the US

## Description of the inflation forecasting exercise

This project consists of an inflation forecasting exercise using a large dataset of macroeconomic variables of monthly frequency. We will use four models: (1) Autoregressive model; (2) Principal Component Regression (PCR) added an autoregressive component; (3) Ridgre Regression; and (4) Lasso.

The data consist of variables from the FED database, FRED-MED. The sample spans from January 1959 to November 2021, or 755 observations. In this exercise, we will only consider variables with all observations available. Eight groups contain all variables, which can be divided into: (1) Output and Income; (2) Labor Market; (3) Housing; (4) Consumption, Orders and Inventories; (5) Money and Credit; (6) Interest and Exchange Rates; (7) Prices; and (8) Stock Market. We will divide the variables among these groups in the second part of the exercise in order to check the importance of each of the eight groups in the forecasts. The sample is contained in the file '2021-12.csv'.

Initially, the variables are transformed to become approximately stationary. Each column of the file '2021-12.csv' corresponds to a variable, and the rows starting from the third correspond to the observations. In the first line of the file '2021-12.csv', there is an identifier of the variable; in the second, there is a number indicating the transformation to be used for the corresponding variable. We will use the 'fbi' package available at https://research.stlouisfed.org/econ/mccracken/fred-databases/ in order to perform the transformations.

For this exercise, the variable used to construct the dependent variable is the CPIAUCSL price index (CPI all items). In this exercise, this variable will not be transformed according to the indication in the csv file, but according to the expression $\pi_t = \frac{\Delta y_t}{y_t}$, where $y_t$ is the price level in period $t$.

We will use four linear models for inflation forecasting: AR, PCR + AR, Ridge and LASSO.

## AR model

Only the variable $\pi$ will be used, with the AR order defined by the BIC criterion.

## AR + PCR model

In the PCR part, the explanatory variables are given by the factor estimates given by the principal component analysis of the whole dataset. The number of factors is defined as the minimum amount able to explain 90% of the variance. In the AR part, we will adopt order equal to 2.

## Ridge & LASSO models

For both models, we will use all variables from the dataset except CPIAUCSL. The penalty term for both models is selected by the BIC criterion.

# Framework

Forecasts using each of the four models are based on a fixed-size rolling-window framework of 492 observations. For each model, the procedure is described according to the three steps:

1) The model is estimated for each window, from observation $a$ to observation $a + 492 - 1$.

2) Then, the estimated model is used to forecast for the observation at position $a + 492$.

3) We set $a = a + 1$, and the above two steps are repeated.

For each model, we compute the squared error for the next observation and plot the cumulative squared errors. We thus check the forecasting performance of each of the models.

# Importance of groups for inflation

Next, for the PCR + AR, Ridge and LASSO models, we compute a measure of importance of groups of variables. For each window, we will compute the importance of each variable in the sample, and then aggregate the importance of each variable in its respective group. In this case, we will create an additional group called 'lag', consisting only of lags of the dependent variable, which is 'inflation'.

For the Ridge and LASSO models, the importance of the variable is measured from the product of the estimated coefficient and the standard deviation of the variable. With respect to the AR part of the PCR + AR model, the relative importance of the "inflation" lags is given in the same way as the previous two models. In the PCR part, each estimated factor is a combination of the original variables, given by $F_{it} = \alpha'_i X_t$, where $X_t$ is the vector containing all variables (except inflation lags) in period $t$. Then, the importance of each variable $X_t$ is given by the product of the corresponding element in the vector $\alpha$ with the parameter estimated in the PCR model.

Finally, the importances are grouped into the nine groups and the results are normalized by setting a value equal to 100 for the highest result in module. Graphs are displayed to compare the distinct group importances between the three models. 
