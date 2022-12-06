
# Modelo de Regresión Random Forest (RF)
# Karen Brigitte Mejía Correal
# Máster en Geoinformática para la Gestión de Recursos Naturales 
# Diciembre 2022
# Desarrollado con R 4.2.0 y RStudio (2022.07.2+576)

############################################################################################################################################################
  
# Eliminar todos los objetos de memoria
rm(list=ls(globalenv()))

# Lista de paquetes necesarios:
.packages = c("raster","sp","rgdal","maptools","rgeos","foreign","corrplot","lattice",
              "ggplot2","mlbench","caret","ggpubr","tidyverse","randomForest","chillR","plyr")

# Instalar los paquetes desde el CRAN packages si no se encuentran instalados en el ordenador:
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Cargar los paquetes eliminando los mensajes advertencia y sin ellos:
suppressMessages(invisible(lapply(.packages, require, character.only=TRUE)))

############################################################################################################################################################


# Cargar librerias
library(raster)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(foreign)
library(corrplot)

library(lattice)
library(ggplot2)
library(mlbench)
library(caret)
library (ggpubr)
library(tidyverse)

library(randomForest)
library(chillR)

############################################################################################################################################################
##------- MODELO RANDOM FOREST (RF)-------------------------------------------------------------------------------------------------------------------------
############################################################################################################################################################

# Establecer el directorio de trabajo
setwd("C:/TFM_SST_Uvas")

# Cargar fichero en formato CSV para cada variedad (Godello, Verdejo, Mencía y Tempranillo)
# Datos crudos (Raw) o SNV. Disponibles en: https://github.com/karenbrigitte/TFM_SST_Uvas.git
datos<-read.csv2(choose.files(default = "C:/TFM_SST_Uvas"))

# Eliminar valores NA (sin valores)
datos<-drop_na(datos[,])


############################################################################################################################################################
##------- Toda la reflectancia (400-2500 nm) ---------------------------------------------------------------------------------------------------------------
############################################################################################################################################################

## 1.PROCESO DE MODELADO -----------------------------------------------------------------------------------------------------------------------------------

# Selección columnas de interés
data<-datos[,c(2,53:2153)]

# Establecer el conjunto de datos
trainData <- subset(data[,])

# Configuración de prueba
# Datos de entrenamiento
myRespName <- 'TSS' # Nombre de la variable dependiente a modelizar
y_trainData <- trainData[,1] # Posición de la variable dependiente a modelizar en el trainData
x_trainData <- trainData[,2:2102] # Conjunto de variables independientes del trainData

# Cálculo estadísticos variable dependiente (SST)
myRespName_mean <- mean(y_trainData)
myRespName_max <- max(y_trainData)
myRespName_min <- min(y_trainData)
myRespName_sd <- sd(y_trainData)

# Asegurar que los resultados sean repetibles
seed <- 7
set.seed(seed)
# Métrica R-cuadrado
metric <- "Rsquared"

# Algoritmo de prueba
# Crear trainControl
rf_fit_control <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=10,
                               search="grid",
                               savePredictions=TRUE)

# Grid Search
rf_fit_tunegrid <- expand.grid(.mtry=c(1:15))
#mtry: número de variables muestreadas aleatoriamente como candidatas en cada división.
#ntree: número de árboles.

set.seed(seed)

# Parámetros del modelo
rf_fit_gridsearch <- train(TSS ~.,
                           data=trainData, na.action=na.exclude, method="rf", metric=metric, ntree=1500, 
                           tuneGrid=rf_fit_tunegrid, trControl=rf_fit_control, allowParallel = TRUE)

# NOTA: en el caso de utilizar una lista de variables reducidas se debe reemplazar el punto en "TSS~.,"
# por las variables a usar, por ejemplo: "TSS~X709+X868+X896+X905,"

print(rf_fit_gridsearch)
plot(rf_fit_gridsearch)

# Importancia de las variables
varImp(rf_fit_gridsearch, scale = FALSE) #  Mean Decrease Gini
varImp(rf_fit_gridsearch)## Importancia en %


## 2. MODELO DE ENTRENAMIENTO (Bondad de ajuste de un modelo estadístico) --------------------------------------------------------------------------------             

# Modelo de entrenamiento: Predichos y Residuos
set.seed(seed)
rf_fit_gridsearch_predicted <- predict(rf_fit_gridsearch)
rf_fit_gridsearch_residuals <- residuals(rf_fit_gridsearch)

# Cálculo RPD
RPD(rf_fit_gridsearch_predicted,y_trainData) 
RPD <- formatC(RPD(rf_fit_gridsearch_predicted,y_trainData), digits = 2, format = "f")

# Cálculo RMSE
rf_fit_RMSE = sqrt(mean(rf_fit_gridsearch_residuals^2))
cat('The root mean square error of the training data is ', round(rf_fit_RMSE,3), '\n')

# R-cuadrado. Cálculo de la suma de los cuadrados
rf_y_fit_mean = mean(y_trainData)
rf_fit_tss =  sum((y_trainData - rf_y_fit_mean)^2)

# Cálculo suma residual de cuadrados
rf_fit_rss =  sum(rf_fit_gridsearch_residuals^2)

# Cálculo R-cuadrado
rf_fit_rsq  =  1 - (rf_fit_rss/rf_fit_tss)
cat('The R-square of the training data is ', round(rf_fit_rsq,3), '\n')

# Cálculo R-cuadrado del modelo
R2 <- formatC(mean(rf_fit_gridsearch$results$Rsquared), digits = 2, format = "f")

# Cálculo RMSE del modelo
RMSE <- formatC(mean(rf_fit_gridsearch$results$RMSE), digits = 2, format = "f")


## 3. VALIDACIÓN CRUZADA (CV) DEL CONJUNTO DE DATOS (Bondad de la prueba de un modelo estadístico)

# Evaluación del modelo
set.seed(seed)
rf_cv_all_gridsearch <- rf_fit_gridsearch$pred
rf_cv_all_gridsearch$Fold <- substr(rf_cv_all_gridsearch$Resample, start=1 , stop=6)
rf_cv_all_gridsearch$Rep <- substr(rf_cv_all_gridsearch$Resample, start=8 , stop=12)
rf_cv_all_gridsearch$Rep_Fold <- paste(rf_cv_all_gridsearch$Rep, rf_cv_all_gridsearch$Fold)
rf_cv_all_gridsearch$residuals <- (rf_cv_all_gridsearch$obs - rf_cv_all_gridsearch$pred)


# Predichos y Observados modelo CV 
# Configuración de los parámetros del mejor modelo
rf_cv_rep_gridsearch <- subset(rf_cv_all_gridsearch, mtry==5) # mtry== mtry obtenido en el mejor modelo

# Cálculo RMSE
rf_cv_rep_gridsearch$residuals2 = rf_cv_rep_gridsearch$residuals^2
rf_cv_rep_RMSE = sqrt(tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, mean))
rf_cv_RMSE_mean = mean(rf_cv_rep_RMSE)
rf_cv_RMSE_sd = sd(rf_cv_rep_RMSE)
cat('The mean of the root mean square error of the testing data is ', round(rf_cv_RMSE_mean,3), '\n')
cat('The sd of the root mean square error of the testing data is ', round(rf_cv_RMSE_sd,3), '\n')
rf_cv_RMSE <- as.data.frame(rf_cv_rep_RMSE)

# Cálculo R-cuadrado
# Cálculo suma residual de cuadrados
rf_cv_rep_rss =  tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, sum)

# Cálculo de la suma de los cuadrados
# Función para calcular la media y la desviación estándar de cada grupo
library(plyr)
function_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
rf_cv_rep_obs_mean <- function_summary(rf_cv_rep_gridsearch, varname="obs", 
                                       groupnames=c("Rep_Fold"))

rf_cv_rep_gridsearch_merge <- merge (rf_cv_rep_gridsearch, rf_cv_rep_obs_mean, by = "Rep_Fold")
rf_cv_rep_tss = tapply((rf_cv_rep_gridsearch_merge$obs - rf_cv_rep_gridsearch_merge$mean)^2,
                       rf_cv_rep_gridsearch_merge$Rep_Fold, sum)
rf_cv_rep_rsq  =  1 - (rf_cv_rep_rss/rf_cv_rep_tss)
rf_cv_rsq_mean = mean(rf_cv_rep_rsq)
rf_cv_rsq_sd = sd(rf_cv_rep_rsq)
cat('The mean R-square of the testing data is ', round(rf_cv_rsq_mean,3), '\n')
cat('The std R-square of the testing data is ', round(rf_cv_rsq_sd,3), '\n')

# Gráfica Predichos vs Observados
options(repr.plot.width=8, repr.plot.height=4)
rf_fit_gridsearch_predvsobs = as.data.frame(cbind(predicted = rf_fit_gridsearch_predicted, 
                                                  observed = y_trainData))
ggplot(rf_fit_gridsearch_predvsobs,aes(predicted, observed)) + 
  geom_point(color = "violetred4", alpha = 0.75) + 
  geom_smooth(method=lm, color="orangered", size=1)+
  ggtitle("Random Forest (RF): Predichos vs Observados") +
  xlab("SST predichos (º Brix)") + ylab("SST observados (º Brix)") + 
  theme(plot.title = element_text(color="black",size=16,hjust = 0.5, family = "serif",face="bold"),
        axis.text.y = element_text(face="bold", size=12, family = "serif"), 
        axis.text.x = element_text(face="bold", size=12,hjust=.5, family = "serif"),
        axis.title.x = element_text(size=14, family = "serif"), 
        axis.title.y = element_text(size=14, family = "serif"),
        axis.line = element_line(size = 0.8, colour = "black"))+
  geom_abline (slope=1, intercept=0, linetype = "dashed", color="blue")


############################################################################################################################################################
##------- Rango VIS (400-700 nm) ---------------------------------------------------------------------------------------------------------------
############################################################################################################################################################

## 1.PROCESO DE MODELADO -----------------------------------------------------------------------------------------------------------------------------------

# Selección columnas de interés
data<-datos[,c(2,53:353)]

# Establecer el conjunto de datos
trainData <- subset(data[,])

# Configuración de prueba
# Datos de entrenamiento
myRespName <- 'TSS' # Nombre de la variable dependiente a modelizar
y_trainData <- trainData[,1] # Posición de la variable dependiente a modelizar en el trainData
x_trainData <- trainData[,2:302] # Conjunto de variables independientes del trainData

# Cálculo estadísticos variable dependiente (SST)
myRespName_mean <- mean(y_trainData)
myRespName_max <- max(y_trainData)
myRespName_min <- min(y_trainData)
myRespName_sd <- sd(y_trainData)

# Asegurar que los resultados sean repetibles
seed <- 7
set.seed(seed)
# Métrica R-cuadrado
metric <- "Rsquared"

# Algoritmo de prueba
# Crear trainControl
rf_fit_control <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=10,
                               search="grid",
                               savePredictions=TRUE)

# Grid Search
rf_fit_tunegrid <- expand.grid(.mtry=c(1:15))
#mtry: número de variables muestreadas aleatoriamente como candidatas en cada división.
#ntree: número de árboles.

set.seed(seed)

# Parámetros del modelo
rf_fit_gridsearch <- train(TSS ~.,
                           data=trainData, na.action=na.exclude, method="rf", metric=metric, ntree=1500, 
                           tuneGrid=rf_fit_tunegrid, trControl=rf_fit_control, allowParallel = TRUE)

print(rf_fit_gridsearch)
plot(rf_fit_gridsearch)

# Importancia de las variables
varImp(rf_fit_gridsearch, scale = FALSE) #  Mean Decrease Gini
varImp(rf_fit_gridsearch)## Importancia en %


## 2. MODELO DE ENTRENAMIENTO (Bondad de ajuste de un modelo estadístico) --------------------------------------------------------------------------------             

# Modelo de entrenamiento: Predichos y Residuos
set.seed(seed)
rf_fit_gridsearch_predicted <- predict(rf_fit_gridsearch)
rf_fit_gridsearch_residuals <- residuals(rf_fit_gridsearch)

# Cálculo RPD
RPD(rf_fit_gridsearch_predicted,y_trainData) 
RPD <- formatC(RPD(rf_fit_gridsearch_predicted,y_trainData), digits = 2, format = "f")

# Cálculo RMSE
rf_fit_RMSE = sqrt(mean(rf_fit_gridsearch_residuals^2))
cat('The root mean square error of the training data is ', round(rf_fit_RMSE,3), '\n')

# R-cuadrado. Cálculo de la suma de los cuadrados
rf_y_fit_mean = mean(y_trainData)
rf_fit_tss =  sum((y_trainData - rf_y_fit_mean)^2)

# Cálculo suma residual de cuadrados
rf_fit_rss =  sum(rf_fit_gridsearch_residuals^2)

# Cálculo R-cuadrado
rf_fit_rsq  =  1 - (rf_fit_rss/rf_fit_tss)
cat('The R-square of the training data is ', round(rf_fit_rsq,3), '\n')

# Cálculo R-cuadrado del modelo
R2 <- formatC(mean(rf_fit_gridsearch$results$Rsquared), digits = 2, format = "f")

# Cálculo RMSE del modelo
RMSE <- formatC(mean(rf_fit_gridsearch$results$RMSE), digits = 2, format = "f")


## 3. VALIDACIÓN CRUZADA (CV) DEL CONJUNTO DE DATOS (Bondad de la prueba de un modelo estadístico)

# Evaluación del modelo
set.seed(seed)
rf_cv_all_gridsearch <- rf_fit_gridsearch$pred
rf_cv_all_gridsearch$Fold <- substr(rf_cv_all_gridsearch$Resample, start=1 , stop=6)
rf_cv_all_gridsearch$Rep <- substr(rf_cv_all_gridsearch$Resample, start=8 , stop=12)
rf_cv_all_gridsearch$Rep_Fold <- paste(rf_cv_all_gridsearch$Rep, rf_cv_all_gridsearch$Fold)
rf_cv_all_gridsearch$residuals <- (rf_cv_all_gridsearch$obs - rf_cv_all_gridsearch$pred)


# Predichos y Observados modelo CV 
# Configuración de los parámetros del mejor modelo
rf_cv_rep_gridsearch <- subset(rf_cv_all_gridsearch, mtry==5) # mtry== mtry obtenido en el mejor modelo

# Cálculo RMSE
rf_cv_rep_gridsearch$residuals2 = rf_cv_rep_gridsearch$residuals^2
rf_cv_rep_RMSE = sqrt(tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, mean))
rf_cv_RMSE_mean = mean(rf_cv_rep_RMSE)
rf_cv_RMSE_sd = sd(rf_cv_rep_RMSE)
cat('The mean of the root mean square error of the testing data is ', round(rf_cv_RMSE_mean,3), '\n')
cat('The sd of the root mean square error of the testing data is ', round(rf_cv_RMSE_sd,3), '\n')
rf_cv_RMSE <- as.data.frame(rf_cv_rep_RMSE)

# Cálculo R-cuadrado
# Cálculo suma residual de cuadrados
rf_cv_rep_rss =  tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, sum)

# Cálculo de la suma de los cuadrados
# Función para calcular la media y la desviación estándar de cada grupo
library(plyr)
function_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
rf_cv_rep_obs_mean <- function_summary(rf_cv_rep_gridsearch, varname="obs", 
                                       groupnames=c("Rep_Fold"))

rf_cv_rep_gridsearch_merge <- merge (rf_cv_rep_gridsearch, rf_cv_rep_obs_mean, by = "Rep_Fold")
rf_cv_rep_tss = tapply((rf_cv_rep_gridsearch_merge$obs - rf_cv_rep_gridsearch_merge$mean)^2,
                       rf_cv_rep_gridsearch_merge$Rep_Fold, sum)
rf_cv_rep_rsq  =  1 - (rf_cv_rep_rss/rf_cv_rep_tss)
rf_cv_rsq_mean = mean(rf_cv_rep_rsq)
rf_cv_rsq_sd = sd(rf_cv_rep_rsq)
cat('The mean R-square of the testing data is ', round(rf_cv_rsq_mean,3), '\n')
cat('The std R-square of the testing data is ', round(rf_cv_rsq_sd,3), '\n')

# Gráfica Predichos vs Observados
options(repr.plot.width=8, repr.plot.height=4)
rf_fit_gridsearch_predvsobs = as.data.frame(cbind(predicted = rf_fit_gridsearch_predicted, 
                                                  observed = y_trainData))
ggplot(rf_fit_gridsearch_predvsobs,aes(predicted, observed)) + 
  geom_point(color = "violetred4", alpha = 0.75) + 
  geom_smooth(method=lm, color="orangered", size=1)+
  ggtitle("Random Forest (RF): Predichos vs Observados") +
  xlab("SST predichos (º Brix)") + ylab("SST observados (º Brix)") + 
  theme(plot.title = element_text(color="black",size=16,hjust = 0.5, family = "serif",face="bold"),
        axis.text.y = element_text(face="bold", size=12, family = "serif"), 
        axis.text.x = element_text(face="bold", size=12,hjust=.5, family = "serif"),
        axis.title.x = element_text(size=14, family = "serif"), 
        axis.title.y = element_text(size=14, family = "serif"),
        axis.line = element_line(size = 0.8, colour = "black"))+
  geom_abline (slope=1, intercept=0, linetype = "dashed", color="blue")


############################################################################################################################################################
##------- Rango NIR (701-1000 nm) ---------------------------------------------------------------------------------------------------------------
############################################################################################################################################################

## 1.PROCESO DE MODELADO -----------------------------------------------------------------------------------------------------------------------------------

# Selección columnas de interés
data<-datos[,c(2,354:653)]

# Establecer el conjunto de datos
trainData <- subset(data[,])

# Configuración de prueba
# Datos de entrenamiento
myRespName <- 'TSS' # Nombre de la variable dependiente a modelizar
y_trainData <- trainData[,1] # Posición de la variable dependiente a modelizar en el trainData
x_trainData <- trainData[,2:301] # Conjunto de variables independientes del trainData

# Cálculo estadísticos variable dependiente (SST)
myRespName_mean <- mean(y_trainData)
myRespName_max <- max(y_trainData)
myRespName_min <- min(y_trainData)
myRespName_sd <- sd(y_trainData)

# Asegurar que los resultados sean repetibles
seed <- 7
set.seed(seed)
# Métrica R-cuadrado
metric <- "Rsquared"

# Algoritmo de prueba
# Crear trainControl
rf_fit_control <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=10,
                               search="grid",
                               savePredictions=TRUE)

# Grid Search
rf_fit_tunegrid <- expand.grid(.mtry=c(1:15))
#mtry: número de variables muestreadas aleatoriamente como candidatas en cada división.
#ntree: número de árboles.

set.seed(seed)

# Parámetros del modelo
rf_fit_gridsearch <- train(TSS ~.,
                           data=trainData, na.action=na.exclude, method="rf", metric=metric, ntree=1500, 
                           tuneGrid=rf_fit_tunegrid, trControl=rf_fit_control, allowParallel = TRUE)

print(rf_fit_gridsearch)
plot(rf_fit_gridsearch)

# Importancia de las variables
varImp(rf_fit_gridsearch, scale = FALSE) #  Mean Decrease Gini
varImp(rf_fit_gridsearch)## Importancia en %


## 2. MODELO DE ENTRENAMIENTO (Bondad de ajuste de un modelo estadístico) --------------------------------------------------------------------------------             

# Modelo de entrenamiento: Predichos y Residuos
set.seed(seed)
rf_fit_gridsearch_predicted <- predict(rf_fit_gridsearch)
rf_fit_gridsearch_residuals <- residuals(rf_fit_gridsearch)

# Cálculo RPD
RPD(rf_fit_gridsearch_predicted,y_trainData) 
RPD <- formatC(RPD(rf_fit_gridsearch_predicted,y_trainData), digits = 2, format = "f")

# Cálculo RMSE
rf_fit_RMSE = sqrt(mean(rf_fit_gridsearch_residuals^2))
cat('The root mean square error of the training data is ', round(rf_fit_RMSE,3), '\n')

# R-cuadrado. Cálculo de la suma de los cuadrados
rf_y_fit_mean = mean(y_trainData)
rf_fit_tss =  sum((y_trainData - rf_y_fit_mean)^2)

# Cálculo suma residual de cuadrados
rf_fit_rss =  sum(rf_fit_gridsearch_residuals^2)

# Cálculo R-cuadrado
rf_fit_rsq  =  1 - (rf_fit_rss/rf_fit_tss)
cat('The R-square of the training data is ', round(rf_fit_rsq,3), '\n')

# Cálculo R-cuadrado del modelo
R2 <- formatC(mean(rf_fit_gridsearch$results$Rsquared), digits = 2, format = "f")

# Cálculo RMSE del modelo
RMSE <- formatC(mean(rf_fit_gridsearch$results$RMSE), digits = 2, format = "f")


## 3. VALIDACIÓN CRUZADA (CV) DEL CONJUNTO DE DATOS (Bondad de la prueba de un modelo estadístico)

# Evaluación del modelo
set.seed(seed)
rf_cv_all_gridsearch <- rf_fit_gridsearch$pred
rf_cv_all_gridsearch$Fold <- substr(rf_cv_all_gridsearch$Resample, start=1 , stop=6)
rf_cv_all_gridsearch$Rep <- substr(rf_cv_all_gridsearch$Resample, start=8 , stop=12)
rf_cv_all_gridsearch$Rep_Fold <- paste(rf_cv_all_gridsearch$Rep, rf_cv_all_gridsearch$Fold)
rf_cv_all_gridsearch$residuals <- (rf_cv_all_gridsearch$obs - rf_cv_all_gridsearch$pred)


# Predichos y Observados modelo CV 
# Configuración de los parámetros del mejor modelo
rf_cv_rep_gridsearch <- subset(rf_cv_all_gridsearch, mtry==5) # mtry== mtry obtenido en el mejor modelo

# Cálculo RMSE
rf_cv_rep_gridsearch$residuals2 = rf_cv_rep_gridsearch$residuals^2
rf_cv_rep_RMSE = sqrt(tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, mean))
rf_cv_RMSE_mean = mean(rf_cv_rep_RMSE)
rf_cv_RMSE_sd = sd(rf_cv_rep_RMSE)
cat('The mean of the root mean square error of the testing data is ', round(rf_cv_RMSE_mean,3), '\n')
cat('The sd of the root mean square error of the testing data is ', round(rf_cv_RMSE_sd,3), '\n')
rf_cv_RMSE <- as.data.frame(rf_cv_rep_RMSE)

# Cálculo R-cuadrado
# Cálculo suma residual de cuadrados
rf_cv_rep_rss =  tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, sum)

# Cálculo de la suma de los cuadrados
# Función para calcular la media y la desviación estándar de cada grupo
library(plyr)
function_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
rf_cv_rep_obs_mean <- function_summary(rf_cv_rep_gridsearch, varname="obs", 
                                       groupnames=c("Rep_Fold"))

rf_cv_rep_gridsearch_merge <- merge (rf_cv_rep_gridsearch, rf_cv_rep_obs_mean, by = "Rep_Fold")
rf_cv_rep_tss = tapply((rf_cv_rep_gridsearch_merge$obs - rf_cv_rep_gridsearch_merge$mean)^2,
                       rf_cv_rep_gridsearch_merge$Rep_Fold, sum)
rf_cv_rep_rsq  =  1 - (rf_cv_rep_rss/rf_cv_rep_tss)
rf_cv_rsq_mean = mean(rf_cv_rep_rsq)
rf_cv_rsq_sd = sd(rf_cv_rep_rsq)
cat('The mean R-square of the testing data is ', round(rf_cv_rsq_mean,3), '\n')
cat('The std R-square of the testing data is ', round(rf_cv_rsq_sd,3), '\n')

# Gráfica Predichos vs Observados
options(repr.plot.width=8, repr.plot.height=4)
rf_fit_gridsearch_predvsobs = as.data.frame(cbind(predicted = rf_fit_gridsearch_predicted, 
                                                  observed = y_trainData))
ggplot(rf_fit_gridsearch_predvsobs,aes(predicted, observed)) + 
  geom_point(color = "violetred4", alpha = 0.75) + 
  geom_smooth(method=lm, color="orangered", size=1)+
  ggtitle("Random Forest (RF): Predichos vs Observados") +
  xlab("SST predichos (º Brix)") + ylab("SST observados (º Brix)") + 
  theme(plot.title = element_text(color="black",size=16,hjust = 0.5, family = "serif",face="bold"),
        axis.text.y = element_text(face="bold", size=12, family = "serif"), 
        axis.text.x = element_text(face="bold", size=12,hjust=.5, family = "serif"),
        axis.title.x = element_text(size=14, family = "serif"), 
        axis.title.y = element_text(size=14, family = "serif"),
        axis.line = element_line(size = 0.8, colour = "black"))+
  geom_abline (slope=1, intercept=0, linetype = "dashed", color="blue")


############################################################################################################################################################
##------- Rango SWIR (1001-2500 nm) ---------------------------------------------------------------------------------------------------------------
############################################################################################################################################################

## 1.PROCESO DE MODELADO -----------------------------------------------------------------------------------------------------------------------------------

# Selección columnas de interés
data<-datos[,c(2,654:2153)]

# Establecer el conjunto de datos
trainData <- subset(data[,])

# Configuración de prueba
# Datos de entrenamiento
myRespName <- 'TSS' # Nombre de la variable dependiente a modelizar
y_trainData <- trainData[,1] # Posición de la variable dependiente a modelizar en el trainData
x_trainData <- trainData[,2:1501] # Conjunto de variables independientes del trainData

# Cálculo estadísticos variable dependiente (SST)
myRespName_mean <- mean(y_trainData)
myRespName_max <- max(y_trainData)
myRespName_min <- min(y_trainData)
myRespName_sd <- sd(y_trainData)

# Asegurar que los resultados sean repetibles
seed <- 7
set.seed(seed)
# Métrica R-cuadrado
metric <- "Rsquared"

# Algoritmo de prueba
# Crear trainControl
rf_fit_control <- trainControl(method="repeatedcv", 
                               number=10, 
                               repeats=10,
                               search="grid",
                               savePredictions=TRUE)

# Grid Search
rf_fit_tunegrid <- expand.grid(.mtry=c(1:15))
#mtry: número de variables muestreadas aleatoriamente como candidatas en cada división.
#ntree: número de árboles.

set.seed(seed)

# Parámetros del modelo
rf_fit_gridsearch <- train(TSS ~.,
                           data=trainData, na.action=na.exclude, method="rf", metric=metric, ntree=1500, 
                           tuneGrid=rf_fit_tunegrid, trControl=rf_fit_control, allowParallel = TRUE)

print(rf_fit_gridsearch)
plot(rf_fit_gridsearch)

# Importancia de las variables
varImp(rf_fit_gridsearch, scale = FALSE) #  Mean Decrease Gini
varImp(rf_fit_gridsearch)## Importancia en %


## 2. MODELO DE ENTRENAMIENTO (Bondad de ajuste de un modelo estadístico) --------------------------------------------------------------------------------             

# Modelo de entrenamiento: Predichos y Residuos
set.seed(seed)
rf_fit_gridsearch_predicted <- predict(rf_fit_gridsearch)
rf_fit_gridsearch_residuals <- residuals(rf_fit_gridsearch)

# Cálculo RPD
RPD(rf_fit_gridsearch_predicted,y_trainData) 
RPD <- formatC(RPD(rf_fit_gridsearch_predicted,y_trainData), digits = 2, format = "f")

# Cálculo RMSE
rf_fit_RMSE = sqrt(mean(rf_fit_gridsearch_residuals^2))
cat('The root mean square error of the training data is ', round(rf_fit_RMSE,3), '\n')

# R-cuadrado. Cálculo de la suma de los cuadrados
rf_y_fit_mean = mean(y_trainData)
rf_fit_tss =  sum((y_trainData - rf_y_fit_mean)^2)

# Cálculo suma residual de cuadrados
rf_fit_rss =  sum(rf_fit_gridsearch_residuals^2)

# Cálculo R-cuadrado
rf_fit_rsq  =  1 - (rf_fit_rss/rf_fit_tss)
cat('The R-square of the training data is ', round(rf_fit_rsq,3), '\n')

# Cálculo R-cuadrado del modelo
R2 <- formatC(mean(rf_fit_gridsearch$results$Rsquared), digits = 2, format = "f")

# Cálculo RMSE del modelo
RMSE <- formatC(mean(rf_fit_gridsearch$results$RMSE), digits = 2, format = "f")


## 3. VALIDACIÓN CRUZADA (CV) DEL CONJUNTO DE DATOS (Bondad de la prueba de un modelo estadístico)

# Evaluación del modelo
set.seed(seed)
rf_cv_all_gridsearch <- rf_fit_gridsearch$pred
rf_cv_all_gridsearch$Fold <- substr(rf_cv_all_gridsearch$Resample, start=1 , stop=6)
rf_cv_all_gridsearch$Rep <- substr(rf_cv_all_gridsearch$Resample, start=8 , stop=12)
rf_cv_all_gridsearch$Rep_Fold <- paste(rf_cv_all_gridsearch$Rep, rf_cv_all_gridsearch$Fold)
rf_cv_all_gridsearch$residuals <- (rf_cv_all_gridsearch$obs - rf_cv_all_gridsearch$pred)


# Predichos y Observados modelo CV 
# Configuración de los parámetros del mejor modelo
rf_cv_rep_gridsearch <- subset(rf_cv_all_gridsearch, mtry==5) # mtry== mtry obtenido en el mejor modelo

# Cálculo RMSE
rf_cv_rep_gridsearch$residuals2 = rf_cv_rep_gridsearch$residuals^2
rf_cv_rep_RMSE = sqrt(tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, mean))
rf_cv_RMSE_mean = mean(rf_cv_rep_RMSE)
rf_cv_RMSE_sd = sd(rf_cv_rep_RMSE)
cat('The mean of the root mean square error of the testing data is ', round(rf_cv_RMSE_mean,3), '\n')
cat('The sd of the root mean square error of the testing data is ', round(rf_cv_RMSE_sd,3), '\n')
rf_cv_RMSE <- as.data.frame(rf_cv_rep_RMSE)

# Cálculo R-cuadrado
# Cálculo suma residual de cuadrados
rf_cv_rep_rss =  tapply(rf_cv_rep_gridsearch$residuals2,rf_cv_rep_gridsearch$Rep_Fold, sum)

# Cálculo de la suma de los cuadrados
# Función para calcular la media y la desviación estándar de cada grupo
library(plyr)
function_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}
rf_cv_rep_obs_mean <- function_summary(rf_cv_rep_gridsearch, varname="obs", 
                                       groupnames=c("Rep_Fold"))

rf_cv_rep_gridsearch_merge <- merge (rf_cv_rep_gridsearch, rf_cv_rep_obs_mean, by = "Rep_Fold")
rf_cv_rep_tss = tapply((rf_cv_rep_gridsearch_merge$obs - rf_cv_rep_gridsearch_merge$mean)^2,
                       rf_cv_rep_gridsearch_merge$Rep_Fold, sum)
rf_cv_rep_rsq  =  1 - (rf_cv_rep_rss/rf_cv_rep_tss)
rf_cv_rsq_mean = mean(rf_cv_rep_rsq)
rf_cv_rsq_sd = sd(rf_cv_rep_rsq)
cat('The mean R-square of the testing data is ', round(rf_cv_rsq_mean,3), '\n')
cat('The std R-square of the testing data is ', round(rf_cv_rsq_sd,3), '\n')

# Gráfica Predichos vs Observados
options(repr.plot.width=8, repr.plot.height=4)
rf_fit_gridsearch_predvsobs = as.data.frame(cbind(predicted = rf_fit_gridsearch_predicted, 
                                                  observed = y_trainData))
ggplot(rf_fit_gridsearch_predvsobs,aes(predicted, observed)) + 
  geom_point(color = "violetred4", alpha = 0.75) + 
  geom_smooth(method=lm, color="orangered", size=1)+
  ggtitle("Random Forest (RF): Predichos vs Observados") +
  xlab("SST predichos (º Brix)") + ylab("SST observados (º Brix)") + 
  theme(plot.title = element_text(color="black",size=16,hjust = 0.5, family = "serif",face="bold"),
        axis.text.y = element_text(face="bold", size=12, family = "serif"), 
        axis.text.x = element_text(face="bold", size=12,hjust=.5, family = "serif"),
        axis.title.x = element_text(size=14, family = "serif"), 
        axis.title.y = element_text(size=14, family = "serif"),
        axis.line = element_line(size = 0.8, colour = "black"))+
  geom_abline (slope=1, intercept=0, linetype = "dashed", color="blue")

