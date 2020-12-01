## ----setup, include=TRUE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache=TRUE, fig.width=7, fig.height=5.5)
options(width=120)

## ----u107432----------------------------------------------------------------------------------------------------------
Sys.setenv(MKL_NUM_THREADS="8")
Sys.setenv(MKL_DYNAMIC="TRUE")
Sys.setenv(OMP_NUM_THREADS=4)
cores=4

## ---------------------------------------------------------------------------------------------------------------------
comienzo <- Sys.time()

## ---------------------------------------------------------------------------------------------------------------------
library(ipv)
data("venta.dep")
data("alquiler.dep")

## ---------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(mgcv)
library(tidyverse, verbose=FALSE)
library(sp)
library(rgdal)
library(rgeos)
library(spgwr)
library(zoo)

## ----PrepDatos--------------------------------------------------------------------------------------------------------
datos <- venta.dep %>% 
  mutate( logpm2 = log(IMP_CREACION / NUM_SUPERF),
          x = as.integer(FEC_CREACION - min(FEC_CREACION)) ) %>%
  filter( !is.na(DX_ETRS89) & !is.na(DY_ETRS89)) %>%
  as.data.frame()

## ----ModeloGlobal-----------------------------------------------------------------------------------------------------
mod1 <- gam( logpm2 ~   SUBREGION + ROOMNUMBER + HASLIFT + 
                HASGARDEN + HASSWMMINGPOOL + HASTERRACE + 
                IND_GARAJE  + s(x,bs="cr", k=16), 
              data=datos)

## ---- ResultadosModeloGlobal------------------------------------------------------------------------------------------
summary(mod1)
anova(mod1)

## ---- ZonasViaSUBREGION-----------------------------------------------------------------------------------------------
ind <- ConsInd(modelo=mod1, basedate="2007-12-31", 
        tit="Índice precio vivienda CAPV.\nEf. espacial: SUBREGION",
        fechas=as.Date(datos$FEC_CREACION),
        conf=0.95,
        plot=TRUE)

