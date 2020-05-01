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
ind <- ConsInd(modelo=mod1, base="2007-12-31", 
        tit="Índice precio vivienda CAPV.\nEf. espacial: SUBREGION",
        fechas=as.Date(datos$FEC_CREACION),
        plot=TRUE)

## ---------------------------------------------------------------------------------------------------------------------
class(ind)
head(ind)

## ---------------------------------------------------------------------------------------------------------------------
ind.prov <- IndZonas(datos, zonas="PROVINCIA", frm=formula(mod1))

## ---------------------------------------------------------------------------------------------------------------------
ind.prov <- IndZonas(datos, zonas="PROVINCIA", 
                     frm=formula( logpm2 ~ SUBREGION + ROOMNUMBER + 
                                    HASLIFT + HASGARDEN + HASSWMMINGPOOL +
                                    HASTERRACE +  IND_GARAJE  + 
                                    s(x,bs="cr", k=16) ))

## ---------------------------------------------------------------------------------------------------------------------
prov <- names(ind.prov)
prov
mggplot(ind.prov, columnas=prov)   

## ---------------------------------------------------------------------------------------------------------------------
mggplot(ind.prov, columnas=prov, tipo="multiple")   

## ---- CreacionSpatialDataFrame----------------------------------------------------------------------------------------
datos$HASLIFT <- as.factor(datos$HASLIFT)
spdatos <- datos[,c("logpm2","SUBREGION","ROOMNUMBER","HASLIFT",
                    "HASGARDEN","HASSWMMINGPOOL","HASTERRACE",
                    "IND_GARAJE","DX_ETRS89","DY_ETRS89","FEC_CREACION")]
spdatos <- spdatos[complete.cases(spdatos),]
coordinates(spdatos) <- ~ DX_ETRS89 + DY_ETRS89

## ---- BackFitting-Global----------------------------------------------------------------------------------------------
frm         <- formula(mod1)
frm.param   <- update(frm, ~ . - s(x,bs="cr",k=16))
smooth.term <- 's(x,bs="cr",k=16)'

sel <- ((1:nrow(datos)) %% 5) == 0                        # One observation out of each 5
sel <- sel & datos$FEC_CREACION > as.Date("2006-12-31")   # Only after this date there is enough data

indice <- BackFitting(frm.param=frm.param,
            smooth.term=smooth.term,
            var.loc="SUBREGION",
            coords=c("DX_ETRS89","DY_ETRS89"),
            datos=datos[sel,],
            cores=cores,
            var.fecha='FEC_CREACION',
            baseday='2008-01-01')

## ---------------------------------------------------------------------------------------------------------------------
plot(indice)

## ---- CreacionPuntosIndLocal------------------------------------------------------------------------------------------
# barrios <- c("Indautxu", "Abando - Albia", "Deusto")
# cal.pts <- matrix(0,length(barrios),2)
# b    <- venta.dep[venta.dep$BARRIO %in% barrios, ]
# dimnames(cal.pts) <- list(barrios, NULL)
# for (i in barrios) {
#   cal.pts[i,] <- apply(b[b$BARRIO==i, c("DX_ETRS89","DY_ETRS89")], 2, FUN="median", na.rm=TRUE)
# }
# cal.pts

xy <- data.frame(ID=c("PlazaEnsanche","GeneralEguia"),
                 X=c(-2.931413,-2.946283),
                 Y=c(43.263538,43.259257))
coordinates(xy) <- c("X","Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
res <- spTransform(xy, CRS("+proj=utm +zone=30 ellps=WGS84"))
res <- as.data.frame(res)
cal.pts <- as.matrix(res[,c("X","Y")])
rownames(cal.pts) <- res[,"ID"]

## ---- BackFittingLocal------------------------------------------------------------------------------------------------
# 
#  Eliminación de observaciones muy lejos de Bilbao
#
sel <- (datos$DX_ETRS89 > 504000) & (datos$DX_ETRS89 < 506000)
sel <- sel & (datos$DY_ETRS89 > 4788000) & (datos$DX_ETRS89 < 4791000)
sel <- sel & datos$FEC_CREACION > as.Date("2010-12-31")   # Sólo donde hay datos suficientes

lista.ind <- BackFittingLocal(frm.param=frm.param,
                        smooth.term='s(x,bs="cr",k=9)',
                        cal.pts=cal.pts,
                        datos=datos[sel,],
                        coords=c("DX_ETRS89","DY_ETRS89"),
                        var.fecha='FEC_CREACION',
                        var.loc='SUBREGION',
                        baseday="2010-12-31",
                        cores=cores,
                        bw=500,
                        bws=200,
                        tol=0.005,
                        plotind=FALSE)

## ---------------------------------------------------------------------------------------------------------------------
names(lista.ind[[2]]) <- rownames(lista.ind[[1]])
p <-mggplot(lista.ind[[2]], columnas=rownames(lista.ind[[1]]), tipo="multiple")
p <- p + ggtitle("Local indices",
              sub=paste("Base: ","31-12-2010"," = 100"))
# pdf(file="localindices.pdf")
print(p)
# scratch <- dev.off()

## ---------------------------------------------------------------------------------------------------------------------
median.houses <- spdatos@data[1:2,]
median.houses$ROOMNUMBER <- 3
median.houses$HASLIFT <- factor("SI", levels=levels(spdatos@data$HASLIFT))
median.houses

## ---------------------------------------------------------------------------------------------------------------------
coordinates(median.houses) <- cal.pts

## ---------------------------------------------------------------------------------------------------------------------
preds <- SlicedIndex(
    cal.pts = median.houses,
    spdatos = spdatos,
    from = as.yearqtr("2016 Q1"),
    to = as.yearqtr("2019 Q1"),
    bw = 1000,
    frm = formula(logpm2 ~ ROOMNUMBER + 
                    HASLIFT + IND_GARAJE),
    date = "FEC_CREACION",
    slice.by=as.yearqtr)
plot(exp(preds))

## ---------------------------------------------------------------------------------------------------------------------
final <- Sys.time()
final - comienzo

