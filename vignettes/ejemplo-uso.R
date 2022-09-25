## ----setup, include=TRUE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache=TRUE, fig.width=7, fig.height=5.5)
options(width=120)

## ----u107432----------------------------------------------------------------------------------------------------------
cores=2

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
datos <- alquiler.dep %>% 
  mutate( logpm2 = log(IMP_CREACION / NUM_SUPERF),
          x = as.integer(FEC_CREACION - min(FEC_CREACION)) ) %>%
  filter( !is.na(DX_ETRS89) & !is.na(DY_ETRS89)) %>%
  as.data.frame()

## ----ModeloGlobal-----------------------------------------------------------------------------------------------------
mod1 <- gam( logpm2 ~   SUBREGION + ROOMNUMBER +
                HASLIFT +  IND_GARAJE +
                NUM_SUPERF + s(x,bs="cr", k=10), 
              data=datos)

## ---- ResultadosModeloGlobal------------------------------------------------------------------------------------------
summary(mod1)
anova(mod1)

## ---- ZonasViaSUBREGION-----------------------------------------------------------------------------------------------
ind <- ConsInd(modelo=mod1, basedate="2007-12-31", 
        tit="Índice precio vivienda CAPV.\nEf. espacial: SUBREGION",
        fechas=as.Date(datos$FEC_CREACION), conf=0.95)

## ---------------------------------------------------------------------------------------------------------------------
class(ind)
head(ind)

## ---------------------------------------------------------------------------------------------------------------------
datos1 <- subset(datos, datos$FEC_CREACION < as.Date('2016-01-01') )
mod1.1 <- gam( logpm2 ~   SUBREGION + ROOMNUMBER +
                HASLIFT +  IND_GARAJE +
                NUM_SUPERF + s(x,bs="cr", k=10), 
              data=datos1)

## ---------------------------------------------------------------------------------------------------------------------
mod1.2 <- gam( logpm2 ~   SUBREGION + ROOMNUMBER +
                HASLIFT +  IND_GARAJE +
                NUM_SUPERF + s(x,bs="cr", k=10), 
              data=datos)

## ---- Mod1.1----------------------------------------------------------------------------------------------------------
ind1.1 <- ConsInd(modelo=mod1.1, basedate="2007-12-31", 
        tit="Índice precio vivienda CAPV.\nEf. espacial: SUBREGION",
        fechas=as.Date(datos1$FEC_CREACION), conf=0.95)
head(ind1.1)
head(ind)

## ---- Mod1.2----------------------------------------------------------------------------------------------------------
ind1.2 <- ConsInd(modelo=mod1.1, basedate="2007-12-31", 
        tit="Índice precio vivienda CAPV.\nEf. espacial: SUBREGION",
        congelado=ind1.1,
        fechas=as.Date(datos$FEC_CREACION), conf=0.95)
head(ind1.2)
tail(ind1.2)
tail(ind)

## ---------------------------------------------------------------------------------------------------------------------
ind.prov <- IndZonas(datos, base="2007-12-31", zonas="PROVINCIA", frm=formula(mod1))

## ----IndicesProvinciales----------------------------------------------------------------------------------------------
sel   <- (datos$FEC_CREACION > as.Date("2007-11-30"))   # Only after this date there is enough data
datos <- datos[sel,]
ind.prov <- IndZonas(datos, zonas="PROVINCIA", base="2007-12-31",
                     frm=mod1$formula)
ind.af   <- IndZonas(datos, zonas="AF", base="2007-12-31",
                     frm=update(mod1$formula, . ~ . - SUBREGION))

## ---------------------------------------------------------------------------------------------------------------------
mggplot(ind.prov, titulo="Indice precios alquiler")   
mggplot(ind.af, titulo="Indice precios alquiler")

## ---------------------------------------------------------------------------------------------------------------------
mggplot(ind.prov, tipo="multiple", titulo="Indice precios alquiler")   
mggplot(ind.af,   tipo="multiple", titulo="Indice precios alquiler")

## ----BackFitting-Global-----------------------------------------------------------------------------------------------
frm         <- formula(mod1)
frm.param   <- update(frm, ~ . - s(x,bs="cr",k=10))
smooth.term <- 's(x,bs="cr",k=10)'

# sel <- ((1:nrow(datos)) %% 5) == 0                        # One observation out of each 5
if (!file.exists("bfg.rda")) {
indice.bfg <- BackFitting(frm.param=frm.param,
            smooth.term=smooth.term,
            var.loc="SUBREGION",
            coords=c("DX_ETRS89","DY_ETRS89"),
            datos=datos[sel,],
            cores=cores,
            var.fecha='FEC_CREACION',
            baseday='2008-01-01')
save(indice.bfg, file="bfg.rda")
} else {
  load("bfg.rda")
}

## ----PlotBackfiffingGlobal, dependson="BackFitting-Global"------------------------------------------------------------
plot(indice.bfg)

## ---- CreacionPuntosIndLocal------------------------------------------------------------------------------------------
barrios <- c("Indautxu", "Abando - Albia", "Deusto")
cal.pts <- matrix(0,length(barrios),2)
b    <- venta.dep[venta.dep$BARRIO %in% barrios, ]
dimnames(cal.pts) <- list(barrios, NULL)
for (i in barrios) {
  cal.pts[i,] <- apply(b[b$BARRIO==i, c("DX_ETRS89","DY_ETRS89")], 2, FUN="median", na.rm=TRUE)
}
cal.pts

xy <- data.frame(ID=c("PlazaEnsanche","GeneralEguia"),
                 X=c(-2.931413,-2.946283),
                 Y=c(43.263538,43.259257))
coordinates(xy) <- c("X","Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84") 
res <- spTransform(xy, CRS("+proj=utm +zone=30 ellps=WGS84"))
res <- as.data.frame(res)
cal.pts <- as.matrix(res[,c("X","Y")])
rownames(cal.pts) <- res[,"ID"]

## ----BackFittingLocal-------------------------------------------------------------------------------------------------
# 
#  Eliminación de observaciones muy lejos de Bilbao
#
sel <- (datos$DX_ETRS89 > 504000) & (datos$DX_ETRS89 < 506000)
sel <- sel & (datos$DY_ETRS89 > 4788000) & (datos$DX_ETRS89 < 4791000)
sel <- sel & datos$FEC_CREACION > as.Date("2010-01-31")   # Sólo donde hay datos suficientes
lista.ind <- BackFittingLocal(frm.param=formula(logpm2 ~ ROOMNUMBER + HASLIFT +
                                                  IND_GARAJE + NUM_SUPERF),
                        smooth.term='s(x,bs="cr",k=10)',
                        cal.pts=cal.pts,
                        datos=datos[sel,],
                        coords=c("DX_ETRS89","DY_ETRS89"),
                        var.fecha='FEC_CREACION',
                        baseday="2010-12-31",
                        cores=cores,
                        bw=500,
                        bws=200,
                        tol=0.005)

## ----BackFittingLocalPlot---------------------------------------------------------------------------------------------
names(lista.ind[[2]]) <- rownames(lista.ind[[1]])
p <-mggplot(lista.ind[[2]], columnas=rownames(lista.ind[[1]]), tipo="multiple")
p <- p + ggtitle("Local indices",
              sub=paste("Base: ","31-12-2010"," = 100"))
# pdf(file="localindices.pdf")
print(p)
# scratch <- dev.off()

## ---------------------------------------------------------------------------------------------------------------------
median.houses <- datos[1000:1001,] 
median.houses$ROOMNUMBER <- 3
median.houses$HASLIFT <- factor("SI", levels=levels(datos$HASLIFT))
median.houses$NUM_SUPERF <- 85
median.houses

## ---------------------------------------------------------------------------------------------------------------------
median.houses <- median.houses[,-match(c("LAT","LNG"),colnames(median.houses))]
coordinates(median.houses) <- c("DX_ETRS89","DY_ETRS89")

## ----SlicedIndex------------------------------------------------------------------------------------------------------
preds <- SlicedIndex(
    cal.pts = median.houses,
    datos = datos,
    coords = c("DX_ETRS89","DY_ETRS89"),
    from = as.yearqtr("2015 Q1"),
    to = as.yearqtr("2019 Q1"),
    bw = 1000,
    frm = formula(logpm2 ~ ROOMNUMBER + NUM_SUPERF +
                    HASLIFT + IND_GARAJE),
    date = "FEC_CREACION",
    slice.by=as.yearqtr)
plot(exp(preds))

## ---------------------------------------------------------------------------------------------------------------------
final <- Sys.time()
final - comienzo

