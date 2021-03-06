#
#  Class definitons for objects of type Indice and derived
#
setOldClass("zoo")
setClassUnion("characterOrNULL",c("character", "NULL"))
setClass("Indice",
         slots=list(ICV="zoo", basedate="ANY", basevalue="numeric"))
setClass("IndiceCB",
         slots=list(lcb="numeric", ucb="numeric", conf="numeric"),
         contains="Indice")
setClass("IndiceLocal",
         slots=list(Indice="Indice",coords="numeric"),
         contains="Indice")
setMethod("initialize",
          signature="Indice",
          definition= function(.Object, ICV, basedate=index(ICV)[1], basevalue=100) {
            .Object@ICV = ICV
            .Object@basedate  = basedate
            .Object@basevalue = basevalue
            return(.Object)
          }
)
setMethod("head",
          signature=c(x="Indice"),
          definition= function(x, ...)  {
            head(x@ICV, ... )
          }
)
setMethod("head",
          signature=c(x="IndiceCB"),
          definition= function(x, ...)  {
            tmp <- zoo(coredata(x@ICV,x@lcb,x@ucb), order.by=index(x@ICV))
            head(tmp, ... )
          }
)
setMethod("initialize",
          signature="IndiceCB",
          definition= function(.Object, Indice, lcb, ucb, conf) {
            .Object@ICV = Indice@ICV
            .Object@basedate = Indice@basedate
            .Object@basevalue = Indice@basevalue
            .Object@lcb = lcb
            .Object@ucb = ucb
            .Object@conf = conf
            return(.Object)
          }
)
setMethod("print",
          c(x="Indice"),
          function(x) {
          print(x@ICV)
            }
          )
setMethod("print",
          c(x="IndiceCB"),
          function(x) {
            tmp <- zoo(cbind(coredata(x@ICV), x@lcb, x@ucb),
                       order.by=index(x@ICV))
            colnames(tmp) <- c("Index", paste(c("lcb ","ucb "),100*x@conf,"%",sep=""))
            print(tmp)
          }
)

setMethod("plot",
          c(x="Indice"),
          function(x,..., tit="characterOrNULL") {
          p <- ggplot(x@ICV, aes(x = index(x@ICV), y = as.numeric(x@ICV)))  +
            geom_line() +
            xlab("Fecha") + ylab("Índice")
          if(!is(x,"IndiceCB")) {
            #
            #  If it is and IndiceCB, its own method will take care of
            #  the title and subtitle, otherwise:
            #
            p <- p + ggtitle(tit,sub=paste("Base: ",x@basedate," = ",x@basevalue))
          }
          return(p)
          }
          )
setMethod("plot",
          c(x="IndiceCB"),
          function(x, ..., tit="characterOrNULL") {
            p <- callNextMethod() +
                 geom_line(aes(x = index(x@ICV), y = x@lcb),
                                   col="blue", linetype="dotted") +
                 geom_line(aes(x = index(x@ICV), y = x@ucb),
                                   col="blue", linetype="dotted") +
                 ggtitle(tit, sub=paste("Base: ",x@basedate," = ",x@basevalue,". Int. conf. ",
                                   100*x@conf,"%"))
            return(p)
          })
setMethod("show",
          c(object="Indice"),
          function(object) {
            print(object)
          }
)

#' Semi-parametric estimation usign backfiffing (global version)
#'
#' This function fits a parametric model whose parametric part is
#' a geographically weigthed regression (GWR) while the non-parametric
#' is a cubic spline global, to all observations
#'
#' @param frm.param Formula for the parametric part  (GWR) of the model
#' @param smooth.term Non-parametric part of the model
#' @param var.loc Variable que en el modelo base se emplea como proxy de la ubicación. Se elimina en el modelo que realiza GWR
#' @param datos Dataframe espacial conteniendo los datos.
#' @param indice0 Valores iniciales del índice. Optativo: si está ausente, se calcula.
#' @param coords Coordenadas espaciales de las observaciones. Pueden ser geográficas (latitud y longitud) en cuyo caso se calcula distancia geodésica, o coordenadas de una proyección plana, en cuyo caso de emplea distancia euclídea ordinaria.
#' @param var.fecha Variable que recoge la fecha de cada observación. Debe formar parte de \code{datos.}
#' @param baseday Día en que se fija la base 100 del índice.
#' @param bw "Bandwidth" o radio del área espacial que recibe ponderación sustancial en la estimación de parámetros para cada punto de calibración.
#' @param cores Número de cores que se desea emplear. Si no se indica nada, todos los de la máquina.
#' @param tol Máxima discrepancia relativa entre los valores de dos splines en iteraciones sucesivas para continuar la iteración.
#' @param plotind Variable lógica: TRUE si se desea la impresión del índice resultante.
#' @param gwrmod Variable lógica, TRUE si queremos los resultados del ajuste espacial
#'
#' @return Un índice de precios, con base en el día especificado, en formato de serie temporal \code{zoo.}
#' @export
#'
#' @examples
#'
BackFitting <- function(frm.param,
                       smooth.term='s(x,bs="cr",k=24)',
                       var.loc=NULL,
                       datos,
                       indice0,
                       coords,
                       var.fecha,
                       baseday="2012-12-31",
                       bw=5000,
                       cores,
                       tol=0.001,
                       plotind=FALSE,
                       gwrmod=FALSE) {
    require(mgcv)
    require(spgwr)
    require(parallel)
    if (missing(cores))
      cl <- makeCluster(detectCores())
    else
      cl <- makeCluster(cores)
    resp  <- as.character(frm.param[[2]])
    datos <- as.data.frame(datos)   #  Tibbles dan problemas
    datos <- cbind(datos,y.orig=datos[,resp])
    #
    #  Modificadores de la fórmula y fragmentos: añadimos a
    #  la parte paramétrica el(los  ) termino(s) suave(s) para
    #  estimar el GAM inicial. Eliminamos las variables de
    #  significado espacial de la fórmula a emplear en la
    #  estimación GWR.
    #
    frm.global <- update(frm.param,
                         eval(parse(text=paste("~ . + ",
                                               smooth.term))))
    if (!is.null(var.loc)) {
      frm.gwr    <- update(frm.param,
                           eval(parse(text=paste("~ . - ",
                                                 var.loc))))
    } else {
      frm.gwr <- frm.param
    }
    frm.tot    <- update(frm.param,
                         eval(parse(text=paste("~ . + y.orig + ",
                                               paste(c(coords,var.fecha),
                                                     collapse="+")))))
    frm.suav <- formula(eval(parse(text=paste(resp," ~  ", smooth.term))))
    #
    #  Creación de una SpatialDataFrame y selección de casos
    #  completos para el backfitting. (gwr no admite NA's)
    #
    spdatos   <- get_all_vars(frm.tot, datos)
    completos <- complete.cases(spdatos)
    spdatos   <- spdatos[completos,]
    coordinates(spdatos) <- eval(parse(text=paste("~ ",
                                                  paste(coords,
                                                        collapse="+"))
    )
    )
    #
    x.base <- min(spdatos@data[,var.fecha])
    spdatos@data <- cbind(spdatos@data,
                          x=as.numeric(spdatos@data[,var.fecha, drop=TRUE] -
                                         x.base
                          )
    )
    x <- spdatos@data$x
    #
    #  Estimación modelo global semi-paramétrico para solución
    #  inicial, si es que no se ha pasado 'indice0' entre los
    #  argumentos.
    #
    if (missing(indice0)) {
      mod.gam <- gam(formula=frm.global, data=spdatos@data)
      lastind <- newind <-  ConsInd(modelo = mod.gam,
                         fechas = spdatos@data[,var.fecha, drop=TRUE],
                         basedate = baseday)@ICV
    }
    #
    #  Comienza ahora la alternancia entre estimaciones de la
    #  parte paramétrica y no paramétrica del modelo, hasta
    #  (esperamos) convergencia
    #
    lastind <- 0 * newind
    iter <-0
    #
    #  Mientras la máxima discrepancia entre dos índices
    #  sucesivos sea pequeña, continúa
    #
    while( max(abs(newind-lastind)) > tol * max(abs(newind)) ) {
      iter <- iter + 1
      lastind <- newind
      #
      #  Deflactamos los datos con el índice provisional
      #
      deflactor <- log(lastind / 100)
      tmp       <- match(spdatos@data[,var.fecha,drop=TRUE],
                         index(deflactor))
      spdatos@data[,resp] <-
        spdatos@data[,"y.orig"] - coredata(deflactor[tmp])
      #
      #  Ajustamos un modelo espacial a los datos deflactados;
      #  los valores ajustados, en mod.gwr$SDF@data$pred
      #
      mod.gwr      <- gwr(frm.gwr, data=spdatos, bandwidth=bw,
                          fit.points=spdatos, hatmatrix=FALSE,
                          gweight=gwr.Gauss,
                          predictions=TRUE,
                          cl=cl)
      #
      #  Los residuos del modelo espacial, a datos[,resp]
      #
      spdatos@data[,resp] <- spdatos@data[,"y.orig"] - mod.gwr$SDF@data$pred
      #
      #  Ajustamos un spline a datos[,resp]
      #
      mod.gam <- gam(frm.suav, data=spdatos@data)
      newind  <- ConsInd(modelo = mod.gam,
                         fechas = spdatos@data[,var.fecha, drop=TRUE],
                         basedate = baseday)@ICV
      cat("Backfitting iteración ",iter,"\n")
    }
    stopCluster(cl)
    #
    #  Finalizado el ajuste, representamos el índice si plotind=TRUE
    #
    ConsInd(modelo = mod.gam,
            fechas = spdatos@data[,var.fecha, drop=TRUE],
            basedate = baseday)
    if (gwrmod) {
      return(list(ind=newind,gwrmod=mod.gwr ))
    } else {
    return(newind)
      }
}

#' Estimación semiparamétrica mediante 'backfitting' local.
#'
#' Realiza la estimación de un modelo semi-paramétrico, cuya parte
#' paramétrica esta formada por una regresión espacialmente ponderada y
#' cuya parte no paramétrica está formada por un spline cúbico para cada
#' ubicación. Esto permite el ajuste de tendencias locales
#'
#' @param frm.param Fórmula de la parte paramétrica (GWR) del modelo.
#' @param smooth.term Término no paramétrico del modelo.
#' @param cal.pts Puntos del espacio para los cuales se calcula el índice
#' @param datos Dataframe espacial conteniendo los datos.
#' @param indice0 Valores iniciales del índice. Optativo: si está ausente, se calcula.
#' @param coords Coordenadas espaciales de las observaciones. Pueden ser geográficas (latitud y longitud) en cuyo caso se calcula distancia geodésica, o coordenadas de una proyección plana, en cuyo caso de emplea distancia euclídea ordinaria.
#' @param var.fecha Variable que recoge la fecha de cada observación. Debe formar parte de \code{datos.}
#' @param var.loc Variable que en el modelo base se emplea como proxy de la ubicación. Se elimina en el modelo que realiza GWR
#' @param baseday Día en que se fija la base 100 del índice.
#' @param bw "Bandwidth" o radio del área espacial que recibe ponderación sustancial en la estimación de parámetros de la GWR para cada punto de calibración.
#' @param bws "Bandwidth" o radio del área espacial que recibe ponderación sustancial en la estimación del spline
#' @param cores Número de cores que se desea emplear. Si no se indica nada, todos los de la máquina.
#' @param tol Máxima discrepancia relativa entre los valores de dos splines en iteraciones sucesivas para continuar la iteración.
#' @param plotind Variable lógica: TRUE si se desea la impresión del índice resultante.
#'
#' @return Una lista de índices de precios, con base en el día especificado, en formato de serie temporal \code{zoo;} uníndice para cada punto especificado en \code{cal.pts}.
#' @export
#'
#' @examples
#'
BackFittingLocal <- function(frm.param,
                             smooth.term='s(x,bs="cr",k=8)',
                             cal.pts=NULL,
                             datos,
                             indice0,
                             coords,
                             var.fecha,
                             var.loc=NULL,
                             baseday="2012-12-31",
                             bw=5000,
                             bws=NULL,
                             cores,
                             tol=0.001) {
  require(mgcv)
  require(spgwr)
  require(parallel)
  if (is.null(cal.pts) & is.null(bws)) {
    global.fit <- TRUE
  } else {
    if (xor(is.null(cal.pts), is.null(bws))) {
      stop("Both cal.pts and is.null must be NULL or non-NULL.")
    } else {
    global.fit <- FALSE
    }
  }
  if (missing(cores)) {
    cl <- makeCluster(detectCores())
  } else {
    cl <- makeCluster(cores)
  }
  resp  <- as.character(frm.param[[2]])
  datos <- as.data.frame(datos)   #  Tibbles dan problemas
  datos <- cbind(datos,y.orig=datos[,resp])
  #
  #  Modificadores de la fórmula y fragmentos: añadimos a
  #  la parte paramétrica el(los) termino(s) suave(s) para
  #  estimar el GAM inicial. Eliminamos las variables de
  #  significado espacial de la fórmula a emplear en la
  #  estimación GWR.
  #

  if (!is.null(var.loc)) {
    frm.param    <- update(frm.param,
                         eval(parse(text=paste("~ . - ",
                                               var.loc))))
  } else {
    frm.param <- frm.param
  }
  frm.global <- update(frm.param,
                       eval(parse(text=paste("~ . + ",
                                             smooth.term))))
  frm.tot    <- update(frm.param,
                       eval(parse(text=paste("~ . + y.orig + ",
                                             paste(c(coords,var.fecha),
                                                   collapse="+")))))
  frm.suav <- formula(eval(parse(text=paste(resp," ~  ", smooth.term))))

  #
  #  Creación de una SpatialDataFrame y selección de casos
  #  completos para el backfitting. (gwr no admite NA's)
  #
  spdatos   <- get_all_vars(frm.tot, datos)
  completos <- complete.cases(spdatos)
  spdatos   <- spdatos[completos, ]
  coordinates(spdatos) <- eval(parse(text=paste("~ ",
                                                paste(coords,
                                                      collapse="+"))))
  #
  x.base <- min(spdatos@data[,var.fecha])
  spdatos@data <- cbind(spdatos@data,
                        x=as.numeric(spdatos@data[,var.fecha, drop=TRUE] -
                                       x.base,
                                     w=rep(1,nrow(spdatos@data))
                        )
  )
  x <- spdatos@data$x
  if (!global.fit) {
    npts    <- nrow(cal.pts)             # Número de puntos
  } else {
    npts <- 1
  }
  indices <- vector("list", npts)      # Lista de índices
  #
  #  La misma iteración que para un índice global único se lleva
  #  a cabo ahora para cada uno de los puntos en cal.pts
  #
  spdatos.completo <- spdatos
  for (p in 1:npts) {
    if (!global.fit) {
      cat("Iniciado índice",p,"\n")
      punto   <- cal.pts[p,]                       # Punto de cálculo
      d2 <- apply(sweep(coordinates(spdatos.completo),
                        2, punto)^2,1,sum)
      # Distancias al cuadrado al punto de cálculo
      spdatos.completo@data$w <- w <- spgwr::gwr.Gauss(d2,bws)  # Pesos
      spdatos <- spdatos.completo[w > 0.01,]
      cat("Num. obs. empleadas =",nrow(spdatos),"\n")
      if (nrow(spdatos) < 50) {
        cat("Insufficient data for point ",punto,"\n")
        break
      }
      mod.gam <- gam(formula=frm.global, data=spdatos@data, weights=w)
      lastind <- newind <- ConsInd(modelo = mod.gam,
                         fechas = spdatos@data[,var.fecha, drop=TRUE],
                         basedate = baseday)@ICV
      #
      #  Comienza ahora la alternancia entre estimaciones de la
      #  parte paramétrica y no paramétrica del modelo, hasta
      #  (esperamos) convergencia
      #
      lastind <- 0 * newind
      iter <-0
      #
      #  Mientras la máxima discrepancia entre dos índices
      #  sucesivos sea pequeña, continúa
      #
      while( max(abs(newind-lastind)) > tol * max(abs(newind)) ) {
        iter <- iter + 1
        lastind <- newind
        #
        #  Deflactamos los datos con el índice provisional
        #
        deflactor <- log(lastind / 100)
        tmp       <- match(spdatos@data[,var.fecha,drop=TRUE],
                           index(deflactor))
        spdatos@data[,resp] <-
          spdatos@data[,"y.orig"] - coredata(deflactor[tmp])
        #
        #  Ajustamos un modelo espacial a los datos deflactados;
        #  los valores ajustados, en mod.gwr$SDF@data$pred
        #
        mod.gwr      <- gwr(frm.param, data=spdatos, bandwidth=bw,
                            hatmatrix=FALSE,
                            gweight=gwr.Gauss,
                            predictions=TRUE,
                            cl=cl)
        #
        #  Los residuos del modelo espacial, a datos[,resp]
        #
        spdatos@data[,resp] <- spdatos@data[,"y.orig"] - mod.gwr$SDF@data$pred
        #
        #  Ajustamos un spline a datos[,resp]
        #
        mod.gam <- gam(formula=frm.suav, data=spdatos@data, weights=w)
        tmp     <- ConsInd(modelo = mod.gam,
                           fechas = spdatos@data[,var.fecha, drop=TRUE],
                           basedate = baseday)
        newind  <- tmp@ICV
        cat("Backfitting iteración ",iter,"\n")
      }
      indices[[p]] <- tmp
    }
  }
    stopCluster(cl)
    if (global.fit)
      indices=indices[[1]]
    else
      indices <- list(cal.pts, indices)
    return(indices)
}


#' This function creates cloropleth maps
#'
#' Given a spatail dataframe, creates a cloropleth map from the valaues of one of the variables in the @data slot
#'
#' @param sp Spatial dataframe, in whose @data slot variable \code{var.sp} defines the aggregation zones.
#' @param var.sp Variable whose values specify the spatail aggregation.
#' @param df Ordinary dataframe with variable \code{var.df} matching \code{var.sp.}
#' @param var.df Variable matching \code{var.sp}; both need to have the same levels.
#' @param var Variable to code with color. Must be a column of \code{df.}
#' @param sub.var Variable in \code{sp}, optionally used to subset \code{sp}.
#' @param sub.set Vector of values of \code{sub.var} defining a subser of rows of \code{sp}.
#' @param legend.title Legend heading
#' @param logcolor FALSE if linear color scale, TRUE if logarithmic.
#' @param graf.title Heading of the graph.
#' @param ... Optional additional parameters to pass to the underlying \code{ggplot} function.
#'
#' @return A map, of ggplot2 format, printed by default.
#' @export
#'
#' @examples
#'
cloromap <- function(sp, var.sp,
                     df, var.df,
                     var,
                     sub.var=NULL, sub.set=NULL,
                     legend.title=NULL,
                     logcolor=FALSE,
                     graf.title=NULL, ...) {
  pm <- pols(base=sp, var.agreg=var.sp, sub.var=sub.var, sub.set=sub.set)
  pm <- cbind(pm, x=df[match(pm$id,df[,var.df]), var])
  mapa <- ggplot(data=pm) +
    geom_polygon(aes(x=long, y=lat, group=group, fill=x),
                 color="black")
  if (logcolor) {
    mapa <- mapa + scale_fill_gradient2(low="Blue",
                                        high="Red", trans="log",
                                        name=legend.title, ...)
  } else {
    mapa <- mapa + scale_fill_gradient2(low="Blue",high="Red",
                                        name=legend.title, ...)
  }
  proy <- proj4string(sp)
  if(length(grep("proj=utm", proy)) > 0) {
    zona <- gsub("^(.*)zone=(\\d*).*","\\2", proy, perl=TRUE)
    mapa <- mapa + coord_equal() +
      xlab(paste("UTM x (Zona ",zona,")",sep="")) + ylab("UTM y")
  } else if(length(grep("proj=longlat", proy)) > 0)  {
    mapa <- mapa + coord_map("stereographic") + xlab("Longitud") +
      ylab("Latitud")
  } else  stop("Proyección desconocida.")
  mapa <- mapa + ggtitle(graf.title)
  print(mapa)
}


#' CompFechas
#'
#' Given an index, defined over a set of incomplete dates, interpolates a new index over the full set of dates.
#'
#' @param indice An index, such as returned by \code{ConsInd.}, of class \code{Indice} or \code{IndiceCB}
#' @param fechas Vector of dates whose extremes define the new date range. Can be omitted, in which case the extremes are taken from the index of \code{indice}
#'
#' @return A dily index over the period from \code{min(fechas)} to \code{max(fechas).}
#' @export
#'
#' @examples
#'
CompFechas <- function(indice,fechas=NULL) {
  default.fechas <- index(indice@ICV)
  if (is.null(fechas))
    fechas <- default.fechas
  fechas <- unique(fechas)
  scratch     <- zoo(0,order.by=seq.Date(from=min(fechas),to=max(fechas)+1,by="day"))
  indice@ICV  <- na.approx(merge(indice@ICV,scratch)[,1], rule=2)
  if (class(indice)=="IndiceCB") {
     indice@lcb <- as.numeric(na.approx(merge(zoo(indice@lcb,fechas),scratch)[,1], rule=2))
     indice@ucb <- as.numeric(na.approx(merge(zoo(indice@ucb,fechas),scratch)[,1], rule=2))
   }
  return(indice)
}

#' ConsInd
#'
#' Dado un modelo semi-paramétrico ya estimado, que incluye un término spline recogiendo la tendencia,
#' extrae dicho término, lo completa para los días en que no haya observación, y lo devuelve como una serie temporal, representándolo gráficamente
#' si se desea. Se supone que la variable respuesta el log(Precio) o log(Precio/m2); el índice devuelto lo es para la variable Precio o Precio/m2, respectivamente, y datos fechados diariamente.
#'
#' @param modelo Modelo semiparamétrico estimado por \code{gam.}
#' @param base Fecha (día) tomada como base 100 del índice.
#' @param conf Confianza del intervalo; si no se especifica, 95\%.
#' @param fechas Vector de fechas de las observaciones; habitualmente se toma de una columna de la dataframe empleada por \code{gams} para estimar \code{modelo}
#' @param tit Encabezamiento del gráfico producido con el índice.
#' @param ylabel Rótulo eje de ordenadas.
#' @return Un objeto de clase 'Ìndice' o 'IndiceCB'.
#' @export
#'
#' @examples
#'
ConsInd <- function(modelo=NULL,
                    conf=NULL,
                    fechas=NULL,
                    tit=tit,
                    ylabel="Índice",
                    basedate,
                    basevalue=100) {
  if (missing(fechas))
    stop("Argument 'fechas' must be provided.")
  plotdata <- plot(modelo, select=0)[[1]]
  startdate <- min(fechas)
  x <- as.Date(plotdata$x + as.Date(startdate))
  y    <- c(plotdata$fit)
  se   <- plotdata$se
  difs <- as.numeric(abs(x-as.Date(basedate)))
  orig <- seq_along(x)[difs==min(difs)]   # la fecha más próxima a la base
  ind <- new("Indice",
             ICV=zoo(100*exp(y - as.numeric(y[orig])),
                      order.by=x),
             basevalue=100,
             basedate=basedate)
  if (!is.null(conf)) {
    lcb <- 100*exp(y - qnorm(conf)*se - as.numeric(y[orig]))
    ucb <- 100*exp(y + qnorm(conf)*se - as.numeric(y[orig]))
    ind <- new("IndiceCB", ind, lcb=lcb, ucb=ucb, conf=conf)
  }
  #
  #  Antes de devolver el índice, vamos a completarlo para TODAS las fechas posibles entre la
  #  primera y la última. Eso garantiza que nos permitirá deflactar cualquier magnitud
  #  definida en el mismo intervalo de fechas, pero no necesariamente sobre las mismas
  #  abscisas temporales. Empleamos la función previamente
  #  definida 'CompFechas'.
  #
  ind <- CompFechas(ind, fechas)
  return(ind)
}

#' Fusion
#'
#' @param X  Second dataframe to merge
#' @param Y  Second dataframe to merge
#' @param varX Variable acting as key in X
#' @param locX Column in \code{canon} against which to match values of column \code{varX} in \code{X}
#' @param varY Variable acting as key in Y
#' @param locY Column in \code{canon} against which to match values of column \code{varY} in \code{Y}
#' @param canon Correspondence between varX, varY and canonical value of the key (defaults to column 'Canon').
#' @param locCanon Column in \code{canon} holding the canonical names. Defaults to column named 'Canon' if one exists and the argument is not set on call.
#' @param transX Transformation to be applied to varX (if any)
#' @param transY Transformation to be applied to varY (if any)
#'
#' @return Data frame or merge of two data frames, with location names replaced by canonical names.
#' @export
#'
#' @examples
#' X     <- data.frame(a=c("Bilbao","San Sebastian","Vitoria","Teruel"),
#'                     b=c("A","B","B","D"), d=c(TRUE,TRUE,FALSE,TRUE))
#' X
#' Y     <- data.frame(e=c("Bilbao","Donostia","Vitoria/Gazteiz","Soria"),
#'                     g=c("F","G","H","J"), h=c(TRUE,FALSE,TRUE,TRUE))
#' Y
#' canon <- data.frame(X=c("Bilbao","San Sebastian","Vitoria"),
#'                     Y=c("Bilbao","Donostia","Vitoria/Gazteiz"),
#'                     Canon=c("Bilbao","Donostia","Vitoria"))
#' canon
#'
#' Fusion(X, varX="a", locX=1, locCanon=3, canon=canon)
#'
#' Fusion(X=Y, varX="e", locX=2, locCanon=3, canon=canon)
#'
#' Fusion(X,Y,varX="a",varY="e",locX=1, locY=2, canon=canon)
#'
#' canon <- data.frame(X=c("Bilbao","San Sebastian","Vitoria","Teruel","Soria"),
#'                     Y=c("Bilbao","Donostia","Vitoria/Gazteiz", "Teruel","Soria"),
#'                     Canon=c("Bilbao","Donostia","Vitoria","Teruel","Soria"))
#'
#' Fusion(X, varX="a", locX=1, locCanon=3, canon=canon)
#'
#' Fusion(X=Y, varX="e", locX=2, locCanon=3, canon=canon)
#'
#' Fusion(X,Y,varX="a",varY="e",locX=1, locY=2, canon=canon)
#'
#' canon <- data.frame(X=c("Bilbao","San Sebastian","Vitoria","Teruel","Soria"),
#'                     Y=c("Bilbao","Donostia","Vitoria/Gazteiz", "Teruel","Soria"),
#'                     Canon=c("CAPV","CAPV","CAPV","NoCAPV","NoCAPV"))
#'
#' Fusion(X, varX="a", locX=1, locCanon=3, canon=canon)
#'
#' Fusion(X=Y, varX="e", locX=2, locCanon=3, canon=canon)
#'
#' # This hardly makes sense:
#'
#' Fusion(X, Y, varX="a", varY="e", locX=1, locY=2, locCanon=3, canon=canon)
#'
Fusion <- function(X,
                   Y=NULL,
                   varX, locX=1,
                   varY, locY=2,
                   canon, locCanon=NULL,
                   transX = NULL,
                   transY = NULL)
{
  #
  #  Use default, if canonical column has not been specified
  #
  if (is.null(locCanon))
    locCanon=match("Canon",colnames(canon))
  #
  #  Transliterate as needed
  #
  if (!is.null(transX))
    X[, varX] <-
      iconv(x = X[, varX], from = transX, to = "utf8")
  if (!is.null(transY))
    Y[, varY] <-
      iconv(x = Y[, varY], from = transY, to = "utf8")
  #
  #  Match key varX and replace by canonical values
  #
  X <- as.data.frame(X)
  tmp <- match(as.character(X[, varX]), as.character(canon[, locX]))
  if (any(is.na(tmp))) {
    cat("\nThe following values of varX could not be matched:\n")
    cat(as.character(X[is.na(tmp),varX]),"\n")
  }
  tX <- ifelse(!is.na(tmp), canon[tmp, locCanon], NA)
  if (is.factor(canon[tmp, locCanon]))
    X[, varX] <- factor(levels(canon[tmp, locCanon])[tX])
  else
    X[,varX] <- tX
  #
  #  Do likewise for key varY, if second argument is given
  #
  if (!is.null(Y)) {
    tmp <- match(Y[, varY], canon[, locY])
    if (any(is.na(tmp))) {
      cat("\nThe following values of varY could not be matched:\n")
      cat(as.character(Y[is.na(tmp),varY],"\n"))
    }
    tY <- ifelse(!is.na(tmp), canon[tmp, locCanon], NA)
    if (is.factor(canon[tmp, locCanon]))
      Y[, varY] <- factor(levels(canon[tmp, locCanon])[tY])
    else
      Y[,varY] <- tY
  }
  #
  #  If two data frames 'x' and 'y' were passed, merge them to return,
  #  else return the first with the key replaced by canonical names
  #
  if (!is.null(Y))
    tmp <- merge(x=X, y=Y, by.x=varX, by.y=varY, all=TRUE)
  else
    tmp <- X
  #
  return(tmp[!is.na(tmp[,varX]),])
}

#' IndZonas
#'
#' Estimation of global indices by areas; merely a wrapper around  \code{ConsInd}.
#' @param datos Dataframe of observations
#' @param zonas Name of the column of \code{datos} defining the areeas for which we want indices.
#' @param frm Formula specifying the base model, in the notation expected by  \code{gam.} The formula is common to all indices computed byu the function.
#' @param base Base date (index=100). It is a day in format yyyy-mm-dd.
#' @param plot \code{FALSE} (default) if no graphs are desired, \code{TRUE} otherwise.
#'
#' @return A list with the computed indices, in the format in which they are returned by \code{ConsInd}.
#' @export
#'
#' @examples
#'
IndZonas <- function(datos,
                     zonas,
                     frm = NULL,
                     base=NULL,
                     plot=FALSE) {
  if (is.null(base))
    stop("Base date must be set in arbument 'base'.`")
  x <- vector("list",0)
  k <- match(zonas, colnames(datos))
  for (i in unique(datos[,k])) {
    sel <- datos[,k] == i
    mod <- gam( frm, data=datos[sel,])
    indice0 <- ConsInd(modelo = mod,
                       fechas = datos[sel,]$FEC_CREACION,
                       basedate = base,
                       conf= 0.95,
                       tit = paste("Índice precio vivienda. Ámbito: ",
                                   i, sep="")
    )
    x[[i]] <- indice0
  }
  return(x)
}

lonlat2UTM = function(lonlat,lat=NULL) {
  if (!is.null(lat)) {       # In case lon and lat given as different numbers
    lonlat <- c(lonlat,lat)
  }
  utm = (floor((lonlat[1]+180) / 6) %% 60) +1
  if(lonlat[2] > 0) {
    utm + 32600
  } else {
    utm + 32700
  }
}



#' mggplot
#'
#' Representar múltiples índices (típicamente devueltos por la función \code{IndZonas}) en sendos paneles o superpuestos en un único panel.
#'
#' @param x Lista de índices tal como los proporciona \code{ConsInd.} Puede ser una dataframe  (caso especial de lista).
#' @param columnas Nombres de los índices a representar; por defecto, los nombres de la lista \code{x} (o columnas de la dataframe \code{x}.)
#' @param tipo Tipo de plot que se desea: "panel" si se desea cada serie en un panel propio, o "multiple", si se desean todas las series en un único panel.
#' @param leyenda Lugar donde situar la leyenda; por defecto, "bottom", pero puede especificarfse "right".
#'
#' @return Llamada por su efecto secundario, consistente en generar los gráficos.
#' @export
#'
#' @examples
#'
mggplot <- function(x, columnas=names(x),
                    tipo=c("panel","multiple"),
                    titulo=NULL,
                    leyenda=c("bottom","right")) {
  if ( !is.list(x))
    stop("Argument 'x' must be a list of indices.")
  basedate <- x[[1]]@basedate
  x <- do.call("merge", lapply(x[columnas], FUN=function(x) x@ICV))
  fecha <- index(x)
  tmp   <- cbind(fecha, as.data.frame(x))
  tmp   <- tmp %>% gather(Ambito,indice,-fecha)
  if (tipo[1]=="panel")
    p <- ggplot(tmp) +
    geom_line(aes(x=fecha, y=indice, colour=Ambito)) +
    xlab("Fecha") + ylab("Índice") +
    facet_wrap(~ Ambito) +
    theme(legend.position=leyenda[1])
  else if (tipo[1]=="multiple")
    p <- ggplot(tmp) +
    geom_line(aes(x=fecha, y=indice, colour=Ambito)) +
    xlab("Fecha") + ylab("Índice") +
    theme(legend.position=leyenda[1])
  if (!is.null(titulo))
  p <- p + labs(title=titulo,
                subtitle=paste("Base: ",basedate," = 100.", sep="")
  )
  return(p)
}

#' pols
#'
#' Generación de poligonos por agrupación de otros. Por ejemplo, podemos tener una
#' dataframe espacial (\code{base}) con polígonos de municipios y querer otra con polígonos de
#' comarcas o áreas funcionales. Entonces, en \code{var.agreg} daríamos el nombre de la
#' columna que nombra las comarcas o áreas funcionales
#' y obtendríamos una dataframe espacial con los contornos de las mismas.
#'
#' Podemos especificar que queremos solo un fragmento indicado una variable (\code{sub.var})
#' y los niveles de la misma a que deseamos limitarnos (\code{sub.set}). Por ejemplo,
#' \code{sub.var} podría ser \code{PROVINCIA} y \code{sub.set} \code{VIZCAYA}. Aunque la dataframe
#' original contuviera los tres territorios de la CAPV, el resltado contendría la agregación
#' por comarcas o áreas funcionales de sólo Bizkaia.
#'
#'
#' @param base  Dataframe espacial
#' @param sub.var Variable a partir de la cual seleccionaremos un subconjunto, si procede.
#' @param sub.set Niveles de \code{sub.var} definiendo el subconjunto que queremos.
#' @param var.agreg Variable de agregación: se agruparán todos los polígonos que tengan mismos nivel de esta variable.
#'
#' @return Una dataframe espacial de polígonos para los ámbitos agregados.
#' @export
#'
#' @examples
#'
pols <- function(base, sub.var=NULL, sub.set=NULL,
                 var.agreg=NULL) {
  if (!is.null(sub.var))
    base <- subset(base,
                   eval(parse(text=paste(sub.var," %in% ",sub.set))),
                   env=base@data)
  pols <- gUnaryUnion(base, id=as.character(eval(parse(text=var.agreg),
                                                 env=base@data)))
  return(invisible(fortify(pols, region=var.agreg)))
}

#' SlicedIndex
#'
#' Index computed in slices
#'
#' @param cal.pts Points at which the index is to be computed. It has the form of a spatial object, containing in the @data slot all required variables in formula `frm`  with values typically set for a "median house", a type of dwelling characteristic around the calibration point. It also needs to contain a date variable to create the time slices.
#' @param spdatos Spatial data frame containing data and the `date` variable
#' @param from Starting point in time for index computation.
#' @param to End point in time for index computation.
#' @param base Base ( = 100) period for index computation.
#' @param frm Formula to be used in GWR. All variables must be present in the @data slot of `spdatos`.
#' @param bw  Bandwidth to be used in GWR.
#' @param area.radius Area to subselect data. Usefull when we want a local index and wnat to discard observations which are hardly relevant. By default, five times the selected GWR bandwidth.
#' @param date  Variable in the @data slot of `spdatos` which gives the date of each observation.
#' @param slice.by  Function to be applied to `date` in order to obtain the time slices. It must produce a time value of the same type as `to` and  `from`.
#'
#' @return
#' @export
#'
#' @examples
#'
SlicedIndex <- function(cal.pts,
           datos = NULL,
           from = NULL,
           to = NULL,
           base = NULL,
           frm = NULL,
           bw = 1000,
           coords = NULL,
           area.radius = 5 * bw,
           date = NULL,
           slice.by = NULL)
  {
    require(zoo)
    if (is.null(date))
      stop("'date' not specified with no default.")
    if (is.null(slice.by))
      stop("'slice.by' not specified with no default.")
    spdatos <- datos
    datos   <- get_all_vars(frm, datos)
    spdatos <- spdatos[complete.cases(datos), ]
    coordinates(spdatos) <- eval(parse(text=paste("~ ",
                                                  paste(coords,
                                                        collapse="+"))))
    #
    # Convert observations to selected time interval
    #
    spdatos@data[, date] <- sapply(spdatos@data[, date], FUN = slice.by)
    spdatos <- spdatos[(spdatos@data[, date] >= from) &
                         (spdatos@data[, date] <= to), ]
    slices  <- sort(unique(spdatos@data[, date]))
    #
    # For each calibration point, take a subset of observations within
    # area.radius distance (projected coordinates are assumed)
    #
    xy      <- coordinates(cal.pts)
    preds   <- as.data.frame(matrix(0, length(slices), nrow(cal.pts)))
    colnames(preds) <- rownames(cal.pts)
    preds   <- zoo(preds, order.by = slices)
    for (i in nrow(cal.pts)) {
      sel <- rep(FALSE, nrow(spdatos@data))
      sel <- sel | (colSums((t(coordinates(spdatos)) - xy[i, ]) ^ 2) <
                      area.radius ^ 2)
      for (j in seq_along(slices)) {
        subsel <- sel & (spdatos@data[, date] == slices[j])
        cat("Processing slice:  ", as.character(slices[j]), "\n")
        cat(sum(subsel), "\n")
        spdatos.sub <- spdatos[subsel,]
        xx <- gwr(formula=frm,
                  data = spdatos.sub,
                  bandwidth = bw,
                  hatmatrix = TRUE)
        preds[j, ] <- gwr(formula = frm,
                          data = spdatos.sub,
                          bandwidth = bw,
                          fit.points = cal.pts,
                          predict = TRUE,
                          fittedGWRobject=xx)$SDF@data$"(Intercept)"
      }
    }
    return(preds)
  }

