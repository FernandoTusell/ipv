#' Estimación semiparamétrica mediante 'backfitting'.
#'
#' Realiza la estimación de un modelo semi-paramétrico, cuya parte
#' paramétrica esta formada por una regresión espacialmente ponderada y
#' cuya parte no paramétrica está formada por un spline cúbico global
#' (único para todas las observaciones)
#'
#' @param frm.param Fórmula de la parte paramétrica (GWR) del modelo.
#' @param smooth.term Término no paramétrico del modelo.
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
#'
#' @return Un índice de precios, con base en el día especificado, en formato de serie temporal \code{zoo.}
#' @export
#'
#' @examples Ver viñeta.
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
                        plotind=FALSE) {
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
  #  la parte paramétrica el(los) termino(s) suave(s) para
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
    indice0 <- ConsInd(modelo = mod.gam,
                       fechas = spdatos@data[,var.fecha, drop=TRUE],
                       base = baseday,
                       plot=FALSE)
  }
  #
  #  Comienza ahora la alternancia entre estimaciones de la
  #  parte paramétrica y no paramétrica del modelo, hasta
  #  (esperamos) convergencia
  #
  newind <- indice0 ; lastind <- 0 * newind ; iter <-0
  #
  #  Mientras la máxima discrepancia entre dos índices
  #  sucesivos sea pequeña, continúav
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
                       base = baseday,
                       plot=plotind,
                       x.base=x.base)
    cat("Backfitting iteración ",iter,"\n")
  }
  stopCluster(cl)
  #
  #  Finalizado el ajuste, representamos el índice si plotind=TRUE
  #
  ConsInd(modelo = mod.gam,
          fechas = spdatos@data[,var.fecha, drop=TRUE],
          base = baseday,
          plot=plotind,
          x.base=x.base)
  return(newind)
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
#' @examples Ver viñeta.
#'
BackFittingLocal <- function(frm.param,
                             smooth.term='s(x,bs="cr",k=24)',
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
                             tol=0.001,
                             plotind=FALSE) {
  require(mgcv)
  require(spgwr)
  require(parallel)
  if (is.null(cal.pts) & is.null(bws))
    global.fit <- TRUE
  else if (xor(is.null(cal.pts), is.null(bws)))
    stop("Both cal.pts and is.null must be NULL or non-NULL.")
  else
    global.fit <- FALSE
  if (missing(cores))
    cl <- makeCluster(detectCores())
  else
    cl <- makeCluster(cores)
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
  frm.global <- update(frm.param,
                       eval(parse(text=paste("~ . + ",
                                             smooth.term))))
  frm.tot    <- update(frm.param,
                       eval(parse(text=paste("~ . + y.orig + ",
                                             paste(c(coords,var.fecha),
                                                   collapse="+")))))
  if (!is.null(var.loc)) {
    frm.gwr    <- update(frm.param,
                         eval(parse(text=paste("~ . - ",
                                               var.loc))))
  } else {
    frm.gwr <- frm.param
  }

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
  for (p in 1:npts) {
    if (!global.fit) {
      cat("Iniciado índice",p,"\n")
      punto   <- cal.pts[p,]                       # Punto de cálculo
      d2 <- apply(sweep(coordinates(spdatos),2, punto)^2,1,sum)
                                                   # Distancias al cuadrado al punto de cálculo
    spdatos@data$w <- spgwr::gwr.Gauss(d2,bws)  # Pesos
    }
    mod.gam <- gam(formula=frm.global, data=spdatos@data, weights=w)
    indice0 <- ConsInd(modelo = mod.gam,
                       fechas = spdatos@data[,var.fecha, drop=TRUE],
                       base = baseday,
                       plot=FALSE)
    #
    #  Comienza ahora la alternancia entre estimaciones de la
    #  parte paramétrica y no paramétrica del modelo, hasta
    #  (esperamos) convergencia
    #
    newind <- indice0 ; lastind <- 0 * newind ; iter <-0
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
      newind  <- ConsInd(modelo = mod.gam,
                         fechas = spdatos@data[,var.fecha, drop=TRUE],
                         base = baseday,
                         plot=plotind,
                         x.base=x.base)
      cat("Backfitting iteración ",iter,"\n")
    }
    indices[[p]] <- newind
  }
  stopCluster(cl)
  if (global.fit)
    indices=indices[[1]]
  else
    indices <- list(cal.pts, indices)
  return(indices)
}

#' cloromap
#'
#' Crea mapas codificando una variable de interés por zona con diferentes colores
#'
#' @param sp Dataframe espacial, en cuyo slot @data está la variable \code{var.sp} definiendo los ámbitos sobre los que se agrega.
#' @param var.sp Variable proporcionando la desagregación espacial.
#' @param df Dataframe ordinaria, conteniendo una variable \code{var.df} emparejable con \code{var.sp.}
#' @param var.df Variable emparejable con \code{var.sp}; ambas han de tener el mismo número de niveles e idénticamente nombrados.
#' @param var Variable a representar. Debe ser una columna de \code{df.}
#' @param sub.var Variable de \code{sp} utilizada para definir un subconjunto de filas de la misma.
#' @param sub.set Vector de valores de \code{sub.var} que definen un subconjunto de filas de \code{sp}.
#' @param legend.title Encabezamiento de la leyenda.
#' @param logcolor FALSE si la escala de color es lineal, TRUE si logarítmica
#' @param graf.title Encabedzamiento del gráfico.
#' @param ... Parámetros adicionales a pasar a la función \code{ggplot} subyacente.
#'
#' @return Un mapa, en formato ggplot2, que por defecto se imprime.
#' @export
#'
#' @examples Ver viñeta
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
#' Dado un índice, definido sobre un conjunto de abscisas temporales con "huecos", define un nuevo índice interpolando los periodos sin observaciones.
#'
#' @param indice Un índice tal como lo devuelve \code{ConsInd.}
#' @param fechas Vector de fechas cuyos extremos definirán el nuevo fechado del índice
#'
#' @return Un índice diario, sobre todos los días que van de \code{min(fechas)} a \code{max(fechas).}
#' @export
#'
#' @examples Ver viñeta.
#'
CompFechas <- function(indice,fechas) {
  scratch <- zoo(0,order.by=seq.Date(from=min(fechas),to=max(fechas)+1,by="day"))
  indice  <- merge(indice,scratch)[,1]
  return(na.approx(indice,rule=2))
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
#' @param ylabel Rótulo eje de ordenadas
#' @param plot Variable lógica: si es \code{FALSE} (default) no se genera el plot, en caso contrario sí.
#' @param x.base Primera fecha en el plot; habitualmente en blanco y obtenida a partir de \code{fechas.}
#'
#' @return Un indice en formato \code{zoo}, que opcionalmente (plot=TRUE) es representado gráficamente.
#' @export
#'
#' @examples Ver viñeta.
#'
ConsInd <- function(modelo=NULL,base="2008-02-01",
                    conf=0.95,
                    fechas=NULL, tit=tit,
                    ylabel="Índice",
                    plot=FALSE,
                    x.base) {
  plotdata <- plot(modelo, select=0)[[1]]
  if (missing(x.base))
    x.base <- min(fechas) - 1
  x <- as.Date(plotdata$x + as.Date(x.base))


  y   <- plotdata$fit
  se  <- plotdata$se
  ICV <- data.frame(y, y - qnorm(conf)*se,
                    y + qnorm(conf)*se)
  difs <- as.numeric(abs(x-as.Date(base)))
  orig <- seq_along(x)[difs==min(difs)]   # la fecha más próxima a la base
  ICV  <- zoo(100*exp(ICV - as.numeric(ICV[orig,1])),
              order.by=x)
  colnames(ICV) <- c("Indice","lcb", "ucb")
  if (plot) {
    p <- ggplot(ICV, aes(x = index(ICV), y = Indice))  +
      geom_line() +
      xlab("Fecha") + ylab("Índice") +
      geom_line(aes(x = index(ICV), y = lcb),
                col="blue", linetype="dotted") +
      geom_line(aes(x = index(ICV), y = ucb),
                col="blue", linetype="dotted") +
      #geom_vline(xintercept=as.Date(base),
      #            color="red", size=0.5) +
      ggtitle(tit,
              sub=paste("Base: ",base," = 100. Int. conf. ",
                        100*conf,"%"))
    print(p)
  }
  #  Antes de devolver el índice, que sólo está calculado para fechas que aparecían en
  #  el vector "fechas", vamos a completarlo para TODAS las fechas posibles entre la
  #  primera y la última. Eso garantiza que nos permitirá deflactar cualquier magnitud
  #  definida en el mismo intervalo de fechas, pero no necesariamente sobre las mismas
  #  abscisas temporales. Empleamos la función previamente
  #  definida 'CompFechas'.
  #
  return( invisible(CompFechas(indice=ICV[,1], fechas=x)) )
}


#' IndZonas
#'
#' Estimación de índices globales por zonas; es unam mera función de utilidad que invoca reiteradamente  \code{ConsInd}.
#' @param datos Dataframe de datos
#' @param zonas Nombre de una columna de \code{datos} que define las zonas para las que se desea calcular los índices.
#' @param frm Fórmula especificando el modelo base que quiere emplearse, en la forma esperada por \code{gam.} La fórmula es común a todos los índices computados por la función.
#' @param base Fecha base (índice=100). Es una fecha diaria en formato yyyy-mm-dd.
#' @param plot \code{FALSE} (default) si no se desean gráficos, \code{TRUE} en caso contrario.
#'
#' @return Una lista con los índices calculados, en el formato en que los devuelve \code{ConsInd}.
#' @export
#'
#' @examples Ver viñeta
#'
IndZonas <- function(datos,
                     zonas,
                     frm = NULL,
                     base="2008-01-01",
                     plot=FALSE) {
  x <- vector("list",0)
  k <- match(zonas, colnames(datos))
  for (i in unique(datos[,k])) {
    sel <- datos[,k] == i
    mod <- gam( frm, data=datos[sel,])
    indice0 <- ConsInd(modelo = mod,
                       fechas = datos$FEC_CREACION,
                       base = base,
                       conf= 0.95,
                       tit = paste("Índice precio vivienda. Ámbito: ",
                                   i, sep="")
    )
    x[[i]] <- indice0
  }
  return(x)
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
#' @examples Ver viñeta.
#'
mggplot <- function(x, columnas=names(x),
                    tipo=c("panel","multiple"),
                    leyenda=c("bottom","right")) {
  if ( is.list(x))
    x <- do.call("merge", x[columnas])
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
  print(p)
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
#' @examples Ver viñeta
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
