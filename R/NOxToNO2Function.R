#' Perform NOx to NO2 Calculation
#'
#' Perform NOx to NO2 Calculation
#' @param df the receptor data frame
#' @return NO2 concentrations from road NOx concentrations
#' @examples
#' ReceptorFile <- data.frame(ID = c("R1", "R2", "R4"), RoadNOx = c(12.1, 14.65, 13.4),
#' No2Backgd = c(19, 24.5, 17.84), LAs = c("Armagh Banbridge and Craigavon"));
#'
#' NOxToNO2(yearselect = 2019, receptors = ReceptorFile);
#' @export

NOxToNO2 <- function(yearselect, receptors) {

  suppressMessages(library(dplyr, warn.conflict = FALSE, quietly = TRUE, include=FALSE))

  RegionalPols <- dplyr::filter(read.csv("https://raw.githubusercontent.com/AndyDBrown/AQTools/main/R/data/RegionalPols.csv", header=T), year == yearselect)
  fno2_data <- dplyr::filter(read.csv("https://raw.githubusercontent.com/AndyDBrown/AQTools/main/R/data/fno2_data.csv", header=T), year == yearselect)
  receptors <- as.data.frame(receptors)
   if (length(which(is.na(match(c("RoadNOx", "No2Backgd", "ID", "LA"), colnames(receptors))))) == 0) # match "RoadNOx", "No2Backgd", "ID", "LA" to receptors col names then count any F values, i.e any incidences where colnames don't match
   {
  receptors <- dplyr::left_join(receptors,RegionalPols, by = "LA")
  receptors <- dplyr::left_join(receptors, fno2_data, by = "year")
  receptors$RoadNO2 <- 0
  ####1.1 Insert NOx to NO2 calc from shiny tool#######

  for (g in 1:nrow(receptors)){

    userbackNO2 <- receptors$No2Backgd[g]
    backNO2 <- userbackNO2 / 1.91
    locVehfNO2 <- receptors$fno2[g]
    rfNO2 <- receptors$fno2[g]
    fractionno2_default <- receptors$fno2[g]
    fractionno2_back <- receptors$fno2[g]
    rO3 <- receptors$rO3[g]
    rNOx <- receptors$rNOx[g]
    rNO2 <- receptors$rNO2[g]
    userroadNox <- receptors$RoadNOx[g]

    O3back <- rO3 / 2
    noxback = rNOx / 1.91
    no2back = rNO2 / 1.91

    O3back2 <- O3back
    noxback2 <- noxback
    no2back2 <- no2back

    h <- 800
    rbrcO3 <- 500
    hb10 <- h/10
    ustar <- 0.4
    z0 <- 1
    afac <- locVehfNO2
    k1 <- 1.42E-20
    j <- 0.002


    if (backNO2 < no2back) {
      no2back2 <- backNO2

      afac <- fractionno2_back
      NO2_site <- backNO2
      NO2_back <- no2back
      NOx_back <- noxback
      O3_back <- O3back


      NOxb <- 2.59378E+16 * NOx_back
      O3b <- 2.59378E+16 * O3_back
      NO2b <- 2.59378E+16 * NO2_back
      NO2 <- 2.59378E+16 * NO2_site

      a <- k1 * afac
      b <- k1 * (O3b + NO2b - afac * NOxb - NO2 - afac * NO2)

      c <- k1 * NO2 * NO2 + k1 * afac * NO2 * NOxb - k1 * NO2 * O3b - k1 * NO2b * NO2 - j * NO2

      NOx <- (-b + sqrt(b * b - 4 * a * c)) / 2 / a

      O3 <- j * NO2 / k1 / (NOx - NO2)

      backNOx <- NOx / 2.59378E+16
      backO3 <- O3 / 2.59378E+16

      noxback2 <- backNOx
      O3back2 <- backO3
      backnox <- noxback2

      RoadNOx <- userroadNox / 1.91
      HighNOx <- backNOx + 2 * RoadNOx

      fractionno2 <- fractionno2_default
      fractionno2_high <- (fractionno2_back * backnox + 2 * RoadNOx * fractionno2) / HighNOx
    } else {

      ###start grad7####
      afac <- fractionno2_back
      NO2_site <- backNO2
      NO2_back <- no2back
      NOx_back <- noxback
      O3_back <- O3back

      gradient7 <- data.frame(col = 1:11)

      gradient7$NOx <- 2.59378E+16 * NOx_back
      gradient7$O3 <- 2.59378E+16 * O3_back
      gradient7$NO2 <- 2.59378E+16 * NO2_back
      gradient7$NO <- 2.59378E+16 * (NOx_back - NO2_back)

      gradient7$NO2[1] <- 2.59378E+16 * NO2_site


      eO3 <- 0
      ku = 0.41 * ustar
      ra = 1 / ku * log(hb10 / z0)
      dr = ra / 10

      alpha = 1 / ku / dr / dr

      eno2 = (gradient7$NO2[1] - gradient7$NO2[11]) / ra
      eno = eno2 / afac * (1 - afac)

      errorout <- NULL
      errorsum <- 1
      while (errorsum >0.0001) {

        iterNO <-  data.frame(col = 1)
        iterNO$r =0
        iterNO$z =0
        iterNO$a = 0
        iterNO$b = -1/dr
        iterNO$c = 1/dr
        iterNO$d = -eno
        iterNO$b2 = iterNO$b
        iterNO$d2 = iterNO$d

        for (i in 1:11){
          gradient7$tempo[i] <- gradient7$NO[i]
        }

        subiterNO <- data.frame(col = 1:9)
        subiterNO$r <- subiterNO$col * dr
        subiterNO$z <- z0 * exp(ku * subiterNO$r)
        subiterNO$a = alpha
        subiterNO$b = 0
        subiterNO$c = alpha
        subiterNO$d = 0
        subiterNO$b2 =0
        subiterNO$d2 =0

        iterNO <- rbind(iterNO,subiterNO)

        for (i in 2:10){
          iterNO$b[i] = -2 * alpha - k1 * iterNO$z[i] * gradient7$O3[i]
          iterNO$d[i] <- (j*-1) * iterNO$z[i] * gradient7$NO2[i]
        }

        for (i in 2:10){
          #iterNO$d[i] <- (j*-1) * iterNO$z[i] * gradient7$NO2[i]
          iterNO$b2[i] <- iterNO$b[i] - iterNO$a[i] * iterNO$c[i-1] / iterNO$b2[i-1]
          iterNO$d2[i] <- iterNO$d[i] - iterNO$a[i] / iterNO$b2[i-1] * iterNO$d2[i-1]
        }

        gradient7$NO[10] <- (iterNO$d2[10] - iterNO$c[10] * gradient7$NO[11]) / iterNO$b2[10]
        if (gradient7$NO[10] < 0) {gradient7$NO[10] <- 0
        }

        for (i in 2:9) {
          n = 11 - i
          gradient7$NO[n] <- iterNO$d[n+1] / iterNO$a[n+1] - iterNO$b[n+1] / iterNO$a[n+1] * gradient7$NO[n+1] - iterNO$c[n+1] / iterNO$a[n+1] * gradient7$NO[n+2]
          if (gradient7$NO[n] < 0) {gradient7$NO[n] <- 0}
        }

        gradient7$NO[1] <- iterNO$d[1] / iterNO$b[1] - iterNO$c[1] / iterNO$b[1] * gradient7$NO[2]
        if (gradient7$NO[1] < 0) {gradient7$NO[1] <- 0}


        #NO2 iteration

        iterNO2 <-  data.frame(col = 1)
        iterNO2$r =0
        iterNO2$z =0
        iterNO2$a = 0
        iterNO2$b = -1/dr
        iterNO2$c = 1/dr
        iterNO2$d = -eno2
        iterNO2$b2 = iterNO2$b
        iterNO2$d2 = iterNO2$d

        subiterNO2 <- data.frame(col = 1:9)
        subiterNO2$r <- subiterNO2$col * dr
        subiterNO2$z <- z0 * exp(ku * subiterNO2$r)
        subiterNO2$a = alpha
        subiterNO2$b = 0
        subiterNO2$c = alpha
        subiterNO2$d = 0
        subiterNO2$b2 =0
        subiterNO2$d2 =0



        iterNO2 <- rbind(iterNO2,subiterNO2)

        for (i in 2:10){
          iterNO2$b[i] <- -2 * alpha - j * iterNO2$z[i]
          iterNO2$d[i] <- (k1*-1) * iterNO2$z[i] * gradient7$O3[i] * gradient7$NO[i]
        }

        for (i in 2:10){
          iterNO2$b2[i] <- iterNO2$b[i] - iterNO2$a[i] * iterNO2$c[i-1] / iterNO2$b2[i-1]
          iterNO2$d2[i] <- iterNO2$d[i] - iterNO2$a[i] / iterNO2$b2[i-1] * iterNO2$d2[i-1]
        }


        gradient7$NO2[10] <- (iterNO2$d2[10] - iterNO2$c[10] * gradient7$NO2[11]) / iterNO2$b2[10]
        if (gradient7$NO2[10] < 0) {gradient7$NO2[10] <- 0}


        for (i in 2:9) {
          n = 11 - i
          gradient7$NO2[n] <- iterNO2$d[n+1] / iterNO2$a[n+1] - iterNO2$b[n+1] / iterNO2$a[n+1] *
            gradient7$NO2[n+1] - iterNO2$c[n+1] / iterNO2$a[n+1] * gradient7$NO2[n+2]
          if (gradient7$NO2[n] < 0) {gradient7$NO2[n] <- 0}
        }

        tempo2 <- iterNO2$d[1] / iterNO2$b[1] - iterNO2$c[1] / iterNO2$b[1] * gradient7$NO2[2]
        if (tempo2 < 0) {tempo2 <- 0}


        # O3 iteration

        iterO3 <-  data.frame(col = 1)
        iterO3$r =0
        iterO3$z = 0
        iterO3$a = 0
        iterO3$b = -1/dr -1/rbrcO3
        iterO3$c = 1/dr
        iterO3$d = 0
        iterO3$b2 = iterO3$b
        iterO3$d2 = iterO3$d

        subiterO3 <- data.frame(col = 1:9)
        subiterO3$r <- subiterO3$col * dr
        subiterO3$z <- z0 * exp(ku * subiterO3$r)
        subiterO3$a = alpha
        subiterO3$b = 0
        subiterO3$c = alpha
        subiterO3$d = 0
        subiterO3$b2 =0
        subiterO3$d2 =0

        iterO3 <- rbind(iterO3,subiterO3)

        for (i in 2:10){
          iterO3$b[i] = -2 * alpha - k1 * iterO3$z[i] * gradient7$NO[i]
          iterO3$d[i] <- (j*-1) * iterO3$z[i] * gradient7$NO2[i]
        }


        for (i in 2:10){
          iterO3$b2[i] <- iterO3$b[i] - iterO3$a[i] * iterO3$c[i-1] / iterO3$b2[i-1]
          iterO3$d2[i] <- iterO3$d[i] - iterO3$a[i] / iterO3$b2[i-1] * iterO3$d2[i-1]
        }

        gradient7$O3[10] <- (iterO3$d2[10] - iterO3$c[10] * gradient7$O3[11]) / iterO3$b2[10]
        if (gradient7$O3[10] < 0) {gradient7$O3[10] <- 0}


        for (i in 2:9) {
          n = 11 - i
          gradient7$O3[n] <- iterO3$d[n+1] / iterO3$a[n+1] - iterO3$b[n+1] / iterO3$a[n+1] * gradient7$O3[n+1] -
            iterO3$c[n+1] / iterO3$a[n+1] * gradient7$O3[n+2]
          if (gradient7$O3[n] < 0) {gradient7$O3[n] <- 0}
        }

        gradient7$O3[1] <- iterO3$d[1] / iterO3$b[1] - iterO3$c[1] / iterO3$b[1] * gradient7$O3[2]
        if (gradient7$O3[1] < 0) {gradient7$O3[1] <- 0}

        eno2 <- (gradient7$NO2[1] - gradient7$NO2[11]) / (tempo2 - gradient7$NO2[11]) * eno2
        eno <- eno2 / afac * (1 - afac)

        #converge Grad7

        errorsum0 <- 0
        for (i in 1:10) {
          errorsum0 <- errorsum0 + abs(gradient7$NO[i] - gradient7$tempo[i]) / 2.59378E+16

        }
        #errorsum <- errorsum0
        y <- NULL
        y <- (gradient7$NO2[1] + gradient7$NO[1]) / 2.59378E+16

        errorsumm <- cbind(errorsum0,y)
        errorout <- rbind(errorout, errorsumm)

        errorsum <- errorsum0
      }
      grad7out <-data.frame(errorout)
      grad7out <- tail(grad7out,1)
      backNOx <- grad7out[1,2]

      noxback2 <- backNOx
      O3back2 <- O3back
      noxback2 <- noxback

      RoadNOx <- userroadNox / 1.91
      HighNOx <- backNOx + 2 * RoadNOx

      fractionno2 <- fractionno2_default
      fractionno2_high <- (fractionno2_back * backNOx + 2 * RoadNOx * fractionno2) / HighNOx


    }
    ###start grad5 windward side####
    afac <- fractionno2_high
    NOx_site <- backNOx
    NO2_back <- no2back2
    NOx_back <- noxback2
    O3_back <- O3back2

    h <- 800
    rbrcO3 <- 500
    hb10 <- h/10
    ustar <- 0.4
    z0 <- 1

    k1 <- 1.42E-20
    j <- 0.002
    almo <- 0


    gradient5 <- data.frame(col = 1:11)
    gradient5$NOx <- 2.59378E+16 * NOx_back
    gradient5$O3 <- 2.59378E+16 * O3_back
    gradient5$NO2 <- 2.59378E+16 * NO2_back
    gradient5$NO <- 2.59378E+16 * (NOx_back - NO2_back)

    gradient5$NOx[1] <- 2.59378E+16 * NOx_site

    eO3 <- 0
    ku = 0.41 * ustar
    ra = 1 / ku * log(hb10 / z0)
    dr = ra / 10

    alpha = 1 / ku / dr / dr

    errorout5 <- NULL
    errorsum <- 0

    iterNOx <-  data.frame(col = 1)
    iterNOx$r =0
    iterNOx$z =0
    iterNOx$a = 0
    iterNOx$b = 1
    iterNOx$c = 0
    iterNOx$d = gradient5$NOx[1]
    iterNOx$b2 = 1
    iterNOx$d2 = iterNOx$d[1]

    subiterNOx <- data.frame(col = 1:9)
    subiterNOx$r <- subiterNOx$col * dr
    subiterNOx$z <- 0
    subiterNOx$a = alpha
    subiterNOx$b = -2 * alpha
    subiterNOx$c = alpha
    subiterNOx$d = 0
    subiterNOx$b2 =0
    subiterNOx$d2 =0

    for (i in 1:9){
      subiterNOx$z[i] <- z0 * exp(ku * subiterNOx$r[i])
    }

    iterNOx <- rbind(iterNOx,subiterNOx)

    for (i in 2:10){
      iterNOx$b2[i] <- iterNOx$b[i] - iterNOx$a[i] * iterNOx$c[i-1] / iterNOx$b2[i-1]
      iterNOx$d2[i] <- iterNOx$d[i] - iterNOx$a[i] / iterNOx$b2[i-1] * iterNOx$d2[i-1]
    }

    gradient5$NOx[10] <- (iterNOx$d2[10] - iterNOx$c[10] * gradient5$NOx[11]) / iterNOx$b2[10]
    if (gradient5$NOx[10] < 0) {gradient5$NOx[10] <- 0
    }

    for (i in 2:9) {
      n = 11 - i
      gradient5$NOx[n] <- iterNOx$d[n+1] / iterNOx$a[n+1] - iterNOx$b[n+1] / iterNOx$a[n+1] * gradient5$NOx[n+1] -
        iterNOx$c[n+1] / iterNOx$a[n+1] * gradient5$NOx[n+2]
      if (gradient5$NOx[n] < 0) {gradient5$NOx[n] <- 0}
    }

    eno2 = afac * (gradient5$NOx[1] - gradient5$NOx[2]) / dr
    eno = eno2 / afac * (1 - afac)

    # NO iter
    errorout <- NULL
    errorsum <- 1

    while (errorsum >0.0001){

      iterNO <-  data.frame(col = 1)
      iterNO$r =0
      iterNO$z =0
      iterNO$a = 0
      iterNO$b = -1/dr
      iterNO$c = 1/dr
      iterNO$d = -eno
      iterNO$b2 = iterNO$b
      iterNO$d2 = iterNO$d

      subiterNO <- data.frame(col = 1:9)
      subiterNO$r <- subiterNO$col * dr
      subiterNO$z <- z0 * exp(ku * subiterNO$r)
      subiterNO$a = alpha
      subiterNO$b = 0
      subiterNO$c = alpha
      subiterNO$d = 0
      subiterNO$b2 =0
      subiterNO$d2 =0

      iterNO <- rbind(iterNO,subiterNO)

      for (i in 2:10){
        iterNO$b[i] = -2 * alpha - k1 * iterNO$z[i] * gradient5$O3[i]
        iterNO$d[i] <- (j*-1) * iterNO$z[i] * gradient5$NO2[i]
      }

      for (i in 2:10){
        iterNO$b2[i] <- iterNO$b[i] - iterNO$a[i] * iterNO$c[i-1] / iterNO$b2[i-1]
        iterNO$d2[i] <- iterNO$d[i] - iterNO$a[i] / iterNO$b2[i-1] * iterNO$d2[i-1]
      }

      gradient5$NO[10] <- (iterNO$d2[10] - iterNO$c[10] * gradient5$NO[11]) / iterNO$b2[10]
      if (gradient5$NO[10] < 0) {gradient5$NO[10] <- 0
      }

      for (i in 2:9) {
        n = 11 - i
        gradient5$NO[n] <- iterNO$d[n+1] / iterNO$a[n+1] - iterNO$b[n+1] / iterNO$a[n+1] * gradient5$NO[n+1] - iterNO$c[n+1] / iterNO$a[n+1] * gradient5$NO[n+2]
        if (gradient5$NO[n] < 0) {gradient5$NO[n] <- 0}
      }

      gradient5$NO[1] <- iterNO$d[1] / iterNO$b[1] - iterNO$c[1] / iterNO$b[1] * gradient5$NO[2]
      if (gradient5$NO[1] < 0) {gradient5$NO[1] <- 0}


      #NO2 iteration

      iterNO2 <-  data.frame(col = 1)
      iterNO2$r =0
      iterNO2$z =0
      iterNO2$a = 0
      iterNO2$b = -1/dr
      iterNO2$c = 1/dr
      iterNO2$d = -eno2
      iterNO2$b2 = iterNO2$b
      iterNO2$d2 = iterNO2$d

      for (i in 1:11){
        gradient5$tempo[i] <- gradient5$NO2[i]
      }

      subiterNO2 <- data.frame(col = 1:9)
      subiterNO2$r <- subiterNO2$col * dr
      subiterNO2$z <- z0 * exp(ku * subiterNO2$r)
      subiterNO2$a = alpha
      subiterNO2$b = 0
      subiterNO2$c = alpha
      subiterNO2$d = 0
      subiterNO2$b2 =0
      subiterNO2$d2 =0

      iterNO2 <- rbind(iterNO2,subiterNO2)

      for (i in 2:10){
        iterNO2$b[i] <- -2 * alpha - j * iterNO2$z[i]
        iterNO2$d[i] <- (k1*-1) * iterNO2$z[i] * gradient5$O3[i] * gradient5$NO[i]
      }

      for (i in 2:10){
        iterNO2$b2[i] <- iterNO2$b[i] - iterNO2$a[i] * iterNO2$c[i-1] / iterNO2$b2[i-1]
        iterNO2$d2[i] <- iterNO2$d[i] - iterNO2$a[i] / iterNO2$b2[i-1] * iterNO2$d2[i-1]
      }


      gradient5$NO2[10] <- (iterNO2$d2[10] - iterNO2$c[10] * gradient5$NO2[11]) / iterNO2$b2[10]
      if (gradient5$NO2[10] < 0) {gradient5$NO2[10] <- 0}


      for (i in 2:9) {
        n = 11 - i
        gradient5$NO2[n] <- iterNO2$d[n+1] / iterNO2$a[n+1] - iterNO2$b[n+1] / iterNO2$a[n+1] * gradient5$NO2[n+1] -
          iterNO2$c[n+1] / iterNO2$a[n+1] * gradient5$NO2[n+2]
        if (gradient5$NO2[n] < 0) {gradient5$NO2[n] <- 0}
      }

      gradient5$NO2[1] <- iterNO2$d[1] / iterNO2$b[1] - iterNO2$c[1] / iterNO2$b[1] * gradient5$NO2[2]

      # O3 iteration

      iterO3 <-  data.frame(col = 1)
      iterO3$r =0
      iterO3$z = 0
      iterO3$a = 0
      iterO3$b = -1/dr -1/rbrcO3
      iterO3$c = 1/dr
      iterO3$d = 0
      iterO3$b2 = iterO3$b
      iterO3$d2 = iterO3$d

      subiterO3 <- data.frame(col = 1:9)
      subiterO3$r <- subiterO3$col * dr
      subiterO3$z <- z0 * exp(ku * subiterO3$r)
      subiterO3$a = alpha
      subiterO3$b = 0
      subiterO3$c = alpha
      subiterO3$d = 0
      subiterO3$b2 =0
      subiterO3$d2 =0

      iterO3 <- rbind(iterO3,subiterO3)

      for (i in 2:10){
        iterO3$b[i] = -2 * alpha - k1 * iterO3$z[i] * gradient5$NO[i]
        iterO3$d[i] <- (j*-1) * iterO3$z[i] * gradient5$NO2[i]
      }

      for (i in 2:10){
        iterO3$b2[i] <- iterO3$b[i] - iterO3$a[i] * iterO3$c[i-1] / iterO3$b2[i-1]
        iterO3$d2[i] <- iterO3$d[i] - iterO3$a[i] / iterO3$b2[i-1] * iterO3$d2[i-1]
      }

      gradient5$O3[10] <- (iterO3$d2[10] - iterO3$c[10] * gradient5$O3[11]) / iterO3$b2[10]
      if (gradient5$O3[10] < 0) {gradient5$O3[10] <- 0}


      for (i in 2:9) {
        n = 11 - i
        gradient5$O3[n] <- iterO3$d[n+1] / iterO3$a[n+1] - iterO3$b[n+1] / iterO3$a[n+1] * gradient5$O3[n+1] -
          iterO3$c[n+1] / iterO3$a[n+1] * gradient5$O3[n+2]
        if (gradient5$O3[n] < 0) {gradient5$O3[n] <- 0}
      }

      gradient5$O3[1] <- iterO3$d[1] / iterO3$b[1] - iterO3$c[1] / iterO3$b[1] * gradient5$O3[2]
      if (gradient5$O3[1] < 0) {gradient5$O3[1] <- 0}

      #converge Grad5 windward

      errorsum0 <- 0
      for (i in 1:10) {
        errorsum0 <- errorsum0 + abs(gradient5$NO2[i] - gradient5$tempo[i]) / 2.59378E+16 *1.91

      }
      y <- NULL
      y <- gradient5$NO2[1] / 2.59378E+16
      errorsum5 <- cbind(errorsum0,y)
      errorout5 <- rbind(errorout5, errorsum5)
      errorsum <- errorsum0
    }
    grad5out <-data.frame(errorout5)
    grad5out <- tail(grad5out,1)
    backNO2ppb <- grad5out[1,2]
    backNO2 <- backNO2ppb*1.91


    ###start grad5 Leeward####
    afac <- fractionno2_high
    NOx_site <- backNOx
    NO2_back <- no2back2
    NOx_back <- noxback2
    O3_back <- O3back2

    h <- 800
    rbrcO3 <- 500
    hb10 <- h/10
    ustar <- 0.4
    z0 <- 1

    k1 <- 1.42E-20
    j <- 0.002
    almo <- 0


    gradient5 <- data.frame(col = 1:11)
    gradient5$NOx <- 2.59378E+16 * NOx_back
    gradient5$O3 <- 2.59378E+16 * O3_back
    gradient5$NO2 <- 2.59378E+16 * NO2_back
    gradient5$NO <- 2.59378E+16 * (NOx_back - NO2_back)

    gradient5$NOx[1] <- 2.59378E+16 * HighNOx

    eO3 <- 0
    ku = 0.41 * ustar
    ra = 1 / ku * log(hb10 / z0)
    dr = ra / 10

    alpha = 1 / ku / dr / dr

    errorout5 <- NULL
    errorsum <- 0

    iterNOx <-  data.frame(col = 1)
    iterNOx$r =0
    iterNOx$z =0
    iterNOx$a = 0
    iterNOx$b = 1
    iterNOx$c = 0
    iterNOx$d = gradient5$NOx[1]
    iterNOx$b2 = 1
    iterNOx$d2 = iterNOx$d[1]

    subiterNOx <- data.frame(col = 1:9)
    subiterNOx$r <- subiterNOx$col * dr
    subiterNOx$z <- 0
    subiterNOx$a = alpha
    subiterNOx$b = -2 * alpha
    subiterNOx$c = alpha
    subiterNOx$d = 0
    subiterNOx$b2 =0
    subiterNOx$d2 =0

    for (i in 1:9){
      subiterNOx$z[i] <- z0 * exp(ku * subiterNOx$r[i])
    }

    iterNOx <- rbind(iterNOx,subiterNOx)

    for (i in 2:10){
      iterNOx$b2[i] <- iterNOx$b[i] - iterNOx$a[i] * iterNOx$c[i-1] / iterNOx$b2[i-1]
      iterNOx$d2[i] <- iterNOx$d[i] - iterNOx$a[i] / iterNOx$b2[i-1] * iterNOx$d2[i-1]
    }

    gradient5$NOx[10] <- (iterNOx$d2[10] - iterNOx$c[10] * gradient5$NOx[11]) / iterNOx$b2[10]
    if (gradient5$NOx[10] < 0) {gradient5$NOx[10] <- 0
    }

    for (i in 2:9) {
      n = 11 - i
      gradient5$NOx[n] <- iterNOx$d[n+1] / iterNOx$a[n+1] - iterNOx$b[n+1] / iterNOx$a[n+1] * gradient5$NOx[n+1] - iterNOx$c[n+1] / iterNOx$a[n+1] * gradient5$NOx[n+2]
      if (gradient5$NOx[n] < 0) {gradient5$NOx[n] <- 0}
    }

    eno2 = afac * (gradient5$NOx[1] - gradient5$NOx[2]) / dr
    eno = eno2 / afac * (1 - afac)

    # NO iter
    errorout <- NULL
    errorsum <- 1
    while (errorsum >0.0001){

      iterNO <-  data.frame(col = 1)
      iterNO$r =0
      iterNO$z =0
      iterNO$a = 0
      iterNO$b = -1/dr
      iterNO$c = 1/dr
      iterNO$d = -eno
      iterNO$b2 = iterNO$b
      iterNO$d2 = iterNO$d

      subiterNO <- data.frame(col = 1:9)
      subiterNO$r <- subiterNO$col * dr
      subiterNO$z <- z0 * exp(ku * subiterNO$r)
      subiterNO$a = alpha
      subiterNO$b = 0
      subiterNO$c = alpha
      subiterNO$d = 0
      subiterNO$b2 =0
      subiterNO$d2 =0

      iterNO <- rbind(iterNO,subiterNO)

      for (i in 2:10){
        iterNO$b[i] = -2 * alpha - k1 * iterNO$z[i] * gradient5$O3[i]
        iterNO$d[i] <- (j*-1) * iterNO$z[i] * gradient5$NO2[i]
      }

      for (i in 2:10){
        iterNO$b2[i] <- iterNO$b[i] - iterNO$a[i] * iterNO$c[i-1] / iterNO$b2[i-1]
        iterNO$d2[i] <- iterNO$d[i] - iterNO$a[i] / iterNO$b2[i-1] * iterNO$d2[i-1]
      }

      gradient5$NO[10] <- (iterNO$d2[10] - iterNO$c[10] * gradient5$NO[11]) / iterNO$b2[10]
      if (gradient5$NO[10] < 0) {gradient5$NO[10] <- 0
      }

      for (i in 2:9) {
        n = 11 - i
        gradient5$NO[n] <- iterNO$d[n+1] / iterNO$a[n+1] - iterNO$b[n+1] / iterNO$a[n+1] * gradient5$NO[n+1] - iterNO$c[n+1] / iterNO$a[n+1] * gradient5$NO[n+2]
        if (gradient5$NO[n] < 0) {gradient5$NO[n] <- 0}
      }

      gradient5$NO[1] <- iterNO$d[1] / iterNO$b[1] - iterNO$c[1] / iterNO$b[1] * gradient5$NO[2]
      if (gradient5$NO[1] < 0) {gradient5$NO[1] <- 0}


      #NO2 iteration

      iterNO2 <-  data.frame(col = 1)
      iterNO2$r =0
      iterNO2$z =0
      iterNO2$a = 0
      iterNO2$b = -1/dr
      iterNO2$c = 1/dr
      iterNO2$d = -eno2
      iterNO2$b2 = iterNO2$b
      iterNO2$d2 = iterNO2$d

      for (i in 1:11){
        gradient5$tempo[i] <- gradient5$NO2[i]
      }

      subiterNO2 <- data.frame(col = 1:9)
      subiterNO2$r <- subiterNO2$col * dr
      subiterNO2$z <- z0 * exp(ku * subiterNO2$r)
      subiterNO2$a = alpha
      subiterNO2$b = 0
      subiterNO2$c = alpha
      subiterNO2$d = 0
      subiterNO2$b2 =0
      subiterNO2$d2 =0

      iterNO2 <- rbind(iterNO2,subiterNO2)

      for (i in 2:10){
        iterNO2$b[i] <- -2 * alpha - j * iterNO2$z[i]
        iterNO2$d[i] <- (k1*-1) * iterNO2$z[i] * gradient5$O3[i] * gradient5$NO[i]
      }

      for (i in 2:10){
        iterNO2$b2[i] <- iterNO2$b[i] - iterNO2$a[i] * iterNO2$c[i-1] / iterNO2$b2[i-1]
        iterNO2$d2[i] <- iterNO2$d[i] - iterNO2$a[i] / iterNO2$b2[i-1] * iterNO2$d2[i-1]
      }


      gradient5$NO2[10] <- (iterNO2$d2[10] - iterNO2$c[10] * gradient5$NO2[11]) / iterNO2$b2[10]
      if (gradient5$NO2[10] < 0) {gradient5$NO2[10] <- 0}

      for (i in 2:9) {
        n = 11 - i
        gradient5$NO2[n] <- iterNO2$d[n+1] / iterNO2$a[n+1] - iterNO2$b[n+1] / iterNO2$a[n+1] * gradient5$NO2[n+1] -
          iterNO2$c[n+1] / iterNO2$a[n+1] * gradient5$NO2[n+2]
        if (gradient5$NO2[n] < 0) {gradient5$NO2[n] <- 0}
      }

      gradient5$NO2[1] <- iterNO2$d[1] / iterNO2$b[1] - iterNO2$c[1] / iterNO2$b[1] * gradient5$NO2[2]

      # O3 iteration

      iterO3 <-  data.frame(col = 1)
      iterO3$r =0
      iterO3$z = 0
      iterO3$a = 0
      iterO3$b = -1/dr -1/rbrcO3
      iterO3$c = 1/dr
      iterO3$d = 0
      iterO3$b2 = iterO3$b
      iterO3$d2 = iterO3$d

      subiterO3 <- data.frame(col = 1:9)
      subiterO3$r <- subiterO3$col * dr
      subiterO3$z <- z0 * exp(ku * subiterO3$r)
      subiterO3$a = alpha
      subiterO3$b = 0
      subiterO3$c = alpha
      subiterO3$d = 0
      subiterO3$b2 =0
      subiterO3$d2 =0

      iterO3 <- rbind(iterO3,subiterO3)

      for (i in 2:10){
        iterO3$b[i] = -2 * alpha - k1 * iterO3$z[i] * gradient5$NO[i]
        iterO3$d[i] <- (j*-1) * iterO3$z[i] * gradient5$NO2[i]
      }

      for (i in 2:10){
        iterO3$b2[i] <- iterO3$b[i] - iterO3$a[i] * iterO3$c[i-1] / iterO3$b2[i-1]
        iterO3$d2[i] <- iterO3$d[i] - iterO3$a[i] / iterO3$b2[i-1] * iterO3$d2[i-1]
      }

      gradient5$O3[10] <- (iterO3$d2[10] - iterO3$c[10] * gradient5$O3[11]) / iterO3$b2[10]
      if (gradient5$O3[10] < 0) {gradient5$O3[10] <- 0}


      for (i in 2:9) {
        n = 11 - i
        gradient5$O3[n] <- iterO3$d[n+1] / iterO3$a[n+1] - iterO3$b[n+1] / iterO3$a[n+1] * gradient5$O3[n+1] - iterO3$c[n+1] / iterO3$a[n+1] * gradient5$O3[n+2]
        if (gradient5$O3[n] < 0) {gradient5$O3[n] <- 0}
      }

      gradient5$O3[1] <- iterO3$d[1] / iterO3$b[1] - iterO3$c[1] / iterO3$b[1] * gradient5$O3[2]
      if (gradient5$O3[1] < 0) {gradient5$O3[1] <- 0}

      #converge Grad5 windward

      errorsum0 <- 0
      for (i in 1:10) {
        errorsum0 <- errorsum0 + abs(gradient5$NO2[i] - gradient5$tempo[i]) / 2.59378E+16 *1.91

      }
      y <- NULL
      y <- gradient5$NO2[1] / 2.59378E+16
      errorsum5 <- cbind(errorsum0,y)
      errorout5 <- rbind(errorout5, errorsum5)
      errorsum <- errorsum0
    }
    grad5out <-data.frame(errorout5)
    grad5out <- tail(grad5out,1)
    highNO2ppb <- grad5out[1,2]
    HighNO2 <- highNO2ppb*1.91

    RoadNO2 <- (HighNO2 - backNO2) / 2
    receptors$RoadNO2[g] <- round(RoadNO2,2)
  }

  ####1.2 End of that#######
  receptors <- receptors[,c(1,2,3,4,12)]
  return(receptors)} else
  {message("ERROR: There is a problem with column names in the Receptor dataset \n ensure Recepts dataframe includes 'ID', 'RoadNOx', 'No2Backgd' and 'LA' column names")}


}

