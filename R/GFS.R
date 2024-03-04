# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



#' GFS sampling
#'
#' Take an unequal probability sample based on GFS method
#'
#' @param pik a vector of inclusion probabilities
#' @param v GCD of all the inclusion probabilities

####################################################################
####################################################################
################ External Functions ################################
####################################################################
####################################################################
library(sampling)
library(bigmemory)
library(ggplot2)
library(ggpubr)
library(intervals)


### decimal function
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

####################################################################
####################################################################
################ Original GFS       ################################
####################################################################
####################################################################

GFSfixedsize<-function( pik              #vector of inclusion probabilities
                       ,v  = NULL  #size of bars
                       ,chaos   = NULL  #increasing chaosy, leads to decreasing zero second orders
                       ,M.chaos  = NULL  #number of times that we exchange the bars for increasing chaosy
                       ,secondor  = NULL  #1:shows results of second order
                       ,plot      = NULL  #1:plot
                       ,counter   = NULL  #1: active
                       ,round     = NULL  #1: round to decimal places of v
)

{
  start.time = Sys.time()
  if(is.null(v))  v<- 0.01
  if(is.null(chaos))  chaos  <- 1
  if(is.null(M.chaos)) M.chaos <- 10*length(pik)*(1/v)

  GFSoptimal(
    X          = runif(100, 10, 50)
    ,n         = 30
    ,v         = 0.01
    ,Lambda     = 5000
    ,loss.f    = 1
    ,M.chaos  = 100
    ,secondor  = 0
    ,plot      = 0
    ,epsilon   = 5
  )
  if(is.null(counter))   counter <- 1
  if(is.null(round))     round   <- 0
  if(is.null(plot))      plot    <- 0
  if(is.null(secondor))  secondor<- 1

   if(counter == 1){
    if(secondor == 0) cat('After finishing the progress bar with "111...111"100%, the sample will be ready')
    if(secondor == 1) cat('After finishing the progress bar with "222...222"100%, the sample will be ready')
  }

  minmini<-unique(min(pik))
  for(i in 1:10){
    if(round(minmini, digits = i) > 0) {
      vo = min(v, 1/10^i)
      break}
  }
  v <- vo

  N<-length(pik)
  if(round == 1) pik = round(pik, decimalplaces(v))
  pik[pik == 0] <- v
  alpa<-1/(v)#all possible sample
  alind<-array(0,alpa)
  xo<-ceiling(pik/v)
  tre <- 0.0000001
  nn<-sum(pik)
  u<-c(1:N)


  ##################################
  ### arranging the all the bars ###
  ##################################
  options(bigmemory.allow.dimnames=TRUE)
  Ma<-matrix(0,ncol=N,nrow=alpa)
  #Ma<-as.big.matrix(Ma)
  upper = c()
  lower = c()
  group = c()



  lower[1] <- 0
  upper[1] <- pik[1]


  ### Arranging the bars, leading to Madow 1949
  group[1] <- 1
  index <- 1
  gro   <- 1

  if(counter == 1) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = N, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "0")   # Character used to create the bar
  for(i in 2:N){
    if(counter == 1) setTxtProgressBar(pb,i)
    index <- index + 1
    if( (upper[index - 1] + pik[i]) <= 1 ){
      group[index] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- upper[index - 1] + pik[i]
    }else if(upper[index - 1] == 1){
      group[index] <- i
      lower[index] <- 0
      upper[index] <- pik[i]
    } else {
      group[c(index, index + 1)] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- 1
      lower[index + 1] <- 0
      upper[index + 1] <- pik[i] - (1 - lower[index])
      index <- index + 1
    }

  }


  df <- data.frame(
    Units = factor(group),
    lower = lower,
    upper = upper,
    c1    = array(0,length(group)),
    c2    = array(0,length(group)),
    remove= array(0,length(group))
  )

  mat5<-matrix(0,ncol=4,nrow=2*N)
  llo <- seq(0,1,v)
  ldf <- length(df[,1])
  con5<-0

  if(counter == 1) pb2 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                         max = ldf, # Maximum value of the progress bar
                                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                         width = 50,   # Progress bar width. Defaults to getOption("width")
                                         char = "1")   # Character used to create the bar

  for(i in 1:ldf){
    if(counter == 1) setTxtProgressBar(pb2,i)
    for(j in 1:(alpa)){
      if(df$upper[i]+tre>=llo[j+1] & df$lower[i]-tre<=llo[j]) {Ma[j,df$Units[i]]<-1}
      if(df$upper[i]+tre<llo[j+1] & df$upper[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],llo[j],df$upper[i],0.5)
      }
      if(df$lower[i]+tre<llo[j+1] & df$lower[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],df$lower[i],llo[j+1],0.5)
      }

    }
  }

  mat5<-mat5[mat5[,4]>0,]
  mat5 <- data.frame(mat5)

  #######################################################################
  ### Decreasing zero second order inclusions with increasing chaosy ###
  #######################################################################
  while(chaos==1){
    if(counter == 1) pb3 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = M.chaos, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "2")   # Character used to create the bar

    for(e in 1:M.chaos){
      if(counter == 1) setTxtProgressBar(pb3,e)
      xx<-sample(N,2)
      yy<-sample(alpa,2)
      if( Ma[yy[1],xx[1]]==1&Ma[yy[2],xx[1]]==0&Ma[yy[1],xx[2]]==0&Ma[yy[2],xx[2]]==1)
        c(Ma[yy[1],xx[1]]<-0,Ma[yy[2],xx[1]]<-1,Ma[yy[1],xx[2]]<-1,Ma[yy[2],xx[2]]<-0)
    }
    break
  }

  while(chaos==2){
    for(e in 1:M.chaos){
      xx<-sample(N,1)
      yy<-sample(alpa,2)
      if(Ma[yy[1],xx]==1&Ma[yy[2],xx]==0)
        c(Ma[yy[1],xx]<-0,Ma[yy[2],xx]<-1)
    }
    break
  }

  ###############################
  ### Get back to lower upper ###
  ###############################
  df0 <- data.frame(
    "Units" = rep(1:N, each = alpa),
    "lower" = rep(seq(0,1-v,v),N),
    "upper" = rep(seq(v,1,v),N),
    "keep"  = c(Ma)
  )
  names(mat5) <- names(df0)
  df0<-df0[df0$keep==1,]
  df0 <- rbind(df0,mat5)
  df <- df0

  ##########################
  ### Final sample ###
  ##########################
  indi<-runif(1,0,1)
  sam <- as.numeric(df[which(df$lower<=indi&df$upper>=indi),"Units"])
  end.time = Sys.time()

  ######################################################
  ### Calculate second order inclusion probabilities ###
  ######################################################

  if(secondor==2){
    SecInc<-matrix(0,N,N)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = N, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar


    for (ww in 1:N) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:N) {
        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == ww]),as.numeric(df$upper[df$Units == ww])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == vv]),as.numeric(df$upper[df$Units == vv])))
        Inters <- interval_intersection(Inter1, Inter2)
        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }

  }


  ## Sample second order
  if(secondor==1){
    legsam <- length(sam)
    SecInc<-matrix(0,legsam,legsam)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = legsam, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar

    cww<-0
    for (ww in 1:legsam) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:legsam) {

        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == sam[ww]]),as.numeric(df$upper[df$Units == sam[ww]])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == sam[vv]]),as.numeric(df$upper[df$Units == sam[vv]])))
        Inters <- interval_intersection(Inter1, Inter2)

        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }
    colnames(SecInc)<-sam
    row.names(SecInc)<-sam

  }


  sample_size<-array(0,alpa)
  for(o in 1:alpa){
    sample_size[o]<-length(u[Ma[o,]==1])
  }

  ##########################
  ### ready data for plot ###
  ##########################
  #row.names(Ma)<-seq(1-1/alpa,1/alpa-v,-v)
  while (plot==1) {
    if(chaos == 0) GGT <- "Madow version of GFS"
    if(chaos == 2) GGT <- "Poisson Sampling GFS"
    if(chaos == 1) GGT <- "Fixed Size, Reduced Zero SIPs GFS"

    BPlot <- ggplot(df, aes(Units, colours = Units, Inclusion.Probabilities = rep(1,length(df[,1])), color = "darkblue")) +
      geom_linerange(aes(ymin = lower, ymax = upper)) +
      geom_hline(yintercept=indi, linetype="dashed", color = "black")
    plot(BPlot + theme_classic() + labs(x = "Population Units", y = "Inclusion Probabilities") +
           ylim(0, 1) + xlim(1,N) + theme(legend.position = "none") + ggtitle(GGT))
    break
  }

  return(list(Df                      = df,
              Second.Order.Inclusions = if(secondor > 0)  SecInc,
              Final.Sample = sam,
              Total.time = end.time - start.time))
}



####################################################################
####################################################################
################ Optimal GFS        ################################
####################################################################
####################################################################

starttime <- Sys.time()

GFSoptimal<-function(    X                        #vector of inclusion probabilities
                        ,n                        #Sample size
                        ,v         = NULL         #size of bars
                        ,Lambda     = NULL         #Number of tested designs,
                        ,loss.f    = NULL         #loss function, 1: variance, 2:minimax
                        ,M.chaos   = NULL         #number of times that we exchange the bars for increasing chaosy
                        ,secondor  = NULL         #1:shows results of second order
                        ,plot      = NULL         #1:plot
                        ,counter   = NULL         #1: active
                        ,round     = NULL         #1: round to decimal places of v
                        ,epsilon     = NULL         #The shift to use just one variable in using variance of NHT in PPS
){
  start.time = Sys.time()
  if(is.null(v))v    <- 0.01
  if(is.null(Lambda))Lambda          <- 100
  if(is.null(loss.f))loss.f        <- 1
  if(is.null(M.chaos))  M.chaos    <- 10*length(pik)*(1/v)
  if(is.null(counter)) counter     <- 0
  if(is.null(round))   round       <- 0
  if(is.null(plot))    plot        <- 0
  if(is.null(secondor))secondor    <- 1
  if(is.null(epsilon))    epsilon      <- 2


   if(counter == 1){
    if(secondor == 0) cat('After finishing the progress bar with "222...222"100%, the sample will be ready')
    if(secondor == 1) cat('After finishing the progress bar with "333...333"100%, the sample will be ready')
   }
  N <- length(X)
  pik <- inclusionprobabilities(X,n)
  X<-X+epsilon
  X[X<=0] <- v
  minmini<-unique(min(pik))

   for(i in 1:10){
     if(round(minmini, digits = i) > 0) {
      vo = min(v, 1/10^i)
      break}
   }

  v <- vo
  if(round == 1) pik = round(pik, decimalplaces(v))
  pik[pik == 0] <- v
  alpa<-1/(v)#all possible sample
  alind<-array(0,alpa)
  xo<-ceiling(pik/v)
  tre <- 0.0000001
  nn<-sum(pik)
  u<-c(1:N)


  ##################################
  ### arranging the all the bars ###
  ##################################
  options(bigmemory.allow.dimnames=TRUE)
  Ma<-matrix(0,ncol=N,nrow=alpa)
  #Ma<-as.big.matrix(Ma)
  upper = c()
  lower = c()
  group = c()

  lower[1] <- 0
  upper[1] <- pik[1]


  ### Arranging the bars, leading to Madow 1949
  group[1] <- 1
  index <- 1
  gro   <- 1

  if(counter == 1) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = N, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "0")   # Character used to create the bar
  for(i in 2:N){
    if(counter == 1) setTxtProgressBar(pb,i)
    index <- index + 1
    if( (upper[index - 1] + pik[i]) <= 1 ){
      group[index] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- upper[index - 1] + pik[i]
    }else if(upper[index - 1] == 1){
      group[index] <- i
      lower[index] <- 0
      upper[index] <- pik[i]
    } else {
      group[c(index, index + 1)] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- 1
      lower[index + 1] <- 0
      upper[index + 1] <- pik[i] - (1 - lower[index])
      index <- index + 1
    }

  }

  df <- data.frame(
    Units = factor(group),
    lower = lower,
    upper = upper,
    c1    = array(0,length(group)),
    c2    = array(0,length(group)),
    remove= array(0,length(group))
  )



  mat5<-matrix(0,ncol=4,nrow=2*N)
  llo <- seq(0,1,v)
  ldf <- length(df[,1])
  con5<-0

  if(counter == 1) pb2 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                         max = ldf, # Maximum value of the progress bar
                                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                         width = 50,   # Progress bar width. Defaults to getOption("width")
                                         char = "1")   # Character used to create the bar

  for(i in 1:ldf){
    if(counter == 1) setTxtProgressBar(pb2,i)
    for(j in 1:(alpa)){
      if(df$upper[i]+tre>=llo[j+1] & df$lower[i]-tre<=llo[j]) {Ma[j,df$Units[i]]<-1}
      if(df$upper[i]+tre<llo[j+1] & df$upper[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],llo[j],df$upper[i],0.5)
      }
      if(df$lower[i]+tre<llo[j+1] & df$lower[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],df$lower[i],llo[j+1],0.5)
      }

    }
  }

  mat5<-mat5[mat5[,4]>0,]
  mat5 <- data.frame(mat5)

  ##########################
  ### Decreasing zero second order inclusions with increasing chaosy ###
  ##########################

  if(counter == 1) pb3 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                         max = Lambda*M.chaos, # Maximum value of the progress bar
                                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                         width = 50,   # Progress bar width. Defaults to getOption("width")
                                         char = "2")   # Character used to create the bar


  ALpopu <- Ma
  #Alpopu = as.big.matrix(ALpopu)
  ALNHT<-array(0,Lambda)

  for(b in 1:Lambda){
    ### Decreasing zero second order inclusions ###
    for(e in 1:M.chaos){
      if(counter == 1) setTxtProgressBar(pb3,e + (b-1)*M.chaos)

      xx<-sample(N,2)
      yy<-sample(alpa,2)
      if(Ma[yy[1],xx[1]]==1&Ma[yy[2],xx[1]]==0&Ma[yy[1],xx[2]]==0&Ma[yy[2],xx[2]]==1)
        c(Ma[yy[1],xx[1]]<-0,Ma[yy[2],xx[1]]<-1,Ma[yy[1],xx[2]]<-1,Ma[yy[2],xx[2]]<-0)
    }

    #################################
    ### all possible x estimator ####
    ##################################
    NHTM<-t(apply(Ma, 1, function(x) x*X/pik))
    NHT<-apply(NHTM,1,sum)
    if(loss.f == 1){ALNHT[b]<-var(NHT)}else{ALNHT[b]<-max(abs(NHT-sum(X)))}
    if(ALNHT[b] == min(ALNHT[1:b]))  ALpopu <- Ma
  }#close Lambda
  Ma = ALpopu




  ###############################
  ### Get back to lower upper ###
  ###############################
  df0 <- data.frame(
    "Units" = rep(1:N, each = alpa),
    "lower" = rep(seq(0,1-v,v),N),
    "upper" = rep(seq(v,1,v),N),
    "keep"  = c(Ma)
  )
  names(mat5) <- names(df0)
  df0<-df0[df0$keep==1,]
  df0 <- rbind(df0,mat5)
  df <- df0

  ##########################
  ### Final sample ###
  ##########################
  indi<-runif(1,0,1)
  sam <- as.numeric(df[which(df$lower<=indi&df$upper>=indi),"Units"])
  end.time = Sys.time()

  ######################################################
  ### Calculate second order inclusion probabilities ###
  ######################################################

  if(secondor==2){
    SecInc<-matrix(0,N,N)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = N, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar


    for (ww in 1:N) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:N) {
        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == ww]),as.numeric(df$upper[df$Units == ww])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == vv]),as.numeric(df$upper[df$Units == vv])))
        Inters <- interval_intersection(Inter1, Inter2)
        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }

  }


  ## Sample second order
  if(secondor==1){
    legsam <- length(sam)
    SecInc<-matrix(0,legsam,legsam)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = legsam, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar

    cww<-0
    for (ww in 1:legsam) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:legsam) {
        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == sam[ww]]),as.numeric(df$upper[df$Units == sam[ww]])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == sam[vv]]),as.numeric(df$upper[df$Units == sam[vv]])))
        Inters <- interval_intersection(Inter1, Inter2)
        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }
    colnames(SecInc)<-sam
    row.names(SecInc)<-sam

  }


  sample_size<-array(0,alpa)
  for(o in 1:alpa){
    sample_size[o]<-length(u[Ma[o,]==1])
  }

  ##########################
  ### ready data for plot ###
  ##########################
  #row.names(Ma)<-seq(1-1/alpa,1/alpa-v,-v)
  while (plot==1) {
    BPlot <- ggplot(df, aes(Units, Inclusion.Probabilities = rep(1,length(df[,1])), color = "red")) +
      geom_linerange(aes(ymin = lower, ymax = upper)) +
      geom_hline(yintercept=indi, linetype="dashed", color = "black") + theme_classic() + labs(x = "Population Units", y = "Inclusion Probabilities") +
      ylim(0, 1) + theme(legend.position = "none")
    plot(BPlot)
    break

  }

  return(list(Best.Design             = df,
              Second.Order.Inclusions = if(secondor > 0)  SecInc,
              Sample.Size             = length(sam),
              Final.Sample            = sam,
              Total.time              = end.time - start.time))
}




####################################################################
####################################################################
################ Stream   GFS       ################################
####################################################################
####################################################################
starttime <- Sys.time()
GFSstream<-function(    pik        #vector of inclusion probabilities
                        ,windows   #size of the windows
                        ,v = NULL  #size of bars
                        ,chaos   = NULL #increasing chaosy, leads to decreasing zero second orders
                        ,M.chaos  = NULL #number of times that we exchange the bars for increasing chaosy
                        ,secondor = NULL #1:shows results of second order
                        ,plot     = NULL #1:plot
                        ,counter  = NULL #1: active
                        ,round    = NULL #1: round to decimal places of v
)

{
  start.time = Sys.time()
  if(is.null(v)) v    <- 0.01
  if(is.null(chaos))   chaos      <- 1
  if(is.null(M.chaos))  M.chaos     <- 100*length(pik)*(1/v)
  if(is.null(counter))  counter     <- 1
  if(is.null(round))    round       <- 0
  if(is.null(plot))     plot        <- 0
  if(is.null(secondor)) secondor    <- 1


  if(counter == 1){
    if(secondor == 0) cat('After finishing the progress bar with "222...222"100%, the sample will be ready')
    if(secondor == 1) cat('After finishing the progress bar with "333...333"100%, the sample will be ready')
  }

  minmini<-unique(min(pik))
  for(i in 1:10){
    if(round(minmini, digits = i) > 0) {
      vo = min(v, 1/10^i)
      break}
  }
  v <- vo

  E<-M.chaos
  H<-windows
  x<-round(pik, decimalplaces(v))
  N<-length(x)
  alpa<-1/(v)#all possible sample
  Ualpa<-c(1:alpa)
  u<-c(1:N)
  nn<-sum(x)
  alind<-array(0,alpa)
  xo<-ceiling(pik/v)
  tre <- 0.0000001

  ########################
  ### cross border units
  ########################
  V0<-0
  cont<-0
  nnH<-nn/H

  lel<-floor(sum(x)/H)
  lel
  if((sum(x)/H)%%1!=0) lel<-lel+1
  abmat<-matrix(0,ncol=lel-1,nrow=3) #Cross-Border units and a and b
  for(j in 1:(N-1)){
    Vm10=V0
    V0=V0+x[j]+0.00000001
    if(((ceiling(V0)-floor(Vm10))==2&(floor(V0)%%H)==0)|(((V0-floor(Vm10))==1)&any((V0%%H)==c(0,H)))){
      cont<-cont+1
      abmat[2,cont]<-ceiling(Vm10)-Vm10
      abmat[3,cont]<-V0-floor(V0)
      abmat[1,cont]<-j
    }

  }
  if(sum(abmat[,length(abmat[1,])]==0)) abmat <- abmat[,-length(abmat[1,])]
  abmat<-cbind(abmat,c(N,x[N],0))
  abmat<-cbind(c(1,0,x[1]),abmat)
  abmat<-round(abmat,2)
  id<-c(abmat[1,])  #cross-borders
  aa<-c(abmat[2,])  #a
  bb<-c(abmat[3,])  #b
  legid<-length(id)

  ##################################
  ### arranging the all the bars ###
  ##################################
  options(bigmemory.allow.dimnames=TRUE)
  Ma<-matrix(0,ncol=N,nrow=alpa)
  #Ma<-as.big.matrix(Ma)
  upper = c()
  lower = c()
  group = c()

  lower[1] <- 0
  upper[1] <- pik[1]


  ### Arranging the bars, leading to Madow 1949
  group[1] <- 1
  index    <- 1
  gro      <- 1

  if(counter == 1) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                        max = N, # Maximum value of the progress bar
                                        style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                        width = 50,   # Progress bar width. Defaults to getOption("width")
                                        char = "0")   # Character used to create the bar
  for(i in 2:N){
    if(counter == 1) setTxtProgressBar(pb,i)
    index <- index + 1
    if( (upper[index - 1] + pik[i]) <= 1 ){
      group[index] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- upper[index - 1] + pik[i]
    }else if(upper[index - 1] == 1){
      group[index] <- i
      lower[index] <- 0
      upper[index] <- pik[i]
    } else {
      group[c(index, index + 1)] <- i
      lower[index] <- upper[index - 1]
      upper[index] <- 1
      lower[index + 1] <- 0
      upper[index + 1] <- pik[i] - (1 - lower[index])
      index <- index + 1
    }

  }

  df <- data.frame(
    Units = factor(group),
    lower = lower,
    upper = upper,
    c1    = array(0,length(group)),
    c2    = array(0,length(group)),
    remove= array(0,length(group))
  )

  mat5<-matrix(0,ncol=4,nrow=2*N)
  llo <- seq(0,1,v)
  ldf <- length(df[,1])
  con5<-0

  if(counter == 1) pb2 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                         max = ldf,    # Maximum value of the progress bar
                                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                         width = 50,   # Progress bar width. Defaults to getOption("width")
                                         char = "1")   # Character used to create the bar

  for(i in 1:ldf){
    if(counter == 1) setTxtProgressBar(pb2,i)
    for(j in 1:(alpa)){
      if(df$upper[i]+tre>=llo[j+1] & df$lower[i]-tre<=llo[j]) {Ma[j,df$Units[i]]<-1}
      if(df$upper[i]+tre<llo[j+1] & df$upper[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],llo[j],df$upper[i],0.5)
      }
      if(df$lower[i]+tre<llo[j+1] & df$lower[i]-tre>llo[j]){
        Ma[j,df$Units[i]]<-.5
        con5 <- con5 + 1
        mat5[con5,] <- c(df$Units[i],df$lower[i],llo[j+1],0.5)
      }


    }
  }
  Ma0 <- list()
  Ma0[[1]] <- cbind(Ma[,c(id[1]:id[2])])

  for(i in 2:(length(id)-1)){

    Ma0[[i]] <- cbind(Ma[,c((id[i]+1):id[i+1])])
  }

  llist<-length(Ma0)

  mat5<-mat5[mat5[,4]>0,]
  mat5 <- data.frame(mat5)

  ##########################
  ### Decreasing zero second order inclusions with increasing chaosy ###
  ##########################



  while(chaos==1){
    if(counter == 1) pb3 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = M.chaos, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "2")   # Character used to create the bar
    for(i in 1:llist){

      for(e in 1:round(M.chaos/llist)){
        if(counter == 1) setTxtProgressBar(pb3,(i-1)*round(M.chaos/llist)+e)
        if(length(Ma0[[i]][1,])<=2) next
        if(i<llist & length(Ma0[[i]][1,])>2)  xx<-sample(length(Ma0[[i]][1,])-1,2)

        if(i==llist)                          xx<-sample(length(Ma0[[i]][1,]),2)

        yy<-sample(length(Ma0[[i]][,1]),2)
        if(Ma0[[i]][yy[1],xx[1]]==1&Ma0[[i]][yy[2],xx[1]]==0&
           Ma0[[i]][yy[1],xx[2]]==0&Ma0[[i]][yy[2],xx[2]]==1)
          c(Ma0[[i]][yy[1],xx[1]]<-0,Ma0[[i]][yy[2],xx[1]]<-1,
            Ma0[[i]][yy[1],xx[2]]<-1,Ma0[[i]][yy[2],xx[2]]<-0)
      }}

    break
  }


  Ma00 <- matrix(unlist(Ma0),ncol=N)
  Ma <- Ma00

  ###############################
  ### Get back to lower upper ###
  ###############################
  df0 <- data.frame(
    "Units" = rep(1:N, each = alpa),
    "lower" = rep(seq(0,1-v,v),N),
    "upper" = rep(seq(v,1,v),N),
    "keep"  = c(Ma)
  )
  df$cros[df$Units %in% id[-c(1,length(id))]] <- 1
  names(mat5) <- names(df0)
  df0<-df0[df0$keep==1,]
  df0 <- rbind(df0,mat5)
  df <- data.frame(
    df0,
    "cros"  = rep(0, length(df0[,1]))
  )
  for(i in 1:(length(id)-1)){
    df$cros[df$Units %in% c(id[i]:id[i+1])]<-i
  }
  df$cros[df$Units %in% id[-c(1,length(id))]] <- 0
  df$cros <- factor(df$cros)
  ##########################
  ### Final sample ###
  ##########################
  indi<-runif(1,0,1)
  sam <- as.numeric(df[which(df$lower<=indi&df$upper>=indi),"Units"])
  end.time = Sys.time()


  ######################################################
  ### Calculate second order inclusion probabilities ###
  ######################################################

  if(secondor==2){
    SecInc<-matrix(0,N,N)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = N, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar


    for (ww in 1:N) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:N) {
        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == ww]),as.numeric(df$upper[df$Units == ww])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == vv]),as.numeric(df$upper[df$Units == vv])))
        Inters <- interval_intersection(Inter1, Inter2)
        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }

  }


  ## Sample second order
  if(secondor==1){
    legsam <- length(sam)
    SecInc<-matrix(0,legsam,legsam)
    if(counter == 1) pb4 <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                           max = legsam, # Maximum value of the progress bar
                                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                           width = 50,   # Progress bar width. Defaults to getOption("width")
                                           char = "3")   # Character used to create the bar

    for (ww in 1:legsam) {
      if(counter == 1) setTxtProgressBar(pb4,ww)
      for (vv in ww:legsam) {

        Inter1 <- Intervals(c(as.numeric(df$lower[df$Units == sam[ww]]),as.numeric(df$upper[df$Units == sam[ww]])))
        Inter2 <- Intervals(c(as.numeric(df$lower[df$Units == sam[vv]]),as.numeric(df$upper[df$Units == sam[vv]])))
        Inters <- interval_intersection(Inter1, Inter2)

        SecInc[ww,vv]<-round(sum(Inters[,2]-Inters[,1]),3)
      }
    }
    colnames(SecInc)<-sam
    row.names(SecInc)<-sam

  }


  sample_size<-array(0,alpa)
  for(o in 1:alpa){
    sample_size[o]<-length(u[Ma[o,]==1])
  }

  ##########################
  ### ready data for plot ###
  ##########################

  while (plot==1) {
    BPlot <- ggplot(df, aes(Units, colours = Units, Inclusion.Probabilities = rep(1,length(df[,1])),group = cros, colour = cros)) +
      geom_linerange(aes(ymin = lower, ymax = upper )) +
      geom_hline(yintercept=indi, linetype="dashed", color = "black")+
      labs(x = "Population Units", y = "Inclusion Probabilities") +
      ylim(0, 1)  + theme_classic() + theme(legend.position = "none") + ggtitle("Stream Sampling")

      plot(BPlot)
      break
  }
  return(list(Final.Popu   = df,
              Second.Order.Inclusions = if(secondor > 0)  SecInc,
              Final.Sample = sam,
              Total.time = end.time - start.time))
}




