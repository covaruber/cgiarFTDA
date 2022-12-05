simulator <- function(fixed, data, cropName, percentNA=0.1){
  # ------------------------------------------------------------------------- #
  # This function was generated as part of the fully-automated wheat pipeline #
  # Developer: Giovanny Covarrubias-Pazaran 
  # ------------------------------------------------------------------------- #
  if (missing(data)) {
    data <- environment(fixed)
    # data2 <- environment(random)
    nodata <- TRUE
  }else {
    nodata = FALSE
  }
  if (!inherits(fixed, "formula")) 
    stop("\nfixed must be a formula")
  if (length(fixed) != 3) 
    stop("\nFixed model formula must be of the form \"resp ~ pred\"")
  
  
  expi <- function(j) {
    gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", 
                                                j))[[1]])
  }
  expi2 <- function(x) {
    gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl = T)
  }
  #ginvcheck <- grep("ginv\\(", random)
  
  data2 <- data
  if (length(dim(data2)[2]) == 0) {
    stop("Please provide the 'data' argument.\n", call. = FALSE)
  }
  
  mf <- try(model.frame(fixed, data = data2, na.action = na.pass), 
            silent = TRUE)
  mfna <- try(model.frame(fixed, data = data, na.action = na.pass), 
              silent = TRUE)
  if (inherits(mf,"try-error")) {
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", 
         call. = FALSE)
  }
  # mf <- eval(mf, parent.frame())
  # mfna <- eval(mfna, parent.frame())
  # yvar <- model.response(mfna)
  # X <- model.matrix(fixed, mf)
  # 
  # mf[1:5,1:5]
  if(cropName == "Cotton") {
    traitsMetadata <- NULL
  } else if(cropName == "OSR") {
    traitsMetadata <- NULL
  } else if(cropName == "Wheat") {
    traitsMetadata <- NULL
  } else if(cropName == "Soybean") {
    traitsMetadata <- NULL
  }
  
  ## define the covariates
  xvar.names <- gsub(" ", "", strsplit(as.character(fixed[3]), 
                                       split = "[+]")[[1]])
  xvar.names
  ## define the responses
  resp <- expi2(as.character(fixed[2]))
  yvar.names <- gsub(" ", "", strsplit(resp, 
                                       split = "[,]")[[1]])
  yvar.names
  
  #intersect(yvar.names,colnames(data2))
  toquitlist <- list()
  ## add effects
  for(i in 1:length(yvar.names)){ # for each trait
    itrait <- yvar.names[i]
    print(itrait)
    limos <- traitsMetadata@AXIS[, itrait] #limits of the value
    meanlimos <- mean(limos)
    (unitto <- abs(max(limos) - min(limos))/length(xvar.names)) # maximum effect  to add to each "x" term
    
    efflist <- list()
    toquit <- character()
    for(j in 1:length(xvar.names)){
      
      intercheck <- grep(":",xvar.names[j])#check to see if is an interaction term
      
      if(length(intercheck)>0){
        varos <- strsplit(xvar.names[j],":")[[1]]
        
        data2[,paste(varos,collapse = ".")] <- apply(data2[,varos],1, function(x){paste(x,collapse = ".")})
        xvar.names[j] <- paste(varos,collapse = ".")
        toquit[j] <- paste(varos,collapse = ".")
        
        use <- data.frame(x=data2[,xvar.names[j]])
      }else{
        use <- data.frame(x=data2[,xvar.names[j]])
      }
      
      
      
      chelev <- length(unique(use))#check the number of levels existing
      if(chelev > 1){ # if more than 1 level in x.var
        Xj <- overlay(use)
        unitto2 <- unitto/dim(Xj)[2] # # maximum effect  to add to each level of the "x" term
        
        # create a matrix of effects
        dummy <- apply(as.data.frame(unitto2 * seq(1,dim(Xj)[2],1)),1,function(x){
          rnorm(dim(Xj)[1],mean=x, sd=unitto2/4)
        })
        Xj2 <- (dummy)*Xj
        
        efflist[[j]] <- as.vector(apply(Xj2,1,sum))
      }else{
        efflist[[j]] <- rnorm(dim(use)[1],mean=unitto, sd=unitto/4)
      }
      
    }## end of for each 'x" term
    presp <- apply(do.call(cbind,efflist),1,sum) # provisional response
    #print(i)
    #print(max(abs(presp[which(presp > limos[2])] - limos[2]), na.rm=TRUE) + .1)
    
    outof <- which(presp > limos[2]) # out of imits
    if(length(outof)>0){ # if there was some values out of the limits correct
      tosubstract <- max(abs(presp[outof] - limos[2]), na.rm=TRUE) + .1
    }else{tosubstract <- 0} #else don't
    
    toadd <- presp - tosubstract
    ## only add information to missing cells
    miso <- which(is.na(data2[,itrait]))
    if(length(miso)>0){
      data2[miso,itrait] <-  toadd[miso]
    }else{
      data2[,itrait] <-  toadd
    }
    
    ## add some missing data
    if(percentNA > 0){
      nana <- round(length(data2[,itrait])*percentNA) #number of dta points to delete
      todelete <- sample(1:length(data2[,itrait]),nana)
      data2[todelete,itrait] <- NA
    }
    toquitlist[[i]] <- toquit
  }## end of for each trait
  remo <- unique(unlist(toquitlist)) # remove this created columns by us
  bad <- which(colnames(data2) %in% remo)
  if(length(bad) > 0){
    data2 <- data2[,-bad]
  }
  
  return(data2)
}