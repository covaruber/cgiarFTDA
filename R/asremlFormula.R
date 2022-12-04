asremlFormula <- function(fixed, random, rcov, dat, minRandomLevels=NULL, minResidualLevels=NULL, exchangeRandomEffects=NULL, 
                          exchangeResidualEffects=NULL,customRandomLevels=NULL, customResidualLevels=NULL,
                          xCoordinate= "ROWf",yCoordinate ="RANGEf", doubleConstraintRandom=c("ROWf","RANGEf"),
                          verbose=FALSE){
  
  # ------------------------------------------------------------------------- #
  # This function was generated as part of the fully-automated wheat pipeline #
  # Developer: Giovanny Covarrubias-Pazaran 
  # ------------------------------------------------------------------------- #
  if (missing(dat)) {
    dat <- environment(fixed)
    dat2 <- environment(random)
    nodat <- TRUE
  }else {
    nodat = FALSE
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
  
  struct.detect <- function(z){
    res1 <- apply(data.frame(z),1, function(x){ 
      step1 <- strsplit(x,split=":")[[1]]
      structs <- apply(data.frame(step1),1,function(y){regmatches(y, regexpr("^[^\\(]+", y))[[1]]})
      #apply(data.frame(step1),1,function(x){})
      step1 <- expi2(step1) # remove parenthesis
      step1[which(step1 == "")] <- structs[which(step1 == "")]
      structs[which(structs == step1)] <- ""
      step2 <- paste(step1,collapse = ":")
      return(list(var=step2, structs=structs))
    })
    newvars <- unlist(lapply(res1,function(x){x$var}))
    newstructs <- lapply(res1, function(x){x$structs}); names(newstructs) <- newvars
    return(list(newvars=newvars, newstructs=newstructs))
  }
  
  ##====================================================
  ## check that model specified can be matched with data
  
  dat2 <- dat
  if (length(dim(dat2)[2]) == 0) {
    stop("Please provide the 'dat' argument.\n", call. = FALSE)
  }
  
  mf <- try(model.frame(fixed, dat = dat2, na.action = na.pass), 
            silent = TRUE)
  mfna <- try(model.frame(fixed, dat = dat, na.action = na.pass), 
              silent = TRUE)
  if (class(mf) == "try-error") {
    stop("Please provide the 'dat' argument for your specified variables.\nYou may be specifying some variables in your model not present in your datset.", 
         call. = FALSE)
  }
  
  ##====================================================
  ## define the fixed covariates
  xvar.names <- gsub(" ", "", strsplit(as.character(fixed[3]), 
                                       split = "[+]")[[1]])
  prov <- struct.detect(xvar.names)
  xvar.names.str <- prov$newstructs
  xvar.names <- prov$newvars
  xvar.names
  ## define the random variables
  zvar.names <- gsub(" ", "", strsplit(as.character(random[2]), 
                                       split = "[+]")[[1]])
  prov <- struct.detect(zvar.names)
  zvar.names.str <- prov$newstructs
  zvar.names <- prov$newvars
  zvar.names
  ## define residuals variables
  rvar.names <- gsub(" ", "", strsplit(as.character(rcov[2]), 
                                       split = "[+]")[[1]])
  prov <- struct.detect(rvar.names)
  rvar.names.str <- prov$newstructs
  rvar.names <- prov$newvars
  
  ## define the responses
  resp <- as.character(fixed[2])
  yvar.names <- gsub(" ", "", strsplit(resp, 
                                       split = "[,]")[[1]])
  yvar.names
  
  perc.na <- function(x){  length(which(is.na(x)))/length(x)}
  
  
  ##############
  ## prepare the formulation for fixed and random
  terms <- c(strsplit(xvar.names,split=":"),strsplit(zvar.names,split=":") )
  terms <- lapply(terms, function(x){as.data.frame(matrix(x,nrow=1))})
  terms22 <- do.call(plyr::rbind.fill,terms)
  terms22 <- terms22[with(terms22, order(V1)), ] # now terms are grouped and ready to be tested
  terms22 <- terms22[which(terms22[,1] != "1"),] # remove intercept
  terms22
  
  picking <- function(terms2,nc=NULL){
    ## terms 2 is a dataframe with all the terms to be searched
    # > terms2
    #          V1       V2    V3
    #   2 ENVDATE    BLOCK  <NA>
    #   3 ENVDATE     Rowf  <NA>
    #   4 ENVDATE   Rangef  <NA>
    #   5 ENVDATE EVENTZYG  <NA>
    #   6     ENV     DATE BLOCK
    
    ## numeric.constraint is a list specifying the name of the term and a numeric value of the minimum expected
    ## otherwise is set to zero and therefore not picked for that specific level
    ## list(Rowf=8, Rangef=10)
    
    terms3 <- terms2
    groups <- rep(NA,nrow(terms3))
    ##############
    ## find the groups for making the at levels
    for(w in ncol(terms2):1){
      start <- which(!is.na(terms3[,w]))
      prov <- terms3[start,]
      terms3 <- terms3[-start,]
      groups[start] <-c(paste(letters[w],1:1000000))[as.factor(apply(data.frame(prov[,1:(w-1)]),1,function(x){paste(x,collapse=" ")}))]
    }
    
    (unterms <- as.character(unique(terms2[,"V1"])))
    ugroups <- unique(groups)
    ## now run the algorithm to make the at 
    picklist.u <- list()
    counters <-0
    for(u in ugroups){ # u <- ugroups[1] ## for each major term
      counters <- counters+1
      subterms2 <- terms2[which(groups %in% u),]
      subterms2 <- as.data.frame(subterms2[,which(apply(subterms2,2,perc.na) < 1)])
      
      dat2 <- dat
      newnames <- character()
      
      (colss <- as.character(unlist(subterms2[1,1:(ncol(subterms2)-1)]))) # columns to be pasted 
      (newname <- paste(colss, collapse = "."))
      dat2[,newname] <- apply( data.frame(dat2[ , colss ]) , 1 , paste , collapse = "####" )
      dat2 <- as.data.frame(dat2)
      
      (extra <- as.character(subterms2[,ncol(subterms2)]))
      
      missing <- setdiff(c(newname,extra),colnames(dat2))
      if(length(missing) >0 ){stop(paste("Columns",paste(missing, collapse=", "),"are missing"), call. = FALSE)}
      (nblo0 <- plyr::ddply(droplevels(dat2[,c(newname,extra)]), as.formula(paste("~",newname)), plyr::colwise(function(x){length(which(table(na.omit(x)) > 0))})))
      
      if(verbose){
        print("Original")
        print(nblo0)
      }
      ######################################################################################
      #### from the very beggining whatever doesn't have more than one level cannot be fitted
      #### by applying a double-effect constraint
      if(!is.null(doubleConstraintRandom)){
        for(ri in 1:nrow(nblo0)){
          oky <- which(names(minRandomLevels) %in% doubleConstraintRandom)
          if(length(oky) > 0){mino <- min(unlist(minRandomLevels[oky]))}else{mino <- 1}
          
          check.issue1 <- which(nblo0[ri,doubleConstraintRandom] < mino) # chek if one of the double constrained has only one level
          if(length(check.issue1) > 0){#if exist one that has a single level given that they are intrinsecally connected we have to delete both
            nblo0[ri,doubleConstraintRandom] <- NA
          }
        }
      }
      nblo0
      ######################################################################################
      
      
      
      if(!is.null(minRandomLevels)){
        for(tt in 1:length(minRandomLevels)){ # tt <- 2
          constrained <- names(minRandomLevels)[tt] # term
          constrain <- minRandomLevels[[tt]] # number of minimum levels
          matching <- which(colnames(nblo0) %in% constrained) # column in nblo0 matching i.e. HARVESTBLOCK
          
          if(length(matching)>0){ # if we have a constrained
            
            toconst <- which(nblo0[,constrained] < constrain) # which fields to contrain
            
            if(!is.null(exchangeRandomEffects)){ # if user has some specific oposite relationship among RE, i.e. row-range apply it
              exxx <- which(names(exchangeRandomEffects) %in% constrained)
              if(length(exxx)){
                constrained <- exchangeRandomEffects[[exxx]]
              }
            }
            nblo0[toconst,constrained] <- NA
          }
        }
      }
      ## now location by location decide what will you fit
      
      (step1 <- lapply(split(data.frame(nblo0[,-1]),seq(NROW(data.frame(nblo0[,-1])))),function(x,y){
        
        y[which(!duplicated(x) & x != 0)]
        
      },y=names(nblo0)[-1]) )
      
      #(step1 <- apply(data.frame(nblo0[,-1]),1,function(x,y){y[which(!duplicated(x) & x != 0)]}, y=names(nblo0)[-1]))
      
      #uinuse <- unterms[which(ugroups %in% u)]
      
      
      if(verbose){print("Transformed")}
      gggg <- nblo0
      gggg[,-1] <- t(apply(data.frame(gggg),1,function(x,y){
        
        levused <- x[1] # level used
        x <- as.numeric(unlist(x[-1])) # new x just numbers
        
        if(perc.na(x) < 1){
          onedup <- which(duplicated(x) | x == 0) # found a duplicated in that level of the random effect
          if(length(onedup)>0){
            #subidat <- dat[which(dat[,y[1]] == levused),] # data only for that levelof that random effect
            dups <- which(x == x[onedup]) + 1 # factors duplicated, plus one because we removed the first column
            namedups <- y[dups]
            # loop to find if indeed they are the same plots or are different plots and just by chance have the same number
            coorcheck <- length(which(c(xCoordinate,yCoordinate) %in% namedups)) # if coordinates are matched these are ignored
            if(coorcheck == 2){
              return(x)
            }else{
              x[which(duplicated(x) | x == 0)] <- NA; return(x)
            }
            
          }else{
            x[which(duplicated(x) | x == 0)] <- NA; return(x)
          }
        }else{return(x)}#end of if perc.na < 1
        
      }, y=colnames(gggg)))
      
      if(verbose){print(gggg)}
      
      (step1 <- lapply(split(data.frame(gggg[,-1]),seq(NROW(data.frame(gggg[,-1])))),function(x,y){y[which(!duplicated(x) & x != 0)]},y=names(gggg)[-1]) )
      
      
      names(step1) <- nblo0[,1]
      picklist <- list()
      for(kk in names(nblo0)[-1]){ # kk <- (names(nblo0)[-1])[1]
        pick <- unlist(lapply(step1,function(x){if(kk%in%x){return(TRUE)}else{return(FALSE)}}))
        picklist[[kk]] <- names(pick[which(pick)])
      }
      
      termspasted <- as.character(unlist(subterms2[1,-ncol(subterms2)]))
      
      ## each of the terms nwithin the pcik list is confounded by the #### we need to do an lapply for each of the terms
      #names(picklist) <- paste(paste(termspasted,collapse=":"),names(picklist), sep=":")
      
      ## 1 picklist may convert in a huge list depending the number of terms
      picklist.upper <- list()
      for(nn in 1:length(termspasted)){ # nn <- 1
        picklist.upper[[nn]] <- lapply(picklist, function(x){ # for each element of picklist
          
          if(length(x) > 0){
            rec2 <- unique(apply(data.frame(x),1,function(y){ # for each element of x
              if(length(y) > 0){
                rec <- strsplit(y,"####")[[1]][nn]
              }else{
                rec <- character()
              }
              return(rec)
            }))
          }else{
            rec2 <- character()
          }
          
          return(rec2)
          
        })
        names(picklist.upper)[nn] <- termspasted[nn]
      }
      
      picklist.u <- c(picklist.u,picklist.upper)
    }
    return(picklist.u)
  }
  
  picklist.u <-  picking(terms2 = terms22)
  picklist.u <- lapply(picklist.u,function(x){lapply(x,function(x2){na.omit(x2)})})
  ##############
  ## prepare the formulation for residuals
  
  termsR <- c(strsplit(rvar.names,split=":") )[[1]]
  rvar.names.str <- rvar.names.str[[1]] # reduce the list to a vector
  
  atfactor <- termsR[grep("at",rvar.names.str)]
  datsplit <- split(dat, dat[,atfactor])
  termsRnoat <- setdiff(termsR,atfactor)
  termsRnoat.str <- rvar.names.str[-grep("at",rvar.names.str)]
  
  termlist <- lapply(vector(mode="list",length(termsRnoat)), function(x,nn){res <- vector(mode="list",nn); names(res) <- c("met","nomet"); return(res)}, nn=2) #list of lists
  names(termlist) <- termsRnoat
  for(k in 1:length(datsplit)){ # for each level of at(), i.e. ENV levels
    datprov <- droplevels(datsplit[[k]])
    
    for(o in 1:length(termsRnoat)){ # for each other factos in the rsidual formula
      
      missing <- setdiff(termsRnoat[o],colnames(datprov))
      if(length(missing) >0 ){stop(paste("Columns",paste(missing, collapse=", "),"are missing"), call. = FALSE)}
      
      if(!is.null(exchangeResidualEffects)){ # if user specified opposite relationship
        condition <- termsRnoat[o] %in% names(exchangeResidualEffects) 
        if(condition){totable <- exchangeResidualEffects[[termsRnoat[o]]]}else{
          totable <- termsRnoat[o]
        } # if such condition exist
      }else{totable <- termsRnoat[o]}
      
      
      tabo <- table(na.omit(datprov[,totable])) # table to check condition
      ifcondexist <- termsRnoat[o] %in% names(minResidualLevels) # if condition was provided
      if(ifcondexist){
        bad <- length(tabo) < minResidualLevels[[termsRnoat[o]]] # is there any that doesn't fill the condition?
        if(bad){
          FALSE
        }else{ 
          termlist[[termsRnoat[o]]][["met"]] <- c(termlist[[termsRnoat[o]]][["met"]], names(datsplit)[k])
          
        }
      }else{
        termlist[[termsRnoat[o]]][["met"]] <- c(termlist[[termsRnoat[o]]][["met"]], names(datsplit)[k])
      }#$$
      
    }# end of: for each other factos in the rsidual formula
  } # end of: for each level of at()
  ## now filll the locations not met
  for(i in 1:length(termlist)){
    termlist[[i]][["nomet"]] <- setdiff(names(datsplit),termlist[[i]][["met"]])
  }
  
  step1 <- lapply(as.list(names(termlist)),function(x){c(paste(x,"met",sep="####"), paste(x,"nomet",sep="####"))})
  names(step1) <- names(termlist)
  gridrcov <- expand.grid(step1)
  
  ## for each row we have all potential residuals to incorporate
  rcov.formula <- character(); counter <- 0
  picklist.r <- list()
  for(u in 1:nrow(gridrcov)){
    ss1 <- gridrcov[u,]
    ss1list <- list()
    ss1formula <- character()
    for(u2 in 1:length(ss1)){ # for each random term
      ss2 <- strsplit(as.character(unlist(ss1[u2])), "####")[[1]]
      levo <- termlist[[ss2[1]]][[ss2[2]]] # levels to extract
      if(is.null(levo)){ss1list[[u2]] <- character()}else{ss1list[[u2]] <- levo}
      metyesno <- grep("nomet",ss2) # asses if this structure is the one that met the condition or not
      if(length(metyesno)==0){ss1formula[[u2]] <- paste(termsRnoat.str[u2],"(",ss2[1],")",sep="")}else{
        ss1formula[[u2]] <- paste("id","(",ss2[1],")",sep="")
      }
    }
    atss <- Reduce(intersect, ss1list) # atts levels for the u row of the gridrcov
    if(!is.null(customResidualLevels)){
      if(length(customResidualLevels$picklist.r[[u]])>0){atss <- customResidualLevels$picklist.r[[u]]}
    }
    picklist.r[[u]] <- atss
    if(length(atss)>0){ # makes sense to include in the model
      counter <- counter + 1
      atss <- paste("'",atss,"'", sep="")
      atss <- paste(atss,collapse ="," )
      atss2 <- paste("at(",atfactor,",c(",atss,") )")
      rcov.formula[counter] <-paste(c(atss2,ss1formula),collapse = ":")
    } # else there's no level where all this conditions met
  }
  
  ##################################
  ## end of best combinations
  ##################################
  
  if(!is.null(customRandomLevels)){ # find the best model
    picklist.u <- customRandomLevels
  }
  
  ## now get the formula
  
  xvar.names2 <- c(strsplit(xvar.names,split=":") )
  zvar.names2 <- c(strsplit(zvar.names,split=":") )
  
  ## =====================================
  ## make the fixed term
  fixedlist <- character()
  for(h in 1:length(xvar.names2)){
    prov <- xvar.names2[[h]]
    if(prov == "1"){fixedlist[h] <- "1"}else{
      important <- prov[length(prov)]
      prov <- prov[-length(prov)]
      levelss0 <- character()
      prov.str <- xvar.names.str[[h]]
      for(kk in 1:length(prov)){
        levelss <- picklist.u[[prov[kk]]][[important]]
        if(prov.str[kk] == ""){ # there's a structure
          levelss0[kk] <- prov[kk]
        }else if(prov.str[kk] == "at"){
          levelss0[kk] <- paste(paste(prov.str[kk],"(",sep=""),prov[kk],", c(", paste(paste("'",levelss,"'",sep=""),collapse = ","),") )")
        }else if(prov.str[kk] == "diag" | prov.str[kk] == "us" | prov.str[kk] == "ar1"){
          levelss0[kk] <- paste(paste(prov.str[kk],"(",sep=""),prov[kk],")")
        }
      }
      fixedlist[h] <- paste(c(levelss0,important),collapse = ":")
    }
  }
  fxed <- paste(yvar.names,"~",paste(fixedlist, collapse = " + "))
  
  ## =====================================
  ## make the random part
  randomlist <- character()
  for(h in 1:length(zvar.names2)){
    prov <- zvar.names2[[h]]
    important <- prov[length(prov)]
    prov <- prov[-length(prov)]
    levelss0 <- character()
    prov.str <- zvar.names.str[[h]]
    
    for(kk in 1:length(prov)){
      levelss <- picklist.u[[prov[kk]]][[important]]
      if(prov.str[kk] == ""){ # there's a structure
        if(length(levelss) >0){ # if there's at least one level fit it
          levelss0[kk] <- prov[kk]
        }
      }else if(prov.str[kk] == "at"){
        if(length(levelss) >0){ # if there's at least one level fit it
          levelss0[kk] <- paste(paste(prov.str[kk],"(",sep=""),prov[kk],", c(", paste(paste("'",levelss,"'",sep=""),collapse = ","),") )")
        }
      }else if(prov.str[kk] == "diag" | prov.str[kk] == "us" | prov.str[kk] == "ar1"){
        if(length(levelss) >0){ # if there's at least one level fit it
          levelss0[kk] <- paste(paste(prov.str[kk],"(",sep=""),prov[kk],")")
        }
      }
    }
    
    if(length(levelss) >0){
      randomlist[h] <- paste(c(levelss0,important),collapse = ":")
    }
    
  }
  present <- which(!is.na(randomlist)) # which random effects did have at least one level to fit according to constraints
  randomlist <- randomlist[present]
  rndom <- paste(paste(randomlist, collapse = " + "))
  
  
  ## =====================================
  ## make the residual
  names(picklist.r) <- paste("row_",1:length(picklist.r),"_gridrcov",sep="")
  
  rsidual <- paste(paste(rcov.formula, collapse = " + "))
  
  return(list(fixed=fxed, random=rndom, rcov=rsidual, used=picklist.u, used.res=list(gridrcov=gridrcov,picklist.r=picklist.r)))
  
}
