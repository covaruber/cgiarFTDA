geneticFormulation <- function(mymodel,noFa=c(1,1),overlayColumnName="half", materialNameColumnName="material_id",
                               femaleMaterialNameColumnName="female_material_id",maleMaterialNameColumnName="male_material_id",
                               pedargName="pedarg"){
  '%!in%' <- function(x,y)!('%in%'(x,y))  
  if(missing(mymodel)){stop("Please provide the model desired",call. = FALSE)}
  if(length(mymodel) != 4){stop("Please provide the 4 terms expected",call. = FALSE)}
  if(mymodel[1] %!in% c("rr","dg","cs")){stop("Please provide one of the three models available; 'rr', 'dg', 'cs' ",call. = FALSE)}
  if(mymodel[2] %!in% c("sib","prt")){stop("Please provide one of the two models available; 'sib', 'prt' ",call. = FALSE)}
  if(mymodel[3] %!in% c("vm","ped","na")){stop("Please provide one of the three models available; 'vm', 'ped', 'na' ",call. = FALSE)}
  if(mymodel[4] %!in% c("ide","na")){stop("Please provide one of the two models available; 'ide', 'na' ",call. = FALSE)}
  noFa <- as.character(noFa)
  noFa <- replaceValues(noFa, Search = c("0","NA"), Replace = c("",""))
  geneticFormulaVector <- character()
  if(mymodel[3] != "avoid"){
    materialStart=list(sib=c(materialNameColumnName),prt=c(femaleMaterialNameColumnName,maleMaterialNameColumnName))
    interactionsStart=list(rr=c("rr(","diag("),cs=c("",""),dg=c("diag("))
    interactionsEnd=list(rr=c(")"),cs=c(""),dg=c(")"))
    
    pedigreeStart <- list()
    pedigreeStart[[mymodel[3]]] <- paste0(mymodel[3],"(")
    pedigreeStart[["na"]] <- ""
    
    if(is.na(pedargName)){mp <- ")"}else{mp <- paste0(", source=",pedargName,")")}
    
    pedigreeEnd <- list()
    pedigreeEnd[[mymodel[3]]] <- mp
    pedigreeEnd[["na"]] <- ""
    
    geneticFormulaBase <- expand.grid(interactionsStart[[mymodel[1]]], materialStart[[mymodel[2]]])
    geneticFormulaBase$Var3 <- replaceValues(geneticFormulaBase$Var1,Search = c("rr(","diag(",materialStart[[mymodel[2]]]),Replace=c(paste0(",",noFa[1]),"",rep("",length(materialStart[[mymodel[2]]]))))
    geneticFormulaBase <- data.frame(V1=geneticFormulaBase[,1],
                                     V2="field_book_id", 
                                     V3=geneticFormulaBase$Var3,
                                     V4=interactionsEnd[[mymodel[1]]],
                                     V5=":",
                                     V6=pedigreeStart[[mymodel[3]]],
                                     V7=geneticFormulaBase[,2],
                                     V8=pedigreeEnd[[mymodel[3]]])
    if(mymodel[2] == "prt"){
      geneticFormulaBase$V9 <- paste0(":",overlayColumnName)
    }
    lastColumnName=paste0("V", ncol(geneticFormulaBase)+1)
    geneticFormulaBase <- geneticFormulaBase[with(geneticFormulaBase, order(V1)), ]
    geneticFormulaBase[,1] <- as.character(geneticFormulaBase[,1])
    if(mymodel[1] == "cs"){ # we add main term if cs
      geneticFormulaBase2 <- geneticFormulaBase
      geneticFormulaBase2[,1] <- " "; geneticFormulaBase2[,2:5] <- ""
      geneticFormulaBase <- rbind(geneticFormulaBase,geneticFormulaBase2)
    }
    geneticFormulaBase[,lastColumnName] <- apply(geneticFormulaBase,1,function(x){paste(x,collapse = "")})
    geneticFormulaBase <- geneticFormulaBase[which(!duplicated(geneticFormulaBase[,lastColumnName])),]
    
    geneticFormulaTerms <- character()
    for(v1 in unique(geneticFormulaBase$V1)){
      point0 <- gsub(paste0(":",overlayColumnName),"",(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName])[1])
      if(mymodel[2] == "prt"){
        geneticFormulaTerms[v1] <- paste0("str(~",paste(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName], collapse = " + and("),")",",~", point0 ,")")
      }else{
        geneticFormulaTerms[v1] <- paste(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName], collapse = " + ")
      }
    }
    geneticFormulaVector[1] <- paste(geneticFormulaTerms, collapse = " + ")
  }
  if(mymodel[4] != "na"){
    materialStart=list(sib=c(materialNameColumnName),prt=c(femaleMaterialNameColumnName,maleMaterialNameColumnName))
    interactionsStart=list(rr=c("rr(","diag("),cs=c("",""),dg=c("diag("))
    interactionsEnd=list(rr=c(")"),cs=c(""),dg=c(")"))
    
    pedigreeStart <- list()
    pedigreeStart[[mymodel[4]]] <- paste0(mymodel[4],"(")
    pedigreeStart[["na"]] <- ""
    pedigreeEnd <- list()
    pedigreeEnd[[mymodel[4]]] <- ")"
    pedigreeEnd[["na"]] <- ""
    
    geneticFormulaBase <- expand.grid(interactionsStart[[mymodel[1]]], materialStart[[mymodel[2]]])
    geneticFormulaBase$Var3 <- replaceValues(geneticFormulaBase$Var1,Search = c("rr(","diag("),Replace=c(paste0(",",noFa[2]),""))
    geneticFormulaBase <- data.frame(V1=geneticFormulaBase[,1],
                                     V2="field_book_id", 
                                     V3=geneticFormulaBase$Var3,
                                     V4=interactionsEnd[[mymodel[1]]],
                                     V5=":",
                                     V6=pedigreeStart[[mymodel[4]]],
                                     V7=geneticFormulaBase[,2],
                                     V8=pedigreeEnd[[mymodel[4]]])
    if(mymodel[2] == "prt"){
      geneticFormulaBase$V9 <- paste0(":",overlayColumnName)
    }
    lastColumnName=paste0("V", ncol(geneticFormulaBase)+1)
    geneticFormulaBase <- geneticFormulaBase[with(geneticFormulaBase, order(V1)), ]
    
    geneticFormulaBase[,1] <- as.character(geneticFormulaBase[,1])
    if(mymodel[1] == "cs"){ # we add main term if cs
      geneticFormulaBase2 <- geneticFormulaBase
      geneticFormulaBase2[,1] <- " "; geneticFormulaBase2[,2:5] <- ""
      geneticFormulaBase <- rbind(geneticFormulaBase,geneticFormulaBase2)
    }
    geneticFormulaBase[,lastColumnName] <- apply(geneticFormulaBase,1,function(x){paste(x,collapse = "")})
    geneticFormulaBase <- geneticFormulaBase[which(!duplicated(geneticFormulaBase[,lastColumnName])),]
    
    geneticFormulaTerms <- character()
    for(v1 in unique(geneticFormulaBase$V1)){
      point0 <- gsub(paste0(":",overlayColumnName),"",(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName])[1])
      if(mymodel[2] == "prt"){
        geneticFormulaTerms[v1] <- paste0("str(~",paste(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName], collapse = " + and("),")",",~", point0 ,")")
      }else{
        geneticFormulaTerms[v1] <- paste(geneticFormulaBase[which(geneticFormulaBase$V1 == v1),lastColumnName], collapse = " + ")
      }
    }
    
    geneticFormulaVector[2] <- paste(geneticFormulaTerms, collapse = " + ")
  }
  
  geneticFormulaResult <- paste("~",paste(geneticFormulaVector, collapse=" + "))
  
  return(geneticFormulaResult)
}