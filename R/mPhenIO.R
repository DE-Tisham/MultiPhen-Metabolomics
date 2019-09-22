read.plink<-function(root, indiv = NULL,opts =mPhen.options("geno.input")){
 genoC = mPhen.readGenotypes(root, indiv = indiv, opts = opts)
 return(genoC$genoData)
}

.mergeLists<-function(genoConnection,data){
   nmesg = names(genoConnection)
   nmesd = names(data)
   nmesin = match(nmesd,nmesg)
  
   indna = which(is.na(nmesin))
   nmesin[indna] = (1:length(nmesin[indna]))+length(nmesg)
     for(k in 1:length(nmesin)){
        genoConnection[[nmesin[[k]]]] = data[[k]]
   }
   if(length(indna)>0)names(genoConnection) = c(nmesg, nmesd[indna])
   genoConnection
}

mPhen.readGenotypes<-function(genoConnection, indiv = NULL,opts =mPhen.options("geno.input" )){
  closeConnection = FALSE
  isconnection = is.list(genoConnection)
  if(isconnection) isconnection = !is.null(genoConnection$conn)
  if(!isconnection){
    genoConnection<- list(conn = .openGenoConnection(genoConnection, opts, indiv = indiv))
    samps = attr(genoConnection$conn,"sampleids")
    attr(genoConnection,"rsidToDo") = opts$mPhen.rsidToDo
    if(opts$mPhen.numGenoPCs>0) attr(genoConnection,"K") = array(0,dim = c(length(samps), length(samps)))
    genoConnection$lastPos = opts$mPhen.starting-1
    attr(genoConnection,"closeConnection") = FALSE
  }
  if(attr(genoConnection,"closeConnection")) return(NULL) ## this is the case the genoConnection was already closed

  opts$mPhen.starting  = genoConnection$lastPos+1
  opts$mPhen.rsidToDo = attr(genoConnection,"rsidToDo")
  data <- .readGenoConnection(genoConnection$conn,opts = opts)
  if(!is.null(data)){
    if(opts$mPhen.numGenoPCs>0){
      attr(genoConnection,"K") = .updateRelatedness(data$genoData, attr(genoConnection,"K"))
    }
    genoConnection = .mergeLists(genoConnection,data)
    lastPos = data$lastPos

    if( data$lastPos >= opts$mPhen.ending ) closeConnection = TRUE
    else if(length(data$rsids)<opts$mPhen.batch) closeConnection = TRUE
    rsidToDo = opts$mPhen.rsidToDo
    if(!is.null(rsidToDo)){
      rsidToDo =  rsidToDo[!(rsidToDo %in% data$rsids)]
      if(length(rsidToDo)==0) closeConnection = TRUE
    }
    attr(genoConnection,"rsidToDo") = rsidToDo
  }else{
    closeConnection = TRUE
    genoConnection$genoData = NULL 
  }
  attr(genoConnection,"closeConnection") = closeConnection
  if(!is.null(genoConnection$genoData)) attr(genoConnection$genoData,"closeConnection") = closeConnection
  if(closeConnection){
    .closeGenoConnection(genoConnection$conn)
    if(opts$mPhen.numGenoPCs>0){
      ev = eigen(attr(genoConnection,"K"))
      inds = which(ev$values>1e-5)
      max = min(length(inds), opts$mPhen.numGenoPCs)
      pcs = ev$vectors[,inds[1:max]]
      dimnames(pcs) = list(dimnames(genoConnection$genoData)[[1]],paste("PC",1:max,sep=""))
      genoConnection$pcs = pcs
    }
    ##invisible(genoConnection$genoFiles)
  }
  invisible(genoConnection)
}


##Writes output to output connection
#outC is obtained from mPhen.openOutputConnection
#results is obtained from mPhen.assoc()
#pos_info is vector of marker information obtained from mPhen.readGenoConnection()$pos_info
# If .writeQC is true then a file of quality control information is written for each sample
# this essentially represents the cumulative amount of signal contributed to association across
# all markers.
## This function modifies outC, which is returned (and should be used to update the variable in the executing script)
mPhen.writeOutput<-function(results,
                            output = getOption("mPhen.resultsName","resultsDir/"),
                            geno = NULL,  #for qc and fingerprint
                            ###	     phenoObject = NULL,  #for recording information on phenotype transformations
                            towrite = list(long.txt = getOption("mPhen.writeLong",TRUE),
                                           qc.txt =  getOption("mPhen.writeQC",FALSE), 
                                           wide.txt = getOption("mPhen.writeWide",TRUE)),
                            toplot = list(.manh = TRUE,
                                          .qq = TRUE,
                                          .heatm = TRUE,
                                          .fprint = !is.null(geno)),
                            opts = mPhen.options("plot")

                            ){
  .writeQC =grep("qc",names(towrite))
  .writeLong = grep("long",names(towrite))
  .writeWide = grep("wide",names(towrite))
  towriteQC = if(length(.writeQC)>0) towrite[[.writeQC]] else FALSE
  towriteLong = if(length(.writeLong)>0)  towrite[[.writeLong]] else FALSE
  towriteWide = if(length(.writeWide)>0)  towrite[[.writeWide]] else FALSE
  closeConnection = TRUE
  if(!is.null(geno)){
    if(!is.null(attr(geno,"closeConnection"))) closeConnection = attr(geno,"closeConnection")
  }
  if(is.null(results)) closeConnection = TRUE
  else if(is.null(results$Results)) closeConnection = TRUE
  first = FALSE
  if(is.null(attr(output,"outDirectory"))){
    first=TRUE
    outputL = list();      
    outDirectory = output
    attr(outputL,"outDirectory") = outDirectory
    dir.create(outDirectory,showWarnings = FALSE)
  }else{
    outDirectory = attr(output,"outDirectory")
    outputL = output
  }
  nme = deparse(substitute(results))
  if(nme=="NULL") stop("cannot have name of object as NULL")
  index = which(names(outputL)==nme)
  if(length(index)==0 ){
    index = length(outputL)+1
    limits = NULL
    #if(!is.null(phenoObject) && first) limits = phenoObject$limit
    outputL[[index]] =  .openOutputConnection(outDirectory, nme, towrite,toplot, limits = limits, regopts = attr(results,"opts"))
    names(outputL)[index] = nme
  }
  outC = outputL[[index]]
  if(!outC$open){
    if(towriteQC & !is.null(geno)) warning("sample qc metrics invalid if re-opening old connection")
    outC = .reopenConnection(outC)
    warning("re-opening same connection")
  }
  if(is.null(results)){
    closeConnection = TRUE
    effPhe = 1.0
  }else{
    effPhe = 1.0; 
    pos_info = dimnames(results$Res)[[2]]
    res1 = results$Results
    indices = 1:dim(res1)[3]
    if(!is.null(opts$mPhen.onlyShowJoint)){
      jind = grep("JointModel",dimnames(res1)[[3]])
      if(length(jind)>0){
        if(opts$mPhen.onlyShowJoint) indices = jind
        else indices = indices[-jind] 
      }
    }
    first = (length(outC$pvs)==0)
    if(towriteLong){
      file=outC$files[[.writeLong]]
      toprint =  .flatten(results, pos_info,"chr_start")
      write.table(toprint,file = file,sep="\t",quote=F,col.names=first,row.names=F,append=!first)
      flush(file)
    }
    if(towriteWide){
      filew = outC$files[[.writeWide]]
      .writeTable(results, filew,pos_info, "%3.2e")
      attr(filew,"first")<-FALSE
      flush(filew)
    }
    for(kk in 1:(dim(res1)[[1]])){
      toappendBeta = as.matrix(res1[kk,,indices,1])
      toappend = as.matrix(res1[kk,,indices,2])
      if(dim(toappend)[[2]]!=length(indices)){
        toappend=t(toappend)
        toappendBeta=t(toappendBeta)
      }
      for(k in 1:(dim(toappendBeta)[1])){
          log10p = getOption("mPhen.log10p",FALSE)
          betapvthresh = opts$mPhen.betaPvThresh
	  if(log10p) betapvthresh = log10(betapvthresh)
          naind = is.na(toappendBeta[k,])
          indi = (toappend[k,] > betapvthresh) 
          toappendBeta[k,indi] = 0
	  toappendBeta[k,naind] = NA
      }	  
	  #print(paste(dimnames(res1)[[1]][kk], dimnames(res1)[[3]][k]))
          #print(paste(kk,k,max(toappendBeta,na.rm=T)))
          #print(min(toappendBeta,na.rm=T))
      dimnames(toappend) = list(dimnames(res1)[[2]],dimnames(res1)[[3]][indices])
      dimnames(toappendBeta) = dimnames(toappend)
      dimnames(toappend)[[1]] = dimnames(res1)[[2]]
      dimnames(toappendBeta)[[1]] = 	 dimnames(toappend)[[1]]
      if(first){
        outC$pvs[[kk]] = toappend
        outC$beta[[kk]] = toappendBeta
      } else{
        outC$pvs[[kk]] = rbind(outC$pvs[[kk]],toappend)
        outC$beta[[kk]] = rbind(outC$pvs[[kk]],toappendBeta)
      }
    }
    names(outC$pvs) = dimnames(res1)[[1]]
    names(outC$beta) = dimnames(res1)[[1]]
    if(towriteQC & !is.null(geno)) outC$sampleQC = .updateSampleQC(geno, results, outC$sampleQC, pvthresh = 0.05,maf_thresh = 0.05)
  }
  if(closeConnection){ 
    sampleQC = outC$sampleQC
    if(!is.null(sampleQC) && towriteQC){
      qcfile = outC$files[[.writeQC]]
      toprint = t(.flattenBetaX(sampleQC))
      write.table(toprint,file=qcfile,col.names=FALSE,row.names=TRUE,sep="\t",append=F,quote=F)
      flush(qcfile)
    }
    for(k in 1:length(outC$files)){
      if(!is.null(outC$files[[k]]))    tryCatch(    close(outC$files[[k]]), error = function(e) print(NULL))
    }
    outC$open = FALSE
    if(length(plot)>0){
       mag = length(outC$pvs)
      if(!getOption("mPhen.log10p",FALSE)) 
        logfun<-function(x) -log10(x)
      else     
        logfun<-function(x) -x
      for(kk in 1:length(toplot)){
        if(toplot[[kk]]){
	  filek = outC$plotFiles[[kk]]
          pdf(filek)   #, width = mag *7, height = 7)
          nme =names(toplot)[[kk]]
          if(nme==".fprint"){
            if(!is.null(geno) & !is.null(results)){
              .fprint(geno,results,opts)
            }
          }else if(nme==".heatm"){
            logfun1<-function(x) x
	    if(!is.null(attr(outputL,"Colv"))) Colv = attr(outputL,"Colv")
	    else Colv = opts$mPhen.Colv
	    if(!is.null(Colv)){
		  if(!is.logical(Colv)) {
		#	print("Colv")
                #        print(Colv)
			len = tryCatch(length(as.hclust(Colv)$labels)	,error = function(e) 0)
			if(len!=dim(outC$pvs[[1]])[2]){
				Colv = getOption("mPhen.Colv",TRUE)
			}
                }
             }
#             print(paste("Colv",Colv))
	     opts$mPhen.Colv = Colv
	     opts$mPhen.Colv = .heatm(outC$pvs, effPhe, logfun, opts,FALSE)

            .heatm(outC$beta, effPhe, logfun1, opts,TRUE,filek)
          }else{
            lapply(list(outC$pvs), nme,effPhe, logfun, opts)
          }
          dev.off()
        }
      }
    }
  }
  outputL[[index]] = outC
  attr(outputL,"Colv") = opts$mPhen.Colv
  invisible( outputL)
}


###.phenoFiles is a list across cohorts, within that is a vector across different files per cohort
mPhen.readPhenoFiles<-function(phenoFiles,
                               limitFile = getOption("mPhen.limitFile","./limit.txt"),
                               excludeFile =getOption("mPhen.excludeFile","./exclude.txt"),
                               opts = mPhen.options("pheno.input")){
   todo=NULL

  .naRowThresh = opts$mPhen.naRowThresh
  .naColThresh = opts$mPhen.naColThresh
  .sdThresh = opts$mPhen.sdThresh
  fillMissingPhens = opts$mPhen.fillMissingPhens
  .quantileThresh=opts$mPhen.quantileThresh
  if(!is.list(phenoFiles)) phenoFiles = list(phenoFiles)
  if(!is.list(excludeFile)) excludeFile = list(excludeFile)
  phenos_all = list()
  betasall = list()
  rmsall = list()
  offsets = list()
  for(j in 1:length(phenoFiles)){
    phenoFiles1 = phenoFiles[[j]] 
    phenos = list()
    for(k in 1:length(phenoFiles1)){
      phen1 = .readPhen(phenoFiles1[k],sep=opts$mPhen.sep.pheno,numHeaderRows = opts$mPhen.numHeaderRows.pheno)
      inclcol  = which((apply(phen1,2,.cntNa)<=.naColThresh * (dim(phen1)[1])) & apply(phen1,2,var,na.rm=TRUE)>.sdThresh)
      phen1 = phen1[,inclcol,drop=F]
      inclind  = which(apply(phen1,1,.cntNa)<=.naRowThresh * (dim(phen1)[2]))
      phen2 = phen1[inclind,,drop=F]
      if(fillMissingPhens  &   length(which( apply(is.na(phen2),1,sum)>0))>0){
        #     phen2 = kNNImpute(phen2,3)$x
        phen2 = .fillMissing(phen2)
      }
       inclcol  = which((apply(phen1,2,.cntNa)<=.naColThresh * (dim(phen1)[1])) & apply(phen1,2,var,na.rm=TRUE)>.sdThresh)
      phen1 = phen1[,inclcol,drop=F]
      phenos[[k]] = phen2
    }


    pheno1 = t(.mergeFiles(.trL(phenos)))

    inclind  = which(apply(pheno1,1,.cntNa)<=.naRowThresh * (dim(pheno1)[2]))
    pheno2 = pheno1[inclind,,drop=F]
    if(!is.null(excludeFile[[j]]) && file.exists(excludeFile[[j]])){
      excl = read.table(excludeFile[[j]],as.is=T,header=F,fill=T,sep="\t")
      if(dim(excl)[[2]]>1){
        thresh = quantile(excl[,2],.quantileThresh)
        top10inds = which(excl[,2]>thresh)
        exclnames = unique(excl[top10inds,1])
      }else exclnames = excl[,1]
      exclinds = (match(exclnames,dimnames(pheno2)[[1]]))
      if(length(exclinds)>0) pheno2 = pheno2[-exclinds,,drop=F]
    }
   if(opts$mPhen.imputeMissingPhens & j==1 & !opts$mPhen.onlyCommonPheno){
    if(!is.null(limitFile)){
        if(file.exists(limitFile)){
          todo = .readLimitFile(limitFile,phenNames = dimnames(pheno2)[[2]])
        }
      }
    }

   if(opts$mPhen.imputeMissingPhens & j>1 & !opts$mPhen.onlyCommonPheno){
     
     pheno_1 = phenos_all[[1]]
     print('imputing missing columns in second dataset from first data set')
     dn2 = dimnames(pheno2)[[2]]
     dn1 = dimnames(pheno_1)[[2]]
     pc1 = c()
     pc2 = c()
     todo1 = c()
     todo2 = c()
     if(!is.null(todo)){
      excls = todo[c(grep("strat",todo[,1]),grep("resid",todo[,1]),grep("proj",todo[,1]),grep("covar",todo[,1]),grep("nopheno",todo[,1])),2] 
      todos = todo[c(grep("^pheno",todo[,1])),2]
      if(todos[1]=="all") todos = dn1
      if(length(todos)==1) todos = strsplit(todos[1],";")[[1]] 
      for(excl in excls){
	  pc1 = c(pc1, which(dn1==excl))
     	  pc2 = c(pc2, which(dn2==excl))
       }
       for(t1 in todos){
	  todo1 = c(todo1,which(dn1==t1))
	  todo2 = c(todo2,which(dn2==t1))
       }
      }
      
     mtch = match(dn1,dn2)
     inds = which(is.na(mtch))
     inds1 = which(!is.na(mtch))
     tokeep = mtch[inds1]

     inds =  inds[!inds %in% pc1 & inds %in% todo1]
     inds1 = inds1[!inds1 %in% pc1]
     inds1.2 =  mtch[inds1]
     rms = rep(NA, length(inds))
     offset = rep(NA, length(inds))
     matn = array(dim = c(dim(pheno2)[1],length(inds)), dimnames = list(dimnames(pheno2)[[1]], dn1[inds]))
     betas = array(dim = c(length(inds1),length(inds)), dimnames = list(dn1[inds1], dn1[inds]))
     print(paste(dimnames(betas)[[1]],collapse=";"))
     print(paste(dimnames(betas)[[2]],collapse=";"))
     for(k in 1:length(inds)){
       k1 = inds[k]
       y = pheno_1[,k1]
       x = pheno_1[,inds1,drop=F]
       naind = is.na(y) | apply(x,1,.cntNa)>0
       y = y[!naind]
       x= x[!naind,,drop=F]
       res = .backwardSelection(y,x,NULL,rep(1,length(y)),family =  "gaussian", totweight = 1)[-(dim(x)[[2]]+1),]
       #res = summary(.getGlm(y,x,"gaussian",rep(1,length(y)), NULL))$coefficients[-1,]
       betas[,k] = res[,1]
       y1 =  x %*% betas[,k]
       offset[k] = mean(y) - mean(y1)
       rms[k] = sum( (y-(y1+offset[k]))^2)/length(y)
  print(paste("rms",dn1[k1],rms))
     }
     betasall[[j]] = betas
     rmsall[[j]] = rms
   
     offsets[[j]] = offset
     phenoJ = pheno2[,inds1.2] %*% betas
     pheno2 = cbind(pheno2[,tokeep], phenoJ)
   }else if(opts$mPhen.onlyCommonPheno){
      if(j==1){
	 common_phens = colnames(pheno2)
      }
      else{
	common_phens = common_phens[common_phens %in% colnames(pheno2)]
      }
   }
    phenos_all[[j]] = pheno2

  }
  

  names(phenos_all) = phenoFiles
  pheno3 = .mergeFiles(phenos_all,markerCol=TRUE)
  if(opts$mPhen.onlyCommonPheno){
     pheno3 = pheno3[,colnames(pheno3) %in% c(common_phens,"Cohort_")]
  }
  todo=NULL
  if(!is.null(limitFile)){
    if(file.exists(limitFile)){
      todo = .readLimitFile(limitFile,phenNames = dimnames(pheno3)[[2]])
    }
  }
  if(is.null(todo)){
    todo = .parseLimit(opts$mPhen.limit ,phenNames = dimnames(pheno3)[[2]])
  }
  result = list(pheno=pheno3, limit=todo, betas = betasall, rms = rmsall, offsets = offsets)
  result
}
