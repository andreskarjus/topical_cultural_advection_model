#### Functions to calculate topical-cultural advection ####
# Includes some COHA-specific functions for basic corpus cleaning
# And also some extras like evaluating the role of advection in lexical innovation, 
# and using an alternative LDA based approach.
#
# Andres Karjus, University of Edinburgh
#


#### functions ####

# coha_parallel >calls> fixcohaperiod >calls> coca_csv
coca_csv = function(filepath, minlen=3, 
                    markS = list(S=T,V=F,A=F,D=F,M=F,P=F),
                    rmstopwords=T, removenp=T,
                    stopgrep = "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to$|^xx$|^z|^y$|null|\"|^nnu|^np[12]*|\\.\\.\\.|^p|^v[^v]"   # stopwords by POS
){
  # for compatability with earlier run-scripts, add missing prefix parameters
  fixmarks = setdiff( c("S","V","A","D","M","P", "X", "G"), names(markS))
  if(length(fixmarks)>0){
    tmp = as.list(rep(F, length(fixmarks))); names(tmp) = fixmarks
    markS = c(markS, tmp)
  }
  
  words = read.csv(filepath, sep="\t", quote="", header=F, skipNul=T, 
                   stringsAsFactors=F, fileEncoding="latin1")
  #print(head(words))
  words[ grep("^(##[0-9]|@@[0-9])", words[,1]), 2] = "<doc>" # #=coca, @=coha  

  stopwords=NULL;nps=NULL;badocr=NULL
  if(rmstopwords){
    # "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to$|^xx$|^z|^y$|null|\"|^nnu|\\.\\.\\.|^p|^v[^v]"
    # "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to$|^xx$|^z|^y$|null|\"|^nnu|^np[12]*|\\.\\.\\.|^p|^v[^v]" # current, np as stopwords
    stopwords = grep(stopgrep, words[,3]) 
  }
  
  if(removenp){
    nps = grep("^[A-Z]", words[,1])   # all capitalized is np, hard filter.--------------
  }
  
  shorts = which(nchar(words[,2]) < minlen) #
  
  # remove some known mistagged ocr errors and other unwanted bits and pieces:
  badocr = grep("doesn$|weren$|couldn$|shouldn$|isn$|wouldn$|aren$|hadn$|wasn$|hasn$|didn$|befor$|^dieir|^diat|^diose|^diis|^diere|^diese|^dius|--|\\+|[%;©®™@.*/#=$]|-$|^-|^p[0-9]+", words[,2])
  
  
  toremove = union( c(stopwords,nps, shorts, badocr), NULL) # remove all these
  if(length(toremove) > 0){
    words = words[-toremove, ]
  }
  
  alldocs = words[,2]
  
  if(any(unlist(markS))){
    
    if(markS$S){
      ispos = grep("^nn[^u]*", words[,3])
      alldocs[ispos] = paste("S:", alldocs[ispos], sep="")
    }
    if(markS$V){
      ispos = grep("^vv", words[,3])
      alldocs[ispos] = paste("V:", alldocs[ispos], sep="")
    }
    if(markS$A){
      ispos = grep("^jj[rt]*$", words[,3])
      alldocs[ispos] = paste("A:", alldocs[ispos], sep="")
    }
    if(markS$D){
      ispos = grep("^d", words[,3])
      alldocs[ispos] = paste("D:", alldocs[ispos], sep="")
    }
    if(markS$M){
      ispos = grep("^mc([^gm]|$)", words[,3])
      alldocs[ispos] = paste("M:", alldocs[ispos], sep="")
    }
    if(markS$P){  # proper nouns
      ispos = grep("^np", words[,3])
      alldocs[ispos] = paste("P:", alldocs[ispos], sep="")
    }
    if(markS$X){  # shall, got etc
      ispos = grep("^vvx$", words[,3])
      alldocs[ispos] = paste("X:", alldocs[ispos], sep="")
    }
    
  }
  return(alldocs)
}
library(compiler)
coca_csv = cmpfun(coca_csv, options=list(optimize=3))



fixcohaperiod = function(period1, 
                         path="",
                         minlen=3, 
                         markS = list(S=F,V=F,A=F,D=F,M=F,P=F), 
                         rmstopwords, removenp, saveyear=T, byyear, FOLDER
){
  gc()
  print(paste(Sys.time(), period1))
  if(byyear){
    pfiles = file.path(path,period1) # vector of filenames with path from main folder
  }
  if(!byyear){
    pfiles = list.files(file.path(path,period1), full.names = T)
  }
  
  period = list();
  for(y in pfiles){
    #print(y)
    tryCatch({
      filepath = y 
      alldocs = coca_csv(filepath, minlen=minlen, markS = markS, rmstopwords=rmstopwords, removenp=removenp) 
      if(length(alldocs)>2){
        s = which(alldocs == "<doc>")
        fi = findInterval(seq_along(alldocs), s)
        fi[s] = max(fi)+1
        periodlist2 = split(alldocs, fi)
        if(length(s) > 0) {periodlist2[length(periodlist2)] = NULL} 
        periodlist2[which(sapply(periodlist2, length) < 2)] = NULL 
        period = c(period, periodlist2)
      }
    }, error = function(e) {print(e)})
  }
  
  # rdata filename
  if(byyear){
    subperiod=paste0(gsub( "[^_]*_([0-9]{4}).*","\\1",gsub("[^/]+/(.*)", "\\1", period1[1])),".RData") # gets year
  }
  if(!byyear){
    subperiod = paste0(ifelse(saveyear, gsub("[^0-9]","", period1), period1),".RData") # year from foldername
  }
  tryCatch({save(period, file=file.path(FOLDER, "periods", subperiod ))}, error=function(e) e )
  print(Sys.time())
  #return(period)   #  parallel, just saves
}

coha_parallel = function(theperiods, nfree=0, markS = list(S=F,V=F,A=F,D=F,M=F,P=F), minlen=3, path="", saveyear=T, rmstopwords=T, removenp=T, byyear=F, FOLDER){
  # byyear: do by year instead of decade; theperiods needs to be list where each
  # element is a vector of year files, with full paths
  library(compiler)          
  setCompilerOptions("suppressAll" = T)
  fixcohaperiod = cmpfun(fixcohaperiod, options=list(optimize=3))
  coca_csv = cmpfun(coca_csv, options=list(optimize=3))
  library(parallel)
  nc = detectCores() - nfree
  cl = makeCluster(nc)
  print(paste(Sys.time(), "start"))
  tryCatch({
    clusterExport(cl, c("fixcohaperiod", "coca_csv", "markS", "minlen", "path", "saveyear","rmstopwords", "removenp", "byyear", "FOLDER"),envir = environment()) 
    #clusterExport(cl, "FOLDER")
    tmp = parLapply(cl, theperiods, fixcohaperiod, 
                    markS = markS, minlen=minlen, path=path, saveyear=saveyear, rmstopwords=rmstopwords, removenp=removenp, byyear=byyear, FOLDER=FOLDER
    )
  }, error=function(e){print(e)},  finally = stopCluster(cl) )
  print(paste(Sys.time(), "done"))
}


dotrends = function(foldr, mint=2, nfree=0){
  files = list.files(foldr, full.names = T)  # [1:2]# DEBUGGING
  print(paste(Sys.time(), "start countmat"))
  
  library(parallel)
  nc = detectCores() - nfree
  cl = makeCluster(nc)
  clusterExport(cl, c("mint", "files"),envir = environment())
  clusterEvalQ(cl, c(library(text2vec),library(fastmatch)) )
  tryCatch({
    thecounts = parLapply(cl, files, function(f){
      load(f)
      #period=period[1:2]# DEBUGGING
      it = itoken(period, progressbar = FALSE) # period list from load()
      rm(periods); gc() # clean periods - iterator envir now contains all data
      vocab2 <- create_vocabulary(it, ngram = c(ngram_min = 1L, ngram_max = 1L))
      vocab2 <- prune_vocabulary(vocab2, term_count_min = mint)
      docdisp = vocab2$doc_count/attr(vocab2,"document_count"); names(docdisp)=vocab2$term
      sm=sum(vocab2$term_count)
      freqs=(vocab2$term_count/sm)*1000000; names(freqs)=vocab2$term
      return(list(freqs, sm, docdisp))
    })
    freqs = lapply(thecounts, function(x)   x[[1]])
    thesums = sapply(thecounts, function(x) x[[2]])
    docdisp = lapply(thecounts, function(x) x[[3]])
    rm(thecounts);gc()
    allwords = unique(unlist(lapply(freqs, function(x) names(x) ),use.names = F ) )
    
    print(paste(Sys.time(), "freqs done"))
    
    clusterExport(cl, c("allwords", "freqs", "docdisp"),envir = environment())
    lexorder = parLapply(cl, 1:length(files), function(x){ 
      return(fmatch(allwords, names(freqs[[x]]) ) )
    })
    clusterExport(cl, c("lexorder"),envir = environment())
    
    print(paste(Sys.time(), "lexorder done, aggregating counts and disps"))
    
    countmat1 = parSapply(cl, 1:length(files), function(x){
      freqvec = freqs[[x]]
      tmp = freqvec[lexorder[[x]] ]; tmp[is.na(tmp)] = 0 # if no occurrence, then 0
      return(tmp)
    } )
    rownames(countmat1) = allwords
    
    print(paste(Sys.time(), "done with count matrix, doing dispersion matrix"))
    
    dispmat = parSapply(cl, 1:length(files), function(x){
      freqvec = docdisp[[x]]
      tmp = freqvec[lexorder[[x]] ] # if no occurrence, then !NA! dispersion, not 0
      return(tmp)
    } )
    rownames(dispmat) = allwords
    
  }, error=function(e){print(e)},  finally = stopCluster(cl) )
  
  countmat=list(countmat=countmat1, thesums=thesums, dispmat=dispmat)
  try({ save("countmat", file=file.path(FOLDER, "countmat.RData")) })
  print(paste(Sys.time(), "done with dotrends"))
  return(countmat)
}

dofreqdif = function(countmat, smooth, lexicons){
  #lexicons = getlexicons(countmat[[1]],countmat[[2]], minraw) # 1 comlex for columns, 2 lexicon for dataframes
  Sys.time()
  freqdifmat = matrix(nrow=length(lexicons), ncol=ncol(countmat)); rownames(freqdifmat)=lexicons
  for(i in 2:ncol(countmat)){
    fr1 = countmat[lexicons,i-1] +smooth
    fr2 = countmat[lexicons,i]   +smooth
    #freqdif1 = (abs(fr1 - fr2))/fr1*100; freqdif1 = ifelse(fr1 > fr2, freqdif1*-1, freqdif1)
    freqdifmat[,i] = log(fr2/fr1)    #freqdif1
  }
  gc()
  return(freqdifmat)
}
library(compiler)
dofreqdif = cmpfun(dofreqdif, options=list(optimize=3))

makeppmi = function(pmat, positive=T){
  library(Matrix, quietly = T)
  pmat = Matrix::t(pmat)
  #pmat = (tcmfornews)
  #set.seed(1)
  #pmat = matrix(sample(c(0,0,0,0,0,0,1,10),5*10,T), 5,10, byrow=T)
  #pmat = Matrix::t(Matrix(pmat, sparse=T))
  
  tcmrs = Matrix::rowSums(pmat)
  tcmcs = Matrix::colSums(pmat)
  N = sum(tcmrs)
  colp = tcmcs/N
  rowp = tcmrs/N
  
  pp = pmat@p+1
  ip = pmat@i+1
  tmpx = rep(0,length(pmat@x))
  for(i in 1:(length(pmat@p)-1) ){ 
    #for(i in 1:100 ){ 
    #not0 = which(pmat[, i] > 0)
    ind = pp[i]:(pp[i+1]-1)
    not0 = ip[ind]
    icol = pmat@x[ind]
    #print(icol)
    #tmp = log( (pmat[not0,i]/N) / (rowp[not0] * colp[i] ))
    tmp = log2( (icol/N) / (rowp[not0] * colp[i] ))
    tmpx[ind] = tmp
    #print(tmp)
    # tmp = ifelse(tmp < 0, 0, tmp)
    #pmat2[not0, i] = tmp
  }
  if(positive){
    tmpx[tmpx<0] = 0
  }
  pmat2 = pmat
  pmat2@x = tmpx
  pmat2 = Matrix::t(pmat2) 
  pmat2 = drop0(pmat2,is.Csparse=T) 
  return(pmat2)
}

make_ppmi_tcm = function(nperiod, rdats, minc, gramwin, comlex, wordclassesr, wordclassesc= ".", interpol=c(0), 
                         simuversion=F, maxperiod=NA,segment_indices=NA, 
                         ppmi=T, windowweights ){   
  #  3 for simu only; simuversion is only for compatibility with ac->spok change simulation
  # NB, period files are up to ~400mb, so concatenating will eat RAM fast...
  library(text2vec)
  library(Matrix)
  if(!simuversion){
    periods=list()
    for(ni in (nperiod+interpol) ){ # 0 is current, -1 previous, 1 future, etc
      if(ni > 0 & !(ni > length(rdats)) ){ # avoid indices beyond list
        load(rdats[ni]) # load period
        periods = append(periods, period)
        #print(c(length(period), length(periods)) ) # debug
      }
    }
    rm(period) # clean lastly loaded text data object
  }
  
  if(simuversion){
    periods=list()
    for(ni in (nperiod+interpol) ){ # 0 is current, -1 previous, 1 future, etc
      if(ni > 0 & !(ni > maxperiod) ){ # avoid indices beyond list
        period = c(rdats$academic[segment_indices[[1]][[ni]] ], rdats$spoken[segment_indices[[2]][[ni]] ])
        periods = append(periods, period)
        #print(c(length(period), length(periods)) ) # debug
      }
    }
    rm(period,rdats);gc()
  }
  
  # do iter,vocab, tcm
  it = itoken(periods, progressbar = FALSE) # period list from load()
  rm(periods); gc() # clean periods - iterator envir now contains all data
  vocab2 <- create_vocabulary(it, ngram = c(ngram_min = 1L, ngram_max = 1L))
  #docmin = vocab2$vocab$doc_counts/vocab2[["document_count"]]
  vocab2 <- prune_vocabulary(vocab2, term_count_min = minc, doc_proportion_max=1, doc_proportion_min=0)# 2/vocab2$document_count) # avoid single-doc words
  # should still exclude single-doc-words? calculate dispersion/burstiness instead
  
  vectorizer <- vocab_vectorizer(vocab2)
  if(windowweights == "uniform"){
    tcm = create_tcm(it, vectorizer, skip_grams_window = gramwin, 
                     weights = 1 / rep(1,gramwin) )
  } else {
    tcm = create_tcm(it, vectorizer, skip_grams_window = gramwin,
                     weights =  1 / seq_len(gramwin))
    
  }
  rm(it, vocab2, vectorizer);gc() # free some memory
  tcm = tcm + Matrix::t(triu(tcm))  # originally triangular, copy to make it work
  # NB: vocab_vectorizer and create_tcm are changed in text2vec 0.5.0, 
  # now params in latter; also has symmetric window and weights options!
  
  Matrix::diag(tcm) = 0
  
  # keep only persistent lexicon as context (and pos if set)
  if(!is.na(comlex[1])){
    if(nchar(wordclassesc) > 1){ comlex = grep(wordclassesc, comlex, value=T)} # if filtered comlex
    # NB: currently filtering wordclass in columns works only with comlex, tweak ifelse if needed
    if(nchar(wordclassesr)>1) {
      tcm = tcm[grep(wordclassesr,rownames(tcm)), comlex ]  }
    else { tcm = tcm[, comlex ] }
  } else { 
    # if NO COMLEX
    if(nchar(wordclassesr)>1) {
      tcm = tcm[grep(wordclassesr,rownames(tcm)),  ]  
    }   # NO else, use full tcm, no filtering for wordclass or comlex
    
  }  
  gc()
  if(ppmi){
    tcm = makeppmi(tcm)               # run ppmi weigthing
  }
  #print(paste(Sys.time(), "tcm done")) # in cluster, nowhere to print
  return(tcm)
}



dorelevants = function(ppmimat,relevancy_threshold){
  library(grr)
  relevants = list()
  print(paste(Sys.time(), "do relevance vectors"))
  relevants = list(); length(relevants) = nrow(ppmimat)
  ppmimat = Matrix::as.matrix(ppmimat)
  for(x in 1:nrow(ppmimat)){
    tmp = ppmimat[x,-which(ppmimat[x,]==0)]
    y=tmp[rev(order2(tmp))] # sort2 is decreasing=F BUT MESSES UP NAMES, USE order2
    y=y[1:min(length(y),relevancy_threshold)] # but this need top ones, so rev
    relevants[[x]] = y; names(relevants[[x]]) = names(y)
  }
  print(paste(Sys.time(),"relevance vectors done"))
  return(relevants)
}
dorelevants = cmpfun(dorelevants, options=list(optimize=3))

aggregaterelevants = function(foldr, pattern = NULL){
  # aggregate relevants, 20-30mb per period
  allrelevants = list()
  for(x in list.files(foldr, full.names = T, pattern=pattern)){
    load(x)
    allrelevants[[gsub("[^0-9]","",x) ]] = relevants 
  }
  gc()
  return(allrelevants)
}

# tcm+ppmi calculations
perioddatas = function(rdats = list.files(file.path(FOLDER,"periods"), full.names = T), comlex, minc, tcmwindow, windowweights="weighted", relevancy_threshold, nfree=0,wordclassesr="", tcmfolder=file.path(FOLDER,"periodtcm/"),relefolder=file.path(FOLDER,"periodrelevants/"), interpol=0, yearnames=seq(1810,2000,10)){
  library(parallel)
  nc = detectCores() - nfree
  print(c(paste(Sys.time(), "start tcm phase")))
  
  # TCMs
  if(length(rdats)>0){
    cl = makeCluster(nc)
    clusterExport(cl, c("rdats","comlex", "minc", "tcmwindow","wordclassesr", "interpol", "yearnames", "tcmfolder", "windowweights"),envir = environment())  # from local
    clusterExport(cl, c("make_ppmi_tcm", "makeppmi")) # from global
    clusterEvalQ(cl, c(library(fastmatch), library(text2vec), library(grr), library(Matrix)))
    
    tryCatch({
      parLapply(cl, 1:length(rdats), function(x){
        # load(x) # load period data -> in function instead; send just period index
        # do and save tcm
        # NB: assumes ordered period files!
        tcm = make_ppmi_tcm(nperiod=x,rdats=rdats, comlex=comlex, minc=minc, gramwin = tcmwindow, wordclassesr=wordclassesr, wordclassesc="", interpol=interpol, windowweights=windowweights)
        save(tcm, file=file.path(tcmfolder, paste0(yearnames[x],".RData" )))
        rm(period);gc() # remove period data
      }
      )
    }, error=function(e){print(e)}, finally = stopCluster(cl))  
    gc()
  }
  print(c(paste(Sys.time(), "TCMs done")))
  
  # get relevance sets
  if(length(tcmfolder)>0){
    tcmdats = list.files(tcmfolder, full.names = T)
    cl = makeCluster(nc)
    clusterEvalQ(cl,  c(library(grr), library(Matrix)))
    clusterExport(cl, c("dorelevants"))
    clusterExport(cl, c("relefolder","relevancy_threshold"),envir = environment())  # from local
    tryCatch({
      parLapply(cl, tcmdats, function(x){
        load(x)  # loads tcm
        relevants = dorelevants(tcm, relevancy_threshold=relevancy_threshold)
        names(relevants) = rownames(tcm)
        save(relevants, file=file.path(relefolder, paste0(gsub("[^0-9]","",x),".RData" )))
        rm(tcm);gc()
      })
    }, error=function(e){print(e)}, finally = stopCluster(cl)) 
    gc()
  print(c(paste(Sys.time(), "relevance vectors done")))
  }
}





# generates matrix for advection
dotopicsloop = function(periods, difmat, countmat=NULL, wordclassesr, threshold,
                        estimate0=F, estimate2=F){ 
  if(length(wordclassesr) == 1){
    words = grep(wordclassesr, rownames(difmat), value = T)  # by regex; eg do only noun words
  } else {
    words = wordclassesr
  } # if len>1 then assume wordlist, use that instead
  
  inertiamat=matrix(NA, nrow=length(words),ncol=ncol(difmat)); rownames(inertiamat) = words
  #inertiamat2 = inertiamat
  estimat = inertiamat
  for(x in 2:length(periods)){
    #print(x)
    relevants=periods[[x]]
    inperiod = which(words %in% names(relevants))
    
    # determine if estimating based on previous or future relevant vectors
    # disable if using interpolation at the TCM level
    canestimate0=NULL
    if(estimate0){
      canestimate0 = setdiff(which(words %in% names(periods[[x-1]])), inperiod)
    }
    canestimate2=NULL
    if(estimate2){
      if(x < length(periods)){  # if not last period
        canestimate2 = setdiff(which(words %in% names(periods[[x+1]])),
                               union(canestimate0, inperiod))
      }
    }
    #estimat[inperiod,x] = F; estimat[union(canestimate0, canestimate2), x] = T # others NA
    
    # DO ADVECTION
    # advection from current values of present relevants
    for(i in inperiod){ 
      nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
      inertiamat[i,x] = weighted.mean(difmat   [ names(relevants[[words[i] ]][nrelevants] ), x],
                                      w = relevants[[ words[i] ]][nrelevants] , na.rm = T) 
      #inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]] ), x],
      #                                w = relevants[[ words[i] ]], na.rm = T) 
    }
    
    # interpolate if missing:
    # advection from current values of relevants from t-1:
    if(length(canestimate0)>0){
      relevants=periods[[x-1]]
      for(i in canestimate0){ 
        nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
        inertiamat[i,x] = weighted.mean(difmat   [ names(relevants[[words[i] ]][nrelevants] ), x],
                                        w = relevants[[ words[i] ]][nrelevants] , na.rm = T)
        # inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]]), x],
        #                                 w = relevants[[ words[i] ]], na.rm = T) 
      }
    }
    # advection from current values of relevants from t+1:
    if(length(canestimate2)>0){ # does not go into if last period (length null)
      relevants=periods[[x+1]]
      for(i in canestimate2){
        nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
        inertiamat[i,x] = weighted.mean(difmat[ names(relevants[[words[i] ]][nrelevants] ), x],
                                        w = relevants[[ words[i] ]][nrelevants], na.rm = T)
        # inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]]), x],
        #                              w = relevants[[ words[i] ]], na.rm = T) 
      }
    }
    rm(relevants);gc()
  }
  #attr(inertiamat, "estimat") = estimat
  return( list( inertiamat=inertiamat, estimat=estimat))
}
library(compiler)
dotopicsloop = cmpfun(dotopicsloop, options=list(optimize=3))





#### advection via LDA 

advectionLDA = function(slex, topic_vocab_distr, word_freq_change_per_topic, words_distr_across_topics){
  library(Matrix)
  ix = match(slex, colnames(topic_vocab_distr))
  if(any(is.na(ix))){stop("LDA lexicon <-> components mismatch")}
  adv = rep(NA, length(slex));names(adv) = slex
  
  # turn into sparse to improve speed
  topic_vocab_distr = Matrix(topic_vocab_distr, sparse=T) # default dgCMatrix
  word_freq_change_per_topic = Matrix(word_freq_change_per_topic, sparse=T)
  
  for(i in 1: length(slex)){
    # need to redo topic change every time to exclude word...
    # topic_freq_change = numeric()
    # for(t in 1:nrow(topic_vocab_distr)){
    #   topic_freq_change[t] = sum(topic_vocab_distr[t, -ix] * word_freq_change_per_topic[t,-ix] )
    # }
    # faster:
    topic_freq_change = rowSums(topic_vocab_distr[, -(ix[i]) ] * word_freq_change_per_topic[,-(ix[i]) ])
    adv[i] = weighted.mean(topic_freq_change, words_distr_across_topics[,ix[i] ])
  }
  gc()
  return(adv)
}
library(compiler)
advectionLDA = cmpfun(advectionLDA, options=list(optimize=3))



doLDAperiod = function(nperiod, rdats, interpol=c(0), 
                       topicsperword=25,topics=500, minc,
                       a=0.1, b=0.1,n_iter=5000, word_freq_change){
  if(length(interpol)==1 & interpol[1] != 0 ){stop("interpol value not 0 and not length>1")} # sanity check
  
  library(text2vec)
  periods=list()
  for(ni in (nperiod+interpol) ){ # 0 is current, -1 previous, 1 future, etc
    if(ni > 0 & !(ni > length(rdats)) ){ # avoid indices beyond list
      load(rdats[ni]) # load period
      periods = append(periods, period)
      #print(c(length(period), length(periods)) ) # debug
    }
  }
  rm(period) # clean last one
  
  it = itoken(periods, progressbar = FALSE)
  rm(periods);gc() # it environment has data now
  v = create_vocabulary(it) 
  
  v= prune_vocabulary(v, term_count_min = minc)
  
  # importantly: if flexible number of topics, dependent on number of unique words in the model (above 100 threshold), nwords/25 by default
  if(!topicsperword){ # if false, then use fixed topics
    n = topics
  } else { # else tie to word type count
    n = nrow(v)/topicsperword  
  }
  
  vectorizer = vocab_vectorizer(v)
  dtm = create_dtm(it, vectorizer, type = "dgTMatrix")
  
  rm(v, it, vectorizer);gc() # cleanup
  
  lda_model = LatentDirichletAllocation$new(n_topics = n,
                                            topic_word_prior=a, doc_topic_prior=b)
  doc_topic_distr = lda_model$fit_transform(dtm, n_iter=n_iter, check_convergence_every_n = 50)
  
  # words_dist_across_topics  # for each word, how often it appears in topic 1, topic 2, topic 3 (normalised counts)
  words_distr_across_topics = t(normalize(t(lda_model$components), 'l1')) # [topics, words]
  lex = colnames(words_distr_across_topics)
  
  topic_vocab_distr  = normalize(lda_model$components,'l1')   # [topics, words]
  
  # Then, for each word: 
  word_freq_change_per_topic = words_distr_across_topics
  for(i in lex){
    word_freq_change_per_topic[,i] = words_distr_across_topics[,i] * word_freq_change[i]
  }
  
  slex = grep("^S:", lex, value=T)
  adv = advectionLDA(slex, topic_vocab_distr, word_freq_change_per_topic, words_distr_across_topics)
  
  dispersion_tops = mean(apply(words_distr_across_topics[,slex], 2, sd),na.rm=T) # do dispersion (standard dev across topics) while at it  -> remove mean if not for debug
  
  #return(list(advection=adv, dispersion_tops=dispersion_tops)) # dispersion currently for debug
  #return(list(advection=adv, dispersion_docs=dispersion_docs,dispersion_tops=dispersion_tops ))
  return(adv)
}

newadvections = function(countmat, newsrelevants, freqdifmat1, relevancy_threshold, xbefore=10){
  library(Rmisc)
  smat = countmat[[1]][unique(names(newsrelevants[[1]])),]
  newslist = list();entry=character(); nalist=list()
  for(i in 17:20){
    xnew = rownames(smat)[apply(smat,1, function(x) x[i]>0 & all(x[1:(i-1)]==0) )]
    #ai = sapply(xnew, function(x) which(!is.na(advection1$inertiamat[x,1:min(i+3, 20)]))[1])
    ai=xnew
    for(x in 1:length(ai)){
      if(!is.na(ai[x])){
        #rel=allrelevants1[[ai[x]]][[xnew[x]]][1:75]
        rel=newsrelevants[[1]][[xnew[x]]][1:relevancy_threshold] # relevants for the advection
        xa = numeric()
        notna=numeric()
        jj=1
        for(j in (i-xbefore):i){
          # set changes to NA if 0 frequency:
          freqdifvec = freqdifmat1[names(rel),j]
          freqdifvec[which(countmat[[1]][names(rel),j]==0)] = NA
          
          xa[jj] = weighted.mean(freqdifvec, rel, na.rm=T )  # does advection
          notna[[jj]] = length(which(!is.na(freqdifvec)))
          jj=jj+1
        }
        #newslist[[paste(xnew[x],i,sep="_") ]] = xa
        newslist[[xnew[x] ]] = xa
        #entry[xnew[x] ] = i
        #nalist[[ xnew[x] ]] = notna
      }
    }
  }
  return(newslist)
}

plotword = function(x, freqdif = freqdifmat, countm = countmat[[1]], adv=advection1$inertiamat, sem=F, ylims=c(1,100), xlims=c(1,20), uselims=F, decades=seq(1810,2000,10), ab=NA,yround=100, returnit=F, red=T ){
  x1 = (which(countm[x,] > 0)[1])
  #   points(x1, countm[x,x1])           # debug
  pseq = (x1+1):ncol(countm)
  nas = rep(0, x1-1)
  logfreq0 = (c(nas, log(countm[x,x1]), freqdif[x,pseq]))
  logfreq = cumsum(logfreq0)
  logfreq[0:length(nas)]=NA
  
  
  x2 = which(!is.na(adv[x,]) )[1]
  #   points(x2, countm[x,x2])           # debug
  nas = rep(0, x2-1)
  pseq = (x2+1):ncol(countm)
  #logadv0 = c(nas, logfreq[x2], adv[x,pseq])
  logadv0 = c(nas, log(mean(countm[x,], na.rm=T)), adv[x,pseq])
  logadv = logadv0
  logadv[is.na(logadv)]=0
  logadv = cumsum(logadv)               
  logadv[0:length(nas)]=NA 
  logadv[is.na(adv[x,])] = NA
  
  fa=c(nas, logfreq[x2], logfreq0[pseq] - logadv0[pseq])
  fa[is.na(fa)]=0 
  
  cfa = cumsum(fa)
  cfa[0:length(nas)]=NA
  cfa[is.na(logadv)]=NA
  
  if(!uselims){
    xlims=c(1,length(logfreq))
    mx=max(exp(c( logfreq,logadv,cfa )), na.rm=T)
    ylims=c(0, ceiling(mx/yround)*yround )
  }
  plot( 0, type="n", ylim=ylims, col="darkgray", ylab="", xlab="", 
        #ylab="permillion frequency", xlab="decades",
        main=gsub("^([A-Z]:){1,2}","",x ),xlim=xlims, xaxt="n")
  axis(1, at = 1:length(logfreq), labels = decades)
  abline(v=ab, col=rgb(0,0,0,0.2))
  
  abline(v=seq(xlims[1],xlims[2], by=5)-3, col=rgb(0,0,0,0.2),lty=3,lwd=1 )# grid
  abline(h=seq(ylims[1],ylims[2], length.out = 4)[c(2,3)], col=rgb(0,0,0,0.2),lty=3,lwd=1 )# grid
  
  logfreq[countm[x,]==0]=NA   
  points( exp(logadv ), type="o", col=rgb(0,0.6,0,0.1), cex=0.3, lty=1, lwd=6)#rgb(1,0,0,0.5)
  points( exp(logfreq ), type="o", col=rgb(0,0,0,1), cex=0.5, lty=1)
  red = ifelse(red==T, rgb(0.7,0,0, 0.7), NA)
  points( exp(  cfa ), type="o", col=red, lwd=2, lty=3, cex=0.7)

  
  if(is.na(exp(logfreq)[1]) | exp(logfreq)[1]<ylims[2]*0.6){
    y = seq(ylims[2]*0.8, ylims[2]*0.9, length.out = 3)
  }  
  else y = seq(ylims[2]/100, ylims[2]/10, length.out = 3) 
  text(1, y, rev(c(" - frequency", " - topic", "... adjusted")), col=rev(c("black", rgb(0,0.6,0,0.6), "darkred")), font=c(1,2,1), cex=0.9, pos=4)
  
  if(returnit){
    return(list(freq=exp(logfreq ), corrected = exp(  cfa )  ))
  }
  
}
