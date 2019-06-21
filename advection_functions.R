# Functions for the calculating topical advection
# Andres Karjus



print(rbind("Package"="status:",cbind(ifelse(!sapply(p<-c("magrittr", "text2vec", "compiler", "grr", "parallel", "Matrix", "fastmatch", "Rmisc", "xtable", "RColorBrewer","entropy", "beeswarm"), function(x) require(x, character.only = T)), "missing, please install!", "loaded"))), quote=F)


coha_csv = function(filepath, minlen=3, 
                    markS=F,rmstopwords=T,removenp,
                    stopliteral=T,
                    fixnums = T,
                    fixhyphens = T,
                    docborders= "^##[0-9]|@@[0-9]|^@|<p>",
                    untouchables = NULL
                    
){
  # docborders:
  # ##=coca, @@=coha  
  # |^@$ = censorship holes, leave out to not split docs. But bad if need to filter single-document words
  # <p> paragraph tag appears in more recent coha, might as well make use of it.
  # if left in, should be filtered by minlen and/or stopwords
  # coca spoken: @!xxx speaker tags; now @ counts; if not: their 2nd column (lemma) is zero-lenght, minlen would filter out
  
  
  # define grep to match stopword pos tags:
  stopgrep = "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to|^xx|^z|^y|null|\"|^nnu|np|\\.\\.\\.|!|^p|^v[^v]|^u"
  # Note: all others rely on the primary tag, except for NP proper noun which is matched is subsequent
  # tags as well (in case of multitag like nn1_np1), since it's such a common error in coha.
  # Note that the tag for some punctuation is the actual punctuation itself: ! ...
  
  # define list of additional stopwords - in case they escape the pos filter (e.g. the weird hyphenated ones that will have hyphens filtered down the line) and the max length filter (current default 3)
  stopliterals = unique( c("about","above","after","again","against","because","been","before","being","below","between","both","could","does","doing","down","each","from","further","have","having","he'd","he'll","he's","here","here's","hers","how's","i'll","i've","into","it's","let's","more","most","once","only","other","ought","ours","over","same","she'd","she'll","she's","should","some","such","than","that","that's","their","theirs","them","then","there","there's","these","they","they'd","they'll","they're","they've","this","those","through","under","until","very","we'd","we'll","we're","we've","were","what","what's","when","when's","where","where's","which","while","who's","whom","why's","with","would","you'd","you'll","you're","you've","your","yours","don't","ain't","can't","the", "but", 
  "be", "is", "was", "we", "me", "his", "her", "you", "our", "us", "my", "she", "he", "they", "or", "of", "to", "are", "has","all", "and", "did", "it", "for", "any", "who", "so", "do","by", "thi", "go", "got", "get", "also", "maybe", "no", "yes", "if", "how", "as", "at", "didn't") )
   
  
  # for compatability with earlier run-scripts, add missing prefix parameters
  fixmarks = setdiff( c("N","V","A","D","M","P"), names(markS))
  if(length(fixmarks)>0){
    tmp = as.list(rep(F, length(fixmarks))); names(tmp) = fixmarks
    markS = c(markS, tmp)
  }
  
  ## load data
  words = read.csv(filepath, sep="\t", quote="", header=F, skipNul=T, 
                   stringsAsFactors=F, fileEncoding="latin1")
  #print(head(words))
  doctags = grep(docborders, words[,1])
  words[doctags, 2] = "<doc>" 
  words[doctags, 3] = "<doc>"  # to avoid filters on borders (some tagged as mc for some reason in corpus)
  # NB: COHA lemma column is always lowercase!
  
  
  stopwords=NULL; stopwords2=NULL; nps=NULL;badocr=NULL
  if(rmstopwords){
    # "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to$|^xx$|^z|^y$|null|\"|^nnu|\\.\\.\\.|^p|^v[^v]"
    # "^ge|^a|^c|^d|^ex$|^i|^r[egpa]|^rrq|^mc|^to$|^xx$|^z|^y$|null|\"|^nnu|^np[12]*|\\.\\.\\.|^p|^v[^v]" # current, np as stopwords
    stopwords = grep(stopgrep, words[,3]) 
    # coha <p> tag is null, " tag is actual ", ! tag is !
    # ... tag is ...
    # nnu is units of measurement eg $1000 but sometimes lemma field empty
    # fo is formula, mostly bad ocr eg &+++ BUT ALSO doc borders, which are neede
    # z is -- dash, zz is letters
    # u = interjection
  }
  
  if(removenp){
    #nps = grep("^np", words[,3]) # also under stopwords
    nps = grep("^[A-Z]", words[,1])   # all capitalized considered np, hard filter.----
  }
  
  shorts = which(nchar(words[,2]) < minlen) #
  
  badocr = grep(
    "^doesn|weren[t]*$|^could[nta]*$|^should[nta]*$|^isn[t]*$|^would[nta]*$|^aren[t]*$|^hadn[t]*$|^wasn[t]*$|^hasn[t]*$|^didn[t]*$|^befor$|^dieir|^diat|^diose|^diis|^diere|^diese|^dius|\\+|[%;©®™@.*/#=$]|^p[0-9]+|^[^a-z]*$|^[01][a-z]+|^[[:punct:]]|[[:punct:]]$", 
    words[,2])
  # ^- -$ --  broken compounds, but maybe still keep? especially now if removing hyphens
  #  ---- , &+++++ ; p[0-9]+ are page numbers
  # doesn weren aren  ->  't short forms! but sometimes tagged as nouns...
  # ^[^a-z]*$ - to catch if no (English) letters present but still tagged as content word class
  # ^[0-9][a-z]+ - mostly bad ocr, single number instead of capital (O -> 0)
  # remove words starting with "-" (or any punctuation), mostly broken hypenation (eg gets "ing" as word)
  # new, for COCA: remove words ending with .!? (unsegmented punctuation)
  
  alldocs = words[,2]  # keep only lemmas after this
  ## add POS prefixes is applicable
  # (need to do this before removing stopwords and such so indices match)
  
  # match for the literal stopwords once before adding tags
  if(stopliteral){
    stopwords21 = which(alldocs %fin% stopliterals) # fastmatch
  }
  
  if(any(unlist(markS))){
    # problem in coha/coca: double tags ( _ ) and @-signs...!
    if(markS$N){
      ispos = grep("^nn([^u]|$)", words[,3])
      alldocs[ispos] = paste0("N:", alldocs[ispos])
    }
    if(markS$V){
      ispos = grep("^vv", words[,3])
      alldocs[ispos] = paste0("V:", alldocs[ispos])
    }
    if(markS$A){
      ispos = grep("^j", words[,3])
      alldocs[ispos] = paste0("A:", alldocs[ispos])
    }
    if(markS$D){
      ispos = grep("^d", words[,3])
      alldocs[ispos] = paste0("D:", alldocs[ispos])
    }
    if(markS$M){
      ispos = grep("^mc([^gm]|$)", words[,3])
      alldocs[ispos] = paste0("M:", alldocs[ispos])
    }
    if(markS$P){  # proper nouns   # see stopwords above
      ispos = grep("^np", words[,3])
      alldocs[ispos] = paste0("P:", alldocs[ispos])
    }
  }
  rm(words) # no need for the full csv anymore, just lemmas
  
  ### only these following two modify the actual words:
  # fix hyphens and nums here before removing, to be able to do newwords-filtering
  # note: this currently allows some actually <minlen words to sneak in
  if(fixhyphens){
    alldocs = gsub("-", "", alldocs)
  }
  if(fixnums){  # replace all numbers with generic 000; numbers separated by ., considered single number
    alldocs = gsub("[0-9]+[0-9.,]*", "000", alldocs)
  } # importantly, numbers like 000-000 becomes one single 000 if hyphens are removed too
  ###
  
  
  # match for the literal stopwords here once more (now that the hyphens are gone)
  if(stopliteral){
    stopwords22 = which(alldocs %fin% stopliterals) # fastmatch
  }
  
  # don't touch targets (likely mistagged if on the remove list)
  # and don't touch doc tags!
  # if targets have tags - tags were already applied to current data above
  if(is.null(untouchables)){
    dontremove=doctags
  } else {
    dontremove = unique( c(doctags,          # shouldn't need unique but anyway
                           which(alldocs %fin% untouchables) ) )  # fastmatch
  }
  
  # concatenate remove indices
  toremove = setdiff( c(stopwords, stopwords21, stopwords22, nps, shorts, badocr ), 
                      dontremove) 
  # but keep if on the list
  if(length(toremove) > 0){
    alldocs = alldocs[ -toremove ]
  }
  
  return(alldocs)
}



# docorpus_parallel > fixcohaperiod > coha_csv
fixcohaperiod = function(period1, 
                         path,
                         minlen=3, 
                         markS, rmstopwords, removenp, saveyear=T, byyear,
                         onlyvocab=F,
                         untouchables,
                         FOLDER,periodfolder
){
  gc()
  print(paste(Sys.time(), period1))
  #Magazine/wlp_mag_1996.txt
  
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
      filepath = y #paste0(path, period1,"/", y) # now uses fullnames that has path
      
      alldocs = coha_csv(filepath, minlen=minlen, markS = markS, rmstopwords=rmstopwords, removenp=removenp, untouchables=untouchables) 
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
    subperiod=gsub( "[^_]*_([0-9]{4}).*","\\1",gsub("[^/]+/(.*)", "\\1", period1[1]) ) # gets year
  }
  if(!byyear){
    subperiod = ifelse(saveyear, gsub("[^0-9]","", period1), period1) # year from foldername
  }
  
  if(!onlyvocab){
    try(save(period, file=file.path(FOLDER, periodfolder, paste0(subperiod,".RData" ) )))
    print(Sys.time())
    return(NULL) #  parallel, just saves into folder
  }
  if(onlyvocab){
    vocab = itoken(period, progressbar = FALSE) %>% 
      create_vocabulary(ngram = c(ngram_min = 1L, ngram_max = 1L))
    tmp=vocab$term
    vocab = vocab$term_count; names(vocab) = tmp 
    #try(save(vocab, file=file.path(FOLDER, periodfolder, paste0(subperiod,".RData" ) )))
    return(vocab)
  }  
}


docorpus_parallel = function(theperiods, 
                             nfree=0, # how many cpu cores to leave unused (increase if RAM issues)
                             markS = list(N=F,V=F,A=F,D=F,M=F,P=F), 
                             minlen=3, 
                             path, 
                             saveyear=T,    # use folder decade as filename
                             rmstopwords=T, # remove stopwords?
                             removenp=F,    # remove proper nouns by removing all Capitalized?
                             byyear=F,      # do by year instead of decade; theperiods needs to be list where each element is a vector of year files, with full paths
                             onlyvocab=F, # only collects lexicon frequencies, does not save period files
                             untouchables=NULL, # optional vector of words to bypass filtering
                             FOLDER,
                             periodfolder="cohaperiods",
                             minc=100,
                             skipperiods=F # debug
){
  require("parallel")
  require("grr")
  require("compiler") # pre-compile functions for faster execution in the cluster:       
  setCompilerOptions("suppressAll" = T)
  fixcohaperiod = cmpfun(fixcohaperiod, options=list(optimize=3))
  coha_csv = cmpfun(coha_csv, options=list(optimize=3))
  try(dir.create(file.path(FOLDER, periodfolder) )) # try to create folder in case doesn't exist
  
  if(!skipperiods){
    # prep cluster
    nc = detectCores() - nfree
    cl = makeCluster(nc)#, outfile=file.path(FOLDER, periodfolder, "log.txt"))
    print(paste(Sys.time(), "start"))
    tryCatch({
      clusterExport(cl, c("fixcohaperiod", "coha_csv", "markS", "minlen", "path", "saveyear","rmstopwords", "removenp", "byyear", "onlyvocab", "FOLDER", "periodfolder"),envir = environment()) # exports params and compiled functions to cluster
      clusterEvalQ(cl, c(library(text2vec), library(magrittr), library(fastmatch) ) ) # exports packages to cluster
      
      tmp = parLapply(cl, theperiods, fixcohaperiod, 
                      markS = markS, minlen=minlen, path=path, saveyear=saveyear, rmstopwords=rmstopwords, removenp=removenp,byyear=byyear, onlyvocab=onlyvocab, untouchables=untouchables,
                      FOLDER=FOLDER,periodfolder=periodfolder
      )
    }, error=function(e){print(e)},  finally = stopCluster(cl) )
    
    # if(onlyvocab){   # not functional in this version
      #countmat = fixvocablist(tmp)
      #countmat = countmat[rowSums(countmat)>1, ] # remove hapaxes (single occurrence in entire corpus)
      #print(paste(Sys.time(), "corpus files done"))
      #return(countmat)
      #return(NULL)
    #} else {
    # }
    #print(paste(Sys.time(), "corpus done, doing counts"))
  }
  countmat = dotrends(foldr=periodfolder, mincount=1, nfree=nfree, rm_global_hapaxes=T, 
                      FOLDER=FOLDER, minc=minc)
  #print(paste(Sys.time(), "corpus and counts done"))
  return(countmat)
}



dotrends = function(foldr, mincount=1, smooth0=F, nfree=0, rm_global_hapaxes=T, FOLDER, minc){
  files = list.files( file.path(FOLDER, foldr), full.names = T, pattern="RData$")
  print(paste(Sys.time(), "Start counting words"))
  
  library(parallel)
  nc = detectCores() - nfree
  cl = makeCluster(nc)
  clusterExport(cl, c( "files", "minc"),envir = environment())
  clusterEvalQ(cl, c(library(text2vec),library(fastmatch)) )
  tryCatch({
    freqs = parLapply(cl, files, function(f){
      load(f)
      #period=period[1:2]# DEBUGGING
      it = itoken(period, progressbar = FALSE) # period list from load()
      rm(periods); gc() # clean periods - iterator envir now contains all data
      vocab2 <- create_vocabulary(it, ngram = c(ngram_min = 1L, ngram_max = 1L))
      vocab2 <- prune_vocabulary(vocab2, term_count_min = mincount)
      #docdisp = vocab2$doc_count/attr(vocab2,"document_count"); names(docdisp)=vocab2$term
      
      # now doing normalization in the end!
      # sm=sum(vocab2$term_count)
      freqs=vocab2$term_count
      names(freqs)=vocab2$term 
      return(freqs  )
    })
    
    allwords = unique(unlist(lapply(freqs, function(x) names(x) ),use.names = F ) )
    #print(paste(Sys.time(), "Period counts done"))
    
    clusterExport(cl, c("allwords", "freqs"),envir = environment())
    lexorder = parLapply(cl, 1:length(files), function(x){ 
      return(fmatch(allwords, names(freqs[[x]]) ) )
    })
    clusterExport(cl, c("lexorder"),envir = environment())
    
    print(paste(Sys.time(), "Aggregating counts"))
    
    countmat1 = parSapply(cl, 1:length(files), function(x){
      freqvec = freqs[[x]]
      tmp = freqvec[lexorder[[x]] ]; tmp[is.na(tmp)] = 0 # if no occurrence, then 0
      return(tmp)
    } )
    rownames(countmat1) = allwords
    
  }, error=function(e){print(e)},  finally = stopCluster(cl) )
  
  
  
  ## Fix the counts matrix, normalize, save info ##

  if(rm_global_hapaxes){ # remove useless words that occur just once per entire corpus
    countmat1 = countmat1[rowSums(countmat1)>1, ]
  }
  
  # save list of words that occur frequently enough to be included in future topic models
  #lexicon = rownames(countmat1)[apply(countmat1, 1, function(x) any(x>=minc) ) ]
  
  # normalize matrix columns to per-million
  thesums = colSums(countmat1)
  for(i in 1:ncol(countmat1)){
    countmat1[,i] =  (countmat1[,i]/thesums[i])*1000000
  }
  ones = (1/thesums)*1000000 # vector of (normalized) values that correspond to 1 occurrence 
  # for smoothing in log() frequency difference calculations
  
  countmat=list(countmat=countmat1, thesums=thesums, ones=ones )
  try({ save("countmat", file=file.path(FOLDER, "countmat.RData")) })
  print(paste(Sys.time(), "Done with counting words."))
  return(countmat)
}




fixvocablist = function(freqs,nfree=0, fromfiles=F, foldr){
  require(grr)
  require(fastmatch)
  require(parallel)
  nc = detectCores() - nfree
  cl = makeCluster(nc)
  
  if(fromfiles){
    files=list.files(foldr, full.names = T, pattern="RData$")
    clusterExport(cl, c("files"),envir = environment())
    clusterEvalQ(cl, c(library(text2vec),library(fastmatch)) )
    tryCatch({
      freqs = parLapply(cl, files, function(f){
        load(f)
        it = itoken(period, progressbar = FALSE) # period list from load()
        rm(period) # clean periods - iterator envir now contains all data
        vocab2 <- create_vocabulary(it, ngram = c(ngram_min = 1L, ngram_max = 1L))
        #vocab2 <- prune_vocabulary(vocab2, term_count_min = mint)
        #docdisp = vocab2$doc_count/attr(vocab2,"document_count"); names(docdisp)=vocab2$term
        #sm=sum(vocab2$term_count)
        freq=(vocab2$term_count) #/sm)*1000000; 
        names(freq)=vocab2$term
        return(freq)
      })
    })
  }
  wordl=lapply(freqs, function(x) names(x) ) 
  allwords = unique(.Internal(unlist(wordl, FALSE, FALSE)) )
  allwords = sort2(allwords) # alphabetical; requires grr
  rm(wordl)
  print(paste(Sys.time(), "lexicon done"))
  
  tryCatch({
    clusterEvalQ(cl, c(library(fastmatch)) )
    clusterExport(cl, c("allwords", "freqs"),envir = environment())
    
    lexorder = parLapply(cl, 1:length(freqs), function(x){ 
      return(fmatch(allwords, names(freqs[[x]]) ) )
    })
    clusterExport(cl, c("lexorder"),envir = environment())
    print(paste(Sys.time(), "lexorder done, aggregating counts"))
    
    countmat1 = parSapply(cl, 1:length(freqs), function(x){
      freqvec = freqs[[x]]
      tmp = freqvec[lexorder[[x]] ]; tmp[is.na(tmp)] = 0 # if no occurrence, then 0
      return(tmp)
    } )
  }, error=function(e){print(e)},  finally = stopCluster(cl) )
  rownames(countmat1) = allwords
  countmat1 = countmat1[rowSums(countmat1)>1, ]
  tryCatch({stopCluster(cl)}, error=function(e){}) # double check cluster closing
  return(countmat1)
}


dofreqdif = function(countmat, ones, uselog=T){
  countmat = round(countmat, 10) # to avoid any chance of floating point errors when looking for 0s
  freqdifmat = countmat
  freqdifmat[] = NA
  for(i in 2:ncol(countmat)){
    fr1 = countmat[,i-1] 
    fr2 = countmat[,i]   
    if(uselog){
      fr1[fr1==0] = ones[i-1]  # smoothing by pm value that equals 1 occurrence in that subcorpus
      fr2[fr2==0] = ones[i]
      x = log(fr2/fr1)  
      # Super important: return changes from 0 to 0 back into zeroes - the smoothing above
      # makes all zeroes non-zero, which would lead to weird numbers, so this fixes it:
      x = ifelse( (fr1+fr2)==0, 0, x)
      #
      freqdifmat[,i] = x
    } else {  # raw freq difference
      freqdifmat[,i] = fr2-fr1
    }
  }
  gc()
  return(freqdifmat)
}; dofreqdif = cmpfun(dofreqdif, options=list(optimize=3))


#### ppmi and advection ####


fullppmi = function(pmat, N, colp,rowp, positive=T){
  library(Matrix, quietly = T)
  # throws error is supplied vectors don't make sense:
  if(length(colp) != nrow(pmat) | length(rowp) != ncol(pmat) ) { 
    stop("mismatching norm vector length(s)")
  } 
  # note colp~nrow comparison - because the matrix is transposed for efficiency below
  
  pmat = Matrix::t(pmat)
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
    tmp = log2( (icol/N) / (rowp[not0] * colp[i] ) )
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
}; fullppmi = cmpfun(fullppmi, options=list(optimize=3))




doppmimat = function(years,winsize=5,minc = 100, ppmi=T, foldr){
  require(text2vec)
  require(Matrix)
  
  toload = years
  files=list.files(foldr, full.names = T, pattern = "RData$")
  plist = list()
  for(f in toload){
    load(files[f]) # load period
    plist = c(plist, period) # could optimize (but need to keep it 1-deep list!)
  }
  rm(period) # remove last
  it = itoken(plist, progressbar = FALSE); rm(plist)
  voc0 = prune_vocabulary(create_vocabulary(it), term_count_min = minc)
  vectorizer = vocab_vectorizer(voc0)
  tcm = create_tcm(it, vectorizer, skip_grams_window = winsize)
  tcm = tcm + Matrix::t(Matrix::triu(tcm)) ; Matrix::diag(tcm) = 0 
  if(ppmi){
    normvoc = voc0$term_count/sum(voc0$term_count)
    N=sum(tcm@x)
    rm(it,vectorizer, voc0) # free some memory
    ptcm = fullppmi(pmat=tcm, N=N, colp=normvoc,rowp=normvoc)
    return(ptcm)
  } else {
    return(tcm)
  }
}

dorelevants = function(ppmimat,relevancy_threshold){
  library(grr)
  relevants = list()
  print(paste(Sys.time(), "do relevance vectors"))
  #excludeself = match(rnames, colnames(ppmimat))  # might not have match in comlex ->NA
  #excludeself[is.na(excludeself)] = ncol(ppmimat)+100  # idea: those with no match get index > length of vector, so they can't be -indexed but won't be NA when -indexed, no need for ifelse
  
  # if(paral){
  # library(parallel)
  # nc = detectCores() - nfree
  # cl = makeCluster(nc)
  # clusterEvalQ(cl, library(Matrix)) 
  # clusterExport(cl, c("ppmimat","relevancy_threshold"),envir = environment()) 
  # print(paste(Sys.time(), "export ready, start"))
  # tryCatch({
  #   relevants = parLapply(cl, 1:nrow(ppmimat), 
  #       function(x){
  #         #y=sort(ppmimat[x,-(union(which(ppmimat[x,]==0), excludeself[x]))], decreasing=T)
  #         # self-ppmi is always zero due to how text2vec handles the tcm
  #         y=sort(ppmimat[x,-which(ppmimat[x,]==0)], decreasing=T) 
  #         return(names(y[1:min(length(y),relevancy_threshold)])) 
  #         }
  #       )
  # }, error=function(e){print(e)}, finally = stopCluster(cl))
  # }
  # 
  # if(!paral){
  relevants = list(); length(relevants) = nrow(ppmimat)
  names(relevants) = rownames(ppmimat)
  ppmimat = Matrix::as.matrix(ppmimat)
  for(x in 1:nrow(ppmimat)){
    #y=sort(ppmimat[x,-(union(which(ppmimat[x,]==0), excludeself[x]))], decreasing=T)
    # self-ppmi is zero, fixed in ppmimat (old: thought text2vec handles itself)
    tmp = ppmimat[x,-which(ppmimat[x,]==0)]
    y=tmp[rev(order2(tmp))] # sort2 is decreasing=F BUT MESSES UP NAMES, USE order2
    y=y[1:min(length(y),relevancy_threshold)] # but this need top ones, so rev
    relevants[[x]] = y; names(relevants[[x]]) = names(y)
  }
  print(paste(Sys.time(),"relevance vectors done"))
  return(relevants)
}; dorelevants = cmpfun(dorelevants, options=list(optimize=3))


doadvection = function(relevants, difvec, threshold){ 
  words = names(relevants)
  adv = rep(NA, length(words)); names(adv)=words
  #diffs = difvec[words]
  # advection from current values of present relevants
  for(i in 1:length(relevants)){ 
    nrelevants = 1:(min(threshold, length(relevants[[i]])))
    adv[i] = weighted.mean(difvec[ names(relevants[[i]][nrelevants] ) ],
                           w = relevants[[ i ]][nrelevants] , na.rm = T) 
  }
  return(adv)
}


do_advections = function(winsize, minc, ppmi=T, relevancy_threshold, smoothing=c(0), freqdifmat, foldr, pos=NULL){
  advs = vector("list", ncol(freqdifmat))
  for(i in 2:ncol(freqdifmat)){
    gc() # force memory management
    ii = i+smoothing  # adding previous period to list if smoothing
    ppmimat = doppmimat(years=ii, winsize=winsize,minc = minc, ppmi=ppmi,
                         foldr=foldr)
    
    if(!is.null(pos)){
      ppmimat = ppmimat[grep(pos, rownames(ppmimat)), , drop=F]
    }
    
    relevants = dorelevants(ppmimat, relevancy_threshold) 
    advs[[i]] = doadvection(relevants, difvec=freqdifmat[,i], threshold=relevancy_threshold) # both same threshold in this version
  }
  rm(ppmimat, relevants)
  allwords = unique(unlist(lapply(advs, function(x) names(x) ),use.names = F ) )
  
  advmat = matrix(NA, nrow=length(allwords), ncol=ncol(freqdifmat), dimnames=list(allwords, NULL))
  for(i in 2:ncol(freqdifmat)){
    advmat[,i] = advs[[i]][allwords]
  } 
  return(advmat)
}


#### New words ####

do_newwords_relevants = function(winsize, minc, ppmi=T, relevancy_threshold, freqdifmat, foldr, pos="^N:", targets=17:20){
  advs = vector("list", ncol(freqdifmat))

  ppmimat = doppmimat(years=targets, winsize=winsize, minc = minc, ppmi=ppmi, foldr=foldr)
  if(!is.null(pos)){
      ppmimat = ppmimat[grep(pos, rownames(ppmimat)), , drop=F]
  }
  relevants = dorelevants(ppmimat, relevancy_threshold) 
  return(relevants)
}

do_newadvections = function(countmat, newrels, freqdifmat1, xbefore=10){
  library(Rmisc)
  cm = round(countmat[[1]],10)
  freqdifmat1[freqdifmat1==0] = NA  # this avoids deflating the historical topic scores
  # (some of the topic words are at zero freq in the past, so their "change" would be 0)
  smat = cm[unique(names(newrels)),]
  entry=numeric(); xnew=character()
  for(i in 17:20){
    x1=rownames(smat)[apply(smat,1, function(x) x[i]>0 & all(x[1:(i-1)]==0) )]
    xnew = c(xnew, x1 )
    entry = c(entry, rep(i, length(x1)))
  }
  names(entry) = xnew
  nrels=newrels[xnew]
  newsmat = matrix(NA, ncol=20, nrow=length(xnew), dimnames=list(xnew,1:20))
  for(i in 2:20){
    newsmat[,i] = doadvection(nrels, freqdifmat1[,i], threshold=75)
  }
  namat = newsmat; namat[] = NA
  for(i in 1:length(xnew)){
    xrange =  (entry[i]-xbefore):(entry[i]-1)
    namat[i,xrange] = colSums(cm[names(nrels[xnew[i]][[1]]), xrange] > 0)
  }
  
  cilist = vector("list", length(xnew));names(cilist)=xnew
  for(i in 1:length(xnew)){
    xrange =  (entry[i]-xbefore):(entry[i]-1)
    cilist[[i]] = CI(newsmat[i, xrange])
  }
  
  return(list(newsmat=newsmat, namat=namat, entry=entry, cilist=cilist))
}


summarize_news = function(newslist){
  n=nrow(newslist$newsmat)
  newsmat=data.frame(tomean=rep(NA,n), val=rep(NA,n),row.names = rownames(newslist$newsmat), stringsAsFactors = F)
  for(i in 1:n){
    x = newslist$newsmat[i, newslist$entry[i] ] 
    newsmat$tomean[i] = x - newslist$cilist[[i]]["mean"]
    val="mid"
    if( x > newslist$cilist[[i]][1]) val = "above"
    if( x < newslist$cilist[[i]][3]) val = "below"
    newsmat$val[i] = val
  }
  return(newsmat[order(newsmat$tomean),])
}






#### LDA version ####


advectionLDA = function(slex, topic_vocab_distr, word_freq_change_per_topic, words_distr_across_topics){
  library(Matrix)
  ix = match(slex, colnames(topic_vocab_distr))
  if(any(is.na(ix))){stop("LDA lexicon <-> components mismatch")}
  adv = rep(NA, length(slex));names(adv) = slex
  
  # turn into sparse to improve speed
  topic_vocab_distr = Matrix(topic_vocab_distr, sparse=T) # default dgCMatrix
  word_freq_change_per_topic = Matrix(word_freq_change_per_topic, sparse=T)
  
  for(i in 1: length(slex)){
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
                       a=0.1, b=0.1,n_iter=5000, word_freq_change, pos){
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
  rm(dtm); gc() # cleanup
  
  # word_freq_change = freqdif
  
  # words_dist_across_topics  # for each word, how often it appears in topic 1, topic 2, topic 3 (normalised counts)
  words_distr_across_topics = t(normalize(t(lda_model$components), 'l1')) # [topics, words]
  lex = colnames(words_distr_across_topics)
  # for each topic, how often word 1, word 2 etc appears in it (technically/ideally this is phi in the LDA notation, where it will be smoothed by the beta parameter; you may be using normalised counts from your topic-word matrix representing the zs (ignoring documents))
  topic_vocab_distr  = normalize(lda_model$components,'l1')   # [topics, words]
  
  # Then, for each word: 
  word_freq_change_per_topic = words_distr_across_topics
  for(i in lex){
    word_freq_change_per_topic[,i] = words_distr_across_topics[,i] * word_freq_change[i]
  }
  
  if(!is.null(pos)){
    slex = grep(pos, lex, value=T) 
  } else { slex = lex}
  
  adv = advectionLDA(slex, topic_vocab_distr, word_freq_change_per_topic, words_distr_across_topics)
  
  #  if(length(interpol)==1){
  #dispersion_tops = mean(apply(words_distr_across_topics[,slex], 2, sd),na.rm=T) # do dispersion (standard dev across topics) while at it  -> remove mean if not for debug
  # } else { dispersion_tops = NA}
  
  return(adv) 
  #return(list(advection=adv, dispersion_docs=dispersion_docs,dispersion_tops=dispersion_tops ))
}



ldatester = function(foldr, 
                     smooth=c(-1, 0), topics=500, minc=100, 
                     freqdif=freqdifmat, 
                     nfree=1, periodsx = 2:20, pos=NULL
){
  rdats = rdats = list.files(foldr,  full.names = T, pattern = "RData$")
  
  library(parallel)
  nc = detectCores()-nfree
  cl = makeCluster(nc)
  print(paste(Sys.time(), "start"))
  tryCatch({
    clusterExport(cl, c("doLDAperiod", "advectionLDA"))  # global
    clusterExport(cl, c("smooth", "freqdif", "rdats", "topics", "minc", "pos"), envir = environment())  # local
    clusterEvalQ(cl, c(library(Matrix), library(text2vec) )) 
    ldalist = parLapply(cl, periodsx, function(x) {
      ldaresults = doLDAperiod(nperiod=x, rdats, interpol=smooth, 
                             topicsperword=F,topics=topics, minc=minc,
                             a=0.1, b=0.1,n_iter=5000, word_freq_change=freqdif[,x], pos=pos )
      return(ldaresults)
      } )
  }, error=function(e){print(e)},  finally = stopCluster(cl) )
  
  allwords = unique(unlist(lapply(ldalist, function(x) names(x) ),use.names = F ) )
  
  ldamat = matrix(NA, nrow=length(allwords), ncol=ncol(freqdif), dimnames=list(allwords, NULL))
  for(i in 2:ncol(freqdif)){
    ldamat[,i] = ldalist[[i-1]][allwords]   # list length is -1, as starts with 2
  } 
  
  print(paste(Sys.time(), "LDA version done"))
  return(ldamat) 
}


#### summary functions ####

do_r2 = function(advs, fm){
  library(beeswarm)
  corrs=vector("list", length(advs))
  for(li in seq_along(advs) ){
    xcorrs=rep(NA,ncol(advs[[li]] ))
    fmx = fm[rownames(advs[[li]]), ]
    for(cl in 2:ncol(advs[[li]])){
      xcorrs[cl] = cor(advs[[li]][,cl], fmx[, cl], use="pairw")^2
    }
    #beeswarm(xcorrs, ylim=c(0,1))
    corrs[[li]] = xcorrs
  }
  return(corrs)
}

do_adv_summary = function(advs){
  nu=c(); nn=c()
  for(li in seq_along(advs) ){
    nu[li]=nrow(advs[[li]][apply(advs[[li]], 1, function(x) any(!is.na(x))) ,])
    nn[li]=length(which(!is.na(c(advs[[li]]))))
  }
  return(c(nu,nn) )
}


do_corpus_ks=function(COHAPATH){
  # counts number of lines (~words) in each genre in each decade
  # and calculates successive Kullback-Leibler divergences
  library(entropy)
  gs = c("^fic", "^mag", "^news", "^nf")
  periods = list.files(COHAPATH)
  res = matrix(0, nrow=length(periods), ncol=length(gs), dimnames=list(periods,gs))
  for(i in seq_along(periods)){
    fs = list.files(file.path(COHAPATH, periods[i]))
    for(g in gs){
      ng = c()
      gfiles = grep(g, fs, value=T)
      for(gf in gfiles){
        ng[gf] = length(readLines(file.path(COHAPATH, periods[i], gf),warn = F))
      }
      res[i, g] = sum(ng)
    }
  }
  ks = c()
  res2 = ifelse(res==0, res+1, res)
  for(i in 2:nrow(res2)){
    x1 = res2[i-1,]/sum(res2[i-1,])
    x2 = res2[i,]/sum(res2[i,])
    ks[i] = KL.plugin(x1,x2)
  }
  return(ks)
}

newwordplots = function(w, newslist, newrels, ylab=F, nice=F){
  library(wordcloud)
  m=0.7#-(0.8/4)
  par( fg="gray40")
  if(ylab) par( mar=c(2,2.5,0.3,0.4), cex.axis=m, cex.lab=m, mgp=c(3, .3, 0))
  if(!ylab) par( mar=c(2,1.5,0.3,0.4), cex.axis=m, cex.lab=m, mgp=c(3, .3, 0))
  cols=brewer.pal(11,"Spectral")
  pcols= c("gray", colorRampPalette(cols)(19))
  
  plot(newslist$newsmat[w,], type="n",
       ylab="", xlab="", yaxt="n",xaxt="n", xlim=c(14,20),ylim=c(-0.1,0.57))
  grid(lty=1, lwd=1,col="gray95")
  ci=newslist$cilist[[w]]
  rect(10, ci[3],newslist$entry[w]-0.3, ci[1], border = NA, col="gray91")   #col = cols[7])
  lines(x=c(10,newslist$entry[w]-0.3), y=rep(ci[2],2), lwd=2, col = "gray70"  )
  points(newslist$entry[w],newslist$newsmat[w, newslist$entry[w]], col=pcols[newslist$entry[w]], pch=16,lwd=1, cex=3)
  lines(newslist$newsmat[w,], cex=1, pch=16, lwd=0.9, lty=3,type="b",
        col="gray35" ) #"black", )
  
  
  axis(side = 2, tck=-0.008,  las=2)
  # axis(side = 1, 10:20,labels=rep("",11),   tck=-0.008)
  axis(side = 1, seq(10,20,1),labels=paste0(seq(1900,2000,10),"s"),   tck=-0.008)
  if(ylab) mtext("advection (log topic change)",2, outer = F, line=1.6,cex=m, adj=0.1 , col="black")
  
  text(14, 0.56, gsub("N:","", w), cex=1.2, font=2, adj=c(0,1),col="black")
  lines(x=c(13,13.9), y=rep(0,2), col="black")
  lines(x=c(20.1,21), y=rep(0,2), col="black")
  
  #n=75
  # do.call(rgb,as.list(c(col2rgb(cols[9])/255)+0.2))
  # do.call(rgb,as.list(c(col2rgb(cols[8])/255)-0.3))
  if(nice){
    par( mar=c(0,0,0,0))
    cols2 = colorRampPalette(c("white", "gray25"))(11)[cut(c(4,newrels[[w]]),11,include.lowest = T)][-1]
    wordcloud(gsub("N:|^0.*","", unique(names(newrels[[w]]))), freq = newrels[[w]], min.freq = 1,scale = c(1.5,0.001), random.order = F,ordered.colors = T,  col=cols2, fixed.asp = F, rot.per = 0 )
  } else {
    
    p=par()$mar; p[1]=0; par(mar=p)
    nr = sort(newrels[[w]], decreasing = T)
    names(nr)=gsub("N:","",names(nr)); nr=names(nr); nr=nr[!grepl("0", nr)]
    
    nc=cumsum(nchar(nr)) %>% cut(., seq(0,max(.),  41)); nrs=5
    plot(NA, type="n", ylim=c(0.5,nrs+0.5),xlim=c(1,7), bty="n", xaxt="n", yaxt="n",ylab="",xlab="")
    if(ylab) mtext("top topic words",2, outer = F, line=1.6,cex=m, at=0.7, adj=0, col="black" )
    for(i in 1:nrs){
      l=which(as.numeric(nc)==i)
      text(par("usr")[1]+0.1, nrs-(i-1), paste0(paste0(nr[l], collapse=", "), 
                                                ifelse(i==nrs,", ...", ",")), cex=1, 
           col=do.call(rgb, as.list(c(0,0,0)+(i/8)) ), adj=c(0,0), font=3)
    }
  }
}
