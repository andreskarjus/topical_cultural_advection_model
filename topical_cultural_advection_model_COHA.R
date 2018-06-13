#### Code to calculate advection values of words in the COHA (or COCA) corpus ####
# And some extras like the role of advection in lexical innovation. 
#
# Andres Karjus, University of Edinburgh
#
#
# Specify the parameters and run the entire script to produce advection values (and comparative frequency change values) for all words or only of a chosen POS
#
# Run this to install any missing packages if needed (recommended to make sure everything works):
install.packages(c("compiler", "parallel", "text2vec", "fastmatch", "grr", "Matrix", "Rmisc"))
# note: text2vec is in beta as of writing this; the following code will work with version 0.5.0, no guarantees for later versions.
#
# If you do not have access to COHA or prefer not to use R, feel free to explore the results using the pre-computed matrices based on COHA data, distributed as part of this repo (coha_advection.zip), which contains the matrix of log frequency changes and the matrix of advection values (using the default parameters specified here, i.e. including only nouns that occur >100 times in at least one of the (smoothed/concatenated) periods), and the list of historic advection values of new nouns from 1970-2009 across 10 decades preceding first occurrence of each noun.
#
# See the paper for more details: Karjus, Blythe, Kirby, Smith 2018: Quantifying the dynamics of topical fluctuations in language, https://arxiv.org/abs/1806.00699
# Also see how to apply this model to various culture datasets: https://andreskarjus.github.io/cultevol_tartu_slides, code here: https://github.com/andreskarjus/cultural_advection_TartuCE


#### Parameters ####

####### Important:
#
FOLDER = getwd() # use default working dir OR define a folder where to put the temporary processing files
# Regardless, make sure to put the advection_functions.R file into wherever FOLDER is so it can be sourced. If you did not change the line above, then it is whatever the working directory is, which is:
getwd()


COHAFOLDER = "COHA/wlp" # define the *full* path to the COHA distribution's wlp folder that contains the decade folders, which in turn contain the POS-tagged text files (3 columns, word-lemma-pos). Note that the same functions can also be used to parse COCA (but note the different folder structure and specify the paths accordingly).
#
##################


# Which decades to use; should correspond to COHA folder names:
theperiods = c('wlp_1810s_ksw', 'wlp_1820s_hsi', 'wlp_1830s_ksq', 'wlp_1840s_nsj', 'wlp_1850s_qis', 'wlp_1860s_los', 'wlp_1870s_hiw', 'wlp_1880s_bss', 'wlp_1890s_zna', 'wlp_1900s_nsi', 'wlp_1910s_hse', 'wlp_1920s_bsq', 'wlp_1930s_ney', 'wlp_1940s_gax', 'wlp_1950s_iie', 'wlp_1960s_vus', 'wlp_1970s_qrp', 'wlp_1980s_bsc', 'wlp_1990s_iuw', 'wlp_2000s_iey') 


print(all(file.exists(file.path(COHAFOLDER, theperiods)))) # if this does not print "TRUE", the COHA folders have not been found (likely misspelled path), so the following code will not run.

minraw = 100             # minimal occurrence threshold to include word
tcmwindow = 5            # +- context size for co-occurrence window
relevancy_threshold = 75 # how many top context words to use to generate a topic (in the PPMI based model)
# which POS to tag with POS prefixes, for the purposes of filtering them later; S corresponds to noun/substantive, V=verb, A=adjective, D=determiner, M=numeral, P=proper noun (conflicts with removenp); will be prefixed to the the word with a ':' separator:
markS=list(S=T,V=F,A=F,D=F,M=F,P=F) 
removenp = TRUE          # remove proper nouns (actually removes all capitalized words, in order to make sure none slip in; many PNs are tagged as Ns)
wordclassesr = "^[S]:"   # which POS to use, regex to define the prefix; a "." or an empty string will match everthing.
# "topic smoothing", or how many periods to concatenate for the purposes of generating topics. -1,0 corresponds to "preceeding t-1 and the current (t0) period". 0 disables:
topicsmoothing = c(-1,0) 
# how many cores to leave free when parallel processing; if you have issues with RAM (~<8GB), increase this value to use less memory; 0 uses all available cores:
nfree=0


# Run one or all of the following models by executing the chuncks of code below:


#### COHA, default decade bins, 20 of them: 1810s-2000s; nouns only, t-1 smoothing ####
#
# creates some folders to store the intermediate steps:
dir.create(FOLDER) # will fail with a warning if already present
sapply(c("cohaperiods", "cohaperiodrelevants", "cohaperiodtcm"), function(x) dir.create(file.path(FOLDER, x)))
source(file.path(FOLDER, "advection_functions.R"))  # load functions; if this fails put the R script in the FOLDER folder


# This parses the entire corpus and stores the data as lists within RData files in a specified folder. This needs to be run only once, unless different parsing behaviour is desired.
coha_parallel(theperiods= theperiods, path=COHAFOLDER, nfree=nfree, markS = markS, removenp = removenp, rmstopwords = T, FOLDER=FOLDER) 
countmat = dotrends(foldr = file.path(FOLDER, "cohaperiods"), mint=2, nfree=0)

# These calls below use the parsed data (above) to construct the PPMI-based topic vectors and calculate advection. These make use of the parsed files stored in cohaperiods folder.
perioddatas(rdats = list.files(file.path(FOLDER, "cohaperiods"), full.names = T), 
            comlex=NA, minc=minraw, tcmwindow = tcmwindow, relevancy_threshold=relevancy_threshold, nfree=nfree, wordclassesr=wordclassesr, interpol=topicsmoothing, tcmfolder=file.path(FOLDER,"cohaperiodtcm/"), relefolder=file.path(FOLDER,"cohaperiodrelevants/"))

allrelevants = aggregaterelevants(file.path(FOLDER, "cohaperiodrelevants"))
freqdifmat1 = dofreqdif(countmat[[1]], smooth=1, rownames(countmat[[1]])) 

# advection calculation:
advection1 = dotopicsloop(periods=allrelevants, difmat = freqdifmat1, countmat=countmat[[1]], wordclassesr=wordclassesr, threshold=relevancy_threshold, estimate0 = F, estimate2 = F)

# correlation:
print(cor(as.vector(advection1$inertiamat), as.vector(freqdifmat1[rownames(advection1$inertiamat),]), use="pairw")^2 ) # R^2

# plot:
m=0.8; par(mar=c(3,3,0.1,0.1), cex.axis=m, cex.lab=m)
plot(as.vector(advection1$inertiamat)~ as.vector(freqdifmat1[rownames(advection1$inertiamat),]), cex=0.3, pch=".", col=rgb(0,0,0,0.5), ylab="", xlab="", yaxt="n",xaxt="n", type="p", xlim=c(-4,4.2),ylim=c(-1,1.2))
axis(side = 2, at=seq(-1,1.2,0.5)); axis(side = 1, at=seq(-4,4.2,1))
mtext("log frequency change", 1, outer = F,line=2,cex=m); mtext("advection (log topic frequency change)",2, outer = F, line=2,cex=m); abline(v=0,h=0,col=rgb(0,0,0,0.3))





#### COHA, model for evaluating advection for lexical innovations in 1970s-2000s ####
#
# Assumes cohaperiods folder is populated with the parsed COHA files; if not, run the coha_parallel() command above.
# Then run the following:
#
# will load or create the counts matrix if needed:
if(exists("countmat")){
  "countmat in workspace"
} else { if(file.exists(file.path(FOLDER, "countmat.RData"))){
  load(file.path(FOLDER, "countmat.RData"))
} else { countmat = dotrends(foldr = file.path(FOLDER, "cohaperiods"), mint=2, nfree=0)  
} 
}

dir.create(file.path(FOLDER,"cohaperiodtcmnews")); dir.create(file.path(FOLDER,"relevantsnews"))
perioddatas(rdats = list.files(file.path(FOLDER, "cohaperiods"), full.names = T)[17:20], 
            comlex=NA, minc=minraw, tcmwindow = tcmwindow, relevancy_threshold=relevancy_threshold, nfree=nfree, wordclassesr=wordclassesr, interpol=c(-3,-2,-1, 0), tcmfolder=file.path(FOLDER,"cohaperiodtcmnews/"), relefolder=file.path(FOLDER,"relevantsnews/"),yearnames=seq(1970,2000,10))

newsrelevants = aggregaterelevants(file.path(FOLDER, "relevantsnews"), pattern="2000.RData") # need only the last one (that smoothes over previous periods)

newslist = newadvections(countmat, newsrelevants, freqdifmat1, relevancy_threshold, 
                         xbefore=10) # how many decades before entry decade to consider
library(Rmisc)
above=length(which(sapply(newslist, function(x) CI(x[1:(length(x)-1)])[1] < x[length(x)] )))
within=length(which(sapply(newslist, function(x) CI(x[1:(length(x)-1)])[1] >= x[length(x)]
                       & CI(x[1:(length(x)-1)])[3] <= x[length(x)] )))
below=length(which(sapply(newslist, function(x) CI(x[1:(length(x)-1)])[3] > x[length(x)] )))

# results:
length(newslist) # number of successful new nouns, out of which, respectively:
above; within; below
"of the confidence interval around the mean of the previous advection values of that noun's topic."

## test: z-scores based on mean and sd of each distribution, then do 1-sample t-test
# to test if the set of advection values at entry decades is significantly different from 0
zeds=rep(NA, length(newslist));names(zeds)=names(newslist)
for(i in 1:length(newslist)){
  ni = newslist[[i]]
  l = length(ni)
  zeds[i] = (ni[l] -  mean(ni[1:(l-1)]) ) / sd(ni[1:(l-1)])
}
hist(zeds,20)
t.test(zeds) # 



#### COHA, time series decomposition example ####

# Run the COHA advection section first to generate the matrices.
wrds = c("S:car", "S:payment", "S:negro", "S:happiness") # example nouns with interesting trajectories
par(mfrow=c(1,4), mar=c(3,3,2,1.5))
for(x in wrds){ 
  plotword(x,freqdif = freqdifmat1[,3:20], countm = countmat[[1]][,3:20], adv=advection1$inertiamat[,3:20], sem=F, ylims=c(1,2000), xlims=c(1,11), uselims=F, decades=seq(1830,2000,10))
  r=round(cor(freqdifmat[x,3:20], advection1$inertiamat[x,3:20], use="pairw"), 2)
  text(17,1, paste0("cor=",r), col=rgb(0,0,0,0.7), cex=0.9)
  par(mar=c(3,1,2,1.5))
}
mtext("decades",1, outer = T, line=-1,cex=0.8); mtext(2, text = "permillion frequency", cex=0.8, line=-1, outer=T)





#### Calculate advection using topics from an LDA model instead ####

# Assumes cohaperiods folder is populated with the parsed COHA files; if not, run the coha_parallel() command above.
# The following example calculates the advection values for period 20 (the 2000s), using t-1 smoothing, k=500 latent topics.
lda_advection_20 = doLDAperiod(nperiod=20, 
                               rdats = list.files(file.path(FOLDER, "cohaperiods")), 
                               interpol=c(-1,0), minc=100,
                               topicsperword=F, topics=500, 
                               a=0.1, b=0.1, n_iter=5000, 
                               word_freq_change=freqdifmat1[,20])

print(cor(lda_advection_20, freqdifmat1[names(lda_advection_20), 20], use="pairw")^2 ) # R^2
