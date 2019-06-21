#### Scripts to accompany the paper "Quantifying the dynamics of topical fluctuations in language" ####
#
# Runs the advection model on COHA; can also be used on COCA 
# (and possibly other CLAWS-tagged corpora in the format of three-column csv, 
# where the columns are word, lemma, pos). 
# For other corpora, a custom corpus reader script would be required.
#
# Andres Karjus, University of Edinburgh
# andreskarjus.github.io

#### Preparations ####

# 1. Define parameters:
MINCOUNT = 100            # min frequency to include word
WINSIZE = 10              # ppmi model context window size
RELEVANCY_THRESHOLD = 75  # how many context words to include in advection model
SCRIPTS_PATH  = "advection_functions.R" # full path to the "advection_functions.R" script file
FREE_CORES = 0            # how many cpu cores to keep free when parallel-processing; change this to N-1 (where N is the number of cores), if you have small amounts of RAM (e.g. <4GB when parsing corpus; <16GB when running LDA)

# 2. Define the following paths and create folders:
COHAPATH = "/COHA/wlp"         # path to the COHA distribution wlp folder
FOLDER = ""                    # path to (possibly new) folder where intermediate files may be stored
try(dir.create(FOLDER))
PERIODFOLDER =  "cohaperiods"  # name of subfolder to be created, for period files

# 3. Load the functions and packages; will complain if packages missing and list them:
source(SCRIPTS_PATH)

# 4. Install the required packages:
# Uncomment and run this to install all required packages in one go:
# install.packages(c("magrittr", "text2vec", "compiler", "grr", "parallel", "Matrix", "fastmatch", "Rmisc", "xtable", "RColorBrewer","entropy","beeswarm"))


#### The following runs the model on COHA, provided all required packages are installed ####

# Parse and save corpus files, return (and also save) counts matrix
countmat = docorpus_parallel(theperiods=list.files(COHAPATH) , # decade folders
                  nfree=FREE_CORES,       # use all cpu cores for parallel
                  markS = list(N=T, V=F,A=F,D=F,M=F,P=F), # optional: POS class to prefix with a tag
                  minlen=3,      # min word length (exclude short ones) 
                  path=COHAPATH, 
                  # decades: decade folder path; years: need full paths to files
                  saveyear=T,    # use folder decade/file year as filename
                  rmstopwords=T, # filter out stopwords?
                  removenp=T,    # filter out proper nouns by removing all Capitalized?
                  byyear=F,      # do by year instead of decade? theperiods needs to be list where each element is a vector of year files, with full paths
                  FOLDER = FOLDER,
                  periodfolder=PERIODFOLDER,
                  minc = MINCOUNT, # frequency threshold for future advection models
                  skipperiods=F    # debug
)

freqdifmat = dofreqdif(countmat=countmat$countmat, 
                       ones=countmat$ones, # permil values that == 1 occurrence in each period
                       uselog=T)           # log dfference or raw difference

# PPMI-based model, without and with "smoothing" (contatenating data from target and preceding subcorpus)
advection0 = do_advections(WINSIZE, MINCOUNT, ppmi=T, RELEVANCY_THRESHOLD, smoothing=c(0), freqdifmat=freqdifmat, foldr=file.path(FOLDER, PERIODFOLDER), pos="^N:")

advection1 = do_advections(WINSIZE, MINCOUNT, ppmi=T, RELEVANCY_THRESHOLD, smoothing=c(-1,0), freqdifmat=freqdifmat, foldr=file.path(FOLDER, PERIODFOLDER), pos="^N:")

# LDA-based model, without and with "smoothing"

advection_lda0 = ldatester(foldr=file.path(FOLDER, PERIODFOLDER), smooth=c(0), topics=500, minc=100,  freqdif=freqdifmat, nfree=FREE_CORES, periodsx = 2:ncol(freqdifmat), pos="^N:")
advection_lda1 = ldatester(foldr=file.path(FOLDER, PERIODFOLDER), smooth=c(-1, 0), topics=500, minc=100,  freqdif=freqdifmat, nfree=FREE_CORES, periodsx = 2:ncol(freqdifmat),pos="^N:")


#### numbers for paper ####

do_adv_summary(list(advection0, advection1))

corrs = do_r2(list(advection0, advection1, advection_lda0, advection_lda1 ), freqdifmat)

