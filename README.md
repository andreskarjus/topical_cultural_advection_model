# The topical-cultural advection model

R scripts to calculate the topical-cultural advection values of words in the Corpus of Historical American English (COHA), as described in: _Andres Karjus, Richard A. Blythe, Simon Kirby, Kenny Smith, 2018. Quantifying the dynamics of topical fluctuations in language. arXiv preprint. https://arxiv.org/abs/1806.00699_ 

To successfully run these scripts, you will need a distribution of COHA and a number of R packages (see topical_cultural_advection_model_COHA.R for details and instructions on how to make everything work). Alternatively, feel free to browse the pre-computed matrices of frequency change and advection values of all nouns or for the ones featured in the lexical innovation analysis (see the zipped file).

If you want to use these scripts on a different corpus/dataset, organize the data in the following way. Each period (year, decade, etc., depending on binning) should be stored as an RData file in the "periods" folder (ordered alphanumerically). Each such file should contain a single R `list` object labelled "period". Each element in that list - a character vector where each element is a tokenized word or lemma - is treated as a "document" in the NLP sense (which may be a literal document, a newspaper article, a book, or in other domains, some self-contained grouping of elements). Importantly, the functions that construct the co-occurrence matrices (which are used to constuct topics) do not cross document boundaries when aggregating co-occurring words/elements (the window size is specified by the `tcmwindow` parameter). So e.g. to avoid co-occurrence windows crossing sentence borders, make each sentence a list element instead. By default the windows are linearly distance-weighted (using text2vec defaults), but uniform weighting (set `perioddatas(windowweights="uniform")`) and a larger window may be desirable if the list element/document is a set (i.e. the order or distance of elements does not matter).

Besides language, the advection model has also been applied to various culture datasets (of boardgame mechanics, movie tags and cookbook recipes), see https://andreskarjus.github.io/cultevol_tartu_slides and the corresponding git repo: https://github.com/andreskarjus/cultural_advection_TartuCE

bib entry of the paper:
```
@ARTICLE{karjus_quantifying_2018,
  author = {{Karjus}, A. and {Blythe}, R.~A. and {Kirby}, S. and {Smith}, K.},
  title = "{Quantifying the dynamics of topical fluctuations in language}",
  journal = {ArXiv e-prints},
  archivePrefix = "arXiv",
  eprint = {1806.00699},
  primaryClass = "cs.CL",
  year = 2018
}
```