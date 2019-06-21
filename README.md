# The topical-cultural advection model

R scripts to calculate the topical-cultural advection values of words in the Corpus of Historical American English (COHA), as described in: _Andres Karjus, Richard A. Blythe, Simon Kirby, Kenny Smith, 2018. Quantifying the dynamics of topical fluctuations in language. arXiv preprint. https://arxiv.org/abs/1806.00699_ (and upcoming in Language Dynamics and Change).

To successfully run these scripts, you will need a distribution of COHA and a number of R packages (listed in the script).

It also works on COCA and possibly other CLAWS-tagged corpora in the format of three-column csv (where the columns are word, lemma, pos). For other corpora, a custom corpus reader/parser script would be required.

If you want to use these scripts on a different corpus/dataset that is already parsed, organize the data in the following way. Each period (year, decade, etc., depending on binning) should be stored as an RData file in the "periods" folder (ordered alphanumerically). Each such file should contain a single R `list` object named "period". Each element in that list - a character vector where each element is a tokenized word or lemma - is treated as a "document" in the NLP sense (which may be a literal document, a newspaper article, a book, or in other domains, some self-contained grouping of elements). 

The advection model is also being used in an upcoming work on lexical competition by the same authors as this paper, as a proxy to communicative need.
Besides language, the advection model has also been applied to various culture datasets (of boardgame mechanics, movie tags and cookbook recipes). These analyses are not written up yet, but there are some slides here: https://andreskarjus.github.io/cultevol_tartu_slides and the corresponding git repo: https://github.com/andreskarjus/cultural_advection_TartuCE 
