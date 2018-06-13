# The topical-cultural advection model

R scripts to calculate the topical-cultural advection values of words in the Corpus of Historical American English (COHA), as described in: _Andres Karjus, Richard A. Blythe, Simon Kirby, Kenny Smith, 2018. Quantifying the dynamics of topical fluctuations in language. arXiv preprint. https://arxiv.org/abs/1806.00699_ 

To run these scripts, you will need a distribution of COHA and a number of R packages (see topical_cultural_advection_model_COHA.R for details and instructions on how to make everything work). Alternatively, feel free to browse the pre-computed matrices of frequency change and advection values of all nouns or for the ones featured in the lexical innovation analysis (see the zipped file).

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