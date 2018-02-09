# Next Generation Sequencing Analaysis

These set of pages aim to walk you through analysis of your Next Generation Sequencing.

Start with raw fastq files and generate publication quality figures.

Briefly this analysis can be broken down into:
1. Merging, trimming and quality filtering of reads - USEARCH
2. Deprelication and clustering to produce an OTU (or ZOTU) table - USEARCH
3. Exploratory diversity analysis in QIIME2
4. Taxonomy assignment in QIIME2
5. Statistical analysis and production of publication quality figures using RStudio and packages including Phyloseq, Microbiome and ggplot.

***
# References

Main software and packages used

* [USEARCH](https://www.drive5.com) Developer: Robert Edgar.
* [QIIME2](https://qiime2.org) Developers: J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight.
* [RStudio](http://www.rstudio.com/.) Developers: RStudio Team
* [Phyloseq](https://joey711.github.io/phyloseq/) Developers: Paul McMurdie and Susan Holmes.
* [Microbiome](http://microbiome.github.io/microbiome/): Developers: Leo Lahti, Sudarshan Shetty, Tineka Blake and Jarkko Salojarvi.
* [Tidyverse (including ggplot2)](https://www.tidyverse.org/): Developer: Hadley Wickham.

I have tried to reference and provide links where I have used other developers code and workflows but please email crpytick.lab@gmail if you feel we have not acknowledged or referenced correctly.
