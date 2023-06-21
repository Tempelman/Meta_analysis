<p style="text-align: center;"> ================================================ </p>
<p style="text-align: center;"> Meta Analysis Workshop</p>
<p style="text-align: center;"> Robert J. Tempelman </p>
<p style="text-align: center;"> Michigan State University </p>
<p style="text-align: center;"> National Animal Nutrition Program Workshop</p>
<p style="text-align: center;"> 2023 American Dairy Science Association Meetings</p>
<p style="text-align: center;"> June 25, 2023</p>
<p style="text-align: center;"> Ottawa, Canada</p>
<p style="text-align: center;"> =============================================== </p>


These are the meta-analysis workshop materials developed for the NANP workshop June 25, 2023

The workshop slides can be found here:  [Slides](https://github.com/Tempelman/Meta_analysis/blob/main/TEMPELMAN_META_ANALYSIS.pdf) 

R and R studio will be used for this workshop.  R can be installed from https://cran.r-project.org/ whereas R studio can be installed from https://posit.co/download/rstudio-desktop/ .  Once you've installed R, be sure to install various packages in the following manner:

**<span style="font-family:courier; font-size:4em;">install.packages(c('tidyverse', 'ggplot2', 'multcomp','car', 'emmeans', 'broom','glmmTMB', 'metafor' ))</span>**

The applications will involve three different sets of codes/outputs

1. Analysis of St-Pierre (2001) data.  Open access to the corresponding publication can be found [here](https://pubmed.ncbi.nlm.nih.gov/11352149/).  The corresponding data in the Appendix can be found [here](https://github.com/Tempelman/Meta_analysis/blob/main/Dataregs2.csv) This data actually doesn't match what was actually analyzed in the paper so the actual corrected data is provided [here](https://github.com/Tempelman/Meta_analysis/blob/main/Dataregscorrected.csv).  The corresponding [regular R codes](https://github.com/Tempelman/Meta_analysis/blob/main/StPierre.R) and [R markdown codes](https://github.com/Tempelman/Meta_analysis/blob/main/StPierre.Rmd) are also provided.  I would suggest that you either run either the provided regular R codes and R markdown codes (the latter might be preferable) or just simply sit back and relax using the [R markdown file](https://rpubs.com/TEMPELMAN/1051562) provided.  You may need to be mindful of the directory name where you stored the input data on your laptop/PC.

2. Analysis of simulated completely randomized designs. The corresponding regular R codes are provided as a [regular script](https://github.com/Tempelman/Meta_analysis/blob/main/CRD_study.R) or a [markdown script](https://github.com/Tempelman/Meta_analysis/blob/main/CRD_study.Rmd).  Or again, you could just sit and back and watch from [here](https://rpubs.com/TEMPELMAN/1054063)

3. Analysis of simulated Latin square designs. The corresponding regular R codes are provided as a [regular script](https://github.com/Tempelman/Meta_analysis/blob/main/Latin_square.R) or a [markdown script](https://github.com/Tempelman/Meta_analysis/blob/main/Latin_square.Rmd).  Or again, you could just sit and back and watch from [here](https://rpubs.com/TEMPELMAN/1055615)  If time allows, this portion will also address network meta-analysis whereby both direct and indirect comparisons between two treatments are used to infer upon their mean differences using information on experiments where those treatments are not directly compared against each other.   

Some of you (as am I) are heavy SAS users.  The SAS mixed model procedures are particularly flexible for meta-analysis with some excellent expositions and SAS code provided by [Madden et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/27111798/), [Piepho et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/25410043/) and [Piepho and Madden (2022)](https://pubmed.ncbi.nlm.nih.gov/35638104/)  

