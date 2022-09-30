# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl
    # NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

#Run latex diff by hand in MS directory
shell(cmd="latexdiff 2021-Evol-ConawayAdams-MS-Revised1.tex 2021-Evol-ConawayAdams-MS-R2.tex > diff.tex")
tinytex::latexmk("diff.tex") 

###################
# not working
#library(latexdiffr)
#latexdiff("2021-Evol-ConawayAdams-MS-3.Rmd", "2021-Evol-ConawayAdams-MS-3 - Copy.Rmd")
