# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl
    # NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

#Run latex diff by hand in MS directory
shell(cmd="latexdiff PristurusMS-2.tex PristurusMS-2-DCA.tex > MS-diff.tex")
tinytex::latexmk("MS-diff.tex") 

###################
# not working
#library(latexdiffr)
#latexdiff("2021-Evol-ConawayAdams-MS-3.Rmd", "2021-Evol-ConawayAdams-MS-3 - Copy.Rmd")
