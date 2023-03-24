# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl
    # NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

#Run latex diff by hand in MS directory
shell(cmd="latexdiff PristurusMS-5.tex PristurusMS-5-HTC-IM.tex > MS-diff.tex")
tinytex::latexmk("MS-diff.tex") 
