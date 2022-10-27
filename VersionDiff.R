# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl
    # NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

#Run latex diff by hand in MS directory
shell(cmd="latexdiff PristurusMS-3.tex PristurusMS-2-DCA.tex > MS-diff-2DCA-3.tex")
tinytex::latexmk("MS-diff-2DCA-3.tex") 
