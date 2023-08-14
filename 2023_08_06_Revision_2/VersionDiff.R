# Script to compare two Rmd files
 #Note, requires installation of 'latexdiff' in your tex components & a version of Perl
    # NOTE: need to install latexdiffr from: https://github.com/hughjonesd/latexdiffr

#Run latex diff by hand in MS directory
shell(cmd="latexdiff PristurusMS-R1-text.tex PristurusMS-R2-DCA-text.tex > TrackChangesText.tex")
tinytex::latexmk("TrackChangesText.tex") 

