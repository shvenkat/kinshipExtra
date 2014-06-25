#!/bin/bash

echo "BUILDING PACKAGE DOCUMENTATION ..."
R --slave -e 'library(roxygen2); roxygenize(clean = TRUE)'
if [[ $(uname -s) =~ "CYGWIN_" ]]; then
    sed -e 's/\xc2\xa0/ /g' -i man/*.Rd
fi
echo "BUILDING PACKAGE TARBALL in parent directory"
R --slave -e 'library(devtools); build()'
echo "DONE. Package can be installed with R CMD INSTALL"
