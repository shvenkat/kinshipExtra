#!/bin/bash

echo "BUILDING PACKAGE DOCUMENTATION ..."
R --slave -e 'library(roxygen2); roxygenize(clean = TRUE)'
echo "BUILDING PACKAGE TARBALL in parent directory"
R --slave -e 'library(devtools); build()'
echo "DONE. Package can be installed with R CMD install"
