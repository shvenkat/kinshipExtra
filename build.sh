#!/bin/bash

if [[ $# -ne 0 ]]; then
    echo "ERROR: arguments are not supported"
    exit 2
fi
if [[ ! -f "DESCRIPTION" ]]; then
    echo "ERROR: Current directory is not the top level directory of an R" \
      "package. DESCRIPTION file not found."
    exit 2
fi

echo -n "Building package documentation ... "
R --slave <<-EOF 1>&2 && echo "done" || { echo "failed"; exit 1; }
    suppressPackageStartupMessages(library(roxygen2))
    roxygenize(clean = TRUE)
EOF
if [[ $(uname -s) =~ "CYGWIN" ]]; then
    echo "Detected Cygwin environment, applying charset fix"
    sed -e 's/\xc2\xa0/ /g' -i man/*.Rd
fi
echo -n "Building package tarball ... "
R --slave <<-EOF 1>&2 && echo "done" || { echo "failed"; exit 1; }
    suppressPackageStartupMessages(library(devtools))
    build()
EOF
echo "See parent directory for package tarball." \
    "Install using R CMD INSTALL."
