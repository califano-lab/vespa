# vespa - Virtual Enrichment-based Signaling Protein-activity Analysis

## Installation
The ``vespa`` R-packages requires installation of several dependencies, including the [``vespa.db``](https://github.com/califano-lab/vespa.db) R-package. In addition, it requires roxygen2, devtools and [pkghelper](https://github.com/markusdumke/pkghelper) to be installed to generate all necessary files for building the documentation.

Within R execute the following commands:
```
# install dependencies from CRAN
install.packages(c('devtools','data.table','pbapply','seqinr','stringr','mixtools','plyr','reshape2','tidyr','BiocManager','caret') ,dependencies=TRUE, repos='http://cran.rstudio.com/')

# install dependencies from BioConductor
BiocManager::install(c('viper','preprocessCore','limma'))

# install vespa
library(devtools)

install_github("califano-lab/vespa.db")

install_github("califano-lab/vespa")
```

## Usage
For further information on usage, refer to the [`vespa.tutorial`](https://github.com/califano-lab/vespa.tutorial) documentation, as well as the R-package help pages for information on individual functions:

```
library(vespa)

help(package="vespa")
```

