# TRACER

This repository facilitates the reproduction of results of the following paper:

Karacosta, L. G., Anchang, B., Ignatiadis, N., Kimmey, S. C., Benson, J. A., Shrager, J. B., Tibshirani, R., Bendall, S.C. & Plevritis, S. K. (2019). Mapping Lung Cancer Epithelial-Mesenchymal Transition States and Trajectories with Single-Cell Resolution. [bioRxiv, 570341](https://www.biorxiv.org/content/10.1101/570341v1)

## Basic requirements to generate plots

To generate the plots you need to have R/RStudio installed together with the following packages:
- `tidyverse`
- `wesanderson`
- `igraph`
- `ggraph`
- `cowplot`

Then you can navigate with RStudio to the vignette folder, set it as your working directory (Session->Set working directory) and knit the vignette "transition_analysis.Rmd". The output has already been precomputed and can be seen in "transition_analysis.html".


## Rerun bootstrap simulations

Rerunning the bootstrap analysis (the results of which has been cached and precomputed in the folder "vignette/boot_samples") requires a few more ingredients.

First, you will need a working [Julia](https://julialang.org/), version 1.2., installation with the packages [JuMP](https://github.com/JuliaOpt/JuMP.jl) and Gurobi installed. See e.g. https://github.com/JuliaOpt/Gurobi.jl for instructions for how to install Gurobi (free academic license available at http://www.gurobi.com/).

With these installed, you can install the SparseTransitions.jl package, either manually from this repository ("SparseTransitions.jl" folder) or using the Julia package manager: Start the Julia console, then press `]` to activate the Package manager mode and finally type: 

```
add https://github.com/nignatiadis/SparseTransitions.jl
```

Moving back to R, install the SparseTransition package, by typing the following in the R console (after installing the `devtools` package):

```
devtools::install_github("nignatiadis/TRACER/SparseTransitions")
```

Once it has been installed, you can load the R package as follows:

```
library(SparseTransitions)
jl_pkg_setup('/Applications/Julia-1.2.app/Contents/Resources/julia/bin')
```

You might have to change the last line from above, to match your own Julia installation. See the corresponding documentation file from `JuliaCall`, e.g. ?jl_pkg_setup for instructions on how to do this.


Finally, run the file "vignette/julia_boot_sampling.R", i.e. make sure your working directory is the "vignette" folder and execute:

```
source("julia_boot_sampling.R")
```

This generates all results from the bootstrap sampling and now you can proceed to regenerate the plots from Step 1 by knitting the vignette "transition_analysis.Rmd"

## Software versions
This has been tested with the following versions, although it should work with other versions too.

Julia 1.2 (Julia packages: Gurobi.jl v0.7.2, JuMP.jl v0.18.6)
R 3.6.1
Gurobi 8.1.0
RStudio 1.2.5001



