# Stabilizing Thompson Sampling with Null Hypothesis Bayesian Response-Adaptive Randomization

This repository contains 

1. `./package` The R package **brar** to conduct Bayesian response-adaptive randomization

2. `./paper` Code and data to reproduce result from the paper: *Pawel, S., Held. L. (2025). Stabilizing Thompson Sampling with Null Hypothesis Bayesian Response-Adaptive Randomization. <https://doi.org/10.48550/arXiv.2510.01734>*

To cite our work, use the following BibTeX reference

```BibTeX
@article{PawelHeld2025,
  year = {2025},
  author = {Samuel Pawel and Leonhard Held},
  title = {Stabilizing {Thompson} Sampling with Null Hypothesis {Bayesian} Response-Adaptive Randomization},
  url = {https://github.com/SamCH93/brar},
  doi = {10.48550/arXiv.2510.01734}
}
```

An interactively explorable simulation results dashboard is available at: <https://samch93.github.io/brar/>

## Reproducing the paper results with a local R installation

Make sure to have a recent version of R installed (only tested with version >
4.5.2). Install the required packages with the following R command

```r
install.packages(c("dplyr", "SimDesign", "brar", "mvtnorm", "ggplot2", "ggh4x",
                   "ggpubr", "Ternary", "RARtrials", "knitr"))
```

Then run the R script `/paper/BFBRAR.R` to reproduce the results from the paper.
Make sure that the working directory is set to `/paper`. This uses pre-saved
simulation results. To reproduce the simulation study from scratch, run the R
script `paper/simulation/simulation.R`, but be aware that the simulation took
several days to run on a server with 100 CPU cores. Run the R script
`paper/simulation/simulation-summaries.R` to pre-process the simulation results
for the paper and/or the Quarto results dashboard
(`paper/simulation/simulation-dashboard.Qmd`).

The R session info (R version, platform, and package versions) used to produce
the results in the paper is recorded at the end of the compiled PDF manuscript
(`/paper/BFBRAR.pdf`).

## Reproducing the paper with Docker

Make sure to have Docker and Make installed, then run `make docker-rstudio` from
the root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(http://localhost:8787), and the R scripts in `/paper`, for example,
`/paper/BFBRAR.R`, which contains all code for the results from the paper), can
be rerun. Make sure to change the working directory to `/paper` inside RStudio
Server before running the R scripts. Running `make docker-paper` produces the
`paper/BFBRAR.tex` file from the `paper/BFBRAR.Rnw` source file (dynamically
inserting numbers and figures) and then compiles it to a PDF (requires a local
LaTeX installation; only tested with TeX Live 2023/Debian).

## License

The code is released under the [GPL-3 License](LICENSE).

## Contact

For questions or issues, please contact Samuel Pawel
([samuel.pawel@uzh.ch](mailto:samuel.pawel@uzh.ch)).
