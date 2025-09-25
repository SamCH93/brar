# Bayesian Response Adaptive Randomization with Point Null Bayesian Hypothesis Testing

This repository contains 

1. `./package` The R package **brar** to conduct Bayesian response adaptive randomization

2. `./paper` Code and data to reproduce result from the paper: *Pawel, S., Held. L. (2025). Bayesian Response Adaptive Randomization with Point Null Bayesian Hypothesis Testing. <https://github.com/SamCH93/brar>*

To cite our work, use the following BibTeX reference

```BibTeX
@article{PawelHeld2025,
  year = {2025},
  author = {Samuel Pawel and Leonhard Held},
  title = {{Bayesian} Response Adaptive Randomization with Point Null {Bayesian} Hypothesis Testing}.
  url = {https://github.com/SamCH93/brar}
}
```

An interactively explorable simulation results dashboard is available at: <https://samch93.github.io/brar/>


## Reproducing the paper with Docker

Make sure to have Docker and Make installed, then run `make docker` from the
root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(http://localhost:8787), and the R scripts in `/paper`, for example,
`/paper/BFBRAR.R`, which contains all code for the results from the paper), can
be rerun. Make sure to change the working directory to `/paper` inside RStudio
Server before running the R scripts.
