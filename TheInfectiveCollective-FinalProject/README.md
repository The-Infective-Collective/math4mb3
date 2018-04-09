# Persistence and Synchrony of Recurrent Epidemics
MATH 4MB3 Final Project, by The Infective Collective
* Aurora Basinski-Ferris (basinsa@mcmaster.ca)
* Michael Chong (chongmy@mcmaster.ca)
* Sang Woo Park (parksw3@mcmaster.ca)
* Daniel Presta (prestad@mcmaster.ca)

# Recipe to recreate the main .pdf

`TheInfectiveCollectiveProject.pdf` is the main output file for this project and contains the main paper.
This file is constructed from `TheInfectiveCollectiveProject.tex`, (which can be compiled using your favourite TeX editor) and depend on figures in the `supplementary` folder, namely:
* `supplementary/bifurcation.pdf`
* `supplementary/probabilitycoherence.pdf`
* `supplementary/stochastic.pdf`
* `supplementary/stochastic_illustrate.pdf`

These above 4 figures can be reproduced using the supplementary information source file, `supplementary/supp.Rnw`.
You can run `supplementary/supp.Rnw` by opening it in RStudio and clicking "Compile PDF". This produces `supplementary/supp.pdf`, which contains all the supplementary material and further explanations to accompany the main paper, but also saves the figures for use by the main `.tex` file.

## Important Note regarding computational time
 We highly recommend that the user only compile `supplementary/supp.Rnw` if the following data files are present:
 * `supplementary/coherence_m0.001-0.01.rda`
 * `supplementary/coherence_m0.1-0.5.rda`
 * `supplementary/coherence_m0.0001.rda`
 * `supplementary/stochastic.rda`
 * `supplementary/bifurcation.rda`

 These data files contain statistics generated from many, many simulations and are therefore extremely time-consuming (in fact, so long that we didn't have time to run them again to figure out how long they took; at least 25 hours on a modern student-budget laptop).  If for some reason you do not have access to these files, please contact any of the authors for a copy.
