![Build status: main](https://img.shields.io/github/workflow/status/SwissTPH/VivaxModelR/R-CMD-check/main?style=flat-square)
![Latest commit](https://img.shields.io/github/last-commit/SwissTPH/VivaxModelR/main?style=flat-square)
![coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/clchampag/f66fda014589a897f852edc85f4059fc/raw/vivax.json)

# VivaxModelR

Simulations and calculations on a compartmental model for *Plasmodium Vivax* dynamics.

This code is associated with the following manuscript:  
Champagne C., Gerhards M. , Lana J., García Espinosa B., Bradley C., González O., Cohen J.,Le Menach A., White M., Pothin E.  **"Using observed incidence to calibrate the transmission level of a mathematical model for *Plasmodium vivax* dynamics including case management and importation"**, Mathematical Biosciences, in press.  
The version used in the manuscript corresponds to released version [v1.0.1](https://github.com/SwissTPH/VivaxModelR/tree/v1.0.1).

The current version includes additional malaria control interventions, as detailed in the vignettes.

The package can be installed using the devtools package:  

```{r}
library(devtools)  
install_github("SwissTPH/VivaxModelR",build_vignettes = TRUE)
```

The scripts necessary to reproduce the manuscript figures are situated under demo/figures_UsingObservedIncidence.  

The vignette presenting how to use the package can be accessed with the following command:
```{r}
browseVignettes("VivaxModelR")
```
