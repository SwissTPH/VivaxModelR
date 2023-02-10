![Build status: main](https://img.shields.io/github/actions/workflow/status/SwissTPH/VivaxModelR/r.yml?branch=main&style=flat-square)
![Latest commit](https://img.shields.io/github/last-commit/SwissTPH/VivaxModelR/main?style=flat-square)
![coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/clchampag/691fea8285290758f43b48ce17806edd/raw/vivax_sto.json)

# VivaxModelR

Simulations and calculations on a compartmental model for *Plasmodium Vivax* dynamics.

This code is associated with the following manuscripts:
- Champagne C., Gerhards M. , Lana J., Le Menach A., Pothin E.  **"Quantifying the impact of interventions against \textit{Plasmodium vivax} malaria: a model for country-specific use"**.  
The version used in this manuscript corresponds to the [main version](https://github.com/SwissTPH/VivaxModelR). It includes various control interventions, as detailed in the user [vignette](https://swisstph.github.io/VivaxModelR/articles/vivaxmodelr_delay.html). The scripts necessary to reproduce the manuscript figures are situated under [demo/figures_SimulatedExample](https://github.com/SwissTPH/VivaxModelR/tree/main/demo/figures_SimulatedExample). 

- [Champagne C., Gerhards M. , Lana J., García Espinosa B., Bradley C., González O., Cohen J.,Le Menach A., White M., Pothin E.  **"Using observed incidence to calibrate the transmission level of a mathematical model for *Plasmodium vivax* dynamics including case management and importation"**, Mathematical Biosciences, 2021](https://www.sciencedirect.com/science/article/pii/S0025556421001541).  
The version used in the manuscript corresponds to released version [v1.0.1](https://github.com/SwissTPH/VivaxModelR/tree/v1.0.1), and only focuses on case management interventions, without delayed access to treatment. The scripts necessary to reproduce the manuscript figures are situated under [demo/figures_UsingObservedIncidence](https://github.com/SwissTPH/VivaxModelR/tree/main/demo/figures_UsingObservedIncidence) (to be used with version [v1.0.1](https://github.com/SwissTPH/VivaxModelR/tree/v1.0.1)).  



The package can be installed using the devtools package:  

```{r}
library(devtools)  
install_github("SwissTPH/VivaxModelR")
```

The vignette presenting how to use the package is available [here](https://swisstph.github.io/VivaxModelR/articles/vivaxmodelr_delay.html).  
A simplified version of the package without delayed access to treatment is available [here](https://swisstph.github.io/VivaxModelR/articles/vivaxmodelr.html).  
A quick description for developers of the scripts in the package is available [here](https://swisstph.github.io/VivaxModelR/articles/vivaxmodelr_dev_doc.html).

