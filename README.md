# The Price of Dark Traits: Strategic Exploitation and Its Limitations in Repeated Interactions

## Overview

This repository contains the data and analysis code for the study investigating how the Dark Factor of Personality (D-factor) influences strategic behavior in repeated trust games.

## Abstract

We investigated how individuals high in the Dark Factor of Personality (D) behave in repeated economic exchanges. Participants (N=183) played as trustees in 25-round trust games against computerized investors implemented as Hidden Markov Models (HMMs) that either behaved predictably ("Human-like") or responded quickly to exploitation ("Responsive"). High-D individuals showed strategic exploitation patterns: they initially cooperated but gradually reduced returns when facing predictable, trusting partners. However, when partners punished selfish behavior swiftly, even high-D individuals restrained themselves. Importantly, this strategic exploitation did not result in higher earnings. Our findings suggest that high-D individuals engage in calculated social behavior that adapts to circumstances rather than indiscriminate selfishness.

## Repository Contents

```
├── article.Rmd          # Main manuscript (R Markdown)
├── article.pdf          # Compiled manuscript
├── supplement.Rmd       # Supplementary materials (R Markdown)
├── supplement.pdf       # Compiled supplement
├── bib/
│   └── DFP2.bib         # Bibliography
├── data/
│   ├── anonymized_final_data.csv      # Main experimental data
│   ├── anonymized_demographics.csv    # Participant demographics
│   ├── mod_return_pct.RDS             # Pre-fitted mixed model
│   └── mod_return_pct_sensitivity.RDS # Sensitivity analysis model
└── plots/
    └── figure1.pdf      # HMM schematic figure
```

## Data Description

### anonymized_final_data.csv

Trial-level data from the repeated trust game experiment:

- `playerId`: Anonymized participant identifier
- `gameNum.f`: Game phase (1 = first game, 2 = second game)
- `roundNum`: Round number (1-25)
- `investment`: Amount invested by HMM (0-20 units)
- `return`: Amount returned by participant
- `d_level`: D-factor group (high_D, low_D)
- `volatile_first`: Whether participant faced Responsive HMM first (TRUE/FALSE)

### anonymized_demographics.csv

Participant demographic information including age, gender, and D-factor scores.

## Reproducing the Analysis

### Requirements

- R (>= 4.0)
- RStudio (recommended)
- Required R packages: tidyverse, papaja, afex, emmeans, lmerTest, depmixS4, ggplot2, kableExtra

### Instructions

1. Clone this repository
2. Open `d_factor.Rproj` in RStudio
3. Run `rmarkdown::render("article.Rmd")` to compile the manuscript
4. Run `rmarkdown::render("supplement.Rmd")` to compile the supplement

Note: The pre-fitted model objects (`.RDS` files) are included to reduce computation time. The full analysis can be reproduced from the raw data by modifying the `eval` parameters in the R Markdown files.

## License

This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0).
