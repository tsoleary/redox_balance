# Heat stress drives opposing redox shifts in temperate versus tropical Drosophila melanogaster embryos

## Authors

Thomas S. O’Leary and Brent L. Lockwood

### Abstract

Redox balance is central to aerobic metabolism, yet acute heat stress can destabilize this balance by increasing metabolic rates and shifting the balance of critical electron carriers such as NADH. In early Drosophila melanogaster embryos, maintaining redox balance is particularly critical as embryos undergo a developmental redox shift and rely on oxidative phosphorylation to power nuclear divisions. Here, we assayed six isofemale D. melanogaster lines from temperate (Vermont, USA; France; Japan) and tropical (St. Kitts; Ghana; India) climates to assess metabolic responses to heat in heat-sensitive versus heat-tolerant embryos. We used untargeted LC–MS to measure 33 metabolites and the major redox couples (NADH/NAD+, NADPH/NADP+, and GSH/GSSG) at 25°C and after a 32°C heat shock. In all embryos, heat shock induced shared shifts in metabolic profiles, with increases in nucleotide monophosphates (e.g., AMP, CMP, and GMP) and amino acids (e.g., alanine, glutamic acid, serine). In contrast, redox metabolites diverged by region: heat-sensitive temperate embryos shifted toward a more oxidized state (46.6% decrease in NADH/NAD+ ratio and 4-fold increase in oxidized glutathione), while heat-tolerant tropical embryos maintained glutathione balance and increased the NADH/NAD+ ratio by 52.9%, indicating a more reduced state. These patterns are consistent with higher NADH oxidation and greater oxidative stress (inferred from oxidized glutathione) in the temperate embryos, versus better maintenance of redox balance in tropical embryos. Together, our results suggest that maintaining redox balance is a key determinant of acute heat tolerance, and healthy development overall, during early embryogenesis.

### Project structure

This project is organized into the following directory structure:

- `data/` - all data used in code
  - `raw/` - all raw data
  - `processed/` - all processed data
- `src/` - source code used for analyses -- files saved with numerical prefix
  - `01_metabolomics/` - source code related to metabolomics analysis
  - `02_figs/` - source code related to figures
- `output/` - all output files 
  - `figs/` - all figures generated
  - `metabolomics/` - all results from metabolomics analyses
- `scratch/` - any temporary or extraneous files and information
  - `writing` - early writing scraps (_e.g._ methods sections) to be compiled later
  
    
