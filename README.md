# Heat stress drives contrasting redox shifts in temperate and tropical Drosophila melanogaster embryos 

## Authors

Thomas S. O’Leary and Brent L. Lockwood

### Abstract

Redox balance is central to aerobic metabolism, yet acute heat stress can destabilize this balance by increasing metabolic rates and shifting the balance of critical electron carriers such as NADH. In early Drosophila melanogaster embryos, maintaining redox balance is particularly critical as embryos undergo a developmental redox shift and rely on oxidative phosphorylation to power nuclear divisions. Here, we assayed six isofemale D. melanogaster lines from temperate (Vermont, USA; France; Japan) and tropical (Mexico; Ghana; India) climates to assess metabolic responses to heat in heat-sensitive versus heat-tolerant embryos. We used targeted LC–MS to measure 33 metabolites and the major redox couples (NADH/NAD+, NADPH/NADP+, and GSH/GSSG) at 25°C and after a 32°C heat shock. In all embryos, heat shock induced shared shifts in metabolite profiles, with increases in nucleotide monophosphates (e.g., AMP, CMP, and GMP) and amino acids (e.g., alanine, glutamic acid, serine). In contrast, there were notable differences among genetic lines in the heat-induced responses of redox metabolites. Heat-sensitive temperate lines generally shifted toward a more oxidized state, as indicated by decreases in the NADH/NAD+ ratio, increases in oxidized glutathione, and decreases in the GSH/GSSG ratio, whereas heat-tolerant tropical lines showed little change in redox-couple ratios or, in some cases, shifted toward a more reduced state. These patterns are consistent with greater oxidative stress in the temperate embryos, versus better maintenance of redox balance in tropical embryos. Together, our results suggest that maintaining redox balance is a key element of acute heat tolerance, and healthy development overall, during early embryogenesis. 

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
  
    
