# SeedAdditionSynthesis
Reducing dispersal limitation via seed addition increases species richness but not aboveground biomass

Code used to generate results, and figures for DOI: 
Data DOI: 10.6084/m9.figshare.12319682

This project uses **community and biomass data**,derived from previously published **seed addition experiments**, this README points to the data used, and this repository provides the **R scripts** to investigate and understand models, process data, and generate results and figures.

### Data
Four datasets associated with this project are found on Figshare DOI listed above.

**Plot Level Data** These data summarise species richness and aboveground level biomass responses to a range of diverse seeding treatments across 12 seed addition and diversity experiments. Filename:  _SeedAdd_Plot_Level.csv_

**Beta Diversity** These data summarise the turnover and nestedness components of Beta diversity responses to a range of diverse seeding treatments across 12 seed addition and diversity experiments. Filename:  _betadf.csv_

**Effective Number of Species** These data summarise the effective number of species responses to a range of diverse seeding treatments across 12 seed addition and diversity experiments. Filename:  _seed.pie.csv_

**Species Level Data** For species level data please contact lead author. Filename:  _SeedAdd_Sp_Level.csv_

### **R Scripts** 
Four scripts are provided in this repository.

**Main_Analysis.R** Includes all code to conduct analyses and produce figures found in the main manuscript. This code produces Figure 1-3 and Figures S1a-f.

**Multivariate_Model.R** Includes code to conduct analyses used to understand the correlation between species richness and biomass.

**raw_plot.R** Includes data and code to produce Figure S1.

**div.R** Includes data, model object and analyses used to understand the effect of seeded richness on the effective number of species. This code produces Figure S3.

**Beta.R** Includes data, model objects and analyses used to understand the effect of seeded richness on the turnover and nestedness components of Beta diversity. This code produces Figure S4.

**Citations** Please see details in the Figshare data description.
