# Amazon Algae  

![alt text](https://i.imgur.com/i68N8V5.jpeg)

## Overview  

This repository contains the R project associated with the study **"On the taxonomic richness, evenness, and divergence of periphytic algae in Amazon streams"** by Palheta et al. (2024). The research investigates how landscape alterations impact the limnological and biological conditions of Amazonian streams, affecting the diversity of periphytic algae communities.  

## Biological Context  

Periphytic algae are key components of freshwater ecosystems, contributing to primary production and serving as bioindicators of water quality. However, human-induced landscape modifications—such as deforestation and land use changes—disrupt stream integrity, altering nutrient concentrations and affecting algal diversity.  

The study tested two main hypotheses:  
1. **Nutrient concentrations are higher in streams with low physical integrity**, likely due to increased soil erosion and runoff from deforested areas.  
2. **High nutrient concentrations lead to increased taxonomic richness but lower evenness and divergence**, favoring opportunistic species.  

### Key Findings  
- **Nutrient overload in degraded streams**: Streams with poor habitat integrity had elevated levels of **total phosphorus, ammonia, and aluminum**, which influenced periphytic community structure.  
- **Changes in algal diversity patterns**: Species richness increased in disturbed streams due to nutrient influx, but its effects on evenness were unclear.  
- **The role of riparian vegetation**: Streams with intact riparian zones maintained more stable conditions, while deforested streams exhibited higher algal richness due to increased nutrient availability.  

These results emphasize the importance of monitoring multiple facets of algal diversity to assess the ecological health of Amazonian streams.  

## Folder Structure  

### `datasets/`  

Contains the datasets used in the study:  

- **`FinalPoints.csv`** – Metadata and coordinates of sampling points.  
- **`HII.csv`** – Habitat Integrity Index (HII), quantifying the physical condition of stream environments.  
- **`comm.csv`** – Community composition of periphytic algae at different sampling sites.  
- **`coordinates.csv`** – Geographic coordinates of study sites.  
- **`env.csv`** – Environmental and limnological data, including nutrient concentrations and physical parameters.  

### `scripts/`  

Contains the R scripts for data analysis:  

- **`00. setup.R`** – Loads necessary libraries and functions.  
- **`01. scale of effect.R`** – Assesses the impact of environmental variables on algal diversity metrics.  
- **`02. statistical analysis.R`** – Performs statistical analyses using generalized linear models (GLMs) and beta regression.  
- **`02. statistical analysis_corr.R`** – Alternative statistical analysis focusing on correlation between variables.  

## Getting Started  

1. Clone the repository:  
   ```bash  
   git clone https://github.com/yourusername/amazon-algae.git  
   ```  
2. Open `AlgaeProj #1.Rproj` in RStudio.  
3. Run `00. setup.R` to set up the analysis environment.  
4. Use the provided scripts to explore and analyze the data.  

## Citation  

If you use this repository, please cite:  

> Palheta, L., Colares, L. F., Junqueira, M. G., & Dunck, B. (2024). On the taxonomic richness, evenness, and divergence of periphytic algae in Amazon streams. *Limnology*. https://doi.org/10.1007/s10201-024-00767-4  

## License  

This project is released under the MIT License. 
