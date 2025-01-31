## Packages ----
# This section ensures that all required packages are installed and loaded.

# Check if the "vegan" package is installed, install if missing, then load it.
if(require(vegan)==FALSE){
  install.packages("vegan")  # Install the package if not already installed
  library(vegan)  # Load the package
}

# Check and install "reshape2" (used for reshaping data).
if(require(reshape2)==FALSE){
  install.packages("reshape2")
  library(reshape2)
}

# Check and install "scales" (used for scaling and formatting plots).
if(require(scales)==FALSE){
  install.packages("scales")
  library(scales)
}

# Check and install "colorspace" (used for advanced color manipulation).
if(require(colorspace)==FALSE){
  install.packages("colorspace")
  library(colorspace)
}

# Check and install "viridis" (provides perceptually uniform color maps).
if(require(viridis)==FALSE){
  install.packages("viridis")
  library(viridis)
}

# Check and install "ggplot2" (a powerful visualization package).
if(require(ggplot2)==FALSE){
  install.packages("ggplot2")
  library(ggplot2)
}

# Check and install "ggfortify" (simplifies statistical visualization with ggplot2).
if(require(ggfortify)==FALSE){
  install.packages("ggfortify")
  library(ggfortify)
}

# Check and install "ggrepel" (helps with non-overlapping text labels in plots).
if(require(ggrepel)==FALSE){
  install.packages("ggrepel")
  library(ggrepel)
}
