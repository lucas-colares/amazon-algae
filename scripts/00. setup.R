## Packages ----

if(require(vegan)==FALSE){
  install.packages("vegan")
  library(vegan)
}

if(require(reshape2)==FALSE){
  install.packages("reshape2")
  library(reshape2)
}

if(require(scales)==FALSE){
  install.packages("scales")
  library(scales)
}

if(require(colorspace)==FALSE){
  install.packages("colorspace")
  library(colorspace)
}

if(require(viridis)==FALSE){
  install.packages("viridis")
  library(viridis)
}

if(require(ggplot2)==FALSE){
  install.packages("ggplot2")
  library(ggplot2)
}

if(require(ggfortify)==FALSE){
  install.packages("ggfortify")
  library(ggfortify)
}

if(require(ggrepel)==FALSE){
  install.packages("ggrepel")
  library(ggrepel)
}
