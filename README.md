Analyze binary species occurrence and classification data with a SAR or a cSAR analysis based on a nested or hierarchichal sampling design.

csarGeo allows an advances SAR- and cSAR-analysis to assess biodiversity changes in structurally diverse landscapes. The cSAR analysis function can account for habitat affinity differences across multiple species groups beyond

The package contains an associated vignette with detailed example usages for both of the available analysis methods. For further information regarding the background of countrysideSAR analyses, the papers of ... (pereira, martins)

# Table of Contents

-   [Package Installation](https://github.com/lillyschwietzer/csarGeo/tree/main#1-installation)

-   [Example Analysis](https://github.com/lillyschwietzer/csarGeo/tree/main#2-example-analysis)

-   [References](https://github.com/lillyschwietzer/csarGeo/tree/main#3-references)

# 1. Installation

## 1.1) csarGeo Package

Install package from GitHub:

```{r}
library(pak)

pak("lillyschwietzer/csarGeo")
library(csarGeo)
```

## 1.2) csarGeo Package Data

The package contains three different default data files. One that contains species occurence data and sampling location coordinate information:

```{r}
data("species_data")
head(species_data)
```

One species classification file:

```{r}
data("classes_clusters")
head(classes_clusters)
```

And one SpatRaster land-use file. The latter is a release of the package and may be loaded using a helper function of the csarGeo package called `load_rasterfile()`:

```{r}
library(csarGeo)
library(terra)

# SpatRaster Data
land_use <- load_rasterfile()
plot(land_use)
```

# 2. Example Analysis

The following example only contains details about one of the two possible analysis pathways of csarGeo. For a more detailed explanation, please consult the vignette.

# 3. References
