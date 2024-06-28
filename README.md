BioIndex
================
Walter Zupa
2024-06-28

# Installation

Install the latest BioIndex release from GitHub:

``` r
library(devtools)
#> Caricamento del pacchetto richiesto: usethis
install_github("COISPA/BioIndex",upgrade = c("never"))
#> Downloading GitHub repo COISPA/BioIndex@HEAD
#>          checking for file 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb\remotes409068524384\COISPA-BioIndex-70617f1/DESCRIPTION' ...     checking for file 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb\remotes409068524384\COISPA-BioIndex-70617f1/DESCRIPTION' ...   v  checking for file 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb\remotes409068524384\COISPA-BioIndex-70617f1/DESCRIPTION'
#>       -  preparing 'BioIndex':
#>    checking DESCRIPTION meta-information ...  v  checking DESCRIPTION meta-information
#>       -  checking for LF line-endings in source and make files and shell scripts
#>   -  checking for empty or unneeded directories
#>   -  building 'BioIndex_0.4.04.tar.gz'
#>      
#> 
#> Warning: il pacchetto 'BioIndex' è in uso e non sarà installato

## uncomment the following code line to install the package buiding the vignettes

# install_github("COISPA/BioIndex",upgrade = c("never"), build_vignettes = TRUE)
```

# Description

BioIndex is an R library to perform analysis of trawl survey data using
MEDITS file format. The functions were previously used in the context of
the Black Sea expert group to analyse data collected in the Black Sea
turbot demersal trawl survey. Then, an updated version (v. 3.1) of the
software was issued and upgraded to elaborate data collected during rapa
whelk scientific survey too. The actual version of the software,
BioIndex v4.04 was developed in R language (version 4.2.1) and presented
to the Black Sea experts. Bioindex library is able to perform analysis
on both MEDITS and MEDITS-like data (Black Sea survey data) after that
the integrity of data tables is checked with RoME package functions. The
software allows to perform the analysis on survey data following the
random stratified sampling, (e.g., MEDITS survey) at GSA level, but in
cases two or more countries are included in the same GSA, the analysis
can be also conducted at country level. The functions included in the
package allows to estimate the time series and the spatial distribution
of a wide pool of population state-indicators for a selected species ().
BioIndex routine also offers the possibilities to perform statistical
trend analysis of both abundance and biomass time series based on the
Hotelling-Pabst test (Cotter, 2009) and to perform a spatial analysis,
placing a list of indices in the spatial dimension, using the 0.5\*0.5°
GFCM grid or bubble plots, to describe, for example, the average
distribution and abundance of the species at local and regional spatial
scale or the localization of sensitive life stage of the population.
BioIndex package was built to work both on stand-alone (on a local
machine) or online, embedded in the RDBFIS environment. A complete
analysis on a selected species can be launched on TA, TB and TC tables
(MEDITS data format) running the BioIndex function which calls the
single functions devoted to specific analyses. Each function can also be
used independently from the others, if a given analysis would be carried
out.

# Use

``` r
BioIndex(ta=TA, tb=TB, tc=TC, sspp="MERLMER",rec_threshold=200,
         spaw_threshold=210,sexes="all", depth=c(10,800), GSA=10, 
         country="all", map_lim=c(13.3,15.2,39.9,41.3),
         depth_lines=c(50,200,800), strata=BioIndex::strata_scheme, 
         stratification_tab = BioIndex::stratification, resolution=1, 
         buffer=0.1, wd=tempdir(), zip=TRUE, save=TRUE, verbose=TRUE)
#> - Merging TA-TB files
#> TA-TB files correctly merged
#> Merge TA-TB files saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/mergeTATB_MERLMER.csv'
#> - Merging TA-TC files
#> TA-TC files correctly merged
#> Merge TA-TC files saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/mergeTATC_MERLMER.csv'
#> 
#> ########################
#> spatial metaDB preparation
#> ########################
#> 
#> Catch metaDB saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/MERLMER - allGSAs_metaDB_catch in GRID.csv 
#> 
#> Biological metaDB saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/MERLMER - allGSAs_metaDB_biological in GRID.csv
#> Analysis conducted for the following country: ITA
#> Bubble plot of Hauls position correctly saved
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

    #> 
    #> ########################
    #> Plot of indices by haul
    #> ########################

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

    #> 
    #> Plot of indices by haul - completed
    #> 
    #> ########################
    #> Time series of indices
    #> ########################

![](README_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

    #> 
    #>  Estimation of abundance indices completed

![](README_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->

    #> 
    #>  Estimation of biomass indices completed 
    #> Time series of indices - completed
    #> 
    #> ########################
    #> Time series of indices
    #> ########################
    #> 
    #> [1] "Select the depth range for the analysis"
    #> Warning: Removed 1 row containing missing values (`geom_line()`).
    #> Warning: Removed 6 rows containing missing values (`geom_point()`).
    #> Warning: Removed 1 row containing missing values (`geom_line()`).
    #> Warning: Removed 6 rows containing missing values (`geom_point()`).
    #> 
    #>  Estimation of MIW completed 
    #> Time series of MIW - completed
    #> 
    #> ########################
    #> Sex-ratio time series
    #> ########################
    #> 
    #> Sex-ratio analysis - completed
    #> 
    #> ############################
    #> Spawners' abundance indices
    #> ############################
    #> 
    #> Spawners' indices analysis - completed
    #> 
    #> ############################
    #> Recruits' abundance indices
    #> ############################

![](README_files/figure-gfm/unnamed-chunk-2-7.png)<!-- -->

    #> 
    #> Recruits' indices analysis - completed
    #> 
    #> ############################
    #> LFD, L0.50 & L0.95
    #> ############################

![](README_files/figure-gfm/unnamed-chunk-2-8.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-10.png)<!-- -->

    #> 
    #> LFD & L0.95 analysis - completed
    #> 
    #> ###########################################
    #> Spearman test of trends on short timeseries
    #> ###########################################
    #> 
    #> 
    #> Spearman test - completed
    #> Abundance indices for statistical squares correctly estimated 
    #> inverse of CV of abundance indices for statistical squares correctly estimated 
    #> file of abundance indices for statistical squares saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/MERLMER - GFCM GRID ABUNDANCE.csv 
    #> Biomass indices for statistical squares correctly estimated 
    #> file of Biomass indices for statistical squares saved in the following folder: 'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/MERLMER - GFCM GRID BIOMASS.csv 
    #> MIW for statistical squares correctly estimated 
    #> inverse of CV of MIW for statistical squares correctly estimated 
    #> file of MIW for statistical squares saved in the following folder:
    #>  'C:\Users\Walter\AppData\Local\Temp\RtmpmuRUDb/output/MERLMER - GFCM GRID MIW.csv

![](README_files/figure-gfm/unnamed-chunk-2-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-2-14.png)<!-- -->

    #> 
    #> #############################
    #> Sex-ratio on GFCM grid
    #> #############################

![](README_files/figure-gfm/unnamed-chunk-2-15.png)<!-- -->

    #> 
    #> Sex-ratio on GFCM grid - completed
    #> 
    #> ################################################
    #> Bubble plots - indices of recruits and spawners
    #> ################################################

![](README_files/figure-gfm/unnamed-chunk-2-16.png)<!-- -->

    #> Bubble plot of recruits correctly saved

![](README_files/figure-gfm/unnamed-chunk-2-17.png)<!-- -->

    #> Bubble plot of spawners correctly saved 
    #> 
    #> Bubble plots - indices of recruits and spawners - completed
    #> 
    #> ###################
    #>  Analysis completed
    #> ###################\n
