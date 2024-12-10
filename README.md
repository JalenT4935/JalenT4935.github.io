
## Black Gum Tree Species Distribution
## Biogeography course final project site
- Work description
- [adv GIS repository] _____link


Final Project
Jalen Taliaferro 2024-12-10

getwd()
## [1] "/Users/jalentaliaferro/Downloads"
iv_data <- readRDS("IV_data.RDS")
tree <- readRDS("FIA_tree_master1.RDS")
library(ggplot2)
library(dplyr)
## 
## Attaching package: 'dplyr'

## The following objects are masked from 'package:stats':
## 
##     filter, lag

## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(sf)
## Linking to GEOS 3.11.0, GDAL 3.5.3, PROJ 9.1.0; sf_use_s2() is TRUE
library(rnaturalearth)
library(classInt)
library(gridExtra)
## 
## Attaching package: 'gridExtra'

## The following object is masked from 'package:dplyr':
## 
##     combine
library(eHOF)
## Loading required package: mgcv

## Loading required package: nlme

## 
## Attaching package: 'nlme'

## The following object is masked from 'package:dplyr':
## 
##     collapse

## This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.

## Loading required package: lattice

## This is eHOF 1.15 - build: 2024-01-28
Introduction:
The black gum tree, scientifically known as the Nyssa Sylvatica, is a deciduous tree in the Nyssaceae family. Black gum is also commonly referred to as black Tupelo, Pepperidge, and sour gum. This species grows mainly in the southeastern United States: from New York and Maine, Michigan, Illinois, southern/eastern Oklahoma to Southern Florida and Eastern Texas. Outside of the United States, this species can also be found in central and southern Mexico. Traditionally, the black gum tree is not predominant in any forests.

It has been observed that the black gum tree species has been migrating northward. This all ties back to Rapoport’s Rule, which is the observed trend of increasing latitudinal range sizes of species in higher latitudes. It was originally supported by evidence of the trend in trees, marine mollusks, freshwater, freshwater/coastal fishes, reptiles/amphibians, mammals, and nonmigratory birds. It is based on the suggestion that high environmental variability at higher latitudes call for organisms with broader environmental tolerances to occur over wider latitudinal ranges.

The climate variability hypothesis states that since higher latitudes have more variable climates, species at higher latitudes are selected for more broad environmental tolerances, meaning that they can persist across a broader range of climate conditions.

I will be analyzing different variables as well (height, diameter, elevation, IV) to see variability. From the northern to souther portion of the black gum’s distribution.

Data:
In terms of data, I used a combination of a few different things. For one, I used reputable sources with numerical data on the black gum species’ distribution, climate variability, favorable conditions, and characteristics. I also used the Forest Inventory and Analysis website for information about the species’ northward migration. In addition, I used the IV data and FIA tree data to produce results for variables. I looked at results for elevation, IV, height, and diameter. I also used Excel to display some of the results.

Methods:
For some data and background information about precipitation, temperature, height, and diameter, I used reputable sources to extract key points. For one, I used the USDA Plant Guide for the black gum species to outline these variables. I extracted the height and diameter range and elevation range. Also, I used sources to capture the temperature and precipitation preferences.

I used this code for iv data to produce plot of importance value with latitude. The importance value explains the dominance of the species. Higher importance value could indicate the species is more widespread, abundant, and plays a bigger role in the ecosystem.

hof_plot <- function(species){
  black_gum_data <- iv_data %>%
    filter(COMMON_NAME == species)
  lat_seq <- seq(from = floor(min(black_gum_data$LAT)), to = ceiling(max(black_gum_data$LAT)), by = 0.5)
  black_gum_bands <- data.frame(
    Lat_band = lat_seq,
    avg_IV = numeric(length(lat_seq))
  )
  for (i in 1:length(lat_seq)){
    lat_min <- lat_seq[i] - 0.5
    lat_max <- lat_seq[i] + 0.5
    band_data <- black_gum_data %>%
      filter(LAT >= lat_min & LAT <= lat_max)
    black_gum_bands$avg_IV[i] <- mean(band_data$IV, na.rm = TRUE)
  }
  black_gum_bands <- black_gum_bands[complete.cases(black_gum_bands$avg_IV),]
  hof_model <- HOF(
    black_gum_bands$avg_IV,
    black_gum_bands$Lat_band,
    modeltypes = c("I","II","III","IV","V"),
    family = gaussian,
    bootstrap = 100,
    test = 'AIC'
  )
  best_model <- pick.model(hof_model, modeltypes = c("I","II","III","IV","V"), test = 'AIC')
  predicted_response <- predict(hof_model, model = best_model, 
                                newdata = seq(min(black_gum_bands$Lat_band), max(black_gum_bands$Lat_band), by = 0.1))
  scaled_response <- predicted_response * max(black_gum_bands$avg_IV)
  ggplot(black_gum_bands, aes(x = Lat_band, y = avg_IV)) +
    geom_point(color = "grey60") + 
    geom_line(data = data.frame(lat_band = seq(min(black_gum_bands$Lat_band), 
                                               max(black_gum_bands$Lat_band), by = 0.1),
                                predicted_IV = scaled_response),
              aes(x = lat_band, y = predicted_IV), color = "black") +
    geom_vline(xintercept = Para(hof_model, model = best_model)$opt, color = "red", linetype = "dashed") +
    labs(title = paste("Latitudinal IV pattern", species), x = "Latitude", y = "IV (Importance Value")
  
}

hof_plot("blackgum")
##   |                                                     |                                             |   0%  |                                                     |                                             |   1%  |                                                     |.                                            |   2%  |                                                     |.                                            |   3%  |                                                     |..                                           |   4%  |                                                     |..                                           |   5%  |                                                     |...                                          |   6%  |                                                     |...                                          |   7%  |                                                     |....                                         |   8%  |                                                     |....                                         |   9%  |                                                     |....                                         |  10%  |                                                     |.....                                        |  11%  |                                                     |.....                                        |  12%  |                                                     |......                                       |  13%  |                                                     |......                                       |  14%  |                                                     |.......                                      |  15%  |                                                     |.......                                      |  16%  |                                                     |........                                     |  17%  |                                                     |........                                     |  18%  |                                                     |.........                                    |  19%  |                                                     |.........                                    |  20%  |                                                     |.........                                    |  21%  |                                                     |..........                                   |  22%  |                                                     |..........                                   |  23%  |                                                     |...........                                  |  24%  |                                                     |...........                                  |  25%  |                                                     |............                                 |  26%  |                                                     |............                                 |  27%  |                                                     |.............                                |  28%  |                                                     |.............                                |  29%  |                                                     |..............                               |  30%  |                                                     |..............                               |  31%  |                                                     |..............                               |  32%  |                                                     |...............                              |  33%  |                                                     |...............                              |  34%  |                                                     |................                             |  35%  |                                                     |................                             |  36%  |                                                     |.................                            |  37%  |                                                     |.................                            |  38%  |                                                     |..................                           |  39%  |                                                     |..................                           |  40%  |                                                     |..................                           |  41%  |                                                     |...................                          |  42%  |                                                     |...................                          |  43%  |                                                     |....................                         |  44%  |                                                     |....................                         |  45%  |                                                     |.....................                        |  46%  |                                                     |.....................                        |  47%  |                                                     |......................                       |  48%  |                                                     |......................                       |  49%  |                                                     |......................                       |  50%  |                                                     |.......................                      |  51%  |                                                     |.......................                      |  52%  |                                                     |........................                     |  53%  |                                                     |........................                     |  54%  |                                                     |.........................                    |  55%  |                                                     |.........................                    |  56%  |                                                     |..........................                   |  57%  |                                                     |..........................                   |  58%  |                                                     |...........................                  |  59%  |                                                     |...........................                  |  60%  |                                                     |...........................                  |  61%  |                                                     |............................                 |  62%  |                                                     |............................                 |  63%  |                                                     |.............................                |  64%  |                                                     |.............................                |  65%  |                                                     |..............................               |  66%  |                                                     |..............................               |  67%  |                                                     |...............................              |  68%  |                                                     |...............................              |  69%  |                                                     |................................             |  70%  |                                                     |................................             |  71%  |                                                     |................................             |  72%  |                                                     |.................................            |  73%  |                                                     |.................................            |  74%  |                                                     |..................................           |  75%  |                                                     |..................................           |  76%  |                                                     |...................................          |  77%  |                                                     |...................................          |  78%  |                                                     |....................................         |  79%  |                                                     |....................................         |  80%  |                                                     |....................................         |  81%  |                                                     |.....................................        |  82%  |                                                     |.....................................        |  83%  |                                                     |......................................       |  84%  |                                                     |......................................       |  85%  |                                                     |.......................................      |  86%  |                                                     |.......................................      |  87%  |                                                     |........................................     |  88%  |                                                     |........................................     |  89%  |                                                     |........................................     |  90%  |                                                     |.........................................    |  91%  |                                                     |.........................................    |  92%  |                                                     |..........................................   |  93%  |                                                     |..........................................   |  94%  |                                                     |...........................................  |  95%  |                                                     |...........................................  |  96%  |                                                     |............................................ |  97%  |                                                     |............................................ |  98%  |                                                     |.............................................|  99%  |                                                     |.............................................| 100%


plot_spcd_IV <- function(species, p_size, low_c, high_c){
  IV_spcd <- iv_data %>%
    filter(COMMON_NAME == species)
  ggplot() +
    geom_point(data = IV_spcd, aes(x = LON, y = LAT, color = IV), size = p_size) +
    scale_color_gradient(low = low_c, high = high_c, name = "iv")
}

plot_spcd_IV("blackgum", 0.9, "blue", "red")


I used these codes to produce data for just black gum species. I analyzed a state from the northern, middle, and southern portion of the black gum’s distribution longitudinal strip to compare. The three states used were Illinois, Tennessee, and Florida, with a range of 2014-2021.

blackgumonly1 <- tree %>%
  filter(COMMON_NAME == "blackgum")
head(blackgumonly1)
##   GRIDID      LAT       LON       TREECN       PLT_CN INVYR ABBR STATECD SPCD
## 1  11359 30.33415 -93.06852 5.341129e+13 2.595467e+14  2011   LA      22  693
## 2  12188 30.59704 -90.66415 5.341138e+13 2.595467e+14  2011   LA      22  693
## 3  12188 30.59704 -90.66415 5.341139e+13 2.595467e+14  2011   LA      22  693
## 4  12188 30.59704 -90.66415 5.341140e+13 2.595467e+14  2011   LA      22  693
## 5  11357 30.34407 -93.50173 5.341331e+13 2.595467e+14  2011   LA      22  693
## 6  11357 30.34407 -93.50173 5.341331e+13 2.595467e+14  2011   LA      22  693
##   SPGRPCD   DIA HT FGROWCFAL FMORTCFAL FREMVCFAL COMMON_NAME GENUS   SPECIES
## 1      35 12.20 59  0.000000         0  0.000000    blackgum Nyssa sylvatica
## 2      35  2.60 28  0.000000         0  0.000000    blackgum Nyssa sylvatica
## 3      35  6.40 50  0.099840         0  0.000000    blackgum Nyssa sylvatica
## 4      35  1.20 14  0.000000         0  0.000000    blackgum Nyssa sylvatica
## 5      35  5.72 46  0.041582         0  2.512991    blackgum Nyssa sylvatica
## 6      35  8.18 51  0.079025         0  6.138223    blackgum Nyssa sylvatica
##   SITECLCD STDORGCD ELEV SLOPE ASPECT RESERVCD COND_STATUS_CD TREECLCD STATUSCD
## 1        5        0   20    20      7        0              1       NA        2
## 2        3        0   40     0      0        0              1        2        1
## 3        3        0   40     0      0        0              1        3        1
## 4        3        0   40     0      0        0              1        2        1
## 5        4        0   30     0      0        0              1       NA        2
## 6        4        0   30     0      0        0              1       NA        3
##   GRID20 PET_CLIM     AET  FRS  Aridity TempMean Temp_cold Temp_warm    ART
## 1  16086  1391.16 1146.58 3.28 10565.33   196.42    110.74    272.94 284.72
## 2  16098  1438.94 1209.05 3.49 10997.56   193.55    110.17    269.25 289.25
## 3  16098  1438.94 1209.05 3.49 10997.56   193.55    110.17    269.25 289.25
## 4  16098  1438.94 1209.05 3.49 10997.56   193.55    110.17    269.25 289.25
## 5  16084  1384.86 1138.27 3.32 10479.18   195.47    108.53    272.52 287.28
## 6  16084  1384.86 1138.27 3.32 10479.18   195.47    108.53    272.52 287.28
##   MeanMonthlyTemp PrecipMean Precip_dry Precip_wet     TSN   PSN   DEM RDEM
## 1          112.39    1470.14     318.65     153.06 6328.76 13.63  9.69   23
## 2          122.15    1581.69     324.40     162.88 6226.57 14.25 20.12   54
## 3          122.15    1581.69     324.40     162.88 6226.57 14.25 20.12   54
## 4          122.15    1581.69     324.40     162.88 6226.57 14.25 20.12   54
## 5          112.17    1450.06     289.57     147.00 6388.22 15.16 11.11   30
## 6          112.17    1450.06     289.57     147.00 6388.22 15.16 11.11   30
##   Isothermal WaterBody mean_turnover mean_nestedness mean_beta  PLT_lat
## 1      38.98   9866.93     0.3028711      0.20330153 0.5061727 30.31430
## 2      41.75    493.89     0.2578887      0.06100818 0.3188968 30.51774
## 3      41.75    493.89     0.2578887      0.06100818 0.3188968 30.51774
## 4      41.75    493.89     0.2578887      0.06100818 0.3188968 30.51774
## 5      38.50  11614.00     0.4205357      0.06085149 0.4813872 30.27606
## 6      38.50  11614.00     0.4205357      0.06085149 0.4813872 30.27606
##     PLT_lon
## 1 -93.02408
## 2 -90.72998
## 3 -90.72998
## 4 -90.72998
## 5 -93.61094
## 6 -93.61094
blackgum_IL <- blackgumonly1 %>%
  filter(ABBR == "IL")

summary(blackgum_IL)
##      GRIDID           LAT             LON             TREECN         
##  Min.   :28142   Min.   :37.19   Min.   :-89.84   Min.   :3.541e+14  
##  1st Qu.:28960   1st Qu.:37.51   1st Qu.:-89.17   1st Qu.:4.336e+14  
##  Median :28963   Median :37.54   Median :-88.68   Median :5.593e+14  
##  Mean   :29137   Mean   :37.62   Mean   :-88.79   Mean   :6.575e+14  
##  3rd Qu.:29373   3rd Qu.:37.75   3rd Qu.:-88.42   3rd Qu.:9.439e+14  
##  Max.   :31420   Max.   :38.56   Max.   :-87.83   Max.   :1.078e+15  
##                                                                      
##      PLT_CN              INVYR          ABBR              STATECD  
##  Min.   :1.711e+14   Min.   :2015   Length:145         Min.   :17  
##  1st Qu.:3.041e+14   1st Qu.:2016   Class :character   1st Qu.:17  
##  Median :4.506e+14   Median :2018   Mode  :character   Median :17  
##  Mean   :4.408e+14   Mean   :2018                      Mean   :17  
##  3rd Qu.:6.044e+14   3rd Qu.:2020                      3rd Qu.:17  
##  Max.   :7.218e+14   Max.   :2021                      Max.   :17  
##                                                                    
##       SPCD        SPGRPCD        DIA               HT          FGROWCFAL       
##  Min.   :693   Min.   :35   Min.   : 1.100   Min.   :15.00   Min.   :-2.87598  
##  1st Qu.:693   1st Qu.:35   1st Qu.: 5.000   1st Qu.:37.00   1st Qu.: 0.00000  
##  Median :693   Median :35   Median : 6.600   Median :48.50   Median : 0.02961  
##  Mean   :693   Mean   :35   Mean   : 7.496   Mean   :48.76   Mean   : 0.04682  
##  3rd Qu.:693   3rd Qu.:35   3rd Qu.: 9.275   3rd Qu.:58.25   3rd Qu.: 0.18233  
##  Max.   :693   Max.   :35   Max.   :24.500   Max.   :95.00   Max.   : 2.31514  
##                             NA's   :7        NA's   :13      NA's   :4         
##    FMORTCFAL         FREMVCFAL COMMON_NAME           GENUS          
##  Min.   : 0.0000   Min.   :0   Length:145         Length:145        
##  1st Qu.: 0.0000   1st Qu.:0   Class :character   Class :character  
##  Median : 0.0000   Median :0   Mode  :character   Mode  :character  
##  Mean   : 0.4686   Mean   :0                                        
##  3rd Qu.: 0.0000   3rd Qu.:0                                        
##  Max.   :21.4589   Max.   :0                                        
##  NA's   :4         NA's   :4                                        
##    SPECIES             SITECLCD        STDORGCD      ELEV           SLOPE      
##  Length:145         Min.   :2.000   Min.   :0   Min.   :330.0   Min.   : 0.00  
##  Class :character   1st Qu.:4.000   1st Qu.:0   1st Qu.:500.0   1st Qu.:10.00  
##  Mode  :character   Median :5.000   Median :0   Median :530.0   Median :16.00  
##                     Mean   :4.545   Mean   :0   Mean   :543.4   Mean   :17.39  
##                     3rd Qu.:5.000   3rd Qu.:0   3rd Qu.:600.0   3rd Qu.:22.00  
##                     Max.   :6.000   Max.   :0   Max.   :740.0   Max.   :71.00  
##                                                                                
##      ASPECT         RESERVCD       COND_STATUS_CD    TREECLCD    
##  Min.   :  0.0   Min.   :0.00000   Min.   :1      Min.   :2.000  
##  1st Qu.: 61.0   1st Qu.:0.00000   1st Qu.:1      1st Qu.:2.000  
##  Median :202.0   Median :0.00000   Median :1      Median :2.000  
##  Mean   :167.1   Mean   :0.08276   Mean   :1      Mean   :2.175  
##  3rd Qu.:275.0   3rd Qu.:0.00000   3rd Qu.:1      3rd Qu.:2.000  
##  Max.   :360.0   Max.   :1.00000   Max.   :1      Max.   :4.000  
##                                                   NA's   :8      
##     STATUSCD         GRID20         PET_CLIM         AET             FRS       
##  Min.   :1.000   Min.   : 9662   Min.   :1160   Min.   :859.0   Min.   :12.65  
##  1st Qu.:1.000   1st Qu.:10258   1st Qu.:1187   1st Qu.:882.6   1st Qu.:12.65  
##  Median :1.000   Median :10411   Median :1191   Median :909.5   Median :13.04  
##  Mean   :1.103   Mean   :10409   Mean   :1194   Mean   :904.6   Mean   :13.00  
##  3rd Qu.:1.000   3rd Qu.:10562   3rd Qu.:1205   3rd Qu.:911.6   3rd Qu.:13.40  
##  Max.   :2.000   Max.   :10712   Max.   :1214   Max.   :934.2   Max.   :13.63  
##                                                                                
##     Aridity         TempMean       Temp_cold       Temp_warm    
##  Min.   : 9206   Min.   :125.9   Min.   : 1.73   Min.   :238.5  
##  1st Qu.: 9420   1st Qu.:130.8   1st Qu.:10.13   1st Qu.:241.8  
##  Median :10067   Median :134.3   Median :13.25   Median :241.8  
##  Mean   : 9877   Mean   :134.1   Mean   :13.78   Mean   :243.6  
##  3rd Qu.:10092   3rd Qu.:134.6   3rd Qu.:15.51   3rd Qu.:246.5  
##  Max.   :10423   Max.   :141.4   Max.   :22.77   Max.   :253.2  
##                                                                 
##       ART        MeanMonthlyTemp   PrecipMean     Precip_dry      Precip_wet   
##  Min.   :363.2   Min.   :118.1   Min.   :1074   Min.   :219.8   Min.   :114.7  
##  1st Qu.:365.9   1st Qu.:121.1   1st Qu.:1117   1st Qu.:233.9   1st Qu.:120.3  
##  Median :369.0   Median :121.1   Median :1190   Median :244.0   Median :125.2  
##  Mean   :370.4   Mean   :121.4   Mean   :1179   Mean   :243.5   Mean   :123.6  
##  3rd Qu.:375.4   3rd Qu.:121.9   3rd Qu.:1202   3rd Qu.:254.5   3rd Qu.:125.5  
##  Max.   :381.9   Max.   :124.7   Max.   :1264   Max.   :259.6   Max.   :130.5  
##                                                                                
##       TSN            PSN             DEM              RDEM      
##  Min.   :8616   Min.   :14.50   Min.   : 98.36   Min.   : 47.0  
##  1st Qu.:8726   1st Qu.:15.95   1st Qu.:128.11   1st Qu.: 87.0  
##  Median :8888   Median :15.95   Median :146.47   Median :132.0  
##  Mean   :8858   Mean   :16.17   Mean   :144.13   Mean   :115.7  
##  3rd Qu.:9015   3rd Qu.:16.89   3rd Qu.:157.03   3rd Qu.:135.0  
##  Max.   :9187   Max.   :18.74   Max.   :157.03   Max.   :135.0  
##                                                                 
##    Isothermal      WaterBody      mean_turnover     mean_nestedness  
##  Min.   :31.00   Min.   :   0.0   Min.   :0.07826   Min.   :0.05147  
##  1st Qu.:32.00   1st Qu.:   0.0   1st Qu.:0.14954   1st Qu.:0.09083  
##  Median :32.28   Median : 248.1   Median :0.19620   Median :0.10131  
##  Mean   :32.30   Mean   : 382.6   Mean   :0.22026   Mean   :0.13710  
##  3rd Qu.:32.67   3rd Qu.: 717.7   3rd Qu.:0.23258   3rd Qu.:0.16855  
##  Max.   :32.87   Max.   :1122.7   Max.   :0.62500   Max.   :0.35291  
##                                                                      
##    mean_beta         PLT_lat         PLT_lon      
##  Min.   :0.2286   Min.   :37.15   Min.   :-89.76  
##  1st Qu.:0.2655   1st Qu.:37.48   1st Qu.:-89.29  
##  Median :0.3325   Median :37.55   Median :-88.64  
##  Mean   :0.3574   Mean   :37.61   Mean   :-88.78  
##  3rd Qu.:0.4279   3rd Qu.:37.69   3rd Qu.:-88.42  
##  Max.   :0.8511   Max.   :38.49   Max.   :-87.77  
## 
blackgum_tn <- blackgumonly1 %>%
  filter(ABBR == "TN")

summary(blackgum_tn)
##      GRIDID           LAT             LON             TREECN         
##  Min.   :22825   Min.   :34.92   Min.   :-90.07   Min.   :2.168e+14  
##  1st Qu.:24469   1st Qu.:35.42   1st Qu.:-87.46   1st Qu.:4.349e+14  
##  Median :25303   Median :35.86   Median :-85.14   Median :6.341e+14  
##  Mean   :25446   Mean   :35.83   Mean   :-85.46   Mean   :6.882e+14  
##  3rd Qu.:26525   3rd Qu.:36.26   3rd Qu.:-84.07   3rd Qu.:9.655e+14  
##  Max.   :28174   Max.   :36.65   Max.   :-81.67   Max.   :1.262e+15  
##                                                                      
##      PLT_CN              INVYR          ABBR              STATECD  
##  Min.   :1.766e+14   Min.   :2012   Length:2819        Min.   :47  
##  1st Qu.:2.747e+14   1st Qu.:2014   Class :character   1st Qu.:47  
##  Median :4.733e+14   Median :2016   Mode  :character   Median :47  
##  Mean   :4.605e+14   Mean   :2016                      Mean   :47  
##  3rd Qu.:6.517e+14   3rd Qu.:2018                      3rd Qu.:47  
##  Max.   :7.342e+14   Max.   :2019                      Max.   :47  
##                                                                    
##       SPCD        SPGRPCD        DIA               HT           FGROWCFAL      
##  Min.   :693   Min.   :35   Min.   : 1.000   Min.   :  7.00   Min.   :-6.8683  
##  1st Qu.:693   1st Qu.:35   1st Qu.: 2.100   1st Qu.: 22.00   1st Qu.: 0.0000  
##  Median :693   Median :35   Median : 5.600   Median : 40.00   Median : 0.0000  
##  Mean   :693   Mean   :35   Mean   : 5.925   Mean   : 41.46   Mean   : 0.1044  
##  3rd Qu.:693   3rd Qu.:35   3rd Qu.: 7.900   3rd Qu.: 56.00   3rd Qu.: 0.1709  
##  Max.   :693   Max.   :35   Max.   :31.300   Max.   :115.00   Max.   : 7.6783  
##                             NA's   :6        NA's   :6        NA's   :41       
##    FMORTCFAL         FREMVCFAL       COMMON_NAME           GENUS          
##  Min.   : 0.0000   Min.   : 0.0000   Length:2819        Length:2819       
##  1st Qu.: 0.0000   1st Qu.: 0.0000   Class :character   Class :character  
##  Median : 0.0000   Median : 0.0000   Mode  :character   Mode  :character  
##  Mean   : 0.1667   Mean   : 0.1712                                        
##  3rd Qu.: 0.0000   3rd Qu.: 0.0000                                        
##  Max.   :48.4417   Max.   :57.0451                                        
##  NA's   :41        NA's   :41                                             
##    SPECIES             SITECLCD        STDORGCD      ELEV          SLOPE      
##  Length:2819        Min.   :1.000   Min.   :0   Min.   : 270   Min.   : 0.00  
##  Class :character   1st Qu.:5.000   1st Qu.:0   1st Qu.: 740   1st Qu.:12.00  
##  Mode  :character   Median :5.000   Median :0   Median :1190   Median :24.00  
##                     Mean   :4.875   Mean   :0   Mean   :1352   Mean   :27.32  
##                     3rd Qu.:5.000   3rd Qu.:0   3rd Qu.:1810   3rd Qu.:40.00  
##                     Max.   :7.000   Max.   :0   Max.   :4890   Max.   :94.00  
##                                                                               
##      ASPECT         RESERVCD       COND_STATUS_CD    TREECLCD    
##  Min.   :  0.0   Min.   :0.00000   Min.   :1      Min.   :2.000  
##  1st Qu.: 85.0   1st Qu.:0.00000   1st Qu.:1      1st Qu.:2.000  
##  Median :170.0   Median :0.00000   Median :1      Median :2.000  
##  Mean   :167.3   Mean   :0.04576   Mean   :1      Mean   :2.159  
##  3rd Qu.:256.5   3rd Qu.:0.00000   3rd Qu.:1      3rd Qu.:2.000  
##  Max.   :360.0   Max.   :1.00000   Max.   :1      Max.   :4.000  
##                                                   NA's   :271    
##     STATUSCD         GRID20         PET_CLIM         AET        
##  Min.   :0.000   Min.   :11311   Min.   :1065   Min.   : 908.7  
##  1st Qu.:1.000   1st Qu.:11639   1st Qu.:1198   1st Qu.: 968.3  
##  Median :1.000   Median :12063   Median :1245   Median : 991.9  
##  Mean   :1.135   Mean   :12025   Mean   :1243   Mean   : 985.1  
##  3rd Qu.:1.000   3rd Qu.:12373   3rd Qu.:1287   3rd Qu.:1012.4  
##  Max.   :3.000   Max.   :12830   Max.   :1345   Max.   :1040.0  
##                                                                 
##       FRS           Aridity         TempMean       Temp_cold    
##  Min.   : 9.85   Min.   : 9271   Min.   :104.1   Min.   : 7.68  
##  1st Qu.:12.73   1st Qu.:10548   1st Qu.:129.9   1st Qu.:25.88  
##  Median :13.09   Median :10950   Median :137.4   Median :29.98  
##  Mean   :13.44   Mean   :11202   Mean   :135.2   Mean   :30.83  
##  3rd Qu.:14.14   3rd Qu.:11792   3rd Qu.:142.2   3rd Qu.:36.30  
##  Max.   :15.85   Max.   :13907   Max.   :160.4   Max.   :52.07  
##                                                                 
##    Temp_warm          ART        MeanMonthlyTemp   PrecipMean     Precip_dry   
##  Min.   :194.3   Min.   :322.1   Min.   :114.1   Min.   :1130   Min.   :231.1  
##  1st Qu.:225.7   1st Qu.:336.9   1st Qu.:129.2   1st Qu.:1314   1st Qu.:266.7  
##  Median :232.6   Median :342.1   Median :132.4   Median :1401   Median :274.8  
##  Mean   :232.1   Mean   :344.7   Mean   :132.0   Mean   :1386   Mean   :279.8  
##  3rd Qu.:241.6   3rd Qu.:350.8   3rd Qu.:136.0   3rd Qu.:1455   3rd Qu.:292.8  
##  Max.   :263.5   Max.   :364.7   Max.   :143.1   Max.   :1617   Max.   :342.4  
##                                                                                
##    Precip_wet         TSN            PSN             DEM        
##  Min.   :119.8   Min.   :7200   Min.   :10.14   Min.   : 61.18  
##  1st Qu.:137.0   1st Qu.:7620   1st Qu.:13.34   1st Qu.:209.06  
##  Median :143.3   Median :7792   Median :14.62   Median :372.53  
##  Mean   :146.0   Mean   :7784   Mean   :14.84   Mean   :386.90  
##  3rd Qu.:156.5   3rd Qu.:7935   3rd Qu.:15.98   3rd Qu.:485.40  
##  Max.   :169.1   Max.   :8532   Max.   :20.98   Max.   :974.54  
##                                                                 
##       RDEM          Isothermal      WaterBody      mean_turnover    
##  Min.   :  27.0   Min.   :32.13   Min.   :   0.0   Min.   :0.03286  
##  1st Qu.: 155.0   1st Qu.:36.99   1st Qu.:   0.0   1st Qu.:0.16756  
##  Median : 408.0   Median :37.81   Median :   0.0   Median :0.19300  
##  Mean   : 491.8   Mean   :37.80   Mean   : 172.6   Mean   :0.19348  
##  3rd Qu.: 841.0   3rd Qu.:38.98   3rd Qu.: 103.3   3rd Qu.:0.21531  
##  Max.   :1479.0   Max.   :40.03   Max.   :1277.4   Max.   :0.51825  
##                                                                     
##  mean_nestedness     mean_beta         PLT_lat         PLT_lon      
##  Min.   :0.02084   Min.   :0.1676   Min.   :34.99   Min.   :-89.95  
##  1st Qu.:0.04890   1st Qu.:0.2417   1st Qu.:35.41   1st Qu.:-87.48  
##  Median :0.06616   Median :0.2617   Median :35.86   Median :-85.10  
##  Mean   :0.07297   Mean   :0.2664   Mean   :35.83   Mean   :-85.46  
##  3rd Qu.:0.08489   3rd Qu.:0.2880   3rd Qu.:36.24   3rd Qu.:-84.10  
##  Max.   :0.47027   Max.   :0.7002   Max.   :36.66   Max.   :-81.73  
## 
blackgum_fl <- blackgumonly1 %>%
  filter(ABBR == "FL")

summary(blackgum_fl)
##      GRIDID           LAT             LON             TREECN         
##  Min.   : 6504   Min.   :27.43   Min.   :-87.42   Min.   :3.411e+14  
##  1st Qu.:12217   1st Qu.:30.04   1st Qu.:-85.81   1st Qu.:5.385e+14  
##  Median :12638   Median :30.40   Median :-84.80   Median :8.086e+14  
##  Mean   :12700   Mean   :30.36   Mean   :-84.59   Mean   :7.743e+14  
##  3rd Qu.:13436   3rd Qu.:30.73   3rd Qu.:-83.49   3rd Qu.:9.300e+14  
##  Max.   :14273   Max.   :31.05   Max.   :-80.91   Max.   :1.152e+15  
##                                                                      
##      PLT_CN              INVYR          ABBR              STATECD  
##  Min.   :2.182e+14   Min.   :2014   Length:430         Min.   :12  
##  1st Qu.:3.749e+14   1st Qu.:2016   Class :character   1st Qu.:12  
##  Median :4.489e+14   Median :2017   Mode  :character   Median :12  
##  Mean   :5.164e+14   Mean   :2017                      Mean   :12  
##  3rd Qu.:6.383e+14   3rd Qu.:2018                      3rd Qu.:12  
##  Max.   :7.186e+14   Max.   :2019                      Max.   :12  
##                                                                    
##       SPCD        SPGRPCD        DIA               HT           FGROWCFAL      
##  Min.   :693   Min.   :35   Min.   : 1.000   Min.   : 10.00   Min.   :-3.7974  
##  1st Qu.:693   1st Qu.:35   1st Qu.: 5.400   1st Qu.: 37.00   1st Qu.: 0.0000  
##  Median :693   Median :35   Median : 6.900   Median : 49.00   Median : 0.1214  
##  Mean   :693   Mean   :35   Mean   : 7.472   Mean   : 48.64   Mean   : 0.1501  
##  3rd Qu.:693   3rd Qu.:35   3rd Qu.: 9.200   3rd Qu.: 58.75   3rd Qu.: 0.3030  
##  Max.   :693   Max.   :35   Max.   :21.600   Max.   :117.00   Max.   : 3.9285  
##                                                               NA's   :14       
##    FMORTCFAL        FREMVCFAL       COMMON_NAME           GENUS          
##  Min.   : 0.000   Min.   : 0.0000   Length:430         Length:430        
##  1st Qu.: 0.000   1st Qu.: 0.0000   Class :character   Class :character  
##  Median : 0.000   Median : 0.0000   Mode  :character   Mode  :character  
##  Mean   : 0.448   Mean   : 0.0694                                        
##  3rd Qu.: 0.000   3rd Qu.: 0.0000                                        
##  Max.   :29.191   Max.   :11.1832                                        
##  NA's   :14       NA's   :14                                             
##    SPECIES             SITECLCD       STDORGCD      ELEV       
##  Length:430         Min.   :1.00   Min.   :0   Min.   :  0.00  
##  Class :character   1st Qu.:4.00   1st Qu.:0   1st Qu.: 30.00  
##  Mode  :character   Median :5.00   Median :0   Median :100.00  
##                     Mean   :4.86   Mean   :0   Mean   : 95.42  
##                     3rd Qu.:5.00   3rd Qu.:0   3rd Qu.:140.00  
##                     Max.   :7.00   Max.   :0   Max.   :310.00  
##                                                                
##      SLOPE            ASPECT          RESERVCD        COND_STATUS_CD
##  Min.   : 0.000   Min.   :  0.00   Min.   :0.000000   Min.   :1     
##  1st Qu.: 0.000   1st Qu.:  0.00   1st Qu.:0.000000   1st Qu.:1     
##  Median : 0.000   Median :  0.00   Median :0.000000   Median :1     
##  Mean   : 1.505   Mean   : 20.37   Mean   :0.006977   Mean   :1     
##  3rd Qu.: 1.750   3rd Qu.:  0.00   3rd Qu.:0.000000   3rd Qu.:1     
##  Max.   :21.000   Max.   :330.00   Max.   :1.000000   Max.   :1     
##                                                                     
##     TREECLCD       STATUSCD         GRID20         PET_CLIM         AET        
##  Min.   :2.00   Min.   :1.000   Min.   :15964   Min.   :1392   Min.   : 112.1  
##  1st Qu.:2.00   1st Qu.:1.000   1st Qu.:16267   1st Qu.:1455   1st Qu.:1102.5  
##  Median :2.00   Median :1.000   Median :16574   Median :1478   Median :1189.9  
##  Mean   :2.38   Mean   :1.098   Mean   :16535   Mean   :1473   Mean   :1153.4  
##  3rd Qu.:3.00   3rd Qu.:1.000   3rd Qu.:16726   3rd Qu.:1499   3rd Qu.:1225.1  
##  Max.   :3.00   Max.   :3.000   Max.   :18990   Max.   :1567   Max.   :1248.8  
##  NA's   :30                                                                    
##       FRS           Aridity         TempMean       Temp_cold    
##  Min.   :2.390   Min.   : 7906   Min.   :189.0   Min.   :105.7  
##  1st Qu.:4.220   1st Qu.: 8963   1st Qu.:191.5   1st Qu.:108.2  
##  Median :4.670   Median :10142   Median :195.3   Median :115.4  
##  Mean   :4.583   Mean   : 9923   Mean   :195.8   Mean   :116.7  
##  3rd Qu.:4.830   3rd Qu.:10398   3rd Qu.:197.3   3rd Qu.:118.8  
##  Max.   :6.220   Max.   :11411   Max.   :224.7   Max.   :167.9  
##                                                                 
##    Temp_warm          ART        MeanMonthlyTemp   PrecipMean     Precip_dry   
##  Min.   :265.8   Min.   :231.2   Min.   :112.4   Min.   :1227   Min.   :160.9  
##  1st Qu.:266.9   1st Qu.:284.8   1st Qu.:123.1   1st Qu.:1352   1st Qu.:218.1  
##  Median :267.8   Median :286.0   Median :128.9   Median :1499   Median :266.5  
##  Mean   :267.6   Mean   :284.4   Mean   :127.7   Mean   :1458   Mean   :257.3  
##  3rd Qu.:268.0   3rd Qu.:289.2   3rd Qu.:132.5   3rd Qu.:1513   3rd Qu.:287.5  
##  Max.   :273.5   Max.   :296.0   Max.   :136.5   Max.   :1596   Max.   :314.6  
##                                                                                
##    Precip_wet         TSN            PSN             DEM             RDEM      
##  Min.   :164.9   Min.   :4163   Min.   :19.42   Min.   : 9.17   Min.   :16.00  
##  1st Qu.:180.5   1st Qu.:5848   1st Qu.:22.84   1st Qu.:35.10   1st Qu.:54.00  
##  Median :184.4   Median :5970   Median :26.51   Median :40.19   Median :71.00  
##  Mean   :186.7   Mean   :5920   Mean   :28.42   Mean   :40.97   Mean   :63.07  
##  3rd Qu.:193.7   3rd Qu.:6215   3rd Qu.:31.54   3rd Qu.:46.76   3rd Qu.:74.00  
##  Max.   :219.4   Max.   :6302   Max.   :60.75   Max.   :73.25   Max.   :80.00  
##                                                                                
##    Isothermal      WaterBody       mean_turnover    mean_nestedness  
##  Min.   :40.74   Min.   :  371.6   Min.   :0.1375   Min.   :0.03061  
##  1st Qu.:43.26   1st Qu.: 8555.9   1st Qu.:0.2251   1st Qu.:0.08637  
##  Median :44.61   Median :11725.4   Median :0.2664   Median :0.12174  
##  Mean   :44.42   Mean   :10466.5   Mean   :0.2725   Mean   :0.13328  
##  3rd Qu.:45.63   3rd Qu.:13192.8   3rd Qu.:0.3236   3rd Qu.:0.17259  
##  Max.   :51.71   Max.   :19450.6   Max.   :0.5980   Max.   :0.32270  
##                                                                      
##    mean_beta         PLT_lat         PLT_lon      
##  Min.   :0.2530   Min.   :27.49   Min.   :-87.43  
##  1st Qu.:0.3596   1st Qu.:30.10   1st Qu.:-85.78  
##  Median :0.3924   Median :30.39   Median :-84.85  
##  Mean   :0.4058   Mean   :30.36   Mean   :-84.60  
##  3rd Qu.:0.4387   3rd Qu.:30.74   3rd Qu.:-83.55  
##  Max.   :0.6979   Max.   :30.99   Max.   :-80.95  
## 
I used Excel to create charts for the maximum height, diameter and elevation for the black gum in Illinois, Tennessee, and Florida.

Results:
  
Discussion:
As seen by the importance value plot, the importance value of the black gum species increased with lattitude. This implies that the abundance and/or role for the black gum species increases with latitude, which can support Rapoport’s rule and the idea that increasing latitude has an impact on the latitudinal range of a species.

In terms of the results from the black gum data in Illinois, Tennessee, and Florida, as latitude increases, maximum height decreases. Even with climate change, as you go northward, there are cooler temperatures which may limit tree growth there. In terms of elevation, the spike in elevation in Tenessee is likely due to the mountain coverage. Out of the three states, the maximum tree diameter was highest in Tennessee, which could be an result of the favorable temperatures and consistent rainfall.

In all, future research using more states within the same longitudinal strip would probably be more useful in analyzing the distribution of the black gum species.
