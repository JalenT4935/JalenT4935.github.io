
## Black Gum Tree Species Distribution
## Biogeography course final project site
- Work description
- [adv GIS repository] _____link




```{r}
getwd()
iv_data <- readRDS("IV_data.RDS")
tree <- readRDS("FIA_tree_master1.RDS")
```


```{r}
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(classInt)
library(gridExtra)
library(eHOF)
```



# Introduction:

The black gum tree, scientifically known as the Nyssa Sylvatica, is a deciduous tree in the Nyssaceae family. Black gum is also commonly referred to as black Tupelo, Pepperidge, and sour gum. This species grows mainly in the southeastern United States: from New York and Maine, Michigan, Illinois, southern/eastern Oklahoma to Southern Florida and Eastern Texas. Outside of the United States, this species can also be found in central and southern Mexico. Traditionally, the black gum tree is not predominant in any forests. 

It has been observed that the black gum tree species has been migrating northward. This all ties back to Rapoport’s Rule, which is the observed trend of increasing latitudinal range sizes of species in higher latitudes. It was originally supported by evidence of the trend in trees, marine mollusks, freshwater, 
freshwater/coastal fishes, reptiles/amphibians, mammals, and nonmigratory birds. It is based on the suggestion that high environmental variability at higher latitudes call for organisms with broader environmental tolerances to occur over wider latitudinal ranges. 

The climate variability hypothesis states that since higher latitudes have more variable climates, species at higher latitudes are selected for more broad environmental tolerances, meaning that they can persist across a broader range of climate conditions. 

I will be analyzing different variables as well  (height, diameter, elevation, IV) to see variability. From the northern to souther portion of the black gum’s distribution.


# Data:

In terms of data, I used a combination of a few different things. For one, I used reputable sources with numerical data on the black gum species’ distribution, climate variability, favorable conditions, and characteristics. I also used the Forest Inventory and Analysis website for information about the species’ northward migration. In addition, I used the IV data and FIA tree data to produce results for variables. I looked at results for elevation, IV, height, and diameter. I also used Excel to display some of the results. 

# Methods:

For some data and background information about precipitation, temperature, height, and diameter, I used reputable sources to extract key points. For one, I used the USDA Plant Guide for the black gum species to outline these variables. I extracted the height and diameter range and elevation range. Also, I used sources to capture the temperature and precipitation preferences. 

I used this code for iv data to produce plot of importance value with latitude. The importance value explains the dominance of the species. Higher importance value could indicate the species is more widespread, abundant, and plays a bigger role in the ecosystem. 

```{r}
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
```
![image](https://github.com/user-attachments/assets/15332a94-c0c1-46bc-aec2-68bd5a089731)


```{r}
plot_spcd_IV <- function(species, p_size, low_c, high_c){
  IV_spcd <- iv_data %>%
    filter(COMMON_NAME == species)
  ggplot() +
    geom_point(data = IV_spcd, aes(x = LON, y = LAT, color = IV), size = p_size) +
    scale_color_gradient(low = low_c, high = high_c, name = "iv")
}

plot_spcd_IV("blackgum", 0.9, "blue", "red")
```
![image](https://github.com/user-attachments/assets/6c9fc325-2d11-476d-8d48-ed5349dd910a)


I used these codes to produce data for just black gum species. I analyzed a state from the northern, middle, and southern portion of the black gum’s distribution longitudinal strip to compare. The three states used were Illinois, Tennessee, and Florida, with a range of 2014-2021.

```{r}
blackgumonly1 <- tree %>%
  filter(COMMON_NAME == "blackgum")
head(blackgumonly1)
```

```{r}
blackgum_IL <- blackgumonly1 %>%
  filter(ABBR == "IL")

summary(blackgum_IL)
```

```{r}
blackgum_tn <- blackgumonly1 %>%
  filter(ABBR == "TN")

summary(blackgum_tn)
```

```{r}
blackgum_fl <- blackgumonly1 %>%
  filter(ABBR == "FL")

summary(blackgum_fl)
```


I used Excel to create charts for the maximum height, diameter and elevation for the black gum in Illinois, Tennessee, and Florida. 

# Results: 

```{r echo=FALSE, out.width="100%"}
knitr::include_graphics("height.png",error = FALSE)
```
![image](https://github.com/user-attachments/assets/42f31658-db8a-4938-a2a0-52bef8053fc2)


```{r echo=FALSE, out.width="100%"}
knitr::include_graphics("diameter.png",error = FALSE)
```
![image](https://github.com/user-attachments/assets/c24e8946-f06c-476e-a848-b6f9b7a98554)


```{r echo=FALSE, out.width="100%"}
knitr::include_graphics("elevation.png",error = FALSE)
```
![image](https://github.com/user-attachments/assets/451a1d8b-aa01-4d0b-8ae6-68ef8d745342)

# Discussion:

As seen by the importance value plot, the importance value of the black gum species increased with lattitude. This implies that the abundance and/or role for the black gum species increases with latitude, which can support Rapoport's rule and the idea that increasing latitude has an impact on the latitudinal range of a species.

In terms of the results from the black gum data in Illinois, Tennessee, and Florida, as latitude increases, maximum height decreases. Even with climate change, as you go northward, there are cooler temperatures which may limit tree growth there. In terms of elevation, the spike in elevation in Tenessee is likely due to the mountain coverage. Out of the three states, the maximum tree diameter was highest in Tennessee, which could be an result of the favorable temperatures and consistent rainfall.

In all, future research using more states within the same longitudinal strip would probably be more useful in analyzing the distribution of the black gum species.
