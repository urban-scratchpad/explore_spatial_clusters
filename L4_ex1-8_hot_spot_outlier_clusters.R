# Bridge to ArcGIS
library(arcgisbinding)
# ArcGIS License conxn
arc.check_product()

# Draws heavily from Heatherlee Leary's 
# https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut
# Redo and focus more on this explanation: 
# Kam Tin Seong' lesson: 
# https://rpubs.com/tskam/IS415-Hands-on_Ex07
# AND Also:
# sfdep (https://sfdep.josiahparry.com/)


library(sf)
library(tmap)
library(sfdep)
library(spdep)
library(tidyverse) # includes dplyr, ggplot2, 


# read US Counties layer from SNAP data

gis_data <- st_read("C:/EsriTraining/PatternDetection_SpaceTime/patterndetection_spacetime.gdb", layer = 'US_Counties')

head(gis_data)
dim(gis_data)

# what data is there for the Wayne counties?
gis_data %>%
  filter(NAME == "Wayne")

# Plot SNAPRate histogram with density curve

ggplot(data = gis_data, aes(x = SNAPRate)) +
  geom_histogram(aes(y=..density..), color="lightgray") + # 
  geom_density(alpha=.4, size=1, color="orangered2") 
  
# Plot County SNAPRate histogram 
ggplot(data = gis_data, aes(x = SNAPRate)) +
  geom_histogram(color="lightgray") + # 
  labs(title = "Proportion of SNAP beneficiaries by county", 
       x = "SNAP Rate",
       y = "Counties (count)")

# Visualize County SNAPRates on a MAP using ggplot
library(RColorBrewer)
# note: ArcGIS appears to bin values into quintiles or similar to display fill gradient

ggplot(gis_data) +
  geom_sf(aes(fill=SNAPRate),  color = "seashell4") + # , linewidth =.05
  scale_fill_distiller(palette = "YlOrBr", trans = "reverse",
                       name = "SNAP Rate")


# Visualize County SNAPRates on a MAP using tmap library
tm_shape(gis_data) + 
  tm_polygons("SNAPRate")


## using tmap qtm() function does the same as above (quick thematic map)
qtm(gis_data, "SNAPRate")


## Check for Empty Neighbor sets 
## creates a Queen contiguity matrix
# lib spdep: Create a neighbor list based on queen contiguity
###queen = single shared boundary point (line or vertex) meets contiguity condition

list_nb <- poly2nb(gis_data, queen = TRUE)

# review the list_nb matrix: what does it contain?
summary(list_nb)

# To plot matrix, need to redo, see: https://idblr.rbind.io/post/neighborhoods-ggplot/
# centroid_gis_data <- st_centroid(st_geometry(gis_data))
# plot(gis_data, border="lightgrey")
# plot(list_nb, st_centroid(st_geometry(gis_data)), pch = 19, cex = 0.6, add = TRUE, col= "red")

# Check for empty neighbor sets
# card() calculates number of neighbors for each polygon in the list
# which() finds polygons with 0 neighbors

empty_nb <- which(card(list_nb) == 0)
empty_nb  # 16 without neighbors

# Which counties?
# Subset 'gis_data' to extract polygons with empty neighbor sets
empty_polygons <- gis_data[empty_nb, ]
# print rows for counties with no neighbors
empty_polygons  
# There's some big ones!
empty_polygons$NAME

# drop counties that do not have neighbors

gis_data_sub <- gis_data[-empty_nb, ]


## Global G Test for spatial autocorrelation
# Now that we removed empty neighbor sets (gis_data_sub)
# Identify neighbors sets with queen contiguity (edge/vertex touching)
gis_data_nb <- poly2nb(gis_data_sub, queen = TRUE)


# Binary weighting assigns a weight of 1 to all neighboring features 
# and a weight of 0 to all other features
gis_data_w_binary <- nb2listw(gis_data_nb, style="B")

# Calculate spatial lag of SNAPRate
gis_data_lag <- lag.listw(gis_data_w_binary, gis_data_sub$SNAPRate)

# Test for global G statistic of TreEqty
globalG.test(gis_data_sub$SNAPRate, gis_data_w_binary)

# Getis-Ord global G statistic
# 
# data:  gis_data_sub$SNAPRate 
# weights: gis_data_w_binary 
# 
# standard deviate = 28.978, p-value < 2.2e-16
### 28.978 standard deviations away from what would be expected under the 
###    null hypothesis of no clustering
# alternative hypothesis: greater
### greater : The alternative hypothesis is “greater,” which means that the analysis is 
###    looking for clusters of high SNAPRate values
# sample estimates:
#   Global G statistic        Expectation           Variance 
# 2.177649e-03       1.907781e-03       8.672789e-11 
### Global G : clustering stat for the entire dataset
### Expectation: expected value of the clustering statistic under the null
###    hypothesis of no clustering.
### Variance: an estimate of the variance of the clustering statistic under the null hypothesis.
### In sum, the Global G test suggests there is statistically significant clustering of 
### SNAPRate values

## Local Gi Test

### Test for local spatial autocorrelation (hotspots)

# Identify neighbors, create weights, calculate spatial lag
gis_data_nbs <- gis_data_sub |> 
  mutate(
    nb = st_contiguity(Shape),        # neighbors share border/vertex
    wt = st_weights(nb),                 # row-standardized weights
    gis_data_lag = st_lag(SNAPRate, nb, wt)    # calculate spatial lag of TreEqty
  ) 

# Calculate the Gi using local_g_perm
gis_data_hot_spots <- gis_data_nbs |> 
  mutate(
    Gi = local_g_perm(SNAPRate, nb, wt, nsim = 999)
    # nsim = number of Monte Carlo simulations (999 is default)
  ) |> 
  # The new 'Gi' column itself contains a dataframe 
  # We can't work with that, so we need to 'unnest' it
  unnest(Gi) 

# Cursory visualization
# Plot looks at gi values for all locations
gis_data_hot_spots |> 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
#  scale_fill_distiller(palette = "RdBu")
  scale_fill_gradient2(trans = "reverse") # makes the value 0 (random) be the middle


## Using tmap 
tm_shape(gis_data_hot_spots) + tm_polygons("gi", )



#  Add p value classification to 'gis_data_hot_spots"
gis_data_hot_spots |> 
  # with the columns 'gi' and 'p_folded_sim"
  # 'p_folded_sim' is the p-value of a folded permutation test
  select(gi, p_folded_sim) |> 
  mutate(
    # Add a new column called "classification"
    classification = case_when(
      # Classify based on the following criteria:
      gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) |> 
  # Visualize the results with ggplot2
  ggplot(aes(fill = classification)) +
  geom_sf(color = "black", lwd = 0.1) +
  scale_fill_brewer(type = "div", palette = 5) +
  theme_void() +
  labs(
    fill = "Hot Spot Classification",
    title = "SNAP Rate Hot Spots in by County"
  )

# Get data from ArcGIS
## Following: https://nbviewer.org/github/R-ArcGIS/R-Bridge-Tutorial-Notebooks/blob/master/notebooks/01-basics/R-bridge-reading-converting-writing-data.ipynb
## Data types: "For GIS data that is currently stored in a shapefile, 
## file geodatabase, table, feature service, or ArcGIS supported raster data 
## type, you will begin by using the arc.open() function to read it into R."
# arc.open()
