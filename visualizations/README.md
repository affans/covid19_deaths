These scripts build figure 2 in the manuscript and also the spatiotemporal movie. 


# Data preprocess

Within the `input_data` directory there exists the reported deaths data `Death.csv` and also the estimated death data we calculate `EDeath.csv`. Rows are dates and columns are states. See `R/deaths_data_setup.R` for more information. 

0. The working directory should be this directory (`visualizations/`).
1. Run `R/deaths_data_setup.R` to combine the input data with spatial data. 

# Figure 2 

2. Run  `figure_2_deaths_spatial_map_create.Rmd` to build figure 2 in the manuscript. 
This RMarkdown document was knit to HTML, which can be viewed [here](https://htmlpreview.github.io/?https://github.com/affans/covid19_deaths/blob/master/visualizations/figure_2_deaths_spatial_map_create.html).  

# Spatiotemporal visualization

All frames stored within `figures/frames`, then combined. 

### Run 3D visualization

3. Run `R/deaths_rayrender_render.R` using the `Rscript` command. 
E.g., at the command line run `Rscript R/deaths_rayrender_render.R [date number]` where `[date number]` is a number 
from 1 to the number of dates in the analysis. I ran this on a cluster, using the script in the `scripts_for_cluster` directory to submit all jobs at once. 

### Run 2D visualizations

4. Run `R/deaths_cumulative_ggplot_line.R` and `R/deaths_cumulative_ggplot_spatial.R`

### Combine images

5. Run `R/deaths_combine_images.RData` to compile all images into frames and build the movie. 
