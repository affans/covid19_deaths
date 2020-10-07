# COVID-19 
### Estimating number of deaths in the United States
Affan Shoukat, 2020  
Center for Infectious Disease Modelling and Analysis  
Yale University, New Haven, USA

### Model details: 
We estimated the reporting rate for COVID-19 related deaths using a Poisson-Logistic Bayesian hierarchical framework applied to each state independently, consisting of a Binomial model for the reported (including confirmed and probable) death counts and a latent Poisson model for the actual number of deaths.
 
Complete details of the model in Supplementary Information of the manuscript.

### Reproducibility
The repository is complete and self-contained. All the files, including data files are included in the repository. 

The model is coded in `R`, `Nimble`, and `Julia` for the whole analysis pipeline. 

#### Main Analysis 

The main analysis file is `main_model_run.R` which runs the Bayesian statistical model. There are two flags. On line 264, use the `arg_st` variable to run for a specific state (or comment it out if running the script on the command line with arguments). The `WRITE_DATA` flag is for writing the results to file (make sure file paths through out the file is specific to your system). 

#### Supplmentary files
`create_plots.R` reads the saved data files and creates a plot for all the states. It also has code to process the saved data files and create a table output. 

`sample_national_level.jl` reads the saved data files to perform analysis at the national level. This code also has plotting at the national level. 

`parallel_bsh_script.sh` was run the script in a Slurm enabled HPC. 

`.csv` are the data files. 

The folder `visualizations` contains R code to create an animation of the results. This code uses state of the art ray-tracing technologies. 

### Manuscript
A link to the manuscript will be provided once published. 

### Contribute and Development
PRs are welcome. Message me to get an explanation of how the model works. If submitting PR, please add relevant tests and make sure all tests pass. 
