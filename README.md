Click to run:  
[Pythonanywhere platform: Compute the Cross-Reactivity between two lineages](https://projects-raharinirina.pythonanywhere.com/vasil/FoldR_PNeut/)

(OBSOLETE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AlexiaNomena/VASIL_Extra/HEAD?urlpath=%2Fvoila%2Frender%2FCross_Demo.ipynb)

[Binder platform: Compute the Cross-Reactivity between two lineages](https://mybinder.org/v2/gh/AlexiaNomena/VASIL_Extra/HEAD?urlpath=%2Fvoila%2Frender%2FCross_Demo.ipynb)

## Dependencies

See [VASIL](https://github.com/KleistLab/VASIL/tree/main) repository 

## Run Simulations With Vaccination Timeline Included

### Step 0: 
First, run the main pipeline VASIL as instructed [here](https://github.com/KleistLab/VASIL/tree/main) if the following necessary files are not present in your working directory (e.g., `working_directory = data/ByCountry/Germany`)
- "working_directory/Daily_SpikeGroups_Freq.csv"
- "working_directory/results/Cross_react_dic_spikegroups_ALL.pck"
- "working_directory/Spikegroups_membership.pck"

### Step 1: Adapt `scripts/Vaccination/Process_Vaccination_Timeline.py` to your problem

- Add vaccination timeline dataset in the working directory ... For example, `data/ByCountry/Germany/germany_vaccinations_timeseries_v3.csv`
  
Check if `scripts/Vaccination/Process_Vaccination_Timeline.py` contains elements in the vaccination timeline dataset and edit if needed.
The name of the variant representing the vaccine is added to the name of each vaccine considered (from the front): The mutation profile of this variant is used to compute a neutralization probability that corresponds to the vaccine, and this variant must be listed in the file `working_directory/Daily_SpikeGroups_Freq.csv` dataset, otherwise there will be a code error


### Step 2: Run snakemake pipeline (run twice to complete)

- Configure parameters `working_directory/vacc_config.yaml` ... For example, `working_directory = data/ByCountry/Germany`
- Run the snakemake pipeline from the terminal
  ```
  snakemake --snakefile VACC --configfile working_directory/vacc_config.yaml -j -d working_directory
  ```
- The results will be found in a newly created or updated directory `working_directory/vaccination`. The RELEVANT results are
  * Directory `working_directory/vaccination/Immunological_Landscape_ALL` contains the immune landscape without vaccination, and the weighted mean is `working_directory/vaccination/Susceptible_weighted_mean_over_spikegroups_all_PK.csv`
  * Directory `working_directory/vaccination/ImL_ALL_cs_Vacc_ver2` contains the immune landscape with vaccination, and the weighted mean is `working_directory/vaccination/Susceptible_weighted_mean_over_spikegroups_vs_Vacc_ver2_all_PK.csv`
  * Directory `working_directory/vaccination/Timeline` contains vaccination timeline (total, and per vaccine)
  * Directory `working_directory/vaccination/Cross_Vacc` contains cross-reactivity between variants in the timeline and vaccines
  * Directory `working_directory/vaccination/plots/ALL_vs_Vacc_ver2` contains plots of the results with vaccination
  * Directory `working_directory/vaccination/plots/relative_groups` contains plots of the relative fitness for simulations including vaccination (as requested in `working_directory/vacc_config.yaml`)
 
    
