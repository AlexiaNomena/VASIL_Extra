Click to run:  
[Pythonanywhere platform: Compute the Cross-Reactivity between two lineages](https://projects-raharinirina.pythonanywhere.com/vasil/FoldR_PNeut/)

(OBSOLETE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AlexiaNomena/VASIL_Extra/HEAD?urlpath=%2Fvoila%2Frender%2FCross_Demo.ipynb)

[Binder platform: Compute the Cross-Reactivity between two lineages](https://mybinder.org/v2/gh/AlexiaNomena/VASIL_Extra/HEAD?urlpath=%2Fvoila%2Frender%2FCross_Demo.ipynb)

## Run Simulations With Vaccination Timeline Included

### Step 1: Adapt `scripts/Vaccination/Process_Vaccination_Timeline.py`to your problem

Check if `scripts/Vaccination/Process_Vaccination_Timeline.py` contains elements in the vaccination timeline dataset and edit if needed.
The name of the variant representing the vaccine is added to the name of each vaccine considered (from the front): The mutation profile of this variant is used to compute a neutralization probability that corresponds to the vaccine, and this variant must be present in the genomic surveillance dataset, otherwise there will be a code error


### Step 2: 
