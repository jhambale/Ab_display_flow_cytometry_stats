set up two subfolders within the root of the cloned repository:
1. data: this is where input .csv and intermediate .csv files will go.
within each experiment, create a subfolder within 'data' titled 'fl_xxx' where 'fl' is a two character string and 'xxx' is a 3 digit iterator value
2. outputs: this is where .png images of graphed cell binding will be stored
the experiment subfolder 'fl_xxx' will automatically be generated in 'outputs'

when using the Attune Flow Cytometer, draw the relevant gates using the software interface. 
Then, export the experiment statistics (as one .csv file) which will contain the relevant gate information for each well in the experiment. 
This file is what should be placed in the ./data/{expt_id}/stats folder.

within the 'stats'.ipynb notebook, set the experiment ID to match that which contains the 'stats/stats.csv' file of interest within.
