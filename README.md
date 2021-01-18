# HiCvisualisation

Some visualisation functions developed for Hi-C data

## HiC-integrationmap

usage: HiC-integrationmap.py [-h] [--int INT] [--bin BIN] [--region REGION] [--chromhmm CHROMHMM] [--out OUT]

generating heatmap with arc plot on top to visualise regional Hi-C interaction with integration with chromHMM states.

optional arguments: 
-h, --help show this help message and exit

--int INT significant interaction file with first three columns as bin1 bin2 and interaction count, the file must conatin a header. 

--bin BIN bed file for bin indexes 

--region REGION chr:start-end 

--chromhmm CHROMHMM path to the chromhmm file for integration 

--out OUT name of the sample to plot

### Output example

<https://github.com/ningbioinfostruggling/HiCvisualisation/blob/main/HiC-integrationmap.pdf>

## Integration-tracks_plot.R

For more information about integration-tracks plot, please see <https://github.com/ningbioinfostruggling/3DFAACTS-SNP>.
