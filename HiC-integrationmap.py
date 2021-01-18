# generating heatmap with arc plot on top to visualise regional Hi-C interaction with integration with chromHMM states
import pandas as pd
import pybedtools
import numpy as np
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pylab as plt
from matplotlib.patches import Arc
from matplotlib import gridspec
import argparse

parser = argparse.ArgumentParser(description='generating heatmap with arc plot on top to visualise regional Hi-C interaction with integration with chromHMM states.')
parser.add_argument("--int", help="significant interaction file with first three columns as bin1 bin2 and interaction count, the file must conatin a header.")
parser.add_argument("--bin", help="bed file for bin indexes")
parser.add_argument("--region", help="chr:start-end")
parser.add_argument("--chromhmm", help="path to the chromhmm file for integration")
#parser.add_argument("--name", help="name of the sample to plot")
parser.add_argument("--out",default="out.pdf", help="name of the sample to plot")

args = vars(parser.parse_args())

int_file = args['int']
region = args['region']
c = region.split(':')[0]
s = region.split(':')[0].split('-')[0]
e = region.split(':')[0].split('-')[1]
bin_file = args['bin']
chromhmm_file = args['chromhmm']
out_file = args['out']

def draw_arc(x,y,max_height,ax,colorpalette,colorindex):
    point1 = [x,0]
    point2 = [y,0]
    center = [point1[0]+(point2[0]-point1[0])/2, 0]
    width = (point2[0]-point1[0])
    height = np.sqrt(width)*2

    pac = Arc(center,width,height/2,angle=0,theta1=0,theta2=180,color=colorpalette[colorindex])
    ax.add_patch(pac)

# bins
bins = pd.read_csv(bin_file,
                   sep = '\t',
                   header = None, usecols = [0,1,2,3])

bins = bins[(bins[0] == c) & (bins[1] >= s) & (bins[2] <= e)]
bins.columns = ['chrom','start','end', 'binid']

# interaction
df = pd.read_csv(int_file,
                 sep = '\t', usecols = [0,1,2])

min_bin = min(bins['binid'])
max_bin = max(bins['binid'])

sub_df = df[(df['bin1ID']>=min_bin) & (df['bin1ID']<=max_bin) & (df['bin2ID']>=min_bin) & (df['bin2ID']<=max_bin)]
coordinates = sub_df.values.tolist()

# load annotation
an = pd.read_csv(chromhmm_file,
            sep='\t', skiprows=1, header=None, usecols = [0,1,2,3])

states = list(set(an[3].values.tolist()))
sub_an = an[(an[0] == c) & (an[1] >= s) & (an[2] <= e) & (an[3]!="15_Quies")]

# intersection
bins_bed = pybedtools.BedTool.from_dataframe(bins)
suban_bed = pybedtools.BedTool.from_dataframe(sub_an)

intersect = bins_bed.intersect(suban_bed, wo=True).to_dataframe().iloc[:,[0,1,2,3,7,8]]
intersect.columns = ['chrom','start','end','binid','annotation','overlapbase']
intersect = intersect.drop(['chrom','start','end'],axis=1)

intersect_spread = intersect.groupby(['binid', 'annotation'])['overlapbase'].sum().unstack('annotation')
intersect_spread_t = intersect_spread.transpose()
intersect_spread_t = intersect_spread_t.fillna(0)
for state in states:
    if state not in list(intersect_spread_t.index):
        intersect_spread_t.loc[state] = 0

intersect_spread_t_plot = intersect_spread_t.apply(zscore)
intersect_spread_t_plot = intersect_spread_t_plot.reindex(["1_TssA", "2_TssAFlnk", "10_TssBiv",
                                "3_TxFlnk", "4_Tx", "5_TxWk",
                                "6_EnhG", "7_Enh", "11_BivFlnk", "12_EnhBiv",
                                "8_ZNF/Rpts", "9_Het", "13_ReprPC", "14_ReprPCWk"])

# draw the plot
min_bin = min([min(i[0:2]) for i in coordinates])
max_bin = max([max(i[0:2]) for i in coordinates])
max_width = max([i[1]-i[0] for i in coordinates])
max_height = np.sqrt(max_width)+1

fig = plt.figure(figsize=(20, 20))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

intensity = list(set([i[2] for i in coordinates]))
color_p = sns.color_palette("rocket_r",len(intensity)*2)[len(intensity):]
#color_p.reverse()

ax1.set_ylim(0,max_height/2)
ax1.set_xlim(min_bin-1, max_bin+1)
for i in coordinates:
    draw_arc(int(i[0]),int(i[1]),max_height, ax1, color_p, intensity.index(i[2]))

#ax.axhline(y=0,color ='black', linewidth=1)

cbar_ax = fig.add_axes([.905, .3, .03, .3])
sns.heatmap(intersect_spread_t_plot, ax = ax2, cbar_ax = cbar_ax, cmap = "YlGnBu",
           cbar_kws={'label': 'Z-score'})
cbar_ax.yaxis.label.set_size(20)
cbar_ax.tick_params(labelsize=15)

plt.subplots_adjust(hspace = 0.001)

ax1.axis('off')

ax2.set(xticklabels=[])
ax2.tick_params(bottom=False,labelsize=15)
ax2.set_xlabel(c+':'+str(s)+'-'+str(e),fontsize=20)
ax2.set_ylabel("ChromHMM states",fontsize=20)

plt.savefig(out_file, format='pdf', bbox_inches='tight', dpi=300)
