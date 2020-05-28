"""
createHeatMap.py

Load and format enrichment values with given amino acid and organism sets.
Create heat maps containing enrichment values and optional statistical annotations

python createHeatMap.py --enrichment model_motif1/position1.p --out m1_p1.png
usage: createHeatMap.py [-h] --enrichment ENRICHMENT [-organisms ORGANISMS]
                        [-annotation ANNOTATION] [-thresh THRESH]
                        [-title TITLE] --out OUT

optional arguments:
  -h, --help            show this help message and exit
  --enrichment ENRICHMENT
                        Pickle file to read enrichment values from.
  -organisms ORGANISMS  Ordered organisms to use in plot.
  -annotation ANNOTATION
                        Pickle file containing the statistically significant
                        values
  -thresh THRESH        Threshold value to for the annotations of the plot.
  -title TITLE          Title for plot.
  --out OUT             png file to write plot to.
"""

import argparse
import pickle
import numpy as np
import matplotlib
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable 
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import cm
import matplotlib.pyplot as plt
from fileOutput import parseFileNames

#negative charge, positive charge, aromatic, aliphatic, polar, other
aminoacids = [["D", "E"],["H", "K", "R"],["F", "W", "Y"],["A", "I", "L", "V"],["N", "Q", "S", "T"],["G", "P", "M", "C"]]

"""
getData
picklefile (str)    pickle file containing all frequency and enrichment values
organismfile (str)  provided by user containing ordered organism names

return list of extracted data and the ordered organisms found in pickle file
"""
def getData(picklefile, organismfile):
    data = pickle.load(open(picklefile, "rb"))

    organisms = sorted(list(data["all_f"].keys()))

    if organismfile:
        organism_filenames = parseFileNames(organismfile)    
        new_orgs = [o.split('.')[0].lower() for o in parseFileNames(organismfile)]
        new_orgs = ['_'.join([org.split('_')[0], org.split('_')[1]]) for org in new_orgs]
        organisms = [o for o in new_orgs if o.split('.')[0].lower() in organisms]

    return data, organisms

"""
getEnrichmentArray
data (array)                frequencies and enrichment values extracted from p file 
organisms (list)            organism names extracted from p file and/or provided
aminoacids (list of lists)  amino acid letters separated by property 

Limit the enrichment values displayed to the -0.2 to 2.2 range
return numpy array containing enrichment values to display
"""
def getEnrichmentArray(data, organisms, aminoacids):

    all_enrichment = []
    for amino in aminoacids:
        curr_enrichment = []
        for o in organisms:
            curr_values = []
            for a in amino:
                enrichmentvalue = float(data["pos_e"][o][a][0])
                enrichmentvalue = 2.2 if enrichmentvalue > 2.2 else enrichmentvalue
                enrichmentvalue = -0.2 if enrichmentvalue < -0.2 else enrichmentvalue
                curr_values.append(enrichmentvalue)
            curr_enrichment.append(curr_values)
        curr_enrichment = np.array(curr_enrichment, dtype=np.float)
        all_enrichment.append(curr_enrichment)
    return all_enrichment

"""
getAnnotationArray
annotation_data (array)     array containing statistically significant values
organisms (list)            organism names extracted from p file and/or provided
aminoacids (list of lists)  amino acid letters separated by property
thresh (float)              threshold value to limit annotations to

return numpy array matching enrichment array containing statistically significant values
"""
def getAnnotationArray(annotation_data, organisms, aminoacids, thresh):
    all_annotations = []
    for amino in aminoacids:
        annotation = []
        for o in organisms:
            curr_annotation = []
            for a in amino:
                if a in annotation_data[o].keys():
                    statistic = annotation_data[o][a]
                    if statistic < thresh:
                        curr_annotation.append("*")
                    else:
                        curr_annotation.append("")
                else:
                    curr_annotation.append("")
            annotation.append(curr_annotation)
        annotation = np.array(annotation)
        all_annotations.append(annotation)
    return all_annotations

def heatmap(data, row_labels, col_labels, ax, cbar_kw={}, cbarlabel="", **kwargs):
    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Show x axis tickmarks (amino acid letters)
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(col_labels)

    # Let the horizontal axes labeling appear on the bottom
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)
    return im

def annotateHeatmap(im, data=None, valfmt="{x:.2f}", **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
    
    # Set default alignment to center
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string was supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel"
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color="white")
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

"""
createPlot
enrichment (array)  enrichment values to display 
aminoacids (list)   list of x axis labels
organisms (list)    list of y axis labels
title (str)         plot title
outfile (str)       png file to store plot to
annot (array)       optional, plot annotations

Create heat map with provided x and y axis labels, enrichment values, and optional annotations
"""
def createPlot(enrichment, aminoacids, organisms, title, outfile, annot=False):
    fig, ax = plt.subplots(1, 6, figsize=(8,3))

    # FIRST AMINO ACID GROUP - include y axis labels
    im = heatmap(enrichment[0], organisms, aminoacids[0], ax[0], aspect='equal', vmin=-0.2, vmax=2.2, cmap="seismic")
    ax[0].set_yticks(np.arange(len(organisms))) 
    ax[0].set_yticklabels(organisms)
    if annot:
        annotate_heatmap(im, data=annot[0], valfmt="{x:s}")

    # MIDDLE AMINO ACID GROUPS - only x axis labels
    idx = 1
    for a in aminoacids[1:-1]:
        im = heatmap(enrichment[idx], organisms, aminoacids[idx], ax[idx], aspect='equal', vmin=-0.2, vmax=2.2, cmap="seismic")
        ax[idx].set_yticks([])
        if annot:
            annotate_heatmap(im, data=annot[idx], valfmt="{x:s}")
        idx += 1

    # LAST AMINO ACID GROUP -- add color bar
    im = heatmap(enrichment[-1], organisms, aminoacids[-1], ax[-1], aspect='equal', vmin=-0.2, vmax=2.2, cmap="seismic")
    ax[-1].set_yticks([])
    if annot:
        annotate_heatmap(im, data=annot[-1], valfmt="{x:s}")    
    cbar = ax[-1].figure.colorbar(im, ax=ax[-1], cmap="seismic", ticks=[-0.2, 2.2, 1])

    plt.suptitle(title)
    plt.tight_layout(pad=1.5)
    plt.savefig(outfile)

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--enrichment", required=True, type=str,
        help="Pickle file to read enrichment values from.")
    parser.add_argument("-organisms", type=str,
        help="Ordered organisms to use in plot.")
    parser.add_argument("-annotation", type=str,
        help="Pickle file containing the statistically significant values")
    parser.add_argument("-thresh", type=float, default=0.05,
        help="Threshold value to for the annotations of the plot if the annotation pickle file is provided [default=0.5].")
    parser.add_argument("-title", type=str, default="Enrichment Values",
        help="Title for plot [default=Enrichment Values].")
    parser.add_argument("--out", required=True, type=str,
        help="png file to write plot to.")
    return parser.parse_args()

if __name__ == "__main__":
    # Get all arguments from commandline
    args = parseArguments()
    
    # Load pickle files and store organisms present
    # If an organism file is supplied, compare organisms in pickle file to provided file 
    data, organisms = getData(args.enrichment, args.organisms)
    
    # Create array of enrichment values for heat map
    enrichments = getEnrichmentArray(data, organisms, aminoacids)
    # If annotation is supplied, get array of annotations
    # Create heat map with enrichment values
    if args.annotation:
        annotation = pickle.load(open(args.annotation, "rb"))
        annotations = getAnnotationArray(annotation, organisms, aminoacids, args.thresh)
        createPlot(enrichments, aminoacids, organisms, args.title, args.out, annot=annotations)
    else:
       
        createPlot(enrichments, aminoacids, organisms, args.title, args.out)

