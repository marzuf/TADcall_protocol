
# python plotHiC_withTADs_withCols.py

import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate
from scipy.sparse import coo_matrix, triu

def to_matrix(data, chrom1, chrom2, 
              resolution, assembly="hg19",
              start1=None, end1=None, 
              start2=None, end2=None):
    
    if all([x is not None for x in [start1, end1, start2, end2]]):
        start_bin1 = start1//resolution
        end_bin1 = end1//resolution
        start_bin2 = start2//resolution
        end_bin2 = end2//resolution
    else:
        chromsize = load_chromsizes(assembly)
        chromsize1 = chromsize[chrom1]
        chromsize2 = chromsize[chrom2]
        start_bin1 = 0
        end_bin1 = chromsize1//resolution
        start_bin2 = 0
        end_bin2 = chromsize2//resolution
    
    n_bins_1 = end_bin1 - start_bin1 + 1#int(chromsize1 / resolution) + 1
    n_bins_2 = end_bin2 - start_bin2 + 1#int(chromsize2 / resolution) + 1
    
    bin1s = data.bin1//resolution - start_bin1#(data.bin1 / resolution).astype(int).values
    bin2s = data.bin2//resolution - start_bin2#(data.bin2 / resolution).astype(int).values
    values = data.value.values
    

    print("values=")
    print(values.shape)


    print("bin1s=")
    print(bin1s.shape)

    print("n_bins_1")
    print(n_bins_1)


    print("bin2s=")
    print(bin2s.shape)

    print("n_bins_2")
    print(n_bins_2)    
    
    
    m = coo_matrix( (values, (bin1s, bin2s)), shape=(n_bins_1, n_bins_2) )
    if chrom1 == chrom2:
        m = triu(m, k=1).T + m
    return m

#expr_ds = None


matrix_file = "GM12878_chr19_25kb_matrix_list.txt"
tad_file = "select_TADs_to_plot_topdom.txt"

assert os.path.exists(matrix_file)
assert os.path.exists(tad_file)

#map_start = 2140000
#map_end = 2960000
map_start = 1000000
map_end = 5000000

resolution = 25000

#myargs = sys.argv
#if(len(myargs) == 5):
#    hic_ds = myargs[1]
#    expr_ds = myargs[2]
#    select_tad = myargs[3]
#    matrix_file = myargs[4]
#elif(len(myargs) == 4):
#    hic_ds = myargs[1]
#    select_tad = myargs[2]
#    matrix_file = myargs[3]
#else:    
#    hic_ds = "ENCSR489OCU_NCI-H460_40kb"
#    expr_ds = "TCGAlusc_norm_lusc"
#    #matrix_file = "mega_ENCSR489OCU_NCI-H460_mat_chr10_40kb_ob.txt"
#    matrix_file = "ENCSR489OCU-NCI-H460_10kb_chr10.txt"
#    select_tad = "chr10_TAD268"
#    #matrix_file = "mega_ENCSR489OCU_NCI-H460_mat_chr11_40kb_ob.txt"
#    #select_tad = "chr11_TAD390"


out_dir = os.path.join("PLOTHIC_WITHTADS_WITHCOLS")
output_file = os.path.join(out_dir, "plot_selected_tads_chr19_25kb"+str(map_start)+"_"+str(map_end)+".pdf")

chrom="chr19"

# SOME PARAMETERS FOR PLOTTING
other_col_tad = "black"
select_col_tad = "green"
#other_col_lab = "black"
#select_col_lab = "green"
nAround_toplot = 2
nSurroundBins = 2
shrink = 1
tad_lwd = 1.2
lab_offset = -0.8
labSize = 10    
labBox = True
addGenes = True
geneSymbPos = "right" 
genesOffset = 0.5
symbolOffset = 0.1
withBox = True
        
if nAround_toplot == 0:    
    geneSymbPos = "above" 
    geneName_size = 8
    geneBar_height = 0.2
    genesSpacing = 0.5

elif nAround_toplot == 2:    
    geneName_size = 8
    geneBar_height = 1
    genesSpacing = 1.2

else :    
    geneName_size = 8
    geneBar_height = 0.2
    
gene_col = "darkblue"
gene_lwd = 1
plot_labs = True
plotTit = "GM12878_chr19_25kb 2140001-2960000 "
#if expr_ds:
#    plotTit += " - " + expr_ds

#if addGenes:
#    assert expr_ds

################################################################################################################################# PARAMETERS TO SET HERE  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
setDir = os.path.join("/media", "electron")
setDir = os.path.join("/")



#################################################################################################################################



if not os.path.exists(out_dir):
    os.makedirs(out_dir)


    
tad_dt = pd.read_csv(tad_file, sep="\t", header=None, names = ['chromo', 'start', 'end', "region"] )

#select_idxs=  tad_dt[tad_dt['region'] == select_tad].index
#select_idxs=  list(range(tad_dt.shape[0]))
##assert len(select_idxs) == 1
#select_idx =  int(tad_dt[tad_dt['region'] == select_tad].index[0]) 

#chrom = tad_dt['chromo'][select_idx]
#assert re.match(chrom, select_tad)

#start_pos_list = list(tad_dt['start'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)])
#start_pos_list = list(tad_dt['start'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)] - 1)
#end_pos_list = list(tad_dt['end'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)])
#end_pos_list = list(tad_dt['end'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)] + 1)

start_pos_list = list(tad_dt['start'] - 1)
end_pos_list = list(tad_dt['end'])


#lab_list = ([""] * nAround_toplot) + [select_tad] + ([""] * nAround_toplot)
lab_list = (list(tad_dt['region']) )
col_list_tad = ([select_col_tad] *len(lab_list))
#col_list_tad = ([other_col_tad] * nAround_toplot) + [select_col_tad] + ([other_col_tad] * nAround_toplot)
#col_list_lab = ([other_col_lab] * nAround_toplot) + [select_col_lab] + ([other_col_lab] * nAround_toplot)


#map_start = start_pos_list[0] - nSurroundBins*resolution
#map_end = end_pos_list[len(end_pos_list)-1] + nSurroundBins*resolution



assert len(start_pos_list) == len(end_pos_list) == len(lab_list) == len(col_list_tad)

tad_toplot = [
    start_pos_list,
    end_pos_list,
    lab_list,
    col_list_tad
]

print(tad_toplot)
    
# LOAD AND SUBSET Hi-C DATA


data = pd.read_csv(matrix_file, sep="\t", header=None, names=['chromo', 'bin1', 'bin2', 'value'])

del data['chromo']

print(data.head())

data['bin1'] *= resolution
data['bin2'] *= resolution

print(data.head())

assert all(list(data['bin1'] <= data['bin2']))

data = data[(data.bin1>=map_start) & (data.bin2<=map_end)]
m = to_matrix(data, chrom, chrom, resolution, start1=map_start, end1=map_end, start2=map_start, end2=map_end)

# retrieve the height of the plot (half-diagonal)
diagonal_height = np.sqrt(m.shape[0]**2 / 2)

my_map = plt.get_cmap("Reds")
my_map.set_under('white')
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
X = np.log10(np.triu(m.todense()*shrink) + 1)
ax.matshow( rotate(X, 45), cmap=my_map, vmin=0.01 )

plt.ylim(diagonal_height)

plt.xticks([])
plt.yticks([])
plt.axis('on')

if not withBox:
    plt.box(on=None)

tad_start=tad_toplot[0][0]
tad_end= tad_toplot[1][0]
tad_lab=tad_toplot[2][0]
tad_col=tad_toplot[3][0]
        
for tad_start, tad_end , tad_lab, tad_col in zip(tad_toplot[0], tad_toplot[1], tad_toplot[2], tad_toplot[3]):

    tad_start_bin = (tad_start - map_start)//resolution  
    tad_end_bin = (tad_end - map_start)//resolution

    conv_start = np.sqrt(2) * tad_start_bin  # distance to bin start: on the side of the square, it is the # of bins but in triangle plot, it is the diagonal
    conv_end = np.sqrt(2) * tad_end_bin

    tad_height = diagonal_height-(conv_end-conv_start)*0.5  # needs the "diagonal_heigh" because the way it is plotted and the orientation of the y-axis
    
    tad_triangle = np.array([
        [conv_start,diagonal_height],
        [0.5*(conv_start+conv_end), tad_height],
        [conv_end, diagonal_height]
    ])

    print(tad_triangle)

    tadTri = plt.Polygon(tad_triangle, facecolor='none', edgecolor=tad_col, linewidth=tad_lwd)
    ax.add_patch(tadTri)

    if labBox:
        plt.text(0.5*(conv_start+conv_end),tad_height + lab_offset,tad_lab, fontsize=labSize, horizontalalignment='center', verticalalignment='bottom',
             bbox=dict(facecolor=tad_col, alpha=0.5), fontweight='bold')
    else:
        plt.text(0.5*(conv_start+conv_end),tad_height + lab_offset,tad_lab, fontsize=labSize, horizontalalignment='center', verticalalignment='bottom', fontweight='bold')
    
    


#     first_start = start_pos_list[0]
#     first_end =  0.5*( start_pos_list[0] + end_pos_list[len(end_pos_list)-1])
#     first_lab = "FOO_FIRST_GENE"
#     conv_gene_start = np.sqrt(2) * (first_start - map_start)/resolution
#     conv_gene_end = np.sqrt(2) * (first_end - map_start)/resolution
#     ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=0.2,
#                            facecolor=gene_col, clip_on=False))
#     plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset, first_lab, fontsize=8, horizontalalignment='center', verticalalignment='center')
#     genePos += genesSpacing

#for i in range(tad_genes_symbols_dt.shape[0]):
#        conv_gene_start = np.sqrt(2) * (tad_genes_symbols_dt['start'][i] - map_start)/resolution
#        conv_gene_end = np.sqrt(2) * (tad_genes_symbols_dt['end'][i] - map_start)/resolution
#        
#        ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=geneBar_height, 
#                               facecolor=gene_col, clip_on=False))
#        
##         plt.hlines(y=genePos, xmin=conv_gene_start, xmax=conv_gene_end)
##         plt.hlines(y=genePos, xmin=map_start, xmax=map_end)
#        
#        ax.axhline(y=genePos, xmin=map_start, xmax=map_end, c = 'red')
#        plt.axhline(y=genePos)
#        plt.axhline(y=genePos+geneBar_height)
#        
#        symb_pos = 0.5*(genePos+geneBar_height + genePos)
#        
#        if geneSymbPos == "above":
#            plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset,tad_genes_symbols_dt['symbol'][i], fontsize=geneName_size, horizontalalignment='center', verticalalignment='center')
#        elif geneSymbPos == "right":
#            plt.text(conv_gene_end, symb_pos,tad_genes_symbols_dt['symbol'][i], fontsize=geneName_size, horizontalalignment='left', verticalalignment='center')
#        elif geneSymbPos == "left":
#            plt.text(conv_gene_start, symb_pos,tad_genes_symbols_dt['symbol'][i], fontsize=geneName_size, horizontalalignment='right', verticalalignment='center')

#        genePos += genesSpacing
        
        
#     last_start = 0.5*( start_pos_list[0] + end_pos_list[len(end_pos_list)-1])
#     last_end =  end_pos_list[len(end_pos_list)-1]
#     conv_gene_start = np.sqrt(2) * (last_start - map_start)/resolution
#     conv_gene_end = np.sqrt(2) * (last_end - map_start)/resolution
#     print("LAST");print(conv_gene_start);print(genePos);print(conv_gene_end-conv_gene_start)
#     ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=0.2,
#                            facecolor=gene_col, clip_on=False))
#     plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset, last_lab, fontsize=8, horizontalalignment='center', verticalalignment='center')
#     genePos += genesSpacing
    
        
plt.title(plotTit, fontweight="bold")    

if output_file:
	plt.savefig(output_file, bbox_inches='tight',transparent=True)#, pad_inches=0)
    
	print("... saved: " + output_file + "\n")
    
