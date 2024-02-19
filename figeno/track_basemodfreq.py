import os
import numpy as np
import pandas as pd
import gzip
import pysam

import matplotlib.pyplot as plt
import matplotlib.patches as patches


#from figeno.genes import draw_genes, draw_genes_transcript
#from figeno.plot_ase import draw_chr_axis
from figeno.utils import correct_region_chr, split_box,draw_bounding_box
from figeno.bam import read_read_groups, decode_read_basemods

class basemodfreq_track:
    def __init__(self,bams=[],modbeds=[],show_legend=False,label="Methylation freq",label_rotate=True,fontscale=1,bounding_box=True,height=25,margin_above=1.5):
        self.bams = bams
        self.modbeds = modbeds

    
        self.label=label
        self.label_rotate=label_rotate
        self.show_legend = show_legend
        self.fontscale=fontscale
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

    def draw(self, regions, box ,hmargin):
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            self.draw_region(regions[i][0],boxes[i])
        self.draw_title(box)

    def draw_region(self,region,box):
        margin_y = (box["top"]-box["bottom"]) * 0.03

        if self.bounding_box:
            draw_bounding_box(box)

        for bam_params in self.bams:
            self.draw_region_bam(region,box,**bam_params)


        # Draw metfreq
        #draw_metfreq(self.samfile,box,margin_y,region,self.colors,self.exchange_haplotypes)

        # Controls
        #for control_bam in self.controls:
        #    control_samfile = pysam.AlignmentFile(control_bam, "rb")
        #    draw_metfreq(control_samfile,box,margin_y,region,["#999999","#999999"],alpha=0.33,linewidth=1,split_gap=True,min_coverage=8)
        
        # Legend
        legend_x=box["left"] + 0.05*(box["right"]-box["left"])
        legend_y=box["bottom"] + 0.3*(box["top"]-box["bottom"])
        legend_height = 0.05*(box["top"]-box["bottom"])
        legend_width = 0.05
        for hp in [0,1]:
            if self.show_legend:
                rect = patches.Rectangle((legend_x,legend_y - hp*legend_height*3.0 -legend_height*0.5),legend_width,legend_height,color=self.colors[hp])
                box["ax"].add_patch(rect)
                box["ax"].text(legend_x+legend_width*1.1,legend_y - hp*legend_height*3.0,"Haplotype "+str(hp+1),horizontalalignment="left",verticalalignment="center")

    def draw_region_bam(self,region,box,file,base,mod,min_coverage=5,linewidth=3,opacity=1,fix_hardclip=False,split_by_haplotype=False,labels=[""],colors=["#27ae60"]):
        min_coverage=int(min_coverage)
        linewidth=float(linewidth)
        opacity=float(opacity)
        margin_y = (box["top"]-box["bottom"]) * 0.03
        samfile = pysam.AlignmentFile(file, "rb")
        region = correct_region_chr(region,samfile.references)

        def transform_coord(pos):
            if region.orientation=="+":
                return box["left"] + (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])
            else:
                return box["right"] - (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])

        # Read methylation
        reads_HP1,reads_HP2,reads_unphased,reads_all = read_read_groups(samfile,region)
        if split_by_haplotype:
            df_met1 = smooth_methylation(create_basemod_table(reads_HP1,base,mod,region.chr,region.start-200,region.end+200,min_coverage=4,samfile=samfile,fix_hardclip=fix_hardclip),w=4,start=region.start,end=region.end)
            df_met2 = smooth_methylation(create_basemod_table(reads_HP2,base,mod,region.chr,region.start-200,region.end+200,min_coverage=4,samfile=samfile,fix_hardclip=fix_hardclip),w=4,start=region.start,end=region.end)
            dfs = [df_met1,df_met2]
        else:
            dfs = [smooth_methylation(create_basemod_table(reads_all,base,mod,region.chr,region.start-200,region.end+200,min_coverage=4,samfile=samfile,fix_hardclip=fix_hardclip),w=4,start=region.start,end=region.end)]

        for j in range(len(dfs)):
            df = dfs[j]
            df = df.loc[df["coverage"]>=min_coverage,:]
            x = [transform_coord(pos) for pos in df["pos"]]
            y = [box["bottom"] + margin_y + val/100 * (box["top"]-box["bottom"] -2 * margin_y) for val in df["smoothed"]]

            boundaries=[0]
            for i in range(len(x)-1):
                if abs(x[i+1]-x[i])> abs(box["right"]-box["left"]) / 20:
                    boundaries.append(i+1)
            boundaries.append(len(x))
            for i in range(len(boundaries)-1):
                a,b = boundaries[i], boundaries[i+1]
                if b-a>1:
                    zorder = 1 if linewidth>1 else 0
                    box["ax"].plot(x[a:b],y[a:b],linewidth=linewidth,color=colors[j],alpha=opacity,zorder=zorder) 



    def draw_title(self,box):
        if (not "projection" in box) or box["projection"]!="polar":
            tick_width=0.5
            margin_y = (box["top"]-box["bottom"]) * 0.03
            box["ax"].plot([box["left"],box["left"]+tick_width],[box["bottom"]+margin_y,box["bottom"]+margin_y],color="black",linewidth=0.8,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,box["bottom"]+margin_y,"0",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            box["ax"].plot([box["left"],box["left"]+tick_width],[(box["bottom"]+box["top"])/2,(box["bottom"]+box["top"])/2],color="black",linewidth=0.8,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,(box["bottom"]+box["top"])/2,"0.5",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            box["ax"].plot([box["left"],box["left"]+tick_width],[box["top"]-margin_y,box["top"]-margin_y],color="black",linewidth=0.8,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,box["top"]-margin_y,"1",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"]-tick_width*12*self.fontscale,(box["bottom"]+box["top"])/2,self.label,rotation=rotation,verticalalignment="center",horizontalalignment="right",fontsize=7*self.fontscale)
        else:
            tick_width=-0.005
            margin_y = (box["top"]-box["bottom"]) * 0.03
            box["ax"].plot([box["left"],box["left"]+tick_width],[box["bottom"]+margin_y,box["bottom"]+margin_y],color="black",linewidth=1,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,box["bottom"]+margin_y,"0",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            box["ax"].plot([box["left"],box["left"]+tick_width],[(box["bottom"]+box["top"])/2,(box["bottom"]+box["top"])/2],color="black",linewidth=1,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,(box["bottom"]+box["top"])/2,"0.5",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            box["ax"].plot([box["left"],box["left"]+tick_width],[box["top"]-margin_y,box["top"]-margin_y],color="black",linewidth=1,zorder=4)
            box["ax"].text(box["left"]-tick_width*1.0,box["top"]-margin_y,"1",horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"]-tick_width*6,(box["bottom"]+box["top"])/2,self.label,rotation=rotation,verticalalignment="center",horizontalalignment="right",fontsize=8*self.fontscale)

def draw_metfreq(samfile,box,margin_y,region,colors,exchange_haplotypes=False,linewidth=3,alpha=1,split_gap=False,min_coverage=5):
    def transform_coord(pos):
        if region.orientation=="+":
            return box["left"] + (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])
        else:
            return box["right"] - (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])

    # Read methylation
    reads_HP1,reads_HP2,reads_unphased,reads_all = read_read_groups(samfile,region)
    if exchange_haplotypes: reads_HP1,reads_HP2 = reads_HP2,reads_HP1
    df_met1 = smooth_methylation(create_methylation_table(reads_HP1,region.chr,region.start-200,region.end+200,min_coverage=min_coverage),w=4,start=region.start,end=region.end)
    df_met2 = smooth_methylation(create_methylation_table(reads_HP2,region.chr,region.start-200,region.end+200,min_coverage=min_coverage),w=4,start=region.start,end=region.end)

    # Plot
    hp=0
    for df in [df_met1,df_met2]:
        x = [transform_coord(pos) for pos in df["pos"]]
        y = [box["bottom"] + margin_y + val/100 * (box["top"]-box["bottom"] -2 * margin_y) for val in df["smoothed"]]
        if split_gap: # In case no CpG for a long while, split the plot. This is to avoid straight lines where there is no data.
            boundaries=[0]
            for i in range(len(x)-1):
                if x[i+1]-x[i]> (box["right"]-box["left"]) / 20:
                    boundaries.append(i+1)
            boundaries.append(len(x))
            for i in range(len(boundaries)-1):
                a,b = boundaries[i], boundaries[i+1]
                if b-a>1:
                    box["ax"].plot(x[a:b],y[a:b],linewidth=1,color=colors[hp],alpha=alpha,zorder=1) 
        else:
            box["ax"].plot(x,y,linewidth=linewidth,color=colors[hp],alpha=alpha,zorder=2)
        hp+=1


def read_metbed(filename,chr,start,end):
    if filename.endswith("gz"):
        f=gzip.open(filename,'rt')
    else:
        f=open(filename,"r")

    d={"chr":[],"pos":[],"metPercentage":[],"coverage":[]}
    for line in f:
        linesplit = line.split("\t")
        if linesplit[0]==chr and int(linesplit[1])>=start and int(linesplit[1])<=end and linesplit[3]=="m":
            details = linesplit[9].split(" ")
            coverage = int(details[0])
            if coverage>5:
                d["chr"].append(chr)
                d["pos"].append(int(linesplit[1]))
                d["metPercentage"].append(float(details[1]))
                d["coverage"].append(coverage)
    df = pd.DataFrame(d)
    f.close()
    return df

def smooth_methylation(df,w=2,start=None,end=None):
    df = df.copy(deep=True)
    smoothed_values = []
    for i in range(w,df.shape[0]-w):
        window_start = i-w
        window_end = i+w
        while abs(df.loc[i,"pos"]-df.loc[window_start,"pos"]) > 100 and window_start<i: window_start+=1
        while abs(df.loc[i,"pos"]-df.loc[window_end,"pos"]) > 100 and window_end>i+1: window_end-=1
        smoothed = np.mean(df.loc[window_start:window_end,"metPercentage"])
        smoothed_values.append(smoothed)
    df = df.loc[w:df.shape[0]-w-1,:]
    df["smoothed"] = smoothed_values
    if start is not None: df = df.loc[df["pos"]>=start,:]
    if end is not None: df = df.loc[df["pos"]<=end,:]
    return df

def create_basemod_table(reads,base,mod,chr,start,end,samfile=None,fix_hardclip=False,min_coverage=5):
    d={"chr":[],"pos":[],"metPercentage":[],"coverage":[]}

    pos2methylated={}
    pos2unmethylated={}
    for x in reads:
        for read in x:
            methyl = decode_read_basemods(read,[(base,mod,1)],samfile,fix_hardclip)
            for pos,end2,state in methyl:
                if base=="C" and (read.flag&16)!=0:
                    pos-=1 # if reverse strand, consider the position of the C in the CpG in the forward orientation.
                if pos>=start and pos<=end:
                    if not pos in pos2methylated:
                        pos2methylated[pos]=0
                        pos2unmethylated[pos]=0
                    if state==1: pos2methylated[pos]+=1
                    else: pos2unmethylated[pos]+=1
    for pos in sorted(pos2methylated.keys()):
        total = pos2unmethylated[pos] + pos2methylated[pos]
        if total>=min_coverage:
            d["chr"].append(chr)
            d["pos"].append(pos)
            d["metPercentage"].append(pos2methylated[pos] / total*100)
            d["coverage"].append(total)
    df = pd.DataFrame(d)
    return df







