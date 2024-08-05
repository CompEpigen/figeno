import os
import numpy as np
import pandas as pd
import gzip
import pysam

import matplotlib.pyplot as plt
import matplotlib.patches as patches


#from figeno.genes import draw_genes, draw_genes_transcript
#from figeno.plot_ase import draw_chr_axis
from figeno.utils import KnownException, correct_region_chr, split_box,draw_bounding_box, chr_to_int
from figeno.bam import read_read_groups, decode_read_basemods

class basemodfreq_track:
    def __init__(self,style="lines",smooth=4,gap_frac=0.1,bams=[],bedmethyls=[],show_legend=False,label="Methylation freq",label_rotate=True,fontscale=1,bounding_box=True,height=25,margin_above=1.5,**kwargs):
        self.style=style
        self.smooth=int(smooth) # if set to a value x>0, will have average at each position the methylation values at this position, and at the next and previous x positions. If set to 0, does not perform any smoothing.
        if self.style=="dots": self.smooth=0
        self.gap_frac=float(gap_frac)
        self.bams = bams
        self.bedmethyls = bedmethyls
        if len(bams)==0 and len(bedmethyls)==0: raise KnownException("Please include at least one file in the basemod_freq track.")

        self.label=label
        self.label_rotate=label_rotate
        self.show_legend = show_legend
        self.fontscale=float(fontscale)
        self.bounding_box=bounding_box
        self.height = float(height)
        self.margin_above= float(margin_above)

        self.is_empty=True
        self.kwargs=kwargs

    def draw(self, regions, box ,hmargin,warnings=[]):
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            self.draw_region(regions[i][0],boxes[i],warnings=warnings)
        self.draw_title(box)
        if self.is_empty:
            warnings.append("No data was found in the displayed regions for the basemod_freq track.")

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the basemod_freq track because it is not one of the accepted parameters.")

    def draw_region(self,region,box,warnings=[]):
        margin_y = (box["top"]-box["bottom"]) * 0.03

        if self.bounding_box:
            draw_bounding_box(box)

        for bam_params in self.bams:
            self.draw_region_bam(region,box,**bam_params,warnings=warnings)

        for bedmethyl_params in self.bedmethyls:
            self.draw_region_bedmethyl(region,box,**bedmethyl_params)

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

    def draw_region_bam(self,region,box,file,base,mod,min_coverage=5,linewidth=3,opacity=1,fix_hardclip=False,split_by_haplotype=False,labels=[""],colors=["#27ae60"],warnings=[]):
        min_coverage=int(min_coverage)
        linewidth=float(linewidth)
        opacity=float(opacity)
        margin_y = (box["top"]-box["bottom"]) * 0.03

        if file=="" or file is None:
            raise KnownException("No bam file was provided in the basemod_freq track.")
        if not os.path.isfile(file):
            raise KnownException("The following bam file does not exist (in basemod_freq track): "+file)
        try:
            samfile =pysam.AlignmentFile(file, "rb")
        except: 
            raise KnownException("Failed to open bam file: "+str(file))
        if not samfile.has_index():
            raise KnownException("Missing index file for bam file: "+str(file)+". Such an index (ending in .bai) can be generated with samtools index.")

        samfile = pysam.AlignmentFile(file, "rb")
        region = correct_region_chr(region,samfile.references)

        def transform_coord(pos):
            if region.orientation=="+":
                return min(box["right"]-0.5,max(box["left"]+0.5,box["left"] + (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])))
            else:
                return min(box["right"]-0.5,max(box["left"]+0.5,box["right"] - (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])))

        # Read methylation
        reads_HP1,reads_HP2,reads_unphased,reads_all = read_read_groups(samfile,region)
        if split_by_haplotype:
            df_met1 = smooth_methylation(create_basemod_table_bam(reads_HP1,base,mod,region.chr,region.start-200,region.end+200,min_coverage=min(4,min_coverage),samfile=samfile,fix_hardclip=fix_hardclip,warnings=warnings),w=self.smooth,start=region.start,end=region.end)
            df_met2 = smooth_methylation(create_basemod_table_bam(reads_HP2,base,mod,region.chr,region.start-200,region.end+200,min_coverage=min(4,min_coverage),samfile=samfile,fix_hardclip=fix_hardclip,warnings=warnings),w=self.smooth,start=region.start,end=region.end)
            dfs = [df_met1,df_met2]
        else:
            dfs = [smooth_methylation(create_basemod_table_bam(reads_all,base,mod,region.chr,region.start-200,region.end+200,min_coverage=min(4,min_coverage),samfile=samfile,fix_hardclip=fix_hardclip,warnings=warnings),w=self.smooth,start=region.start,end=region.end)]

        for j in range(len(dfs)):
            df = dfs[j]
            df = df.loc[df["coverage"]>=min_coverage,:]
            x = [transform_coord(pos) for pos in df["pos"]]
            y = [box["bottom"] + margin_y + val/100 * (box["top"]-box["bottom"] -2 * margin_y) for val in df["smoothed"]]
            if df.shape[0]>0: self.is_empty=False
            boundaries=[0]
            for i in range(len(x)-1):
                if abs(x[i+1]-x[i])> abs(box["right"]-box["left"]) / 10:
                    boundaries.append(i+1)
            boundaries.append(len(x))

            if self.style=="lines":
                for i in range(len(boundaries)-1):
                    a,b = boundaries[i], boundaries[i+1]
                    if b-a>1:
                        box["ax"].plot(x[a:b],y[a:b],linewidth=linewidth,color=colors[j],alpha=opacity,zorder=1) 
            else:
                box["ax"].plot(x,y,"o",color=colors[j],alpha=opacity,zorder=1,markersize=linewidth)

            

    def draw_region_bedmethyl(self,region,box,file,mod,min_coverage=5,linewidth=3,opacity=1,label="",color="#27ae60"):
        # bedmethyl or tsv format. tsv format is chr pos percentage

        if file=="" or file is None:
            raise KnownException("No file was provided in the basemod_freq track.")
        if not os.path.isfile(file):
            raise KnownException("The following file does not exist (in basemod_freq track): "+file)
        
        min_coverage=int(min_coverage)
        linewidth=float(linewidth)
        opacity=float(opacity)
        margin_y = (box["top"]-box["bottom"]) * 0.03
        def transform_coord(pos):
            if region.orientation=="+":
                return min(box["right"]-0.5,max(box["left"]+0.5,box["left"] + (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])))
            else:
                return min(box["right"]-0.5,max(box["left"]+0.5,box["right"] - (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])))
            
        # Read first line to determine format
        if file.endswith(".gz"):
            with gzip.open(file,"rt") as f:
                first_line = f.readline()
        else:
            with open(file) as f:
                first_line = f.readline()
            
        firstlinesplit= first_line.split('\t')

        if len(firstlinesplit)>4: # bedmethyl file
            df = smooth_methylation(create_basemod_table_bedmethyl(mod,region,min_coverage=min(4,min_coverage),file=file),
                                    w=self.smooth,start=region.start,end=region.end)
        else: # bedgraph or tsv with 3 columns
            if len(firstlinesplit)==4:
                df = pd.read_csv(file,sep="\t",header=None,names=["chr","pos","end","metPercentage"],dtype={"chr":str})
            elif len(firstlinesplit)==3:
                df = pd.read_csv(file,sep="\t",header=None,names=["chr","pos","metPercentage"],dtype={"chr":str})
            else: raise KnownException("File "+file+" (in basemod_freq track) has only "+str(len(firstlinesplit))+" columns, but at least 3 columns are required (without header): chr pos BaseModPercentage.")
            contigs=set([str(x) for x in df["chr"]])
            region = correct_region_chr(region,contigs)
            df = df.loc[(df["chr"]==region.chr) & (df["pos"]>=region.start) & (df["pos"]<=region.end),:]
            df = smooth_methylation(df,w=self.smooth,start=region.start,end=region.end)

        if df.shape[0]>0: self.is_empty=False
        x = [transform_coord(pos) for pos in df["pos"]]
        y = [box["bottom"] + margin_y + val/100 * (box["top"]-box["bottom"] -2 * margin_y) for val in df["smoothed"]]

        boundaries=[0]
        for i in range(len(x)-1):
            if abs(x[i+1]-x[i])> abs(box["right"]-box["left"]) *self.gap_frac:
                boundaries.append(i+1)
        boundaries.append(len(x))
        if self.style=="lines":
            for i in range(len(boundaries)-1):
                a,b = boundaries[i], boundaries[i+1]
                if b-a>1:
                    zorder = 1 if linewidth>1 else 0
                    box["ax"].plot(x[a:b],y[a:b],linewidth=linewidth,color=color,alpha=opacity,zorder=zorder) 
        else:
            box["ax"].plot(x,y,"o",color=color,alpha=opacity,zorder=1,markersize=linewidth)



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

def smooth_methylation(df,w=2,start=None,end=None):
    df = df.copy(deep=True)
    smoothed_values = []
    for i in range(0,df.shape[0]):
        window_start = max(i-w,0)
        window_end = min(i+w,df.shape[0]-1)
        while abs(df.loc[i,"pos"]-df.loc[window_start,"pos"]) > 100 and window_start<i: window_start+=1
        while abs(df.loc[i,"pos"]-df.loc[window_end,"pos"]) > 100 and window_end>i+1: window_end-=1
        smoothed = np.mean(df.loc[window_start:window_end,"metPercentage"])
        smoothed_values.append(smoothed)
    #df = df.loc[w:df.shape[0]-w-1,:]
    df["smoothed"] = smoothed_values
    if start is not None: df = df.loc[df["pos"]>=start,:]
    if end is not None: df = df.loc[df["pos"]<=end,:]
    return df

def create_basemod_table_bam(reads,base,mod,chr,start,end,samfile=None,fix_hardclip=False,min_coverage=5,warnings=[]):
    d={"chr":[],"pos":[],"metPercentage":[],"coverage":[]}

    pos2methylated={}
    pos2unmethylated={}
    for x in reads:
        for read in x:
            methyl = decode_read_basemods(read,[(base,mod,1)],samfile,fix_hardclip,warnings=warnings)
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

def add_entry_bedmethyl(entry,d,mod,chr,start,end,min_coverage,filename=""):
    entry_split = entry.split("\t")
    if len(entry_split)<9: raise KnownException("Wrong format for bedmethyl file "+filename+" (expected at least 10 columns; or 4 columns for bedgraph format, or 3 for simple tsv.).")
    try: pos = int(entry_split[1])
    except: raise KnownException("Wrong format for bedmethyl file "+filename+". In entry "+str(entry)+", the second value "+entry_split[1]+" should be an integer (position).")
    try: cov=int(entry_split[4])
    except: raise KnownException("Wrong format for bedmethyl file "+filename+". In entry "+str(entry)+", the fifth value "+entry_split[4]+" should be an integer (coverage).")
    if entry_split[0].lstrip("chr")==chr and pos>=start and pos<=end and entry_split[3]==mod and cov>=min_coverage:
        d["chr"].append(chr)
        d["pos"].append(pos)
        if len(entry_split)==10:
            values = entry_split[9].split(" ")
            if len(values)>1:
                try: d["metPercentage"].append(float(values[1]))
                except: raise KnownException("Wrong format for bedmethyl file "+filename+". The methylation percentage could not be converted to float ("+values[1]+").")
            else:
                raise KnownException("Wrong format for bedmethyl file "+filename+".")
        else:
            if len(entry_split)>=11:
                try: d["metPercentage"].append(float(entry_split[10]))
                except: raise KnownException("Wrong format for bedmethyl file "+filename+". The methylation percentage could not be converted to float ("+entry_split[10]+").")
            else:
                raise KnownException("Wrong format for bedmethyl file "+filename+".")

        d["coverage"].append(cov)

def create_basemod_table_bedmethyl(mod,region,file,min_coverage=5):
    d={"chr":[],"pos":[],"metPercentage":[],"coverage":[]}

    if os.path.isfile(file+".tbi") or os.path.isfile(file+".csi"):
        try:
            tabixfile = pysam.TabixFile(file)
        except: raise KnownException("Failed to open bedmethyl file "+str(file))
        region = correct_region_chr(region,tabixfile.contigs)
        for entry in tabixfile.fetch(region.chr,region.start,region.end):
            add_entry_bedmethyl(entry,d,mod,region.chr,region.start,region.end,min_coverage,filename=file)
    else:
        print("Warning: bedmethyl file "+file+ " is not indexed. Index it with tabix to speed up the figure generation.")
        if file.endswith(".gz"): f=gzip.open(file,"rt")
        else: f=open(file)
        for line in f:
            add_entry_bedmethyl(line.rstrip("\n"),d,mod,region.chr,region.start,region.end,min_coverage,filename=file)
        f.close()

    df = pd.DataFrame(d)
    return df