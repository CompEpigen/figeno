import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from figeno.utils import split_box, draw_bounding_box
from figeno.genes import read_transcripts
from figeno.ase import read_SNPs_RNA, read_SNPs_DNA

import importlib_resources as resources
import figeno.data

class ase_track:
    def __init__(self,reference="hg19",genes_file=None,ase_file=None,vcf_DNA=None,min_depth=6,color1="#e55039",color2="#4a69bd",only_exonic=False,grid=False,labels=["RNA","DNA"],labels_rotate=True,fontscale=1,
                bounding_box=False,height=10,margin_above=1.5):
        self.reference=reference
        self.genes_file = genes_file
        self.ase_file=ase_file
        self.vcf_DNA = vcf_DNA
        self.min_depth=min_depth
        self.color1=color1
        self.color2=color2
        self.only_exonic=only_exonic
        self.grid=grid
        self.labels=labels
        self.labels_rotate=labels_rotate
        self.fontscale=fontscale
    
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

        self.box_DNA=None
        self.box_RNA=None

    def draw(self, regions, box ,hmargin):
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            res=self.draw_region(regions[i][0],boxes[i])
            if not res: 
                return False
        self.draw_title(self.box_RNA,self.labels[0])
        if self.use_DNA:
            self.draw_title(self.box_DNA,self.labels[1])

    def draw_region(self,region,box):
        if self.bounding_box: draw_bounding_box(box)

        if self.genes_file is None or self.genes_file=="":
            if self.reference in ["hg19","hg38"]:
                with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                    transcripts = read_transcripts(infile,region.chr,region.start,region.end)
            else:
                raise Exception("Must provide a gene file.")
        else:
            transcripts = read_transcripts(self.genes_file,region.chr,region.start,region.end)
        if len(transcripts)==0: return False

        exons=[]
        for transcript in transcripts:
            exons = exons+transcript.exons
        
        df_ase = read_SNPs_RNA(self.ase_file,chr=region.chr,start=region.start,end=region.end,exons=exons,min_depth=self.min_depth)
        if self.only_exonic: df_ase = df_ase[df_ase["exonic"]]
        if df_ase.shape[0]==0: return False
        if self.vcf_DNA is not None:
            df_ase = read_SNPs_DNA(self.vcf_DNA,df_ase)

        self.use_DNA = "VAF_DNA" in df_ase.columns

        if self.use_DNA:
            box_DNA = {"ax":box["ax"],"left":box["left"]+(box["right"]-box["left"])*0.15,"right":box["left"]+(box["right"]-box["left"])*0.95,"bottom":box["bottom"] + (box["top"]-box["bottom"]) * 0.65, "top":box["top"]}
            box_RNA = {"ax":box["ax"],"left":box["left"]+(box["right"]-box["left"])*0.15,"right":box["left"]+(box["right"]-box["left"])*0.95,"bottom":box["bottom"] + (box["top"]-box["bottom"]) * 0.2, "top":box["bottom"] + (box["top"]-box["bottom"]) * 0.55}
        else:
            box_RNA = {"ax":box["ax"],"left":box["left"]+(box["right"]-box["left"])*0.15,"right":box["left"]+(box["right"]-box["left"])*0.95,"bottom":box["bottom"] + (box["top"]-box["bottom"]) * 0.2, "top":box["top"]}
        x_positions = draw_hist_VAF(box_RNA,list(df_ase["VAF"]),color1=self.color1,color2=self.color2)
        if self.grid: draw_grid(box_RNA)
        if self.box_RNA is None: self.box_RNA = box_RNA

        if self.use_DNA:
            tmp = draw_hist_VAF(box_DNA,list(df_ase["VAF_DNA"]),color1=self.color1,color2=self.color2)
            if self.grid: draw_grid(box_DNA)
            if self.box_DNA is None: self.box_DNA = box_DNA

    

        box_mapping = {"ax":box["ax"],"left":box["left"],"right":box["right"],"bottom":box["bottom"], "top":box["bottom"] + (box["top"]-box["bottom"]) * 0.2}
        draw_mapping_SNP(box_mapping,region,x_positions,list(df_ase["position"]))
        return True

        
       
            
    def draw_title(self,box,label):

        tick_width=1.0
        lw=1

        box["ax"].plot([box["left"],box["left"]],[box["bottom"],box["top"]],color="black",linewidth=lw,zorder=4)

        box["ax"].plot([box["left"],box["left"]-tick_width],[box["bottom"],box["bottom"]],color="black",linewidth=lw,zorder=4)
        box["ax"].text(box["left"]-tick_width*1.3,box["bottom"],"0",horizontalalignment="right",verticalalignment="bottom",fontsize=9*self.fontscale)

        box["ax"].plot([box["left"],box["left"]-tick_width],[box["bottom"]+(box["top"]-box["bottom"])*0.25,box["bottom"]+(box["top"]-box["bottom"])*0.25],color="black",linewidth=lw,zorder=4)

        box["ax"].plot([box["left"],box["left"]-tick_width],[(box["bottom"]+box["top"])/2,(box["bottom"]+box["top"])/2],color="black",linewidth=lw,zorder=4)
        box["ax"].text(box["left"]-tick_width*1.3,(box["bottom"]+box["top"])/2,"0.5",horizontalalignment="right",verticalalignment="center",fontsize=9*self.fontscale)

        box["ax"].plot([box["left"],box["left"]-tick_width],[box["bottom"]+(box["top"]-box["bottom"])*0.75,box["bottom"]+(box["top"]-box["bottom"])*0.75],color="black",linewidth=lw,zorder=4)

        box["ax"].plot([box["left"],box["left"]-tick_width],[box["top"],box["top"]],color="black",linewidth=lw,zorder=4)
        box["ax"].text(box["left"]-tick_width*1.3,box["top"],"1",horizontalalignment="right",verticalalignment="top",fontsize=9*self.fontscale)

        if len(label)>0:
            rotation = 90 if self.labels_rotate else 0
            box["ax"].text(box["left"] - 7.0,(box["top"]+box["bottom"])/2,
                        label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=12*self.fontscale)
        

def draw_hist_VAF(box,VAFs,color1="#e55039",color2="#4a69bd"):
    bar_width=(box["right"]-box["left"]) * 0.7 / len(VAFs)
    x_positions = []
    for i in range(len(VAFs)):
        x = box["left"] + (box["right"]-box["left"]) * (i+0.5)/len(VAFs)
        x_positions.append(x)
        rect = patches.Rectangle((x-bar_width/2,box["bottom"]),width=bar_width, height=(box["top"]-box["bottom"])*VAFs[i],color=color1,lw=0)
        box["ax"].add_patch(rect)
        rect = patches.Rectangle((x-bar_width/2,box["bottom"]+(box["top"]-box["bottom"])*VAFs[i]),width=bar_width,height=(box["top"]-box["bottom"])*(1-VAFs[i]),color=color2,lw=0)
        box["ax"].add_patch(rect)
    return x_positions

def draw_grid(box):
    for y in [0.25,0.5,0.75]:
        y2 = box["bottom"] + y *(box["top"]-box["bottom"])
        box["ax"].plot([box["left"],box["right"]],[y2,y2],linewidth=1,linestyle="dashed",color="gray")

def draw_mapping_SNP(box,region,x_positions,SNP_coords):
    tip_height = (box["top"] - box["bottom"]) * 0.1
    rect_width=0.15
    x_realpositions = [box["left"] + (box["right"]-box["left"]) * (coord-region.start) / (region.end-region.start) for coord in SNP_coords]
    for i in range(len(x_positions)):
        if x_positions[i]<x_realpositions[i]:
            vertices =[ (x_positions[i]-rect_width/2 , box["top"]) , (x_positions[i]-rect_width/2, box["top"]-tip_height-rect_width/2),
                        (x_realpositions[i] - rect_width/2 , box["bottom"]+tip_height-rect_width/2) ,  (x_realpositions[i] - rect_width/2 , box["bottom"]),
                        (x_realpositions[i] + rect_width/2 , box["bottom"]) , (x_realpositions[i] + rect_width/2 , box["bottom"]+tip_height),
                        (x_positions[i]+rect_width/2, box["top"]-tip_height), (x_positions[i]+rect_width/2 , box["top"])]
        else:
            vertices =[ (x_positions[i]-rect_width/2 , box["top"]) , (x_positions[i]-rect_width/2, box["top"]-tip_height),
                        (x_realpositions[i] - rect_width/2 , box["bottom"]+tip_height) ,  (x_realpositions[i] - rect_width/2 , box["bottom"]),
                        (x_realpositions[i] + rect_width/2 , box["bottom"]) , (x_realpositions[i] + rect_width/2 , box["bottom"]+tip_height-rect_width/2),
                        (x_positions[i]+rect_width/2, box["top"]-tip_height-rect_width/2), (x_positions[i]+rect_width/2 , box["top"])]
        box["ax"].add_patch(patches.Polygon(vertices,lw=0,color="black"))
        #box["ax"].plot([x_positions[i],x_positions[i]],[box["top"],box["top"]-tip_height],color="black",linestyle="-",linewidth=0.8)
        #box["ax"].plot([x_positions[i],x_realpositions[i]],[box["top"]-tip_height,box["bottom"]+tip_height],color="black",linestyle="-",linewidth=0.8)
        #box["ax"].plot([x_realpositions[i],x_realpositions[i]],[box["bottom"]+tip_height,box["bottom"]],color="black",linestyle="-",linewidth=0.8)

