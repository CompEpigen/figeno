import pyBigWig
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import pandas as pd
import importlib_resources as resources
import figeno.data
from figeno.genes import read_transcripts
from figeno.utils import split_box,draw_bounding_box , polar2cartesian, cartesian2polar, interpolate_polar_vertices

class copynumber_track:
    def __init__(self,freec_ratios=None,freec_CNAs=None,CNAs=None,purple_cn=None,ploidy=2,grid=True,grid_major=True,grid_minor=True,grid_cn=True, min_cn=None,max_cn=None,round_cn=False,
                 marker_size=0.7,color_normal="#000000",color_loss="#4a69bd",color_gain="#e55039",color_cnloh="#f6b93b", genes=[],reference="hg19",genes_file="",chr_lengths={},
                 label="CN",label_rotate=True,fontscale=1,bounding_box=True,height=20,margin_above=1.5):
        self.freec_ratios = freec_ratios
        if self.freec_ratios=="": self.freec_ratios=None
        self.freec_CNAs = freec_CNAs
        if self.freec_CNAs=="": self.freec_CNAs=None
        self.CNAs=CNAs # Already provide a dict: chr-> list of CNAs , instead of providing freec_CNAs
        self.ploidy= float(ploidy)
        self.purple_cn=purple_cn
        if self.purple_cn=="": self.purple_cn=None
        self.grid=grid # True for showing a grid for the axes
        self.grid_major = grid_major # vertical lines for major ticks
        self.grid_minor=grid_minor # vertical line for minor ticks
        self.grid_cn = grid_cn # horizontal lines for copy numbers
        self.min_cn=min_cn
        if self.min_cn=="": self.min_cn = None
        if self.min_cn is not None: self.min_cn=float(self.min_cn)
        self.max_cn = max_cn
        if self.max_cn=="": self.max_cn = None
        if self.max_cn is not None: self.max_cn=float(self.max_cn)
        self.round_cn=round_cn
        self.marker_size=float(marker_size)
        self.color_normal=color_normal
        self.color_loss=color_loss
        self.color_gain=color_gain
        self.color_cnloh=color_cnloh
        self.genes = genes
        if "," in self.genes: self.genes = self.genes.split(",")
        self.reference=reference
        self.genes_file=genes_file
        self.chr_lengths=chr_lengths
        self.label=label
        self.label_rotate=label_rotate
        self.fontscale=fontscale
        self.height = height
        self.bounding_box=bounding_box
        self.margin_above=margin_above

        self.df_ratios=None
        if self.freec_ratios is not None:
            self.df_ratios = pd.read_csv(self.freec_ratios,sep="\t",dtype={"Chromosome":str,"Start":float})
            self.df_ratios["Chromosome"] = [x.lstrip("chr") for x in self.df_ratios["Chromosome"]]
            self.df_ratios = self.df_ratios.loc[self.df_ratios["Ratio"]>=0,:]
        if self.CNAs is None and (self.freec_CNAs is not None):
            self.CNAs = read_cna_freec(self.freec_CNAs)
        if self.purple_cn is not None:
            self.df_segments = read_cnsegments_purple(self.purple_cn)
        if self.df_ratios is None and self.purple_cn is None:
            if self.CNAs is None:
                raise Exception("Must provide either copy number ratios, CN segments or CNAs.")
            else:
                self.df_segments = read_cnsegments_CNAs(self.CNAs,round_cn=self.round_cn, chr_lengths=self.chr_lengths, ploidy=self.ploidy)
        

    def draw(self, regions, box ,hmargin):

        # find min and max cn across all regions
        if self.min_cn is None or self.max_cn is None:
            min_cn,max_cn = self.compute_min_max_cn(regions)
            if self.min_cn is None: self.min_cn = min_cn
            if self.max_cn is None: self.max_cn = max_cn

        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            if self.df_ratios is not None:
                self.draw_region_ratios(regions[i][0],boxes[i])
            else:
                self.draw_region_segments(regions[i][0],boxes[i])
        self.draw_title(box)
        if "bottom2" in box:
            self.draw_title({"ax":box["ax"],"left":box["left"],"right":box["right"],"top":box["top2"],"bottom":box["bottom2"]})

        #import matplotlib.path as mpath
        #Path = mpath.Path
        #pp1 = patches.PathPatch(mpath.Path([(boxes[0]["left"], boxes[0]["bottom"]),(0, 0),  (boxes[1]["right"], boxes[0]["bottom"]) ],[Path.MOVETO, Path.CURVE3, Path.CURVE3]),facecolor="none",edgecolor="red")
        #boxes[0]["ax"].add_patch(pp1)

    def draw_region_ratios(self,region,box):
        if self.bounding_box: draw_bounding_box(box)
        
        draw_grid_box(box,region,self.min_cn,self.max_cn,self.fontscale,vertical_lines=self.grid,major=self.grid_major,minor=self.grid_minor,cn=self.grid_cn)

        #elif purple_cn_filename is not None: # Purple
        #    df_segments_cn = read_cnsegments_purple(purple_cn_filename=purple_cn_filename,scale=1e6)
        #elif CNAs is not None and chr_arms is not None: # Only CNAs (!=normal cn) are provided
        #    df_segments_cn = read_cnsegments_CNAs(CNAs,chromosomes,round_cn=round_cn,chr_lengths=chr_lengths,scale=1e6)
        #elif freec_CNAs_filename is not None and chr_arms is not None:
        #    CNAs = read_cna_freec(cna_freec_filename=freec_CNAs_filename,return_type=False)
        #    df_segments_cn = read_cnsegments_CNAs(CNAs,chromosomes,round_cn=round_cn,chr_lengths=chr_lengths,scale=1e6)
        df_ratios_reg = self.df_ratios.loc[(self.df_ratios["Chromosome"]==region.chr) & (self.df_ratios["Start"]>=region.start) & (self.df_ratios["Start"]<=region.end),:]
        colors=[]
        x_converted=[]
        y_converted=[]
        for i in df_ratios_reg.index:
            cn = bin2CNV_freec(self.CNAs,region.chr,df_ratios_reg.loc[i,"Start"],ploidy=self.ploidy)
            if cn>self.ploidy:
                colors.append(self.color_gain)
            elif cn<self.ploidy:
                colors.append(self.color_loss)
            else:
                colors.append(self.color_normal)
            if region.orientation=="+":
                x_converted.append(box["left"]+(box["right"]-box["left"]) / (region.end-region.start) * (df_ratios_reg.loc[i,"Start"]-region.start))
            else:
                x_converted.append(box["right"]-(box["right"]-box["left"]) / (region.end-region.start) * (df_ratios_reg.loc[i,"Start"]-region.start))
            y=box["bottom"] + (box["top"]-box["bottom"]) * (self.ploidy*df_ratios_reg.loc[i,"Ratio"]-self.min_cn) / (self.max_cn-self.min_cn)
            y_converted.append(max(min(y,box["top"]),box["bottom"]))
        if box["right"]>box["left"]: 
            x_converted = [min(box["right"]-0.2,max(box["left"]+0.2,x)) for x in x_converted]
        else:
            x_converted = [min(box["left"],max(box["right"],x)) for x in x_converted]
        box["ax"].scatter(x_converted,y_converted,c=colors,s=self.marker_size,marker="o",rasterized=True)

        # Highlight genes
        if len(self.genes)>0:
            if self.genes_file is None or self.genes_file=="":
                if self.reference in ["hg19","hg38"]:
                    with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                        transcripts = read_transcripts(infile,gene_names=self.genes)
                else:
                    raise Exception("Must provide a gene file.")
            else:
                transcripts = read_transcripts(self.genes_file,gene_names=self.genes)
            for transcript in transcripts:
                if transcript.chr ==region.chr and transcript.start<=region.end and transcript.end>=region.start:
                    pos = (transcript.start + transcript.end) / 2
                    pos=min(region.end,max(region.start,pos))
                    cn=self.ploidy
                    cn_end = self.ploidy
                    tx_end = transcript.end if transcript.strand=="+" else transcript.start
                    for i in df_ratios_reg.index:
                        if (tx_end>=df_ratios_reg.loc[i,"Start"]) and (tx_end<=df_ratios_reg.loc[i,"Start"] + 100000):
                            cn_end = self.ploidy * df_ratios_reg.loc[i,"Ratio"]
                        if (pos>=df_ratios_reg.loc[i,"Start"]) and (pos<=df_ratios_reg.loc[i,"Start"] + 100000):
                            cn = self.ploidy * df_ratios_reg.loc[i,"Ratio"]
                    cn = max(cn,cn_end)
                    x_gene = box["left"] + (box["right"]-box["left"]) * (pos-region.start) / (region.end-region.start)
                    y_gene = box["bottom"] + (box["top"]-box["bottom"]) * (cn-self.min_cn) / (self.max_cn - self.min_cn)
                    box["ax"].plot([x_gene],[y_gene],marker="x", markersize=5, markeredgecolor="red")
                    box["ax"].text(x_gene,y_gene + 2,transcript.name,verticalalignment="bottom",horizontalalignment="center",style='italic',fontsize=self.fontscale*8,color="#222222")

    def draw_region_segments(self,region,box):
        if self.bounding_box: draw_bounding_box(box)
        
        draw_grid_box(box,region,self.min_cn,self.max_cn,self.fontscale,vertical_lines=self.grid)
        
        df_segments_reg = self.df_segments.loc[(self.df_segments["chromosome"]==region.chr) & (self.df_segments["end"]>=region.start) & (self.df_segments["start"]<=region.end),:]
        
        for x in df_segments_reg.index:
            expected_CN=self.ploidy
            if df_segments_reg.loc[x,"copyNumber"]<0.8*expected_CN: color = self.color_loss
            elif df_segments_reg.loc[x,"copyNumber"]>1.2*expected_CN: color = self.color_gain
            elif "bafCount" in df_segments_reg.columns and df_segments_reg.loc[x,"bafCount"]>30 and df_segments_reg.loc[x,"baf"]>0.88: color = self.color_cnloh
            else: color = self.color_normal
            if region.orientation=="+":
                start = box["left"] + (box["right"]-box["left"]) * (df_segments_reg.loc[x,"start"]-region.start) / (region.end-region.start)
                end = box["left"] + (box["right"]-box["left"]) * (df_segments_reg.loc[x,"end"]-region.start) / (region.end-region.start)
            else:
                start = box["right"] - (box["right"]-box["left"]) * (df_segments_reg.loc[x,"start"]-region.start) / (region.end-region.start)
                end = box["right"] - (box["right"]-box["left"]) * (df_segments_reg.loc[x,"end"]-region.start) / (region.end-region.start)
            if end-start<0.5:
                s = 0.5-end+start
                start = start-s/2
                end = end+s/2
            if box["left"]<box["right"]:
                start = min(box["right"],max(box["left"],start))
                end = min(box["right"],max(box["left"],end))
            else:
                start = max(box["right"],min(box["left"],start))
                end = max(box["right"],min(box["left"],end))
            cn = df_segments_reg.loc[x,"copyNumber"]
            y = box["bottom"] + (box["top"]-box["bottom"]) * (cn-self.min_cn) / (self.max_cn-self.min_cn)
            rect_width = min(10,(box["top"]-box["bottom"]) /10 /  (self.max_cn-self.min_cn))
            if self.round_cn: cn = round(cn)
            vertices = [(start,y-rect_width) , (start,y+rect_width) , (end,y+rect_width) , (end, y-rect_width)]
            if "projection" in box and box["projection"]=="polar":
                vertices = interpolate_polar_vertices(vertices)
            vertices = [(xv,min(box["top"],max(box["bottom"],yv))) for (xv,yv) in vertices]
            box['ax'].add_patch(patches.Polygon(vertices,color=color,lw=0))

        # Highlight genes
        if len(self.genes)>0:
            if self.genes_file=="":
                if self.reference in ["hg19","hg38"]:
                    with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                        transcripts = read_transcripts(infile,gene_names=self.genes)
                else:
                    raise Exception("Must provide a gene file.")
            else:
                transcripts = read_transcripts(self.genes_file,gene_names=self.genes)
            for transcript in transcripts:
                if transcript==False: continue
                if transcript.chr ==region.chr and transcript.start<=region.end and transcript.end>=region.start:
                    pos = (transcript.start + transcript.end) / 2
                    pos=min(region.end,max(region.start,pos))
                    tx_end = transcript.end if transcript.strand=="+" else transcript.start
                    cn_end = self.ploidy
                    cn=self.ploidy
                    for i in df_segments_reg.index:
                        if (tx_end>=df_segments_reg.loc[i,"start"]) and (tx_end<=df_segments_reg.loc[i,"end"]):
                            cn_end = df_segments_reg.loc[i,"copyNumber"]
                    for i in df_segments_reg.index:
                        if (pos>=df_segments_reg.loc[i,"start"]) and (pos<=df_segments_reg.loc[i,"end"]):
                            cn = df_segments_reg.loc[i,"copyNumber"]
                    cn = max(cn,cn_end)
                    x_gene = box["left"] + (box["right"]-box["left"]) * (pos-region.start) / (region.end-region.start)
                    y_gene = box["bottom"] + (box["top"]-box["bottom"]) * (cn-self.min_cn) / (self.max_cn - self.min_cn)
                    box["ax"].plot([x_gene],[y_gene],marker="x", markersize=5, markeredgecolor="red")
                    box["ax"].text(x_gene,y_gene + 2,transcript.name,verticalalignment="bottom",horizontalalignment="center",style='italic',fontsize=self.fontscale*10,color="#222222")





    def draw_title(self,box):
        if len(self.label)>0:
            if "projection" in box and box["projection"]=="polar":
                x,y= polar2cartesian((box["left"],(box["top"]+box["bottom"])/2))
                x-= 3
                theta,r = cartesian2polar((x,y))
                box["ax"].text(theta,r,self.label,rotation=90,
                               horizontalalignment="right",verticalalignment="center",fontsize=12*self.fontscale)
                
            else:
                rotation = 90 if self.label_rotate else 0
                box["ax"].text(box["left"] - 3.0,(box["top"]+box["bottom"])/2,
                            self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
        
        if "projection" in box and box["projection"]=="polar":
            for y in range(round(self.min_cn),round(self.max_cn)+1):
                if y>=self.min_cn and y<=self.max_cn:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (self.max_cn-self.min_cn) * (y-self.min_cn)
                    xl,yl = polar2cartesian((box["left"],y_converted))
                    xl-=1
                    thetal,rl = cartesian2polar((xl,yl))
                    box["ax"].text(thetal,rl,str(y),horizontalalignment="right",verticalalignment="center",fontsize=self.fontscale*12)
        else:
            for y in range(round(self.min_cn),round(self.max_cn)+1):
                if y>=self.min_cn and y<=self.max_cn:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (self.max_cn-self.min_cn) * (y-self.min_cn)
                    box["ax"].add_patch(patches.Rectangle((box["left"]-0.8,y_converted-0.1),width=0.8,height=0.2,lw=0,color="black"))
                    #box["ax"].plot([box["left"],box["left"]-0.5],[y_converted,y_converted],color="black",linewidth=0.8,zorder=4)
                    box["ax"].text(box["left"]-1.2,y_converted,str(y),horizontalalignment="right",verticalalignment="center",fontsize=self.fontscale*6)

                    
    def compute_min_max_cn(self,regions):
        min_cn=1.2
        max_cn=3.5
        for reg in regions:
            region = reg[0]
            if self.df_ratios is not None:
                df_ratios_reg = self.df_ratios.loc[(self.df_ratios["Chromosome"]==region.chr) & (self.df_ratios["Start"]>=region.start) & (self.df_ratios["Start"]<=region.end),:]
                min_cn = min(min_cn,np.min(df_ratios_reg["Ratio"]*self.ploidy))
                max_cn = max(max_cn,np.max(df_ratios_reg["Ratio"]*self.ploidy))
            else:
                df_segments_reg =  self.df_segments.loc[(self.df_segments["chromosome"]==region.chr) & (self.df_segments["end"]>=region.start) & (self.df_segments["start"]<=region.end),:]
                min_cn = min(min_cn,np.min(df_segments_reg["copyNumber"]))
                max_cn = max(max_cn,np.max(df_segments_reg["copyNumber"]))
        return (min_cn-0.2,max_cn+0.2)

def read_cna_freec(cna_freec_filename,return_type=False):
    CNAs={}
    with open(cna_freec_filename,"r") as infile:
        for line in infile:
            linesplit = line.rstrip("\n").split("\t")
            chr = linesplit[0].lstrip("chr")
            if not chr in CNAs: CNAs[chr] = []
            if return_type:
                CNAs[chr].append((int(linesplit[1]),int(linesplit[2]),linesplit[4]))
            else:
                CNAs[chr].append((int(linesplit[1]),int(linesplit[2]),int(linesplit[3])))
    return CNAs

def bin2CNV_freec(CNAs,chr,start,ploidy):
    result=ploidy
    if CNAs is not None:
        if chr in CNAs:
            for CNA in CNAs[chr]:
                if start>=CNA[0] and start<=CNA[1]:
                    result = CNA[2]
    return result

def read_cnsegments_purple(purple_cn_filename):
    df_segments_cn = pd.read_csv(purple_cn_filename,sep="\t",dtype={"chromosome":str})
    df_segments_cn["chromosome"] = [x.lstrip("chr") for x in df_segments_cn["chromosome"] ]
    return df_segments_cn

def read_cnsegments_CNAs(CNAs,round_cn=False,chr_lengths={},ploidy=2):
    d={"chromosome":[],"start":[],"end":[],"copyNumber":[]}
    for chr in [str(x) for x in range(1,23)] + ["X","Y"]:
        last_x=0
        if not chr in CNAs: CNAs_chr = []
        else: CNAs_chr = CNAs[chr]
        for (start,end,cn) in CNAs_chr:
            if end-start>10000:
                if abs(start-last_x)>100000:
                    d["chromosome"].append(chr)
                    d["start"].append(last_x)
                    d["end"].append(start)
                    d["copyNumber"].append(ploidy)
                if round_cn: cn = round(cn)
                d["chromosome"].append(chr)
                d["start"].append(start)
                d["end"].append(end)
                d["copyNumber"].append(cn)
                last_x=end
        if chr in chr_lengths and abs(last_x-chr_lengths[chr])>200000: 
            d["chromosome"].append(chr)
            d["start"].append(last_x)
            d["end"].append(chr_lengths[chr])
            d["copyNumber"].append(ploidy)
    return pd.DataFrame(d)


def draw_grid_box(box,region,ymin,ymax,fontscale=1,vertical_lines=True,major=True,minor=True,cn=True):
    # Vertical lines
    if vertical_lines:
        if "projection" in box and box["projection"]=="polar":
            if (region.end-region.start) / abs(box["right"]-box["left"]) < 200000000:
                minor_scale = 1e6
                major_scale = 1e7
            else:
                minor_scale = 1e7
                major_scale = 1e8
        else:
            if (region.end-region.start) / abs(box["right"]-box["left"]) < 4000000:
                minor_scale = 1e6
                major_scale = 1e7
            else:
                minor_scale = 1e7
                major_scale = 1e8

        
        pos = (region.start//minor_scale) * minor_scale
        if pos< region.start: pos+=minor_scale
        if pos==region.start: pos+=minor_scale
        while pos<region.end:
            x = box["left"] + (box["right"]-box["left"]) * (pos-region.start) / (region.end-region.start)
            if pos%major_scale==0 and major: box["ax"].plot([x,x],[box["top"],box["bottom"]],zorder=0,linewidth=0.5,color="#AAAAAA")
            else:
                if minor:
                    box["ax"].plot([x,x],[box["top"],box["bottom"]],zorder=0,linewidth=0.2,color="#AAAAAA",linestyle="dashed")
            pos+=minor_scale

    # Horizontal lines
    if cn:
        if "projection" in box and box["projection"]=="polar":
            for y in range(round(ymin),round(ymax)+1):
                if y>=ymin and y<=ymax:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (ymax-ymin) * (y-ymin)
                    n_points = round(abs(box["right"]-box["left"])*100)
                    x_list = np.linspace(box["left"],box["right"],n_points)
                    y_list = [y_converted]*n_points
                    box["ax"].plot(x_list,y_list,zorder=0,linewidth=0.5,color="#AAAAAA")
        else:
            for y in range(round(ymin),round(ymax)+1):
                if y>=ymin and y<=ymax:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (ymax-ymin) * (y-ymin)
                    box["ax"].plot([box["left"],box["right"]],[y_converted,y_converted],zorder=0,linewidth=0.5,color="#AAAAAA")
               