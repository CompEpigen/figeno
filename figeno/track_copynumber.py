import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import pandas as pd
import importlib_resources as resources
import figeno.data
from figeno.genes import read_transcripts
from figeno.utils import KnownException, split_box,draw_bounding_box , polar2cartesian, cartesian2polar, interpolate_polar_vertices

class copynumber_track:
    def __init__(self,input_type=None,freec_ratios=None,freec_CNAs=None,CNAs=None,purple_cn=None,delly_cn=None,delly_CNAs=None,
                 ploidy=2,grid=True,grid_major=True,grid_minor=True,grid_cn=True, min_cn=None,max_cn=None,round_cn=False,
                 marker_size=0.7,color_normal="#000000",color_loss="#4a69bd",color_gain="#e55039",color_cnloh="#f6b93b", genes=[],reference="hg19",genes_file="",chr_lengths={},
                 label="CN",label_rotate=True,fontscale=1,bounding_box=True,height=20,margin_above=1.5, **kwargs):
        self.input_type=input_type
        self.freec_ratios = freec_ratios
        if self.freec_ratios=="": self.freec_ratios=None
        self.freec_CNAs = freec_CNAs
        if self.freec_CNAs=="": self.freec_CNAs=None
        self.CNAs=CNAs # Already provide a dict: chr-> list of CNAs , instead of providing freec_CNAs
        self.purple_cn=purple_cn
        if self.purple_cn=="": self.purple_cn=None
        self.delly_cn=delly_cn
        if self.delly_cn=="": self.delly_cn=None
        self.delly_CNAs=delly_CNAs
        if self.delly_CNAs=="": self.delly_CNAs=None

        try:
            self.ploidy= float(ploidy)
        except: raise KnownException("Ploidy must be a number (for copynumber track).")
        if self.ploidy<=0: raise KnownException("Ploidy must be >0.")
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
        self.reference=reference
        self.genes_file=genes_file
        self.chr_lengths=chr_lengths
        self.label=label
        self.label_rotate=label_rotate
        self.fontscale=float(fontscale)
        self.height = float(height)
        self.margin_above= float(margin_above)
        self.bounding_box=bounding_box

        if self.input_type is None:
            if self.freec_ratios is not None: self.input_type="freec"
            elif self.purple_cn is not None: self.input_type="purple"
            elif self.freec_CNAs is not None: self.input_type="freec"
            elif self.CNAs is not None: self.input_type="freec"
            elif self.delly_cn is not None: self.input_type="delly"
            elif self.delly_CNAs is not None: self.input_type="delly"
            else: raise KnownException("Please provide an input file for the copynumber track.")

        self.df_ratios=None
        if self.input_type=="freec":
            if self.freec_ratios is not None:
                try:
                    self.df_ratios = pd.read_csv(self.freec_ratios,sep="\t",dtype={"Chromosome":str,"Start":float})
                except: raise KnownException("Failed to open file "+str(self.freec_ratios))
                if (not "Chromosome" in self.df_ratios.columns) or (not "Start" in self.df_ratios.columns) or (not "Ratio" in self.df_ratios.columns):
                    raise KnownException("Invalid format for the freec ratios file. Must be tab-separated with 3 columns (with header): "\
                                        "Chromosome, Start, Ratio.")
                self.df_ratios["Chromosome"] = [x.lstrip("chr") for x in self.df_ratios["Chromosome"]]
                self.df_ratios = self.df_ratios.loc[self.df_ratios["Ratio"]>=0,:]

            if self.CNAs is None and (self.freec_CNAs is not None):
                self.CNAs = read_cna_freec(self.freec_CNAs)

            if self.df_ratios is None:
                if self.CNAs is None:
                    raise KnownException("Please provide input files for the copynumber track.") 
                else:
                    self.df_segments = read_cnsegments_CNAs(self.CNAs,round_cn=self.round_cn, chr_lengths=self.chr_lengths, ploidy=self.ploidy)
        
        elif self.input_type=="purple":
            if self.purple_cn is not None:
                self.df_segments = read_cnsegments_purple(self.purple_cn)
            else:
                raise KnownException("Please provide a copy number file for the copynumber track.")
            
        elif self.input_type=="delly":
            if self.delly_cn is not None:
                try:
                    self.df_ratios = pd.read_csv(self.delly_cn,sep="\t",dtype={"chr":str,"start":int,"end":int})
                except:
                    raise KnownException("Failed to open file "+str(self.delly_cn))
                if len(self.df_ratios.columns)!=6: raise KnownException("Expected 6 columns for the delly copy number file, but found "+str(len(self.df_ratios.columns))+".")
                self.df_ratios=self.df_ratios.iloc[:,[0,1,5]]
                self.df_ratios.columns=["Chromosome","Start","Ratio"]
                self.df_ratios["Chromosome"] = [x.lstrip("chr") for x in self.df_ratios["Chromosome"]]
                self.df_ratios["Ratio"]= self.df_ratios["Ratio"]/self.ploidy
            if self.delly_CNAs is not None:
                self.CNAs=read_cna_delly(self.delly_CNAs)
            
            if self.df_ratios is None:
                if self.CNAs is None:
                    raise KnownException("Please provide input files for the copynumber track.")
                else:
                    self.df_segments = read_cnsegments_CNAs(self.CNAs,round_cn=self.round_cn, chr_lengths=self.chr_lengths, ploidy=self.ploidy)

        self.kwargs=kwargs
        

    def draw(self, regions, box ,hmargin,warnings=[]):

        # find min and max cn across all regions
        if self.min_cn is None or self.max_cn is None:
            min_cn,max_cn = self.compute_min_max_cn(regions)
            if self.min_cn is None: self.min_cn = min_cn
            if self.max_cn is None: self.max_cn = max_cn
        self.yticks_freq=1 # by default, put one vertical tick per integer copy number. If too many copy numbers, set a higher interval between ticks.
        cn_amplitude=self.max_cn-self.min_cn
        while cn_amplitude/self.yticks_freq>8:
            if int(str(self.yticks_freq)[0])==2: self.yticks_freq*=2.5
            else: self.yticks_freq*=2
        self.yticks_freq=int(self.yticks_freq)

        if self.df_ratios is None:
            if self.input_type=="freec":
                warnings.append("Only a CNA file was provided for the copynumber track. Figeno will plot segments corresponding to CNAs, and assume "\
                                "that all other positions have a copy number equal to the ploidy. It is recommended to also provide a ratios file, "\
                                "which indicates the copy number of each bin, whether it is in a CNA or not.")
            elif self.input_type=="delly":
                warnings.append("Only a CNA file was provided for the copynumber track. Figeno will plot segments corresponding to CNAs, and assume "\
                                "that all other positions have a copy number equal to the ploidy. It is recommended to also provide a copy number file, "\
                                "which indicates the copy number of each bin, whether it is in a CNA or not.")

        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            if self.df_ratios is not None:
                self.draw_region_ratios(regions[i][0],boxes[i])
            else:
                self.draw_region_segments(regions[i][0],boxes[i])
        self.draw_title(box)
        if "bottom2" in box:
            self.draw_title({"ax":box["ax"],"left":box["left"],"right":box["right"],"top":box["top2"],"bottom":box["bottom2"]})

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the copynumber track because it is not one of the accepted parameters.")

    def draw_region_ratios(self,region,box):
        if self.bounding_box: draw_bounding_box(box)
        
        draw_grid_box(box,region,self.min_cn,self.max_cn,self.fontscale,vertical_lines=self.grid,major=self.grid_major,minor=self.grid_minor,cn=self.grid_cn,yticks_freq=self.yticks_freq)

        df_ratios_reg = self.df_ratios.loc[(self.df_ratios["Chromosome"]==region.chr) & (self.df_ratios["Start"]>=region.start) & (self.df_ratios["Start"]<=region.end),:]
        colors=[]
        x_converted=[]
        y_converted=[]
        for i in df_ratios_reg.index:
            if self.CNAs is not None: cn = bin2CNV_freec(self.CNAs,region.chr,df_ratios_reg.loc[i,"Start"],ploidy=self.ploidy)
            else: cn=df_ratios_reg.loc[i,"Ratio"]*self.ploidy
            if cn>=self.ploidy+0.5:
                colors.append(self.color_gain)
            elif cn<=self.ploidy-0.5:
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
                if self.reference in ["hg19","hg38","mm10"]:
                    with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                        transcripts = read_transcripts(infile,chr=region.chr,start=region.start,end=region.end,gene_names=self.genes)
                else:
                    raise KnownException("Please provide a genes file, if you want to higlight genes in the copy number track (this genes file is only required if you use a custom reference).")
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
        
        draw_grid_box(box,region,self.min_cn,self.max_cn,self.fontscale,vertical_lines=self.grid,yticks_freq=self.yticks_freq)
        
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

            # Ensure the segments are not too narrow
            width=abs(end-start)
            if "projection" in box and box["projection"]=="polar":
                min_width=3.14/300
                if width<min_width:
                    s=min_width-width
                    start= start - s*(1-2*(region.orientation=="-"))
                    end= end+s*(1-2*(region.orientation=="-"))
            else:
                min_width=0.5
                if width<min_width:
                    s=min_width-width
                    start= start - s*(1-2*(region.orientation=="-"))
                    end= end+s*(1-2*(region.orientation=="-"))

            if box["left"]<box["right"]:
                start = min(box["right"],max(box["left"],start))
                end = min(box["right"],max(box["left"],end))
            else:
                start = max(box["right"],min(box["left"],start))
                end = max(box["right"],min(box["left"],end))
            cn = df_segments_reg.loc[x,"copyNumber"]
            y = box["bottom"] + (box["top"]-box["bottom"]) * (cn-self.min_cn) / (self.max_cn-self.min_cn)
            rect_width = min(10,(box["top"]-box["bottom"]) * self.marker_size/0.7 /10 /  (self.max_cn-self.min_cn))
            if self.round_cn: cn = round(cn)
            vertices = [(start,y-rect_width) , (start,y+rect_width) , (end,y+rect_width) , (end, y-rect_width)]
            if "projection" in box and box["projection"]=="polar":
                vertices = interpolate_polar_vertices(vertices)
            vertices = [(xv,min(box["top"],max(box["bottom"],yv))) for (xv,yv) in vertices]
            box['ax'].add_patch(patches.Polygon(vertices,color=color,lw=0))

        # Highlight genes
        if len(self.genes)>0:
            if self.genes_file=="":
                if self.reference in ["hg19","hg38","mm10"]:
                    with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                        transcripts = read_transcripts(infile,gene_names=self.genes)
                else:
                    raise KnownException("Please provide a genes file, if you want to higlight genes in the copy number track (this genes file is only required if you use a custom reference).")
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
                if y>=self.min_cn and y<=self.max_cn and y%self.yticks_freq==0:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (self.max_cn-self.min_cn) * (y-self.min_cn)
                    xl,yl = polar2cartesian((box["left"],y_converted))
                    xl-=1
                    thetal,rl = cartesian2polar((xl,yl))
                    box["ax"].text(thetal,rl,str(y),horizontalalignment="right",verticalalignment="center",fontsize=self.fontscale*12)
        else:
            for y in range(round(self.min_cn),round(self.max_cn)+1):
                if y>=self.min_cn and y<=self.max_cn and y%self.yticks_freq==0:
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
        return (min_cn-0.2,max_cn*1.05+0.2)

def read_cna_freec(cna_freec_filename):
    CNAs={}
    if not os.path.exists(cna_freec_filename): raise KnownException("CNA file does not exist: "+str(cna_freec_filename))
    with open(cna_freec_filename,"r") as infile:
        for line in infile:
            linesplit = line.rstrip("\n").split("\t")
            if len(linesplit)<4: raise KnownException("Invalid format for freec CNA file "+cna_freec_filename+". "\
                                                      "Must be tab-separated with at least 4 columns (without header): "\
                                                      "chr, start, end, copy_number.")
            chr = linesplit[0].lstrip("chr")
            if not chr in CNAs: CNAs[chr] = []
            CNAs[chr].append((int(linesplit[1]),int(linesplit[2]),float(linesplit[3])))
    return CNAs

def read_cna_delly(cna_delly_filename):
    CNAs={}
    if not os.path.exists(cna_delly_filename): raise KnownException("CNA file does not exist: "+str(cna_delly_filename))
    with open(cna_delly_filename,"r") as infile:
        for line in infile:
            linesplit = line.rstrip("\n").split("\t")
            if len(linesplit)<4: raise KnownException("Invalid format for delly CNA file "+cna_delly_filename+". "\
                                                      "Must be tab-separated with 4 or 5 columns (without header): "\
                                                      "chr, start, end, ID (optional), copy_number.")
            chr = linesplit[0].lstrip("chr")
            if not chr in CNAs: CNAs[chr] = []
            if len(linesplit)==4: cn=float(linesplit[3])
            else: cn=float(linesplit[4])
            CNAs[chr].append((int(linesplit[1]),int(linesplit[2]),cn))
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
    if (not "chromosome" in df_segments_cn.columns) or (not "start" in df_segments_cn.columns) or (not "end" in df_segments_cn.columns) \
    or (not "copyNumber" in df_segments_cn.columns):
        raise KnownException("Invalid format for purple copy number file: "+str(purple_cn_filename)+". Must be a tsv file, with at least the following columns: "\
                             "chromosome, start, end, copyNumber (and optionally baf and bafCount to show CNLOH).")
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


def draw_grid_box(box,region,ymin,ymax,fontscale=1,vertical_lines=True,major=True,minor=True,cn=True,yticks_freq=1):
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
                if y>=ymin and y<=ymax and y%yticks_freq==0:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (ymax-ymin) * (y-ymin)
                    n_points = round(abs(box["right"]-box["left"])*100)
                    x_list = np.linspace(box["left"],box["right"],n_points)
                    y_list = [y_converted]*n_points
                    box["ax"].plot(x_list,y_list,zorder=0,linewidth=0.5,color="#AAAAAA")
        else:
            for y in range(round(ymin),round(ymax)+1):
                if y>=ymin and y<=ymax and y%yticks_freq==0:
                    y_converted = box["bottom"] + (box["top"]-box["bottom"]) / (ymax-ymin) * (y-ymin)
                    box["ax"].plot([box["left"],box["right"]],[y_converted,y_converted],zorder=0,linewidth=0.5,color="#AAAAAA")
               