import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath

Path = mpath.Path
import numpy as np
import pandas as pd
import vcfpy
from figeno.utils import correct_region_chr, split_box,draw_bounding_box, chr_to_int

class sv_track:
    def __init__(self,file=None,df_SVs=None,sv_across_regions=True,upside_down=False,label="BP",label_rotate=True,
                 color_del="#4a69bd",color_dup="#e55039",color_h2h="#8e44ad",color_t2t="#8e44ad",color_trans="#27ae60",lw=0.6,
                 fontscale=1,bounding_box=True,height=10,margin_above=1.5):
        self.file=file
        self.df_SVs = df_SVs # Can directly provide a dataframe of SVs instead of providing a vcf.
        self.sv_across_regions=sv_across_regions
        self.upside_down=upside_down
        self.label=label
        self.label_rotate=label_rotate
        self.color_del=color_del
        self.color_dup=color_dup
        self.color_h2h=color_h2h
        self.color_t2t=color_t2t
        self.color_trans=color_trans
        self.lw=lw

        self.fontscale=fontscale
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

        if self.df_SVs is None:
            if self.file is None:
                raise Exception("SVs or VCF must be provided.")
            elif self.file.endswith(".tsv"):
                self.df_SVs = self.read_SVs_tsv()
            else:
                self.df_SVs = self.read_SVs_vcf()
        if not "color" in self.df_SVs.columns:
            self.add_SV_color()

    def draw(self, regions, box ,hmargin):
        if "projection" in box and box["projection"]=="polar":
            self.draw_circle(regions,box,hmargin)
            return
        boxes = split_box(box,regions,hmargin)
        full_box = box
        if "bottom2" in box: full_box = {"ax":box["ax"],"left":box["left"],"right":box["right"],"top":box["top"],"bottom":box["bottom2"]}
        if self.bounding_box:
            if self.sv_across_regions: 
                draw_bounding_box(full_box)
            else:
                for b in boxes:
                    draw_bounding_box(b)
        
        max_dist=10
        arcs=[]
        for i in range(len(regions)):
            for j in range(i,len(regions)):
                if i!=j and not self.sv_across_regions: continue
                df_SV_regions = select_SVs_regions(self.df_SVs,regions[i][0],regions[j][0])
                for a in df_SV_regions.index:
                    if regions[i][0].orientation=="+":
                        x1 = boxes[i]["left"] + (boxes[i]["right"]-boxes[i]["left"]) / (regions[i][0].end-regions[i][0].start) * (df_SV_regions.loc[a,"pos1"]-regions[i][0].start)
                    else:
                        x1 = boxes[i]["right"] - (boxes[i]["right"]-boxes[i]["left"]) / (regions[i][0].end-regions[i][0].start) * (df_SV_regions.loc[a,"pos1"]-regions[i][0].start)
                    if regions[j][0].orientation=="+":                    
                        x2 = boxes[j]["left"] + (boxes[j]["right"]-boxes[j]["left"]) / (regions[j][0].end-regions[j][0].start) * (df_SV_regions.loc[a,"pos2"]-regions[j][0].start)
                    else:
                        x2 = boxes[j]["right"] - (boxes[j]["right"]-boxes[j]["left"]) / (regions[j][0].end-regions[j][0].start) * (df_SV_regions.loc[a,"pos2"]-regions[j][0].start)
                    y1 = boxes[i]["top"] if self.upside_down or ("upside_down" in boxes[i]) else boxes[i]["bottom"]
                    y2 = boxes[j]["top"] if self.upside_down or ("upside_down" in boxes[j]) else boxes[j]["bottom"]
                    arcs.append({"ax":box["ax"],"x1":x1,"x2":x2,"color":df_SV_regions.loc[a,"color"],"y1":y1,"y2":y2,"upside_down":self.upside_down or ("upside_down" in boxes[i] and "upside_down" in boxes[j])})
                    max_dist = max(max_dist,abs(x1-x2))
        for arc in arcs:
            arc["height"] = (box["top"]-box["bottom"]) * 1.8 * np.log(10+abs(arc["x1"]-arc["x2"])) / np.log(10+max_dist)
            add_arc(**arc,lw=self.lw)

        # Show SVs leading to other regions
        for i in range(len(regions)):
            if self.sv_across_regions:
                df_SV_others = select_SVs_otherregions(self.df_SVs,regions[i][0],regions)
            else:
                df_SV_others = select_SVs_otherregions(self.df_SVs,regions[i][0],[regions[i]])
            for a in df_SV_others.index:
                x=boxes[i]["left"] + (boxes[i]["right"]-boxes[i]["left"]) / (regions[i][0].end-regions[i][0].start) * (df_SV_others.loc[a,"pos1"]-regions[i][0].start)
                if self.upside_down or "upside_down" in boxes[i]:
                    y1 = boxes[i]["top"]
                    y2=boxes[i]["top"] - (boxes[i]["top"]-boxes[i]["bottom"]) * chr_to_int(df_SV_others.loc[a,"chr2"]) / 30
                    boxes[i]["ax"].text(x=x,y=y2-0.2,s=df_SV_others.loc[a,"chr2"],horizontalalignment="center",verticalalignment="top",color="black",fontsize=6*self.fontscale)
                else:
                    y1 = boxes[i]["bottom"]
                    y2=boxes[i]["bottom"] + (boxes[i]["top"]-boxes[i]["bottom"]) * chr_to_int(df_SV_others.loc[a,"chr2"]) / 30
                    boxes[i]["ax"].text(x=x,y=y2,s=df_SV_others.loc[a,"chr2"],horizontalalignment="center",verticalalignment="bottom",color="black",fontsize=6*self.fontscale)
                boxes[i]["ax"].plot([x,x],[y1,y2],color=df_SV_others.loc[a,"color"],linewidth=0.5)

        self.draw_title(full_box)


    def draw_circle(self, regions, box ,hmargin):
        boxes = split_box(box,regions,hmargin)
        if self.bounding_box:
            if self.sv_across_regions: draw_bounding_box(box)
            else:
                for b in boxes:
                    draw_bounding_box(b)

    
        #df_SVs = self.read_SVs_vcf()
        max_dist=10
        arcs=[]
        for i in range(len(regions)):
            for j in range(i,len(regions)):
                if i!=j and not self.sv_across_regions: continue
                df_SV_regions = select_SVs_regions(self.df_SVs,regions[i][0],regions[j][0])
                for a in df_SV_regions.index:
                    if regions[i][0].orientation=="+":
                        x1 = boxes[i]["left"] + (boxes[i]["right"]-boxes[i]["left"]) / (regions[i][0].end-regions[i][0].start) * (df_SV_regions.loc[a,"pos1"]-regions[i][0].start)
                    else:
                        x1 = boxes[i]["right"] - (boxes[i]["right"]-boxes[i]["left"]) / (regions[i][0].end-regions[i][0].start) * (df_SV_regions.loc[a,"pos1"]-regions[i][0].start)
                    if regions[j][0].orientation=="+":
                        x2 = boxes[j]["left"] + (boxes[j]["right"]-boxes[j]["left"]) / (regions[j][0].end-regions[j][0].start) * (df_SV_regions.loc[a,"pos2"]-regions[j][0].start)
                    else:
                        x2 = boxes[j]["right"] - (boxes[j]["right"]-boxes[j]["left"]) / (regions[j][0].end-regions[j][0].start) * (df_SV_regions.loc[a,"pos2"]-regions[j][0].start)
                    y = box["top"] if self.upside_down else box["bottom"] - 0.2
                    control_theta = (x1+x2)/2
                    control_r = 0
                    if i==j: control_r = y * max(0,0.9-abs(df_SV_regions.loc[a,"pos1"]-df_SV_regions.loc[a,"pos2"]) / 400000000)

                    pp1 = patches.PathPatch(mpath.Path([(x1, y), (control_theta, control_r), (x2, y)],[Path.MOVETO, Path.CURVE3, Path.CURVE3]),facecolor="none",edgecolor=df_SV_regions.loc[a,"color"])
                    boxes[0]["ax"].add_patch(pp1)



        self.draw_title(box)



    def draw_title(self,box):
        if len(self.label)>0 and ((not "projection" in box) or box["projection"]!="polar"):
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 3.0,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
    def read_SVs_tsv(self):
        df_SVs = pd.read_csv(self.file,sep="\t",dtype={"chr1":str,"chr2":str})
        if not "color" in df_SVs.columns:
            if "strand1" in df_SVs.columns and "strand2" in df_SVs.columns:
                colors=[]
                for i in df_SVs.index:
                    chr1,chr2= df_SVs.loc[i,"chr1"],df_SVs.loc[i,"chr2"]
                    if chr1!=chr2: 
                        colors.append(self.color_trans)
                    else:
                        pos1,pos2=df_SVs.loc[i,"pos1"],df_SVs.loc[i,"pos2"]
                        strand1,strand2=df_SVs.loc[i,"strand1"],df_SVs.loc[i,"strand2"]
                        if pos2<pos1: strand1,strand2=strand2,strand1
                        if strand1=="-" and strand2=="+": colors.append(self.color_del)
                        elif strand1=="+" and strand2=="-": colors.append(self.color_dup)
                        elif strand1=="+" and strand2=="-": colors.append(self.color_h2h)
                        else: colors.append(self.color_t2t)
                df_SVs["color"] = colors
            else:
                df_SVs["color"]=["black" for i in df_SVs.index]

        return df_SVs


    def read_SVs_vcf(self):
        SVs=[]
        reader = vcfpy.Reader.from_path(self.file)
        for record in reader:
            chr1 = record.CHROM.lstrip("chr")
            pos1 = record.POS
            if record.INFO["SVTYPE"] in ["BND","TRA"]:
                chr2 = record.ALT[0].mate_chrom.lstrip("chr")
                pos2 = record.ALT[0].mate_pos 
            else:
                chr2 = chr1
                pos2 = record.INFO["END"]
            
            if record.INFO["SVTYPE"] == "DEL":
                color= self.color_del
            elif record.INFO["SVTYPE"] == "DUP":
                color = self.color_dup
            elif record.INFO["SVTYPE"] =="INV":
                color = self.color_h2h # TODO
            else:  #BND
                if record.CHROM != record.ALT[0].mate_chrom:
                    color = self.color_trans
                else:
                    if record.POS > record.ALT[0].mate_pos:
                        continue
                    if record.ALT[0].orientation =="-" and record.ALT[0].mate_orientation == "+":
                        color = self.color_del
                    elif record.ALT[0].orientation =="+" and record.ALT[0].mate_orientation == "-":
                        color = self.color_dup
                    else:
                        color = self.color_h2h # TODO
            #if chr1==chr2 and abs(pos2-pos1)>max_SV_size: max_SV_size =  abs(pos2-pos1)
            if not (chr1,pos1,chr2,pos2,color) in SVs:
                SVs.append((chr1,pos1,chr2,pos2,color))
            if not (chr2,pos2,chr1,pos1,color) in SVs:
                SVs.append((chr2,pos2,chr1,pos1,color))
        df_SVs = pd.DataFrame(SVs,columns=["chr1","pos1","chr2","pos2","color"])
        return df_SVs
    

def add_SV_color(self):
    colors=[]
    for x in self.df_SVs.index:
        if self.df_SVs.loc[x,"chr1"]!=self.df_SVs.loc[x,"chr2"]: colors.append(self.color_trans)
        elif self.df_SVs.loc[x,"orientation1"]=="-" and self.df_SVs.loc[x,"orientation2"]=="+": colors.append(self.color_del)
        elif self.df_SVs.loc[x,"orientation1"]=="+" and self.df_SVs.loc[x,"orientation2"]=="": colors.append(self.color_dup)
        elif self.df_SVs.loc[x,"orientation1"]=="-" and self.df_SVs.loc[x,"orientation2"]=="-": colors.append(self.color_t2t)
        elif self.df_SVs.loc[x,"orientation1"]=="+" and self.df_SVs.loc[x,"orientation2"]=="+": colors.append(self.color_h2h)
    self.df_SVs["color"] = colors
   

def select_SVs_regions(df_SVs,region1,region2):
    indices=[]
    for i in df_SVs.index:
        if df_SVs.loc[i,"chr1"] == region1.chr and df_SVs.loc[i,"chr2"] == region2.chr:
            if df_SVs.loc[i,"pos1"]>=region1.start and df_SVs.loc[i,"pos1"]<=region1.end:
                if df_SVs.loc[i,"pos2"]>=region2.start and df_SVs.loc[i,"pos2"]<=region2.end:
                    indices.append(i)
    df_SVs2 = df_SVs.loc[indices,:].copy(deep=True)
    return df_SVs2

def select_SVs_otherregions(df_SVs,region1,otherregions):
    """Select SVs from region1 which do not end up in any other regions"""
    indices=[]
    for i in df_SVs.index:
        keep=False
        if df_SVs.loc[i,"chr1"] == region1.chr and df_SVs.loc[i,"pos1"]>=region1.start and df_SVs.loc[i,"pos1"]<=region1.end:
            keep=True
            for r in otherregions:
                region2 = r[0]
                if df_SVs.loc[i,"chr2"]==region2.chr and df_SVs.loc[i,"pos2"]>=region2.start and df_SVs.loc[i,"pos2"]<=region2.end:
                    keep=False
        if keep: indices.append(i)
    df_SVs2 = df_SVs.loc[indices,:].copy(deep=True)
    return df_SVs2

def add_arc(ax,x1,x2,y1,y2,height,color,lw,upside_down=False):
    if y1==y2:
        if upside_down:
            arc = patches.Arc(xy=((x1+x2)/2,y1),width=abs(x2-x1),height=height,theta1=180,theta2=360,fill=False,color=color,lw=lw)
        else:
            arc = patches.Arc(xy=((x1+x2)/2,y1),width=abs(x2-x1),height=height,theta1=0,theta2=180,fill=False,color=color,lw=lw)
        ax.add_patch(arc)
    else:
        pp1 = patches.PathPatch(mpath.Path([(x1, y1), (x1, (y1+y2)/2), (x2, (y1+y2)/2),(x2,y2)],[Path.MOVETO, Path.CURVE4, Path.CURVE4,Path.CURVE4]),facecolor="none",edgecolor=color,lw=lw)
        ax.add_patch(pp1)
   