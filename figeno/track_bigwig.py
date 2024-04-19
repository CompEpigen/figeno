import os
import pyBigWig
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd 
from figeno.utils import correct_region_chr, split_box,draw_bounding_box, interpolate_polar_vertices, polar2cartesian, cartesian2polar

class bigwig_track:
    def __init__(self,file,n_bins=500,scale="auto",scale_max=None,scale_pos="corner",color="gray",upside_down=False,label="",label_rotate=False,fontscale=1,
                 bounding_box=False,height=10,margin_above=1.5):
        if not os.path.exists(file): raise Exception("The following file does not exist (in bigwig track): "+file)
        try:
            self.bw = pyBigWig.open(file)
        except:
            raise Exception("Error when opening the bigwig file: "+file)
        self.n_bins = int(n_bins)
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.scale = scale
        self.scale_max= scale_max
        self.upside_down=upside_down
        if self.scale_max is not None: self.scale_max=float(self.scale_max)
        self.scale_pos = scale_pos
        if scale=="auto per region": self.scale_pos = "corner all" # If each region has its own scale, we cannot use one global label for the whole track
        self.fontscale=fontscale
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

    def draw(self, regions, box ,hmargin):
        # Assign bins to regions depending on their sizes
        total_length_regions = np.sum([abs(reg[0].end-reg[0].start) for reg in regions])
        bins_regions = [max(1,int(self.n_bins/total_length_regions * abs(reg[0].end-reg[0].start))) for reg in regions]
        # Autoscale across all regions
        if self.scale=="auto": self.scale_max = self.compute_max_regions(regions,bins_regions)

        if self.scale=="custom" and isinstance(self.scale_max,str) and "," in self.scale_max:
            scale_max_regions = [float(x) for x in self.scale_max.split(",")]
        else:
            scale_max_regions=None

        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            show_scale_inside = self.scale_pos=="corner all" or (self.scale_pos=="corner" and i==0)
            if scale_max_regions is not None:
                self.scale_max = scale_max_regions[i]
            self.draw_region(regions[i][0],boxes[i],nbins=bins_regions[i],show_scale_inside=show_scale_inside)
        self.draw_title(box)

    def draw_region(self,region,box,nbins,show_scale_inside):
        if self.bounding_box: draw_bounding_box(box)

        region = correct_region_chr(region,self.bw.chroms())
        values_binned=self.bw.stats(region.chr,region.start,region.end,nBins=nbins,exact=abs(region.end-region.start)<1000000,type="mean")
        values_binned = [x if x is not None else 0 for x in values_binned]
        if region.orientation=="-": values_binned = values_binned[::-1]
        n_bases_per_bin = (region.end-region.start) // nbins
        
        rect_width = (box["right"] - box["left"]) / nbins

        if self.scale=="auto per region": self.scale_max = np.max(values_binned) * 1.1
        if self.scale_max<=0: self.scale_max=0.001
        values_binned = [max(min(self.scale_max,x),0) for x in values_binned]
        
        if not self.upside_down:
            polygon_vertices=[(box["right"],box["bottom"]) , (box["left"],box["bottom"])]
        else:
            polygon_vertices=[(box["right"],box["top"]) , (box["left"],box["top"])]
        for i in range(nbins):
            x=box["left"] + i*n_bases_per_bin/(region.end-region.start) * (box["right"] - box["left"])
            if not self.upside_down:
                y= box["bottom"] + values_binned[i]/self.scale_max * (box["top"]-box["bottom"])
            else:
                y= box["top"] - values_binned[i]/self.scale_max * (box["top"]-box["bottom"])
            polygon_vertices.append((x,y))
        if "projection" in box and box["projection"]=="polar":
            polygon_vertices = interpolate_polar_vertices(polygon_vertices)
        polygon = patches.Polygon(polygon_vertices,lw=0.0,color=self.color)
        box["ax"].add_patch(polygon)


        if show_scale_inside and ((not "projection" in box) or box["projection"]!="polar"):
            upperlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
            lowerlimit_string = "0" 
            label="["+lowerlimit_string+"-"+upperlimit_string+"]"
            box["ax"].text(box["left"]+0.05,box["top"]-0.02,label,horizontalalignment="left",verticalalignment="top",fontsize=6*self.fontscale)
       
    def draw_title(self,box):
        if "projection" in box and box["projection"]=="polar": 
            if len(self.label)>0:
                self.label = self.label.replace("\\n","\n")
                rotation = 90 if self.label_rotate else 0
                x,y= polar2cartesian((box["left"],(box["top"]+box["bottom"])/2))
                theta,r = cartesian2polar((x-1,y))

                box["ax"].text(theta,r,
                            self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=10*self.fontscale)
            if self.scale_pos=="left":
                if not self.upside_down:
                    upperlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
                    lowerlimit_string = "0" 
                else:
                    lowerlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
                    upperlimit_string = "0" 
                x,y= polar2cartesian((box["left"],box["top"]))
                theta,r = cartesian2polar((x-0.2,y))
                box["ax"].text(theta,r,
                            upperlimit_string,horizontalalignment="right",verticalalignment="top",fontsize=7*self.fontscale)
                x,y= polar2cartesian((box["left"],box["bottom"]))
                theta,r = cartesian2polar((x-0.2,y))
                box["ax"].text(theta,r,
                            lowerlimit_string,horizontalalignment="right",verticalalignment="bottom",fontsize=7*self.fontscale)
        else:
            if len(self.label)>0:
                self.label = self.label.replace("\\n","\n")
                rotation = 90 if self.label_rotate else 0

                box["ax"].text(box["left"] - 1.0,(box["top"]+box["bottom"])/2,
                            self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            if self.scale_pos=="left":
                upperlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
                lowerlimit_string = "0" 
                if not self.upside_down:
                    box["ax"].text(box["left"] - 0.5,box["top"],
                                upperlimit_string,horizontalalignment="right",verticalalignment="top",fontsize=6*self.fontscale)
                    box["ax"].text(box["left"] - 0.5,box["bottom"],
                                lowerlimit_string,horizontalalignment="right",verticalalignment="bottom",fontsize=6*self.fontscale)
                else:
                    box["ax"].text(box["left"] - 0.5,box["top"],
                                lowerlimit_string,horizontalalignment="right",verticalalignment="top",fontsize=6*self.fontscale)
                    box["ax"].text(box["left"] - 0.5,box["bottom"],
                                upperlimit_string,horizontalalignment="right",verticalalignment="bottom",fontsize=6*self.fontscale)
    
    def compute_max_regions(self,regions,bins):
        m=0
        for i in range(len(regions)):
            reg = regions[i][0]
            region = correct_region_chr(reg,self.bw.chroms())
            values_binned=self.bw.stats(region.chr,region.start,region.end,nBins=bins[i],exact=True,type="mean")
            values_binned = [x if x is not None else 0 for x in values_binned]
            m=max(m,np.max(values_binned)*1.1)
        return m