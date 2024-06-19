import os
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from figeno.utils import KnownException, correct_region_chr, split_box, draw_bounding_box

class coverage_track:
    def __init__(self,file,n_bins=500,color="gray",scale="auto",scale_max=None,scale_pos="corner",group=None,upside_down=False,label="",label_rotate=False,fontscale=1,
                           vcf=None,SNP_colors="auto",exchange_haplotypes=False,bounding_box=False,height=10,margin_above=1.5,**kwargs):
        if file=="" or file is None:
            raise KnownException("Please provide a bam file for the coverage track.")
        if not os.path.isfile(file):
            raise KnownException("The following bam file does not exist (in coverage track): "+file)
        try:
            self.samfile =pysam.AlignmentFile(file, "rb")
        except: 
            raise KnownException("Failed to open bam file (in coverage track): "+str(file))
        if not self.samfile.has_index():
            raise KnownException("Missing index file for bam file (in coverage track): "+str(file)+". Such an index (ending in .bai) is required and can be generated with samtools index.")
        self.filename=file
        self.n_bins = int(n_bins)
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.fontscale=float(fontscale)
        self.scale = scale
        self.scale_max = scale_max
        self.group=group
        self.scale_pos = scale_pos
        self.upside_down= upside_down
        self.vcf = vcf
        self.SNP_colors=SNP_colors
        self.exchange_haplotypes = exchange_haplotypes
        self.bounding_box=bounding_box
        self.height = float(height)
        self.margin_above= float(margin_above)
        self.kwargs=kwargs

        self.atleast_one_region_has_coverage=False
        

    def draw(self, regions, box ,hmargin,warnings=[]):
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            show_scale_inside = self.scale_pos=="corner all" or (self.scale_pos=="corner" and i==0)
            if i<len(self.scale_max): scale_max_region=self.scale_max[i]
            else: 
                scale_max_region=self.scale_max[0]
                if len(self.scale_max)>1: warnings.append("You provided only "+str(len(self.scale_max))+" values for scale_max, even though "+str(len(regions))+" were used. "\
                                                            "The first scale_max parameter will be used for all regions for which no scale_max value was provided.")
            self.draw_region(regions[i][0],boxes[i],scale_max=scale_max_region,show_scale_inside=show_scale_inside)
        self.draw_title(box)

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the coverage track because it is not one of the accepted parameters.")

        if not self.atleast_one_region_has_coverage: warnings.append("The coverage was null for all displayed regions in the bam file "+self.filename+".")

    def draw_region(self,region,box,scale_max,show_scale_inside):
        region = correct_region_chr(region,self.samfile.references)
        if self.bounding_box: draw_bounding_box(box)
        coverage=np.sum(self.samfile.count_coverage(region.chr,region.start,region.end),axis=0)
        coverage=np.nan_to_num(coverage,0)
        n_bins = min(self.n_bins,len(coverage))
        n_bases_per_bin = (region.end-region.start) / n_bins
        coverage_bin = [np.mean(coverage[int(i*n_bases_per_bin):int((i+1)*n_bases_per_bin)]) for i in range(n_bins)]
        
        max_coverage = np.max(coverage_bin) * 1.1
        if max_coverage>=1: self.atleast_one_region_has_coverage=True
        rect_width = (box["right"] - box["left"]) / len(coverage_bin)


        # Draw rectangles for each bin
        if not self.upside_down:
            vertices=[(box["right"],box["bottom"]),(box["left"],box["bottom"])]
        else:
            vertices=[(box["right"],box["top"]),(box["left"],box["top"])]
        for i in range(len(coverage_bin)):
            rect_left = box["left"] + i*rect_width

            if region.orientation=="+":
                rect_height = coverage_bin[i]/scale_max * (box["top"]-box["bottom"])
            else:
                rect_height = coverage_bin[len(coverage_bin)-1-i]/scale_max * (box["top"]-box["bottom"])
            #rect = patches.Rectangle((rect_left,box["bottom"]),
            #                        rect_width,rect_height,color=self.color,lw=0.2,zorder=1)  
            #box["ax"].add_patch(rect)
            if not self.upside_down:
                vertices.append((rect_left,box["bottom"]+rect_height))
                vertices.append((rect_left+rect_width,box["bottom"]+rect_height))
            else:
                vertices.append((rect_left,box["top"]-rect_height))
                vertices.append((rect_left+rect_width,box["top"]-rect_height))
            
            #if i in bin2vaf: # if there is a SNP in the bin
            #    rect = patches.Rectangle((rect_left,box["bottom"]),
            #                        rect_width,rect_height*bin2vaf[i][0],color=self.SNP_colors[0],lw=0.2,zorder=1.1)  
            #    box["ax"].add_patch(rect)
            #   rect = patches.Rectangle((rect_left,box["bottom"] + rect_height*bin2vaf[i][0]),
            #                        rect_width,rect_height*bin2vaf[i][1],color=self.SNP_colors[1],lw=0.2,zorder=1.1)  
            #    box["ax"].add_patch(rect)
        polygon = patches.Polygon(vertices,lw=0,zorder=1,color=self.color)
        box["ax"].add_patch(polygon)

        if show_scale_inside:
            upperlimit_string = "{:d}".format(int(scale_max))
            lowerlimit_string = "0"
            label="["+lowerlimit_string+"-"+upperlimit_string+"]"
            box["ax"].text(box["left"]+0.05,box["top"]-0.02,label,horizontalalignment="left",verticalalignment="top",fontsize=6*self.fontscale)
                
            
    def draw_title(self,box):
        if len(self.label)>0:
            self.label = self.label.replace("\\n","\n")
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 1.0,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
        if self.scale_pos=="left":
            if not self.upside_down:
                upperlimit_string = "{:d}".format(int(self.scale_max[0]))
                lowerlimit_string = "0" 
            else:
                lowerlimit_string = "{:d}".format(int(self.scale_max[0]))
                upperlimit_string = "0" 
            box["ax"].text(box["left"] - 0.5,box["top"],
                        upperlimit_string,horizontalalignment="right",verticalalignment="top",fontsize=6*self.fontscale)
            box["ax"].text(box["left"] - 0.5,box["bottom"],
                        lowerlimit_string,horizontalalignment="right",verticalalignment="bottom",fontsize=6*self.fontscale)


    def compute_max(self,regions,per_region=False):
        l=[]
        max_cov=1
        for i in range(len(regions)):
            if per_region: max_cov=1
            reg = regions[i][0]
            region = correct_region_chr(reg,self.samfile.references)
            coverage=np.sum(self.samfile.count_coverage(region.chr,region.start,region.end),axis=0)
            coverage=np.nan_to_num(coverage,nan=0)
            n_bins = min(self.n_bins,len(coverage))
            n_bases_per_bin = (region.end-region.start) // n_bins
            coverage_bin = [np.mean(coverage[i*n_bases_per_bin:(i+1)*n_bases_per_bin]) for i in range(n_bins)]
            if (region.end-region.start) % n_bins > 0:
                new_value = np.mean(coverage[n_bases_per_bin*(n_bins-1):])
                if new_value==new_value:coverage_bin.append(new_value)
                coverage_bin+= []
            max_cov=max(max_cov,np.max(coverage_bin) * 1.1)
            if max_cov<1 or max_cov!=max_cov: max_cov=1
            if per_region: l.append(round(max_cov))
        if not per_region: l.append(round(max_cov))
        return l


