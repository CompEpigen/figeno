import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from figeno.utils import correct_region_chr, split_box, draw_bounding_box
from figeno.vcf import read_phased_vcf

class coverage_track:
    def __init__(self,file,n_bins=500,color="gray",scale="auto",scale_max=None,scale_pos="corner",upside_down=False,label="",label_rotate=False,fontscale=1,
                           vcf=None,SNP_colors="auto",exchange_haplotypes=False,bounding_box=False,height=10,margin_above=1.5):
        self.samfile = pysam.AlignmentFile(file, "rb")
        self.n_bins = int(n_bins)
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.fontscale=fontscale
        self.scale = scale
        self.scale_max = scale_max
        if self.scale_max is not None: self.scale_max=float(self.scale_max)
        self.scale_pos = scale_pos
        if scale=="auto per region": self.scale_pos = "corner all" # If each region has its own scale, we cannot use one global label for the whole track
        self.upside_down= upside_down
        self.vcf = vcf
        self.SNP_colors=SNP_colors
        self.exchange_haplotypes = exchange_haplotypes
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

    def draw(self, regions, box ,hmargin):
        if self.scale=="auto": self.scale_max = self.compute_max_regions(regions)
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            show_scale_inside = self.scale_pos=="corner all" or (self.scale_pos=="corner" and i==0)
            self.draw_region(regions[i][0],boxes[i],show_scale_inside=show_scale_inside)
        self.draw_title(box)

    def draw_region(self,region,box,show_scale_inside):
        region = correct_region_chr(region,self.samfile.references)
        if self.bounding_box: draw_bounding_box(box)
        #if region.end-region.start<5000000:
        coverage=np.sum(self.samfile.count_coverage(region.chr,region.start,region.end),axis=0)
        coverage=np.nan_to_num(coverage,0)
        n_bins = min(self.n_bins,len(coverage))
        n_bases_per_bin = (region.end-region.start) / n_bins
        coverage_bin = [np.mean(coverage[int(i*n_bases_per_bin):int((i+1)*n_bases_per_bin)]) for i in range(n_bins)]

        #else:
        #    n_bins = self.n_bins
        #    coverage_bin=[]
        #    for i in range(n_bins):
        #        start = region.start + i* (region.end-region.start)/n_bins
        #        end = start + (region.end-region.start)/n_bins
        #        cov = self.samfile.count(region.chr,start,end)
        #        coverage_bin.append(cov)


        
        max_coverage = np.max(coverage_bin) * 1.1
        if max_coverage<1 or max_coverage!=max_coverage: max_coverage=1
        if self.scale=="auto per region": self.scale_max=max_coverage
        rect_width = (box["right"] - box["left"]) / len(coverage_bin)

        #SNPs=[]
        #if self.vcf is not None:
        #    SNPs = read_phased_vcf(self.vcf,region.chr,region.start,region.end)
        #if self.SNP_colors=="auto": self.SNP_colors = ["#1a7242","#e67e22"]
        #bin2vaf={}
        #nuc_to_index = {"A":0,"C":1,"G":2,"T":3}
        #for i in range(len(SNPs)):
        #    pos,nuc1,nuc2 = SNPs[i]
        #    if nuc1 in nuc_to_index and nuc2 in nuc_to_index:
        #        if self.exchange_haplotypes: nuc1, nuc2 = nuc2, nuc1
        #        coverage_SNP=self.samfile.count_coverage(region.chr,pos,pos+1)
        #        cov1=coverage_SNP[nuc_to_index[nuc1]][0]
        #        cov2=coverage_SNP[nuc_to_index[nuc2]][0]
        #        total_cov = np.sum(coverage_SNP) 
        #       if total_cov>=5:
        #           bin2vaf[(pos-region.start)//n_bases_per_bin] = (cov1/total_cov,cov2/total_cov)

        # Draw rectangles for each bin
        if not self.upside_down:
            vertices=[(box["right"],box["bottom"]),(box["left"],box["bottom"])]
        else:
            vertices=[(box["right"],box["top"]),(box["left"],box["top"])]
        for i in range(len(coverage_bin)):
            rect_left = box["left"] + i*rect_width

            if region.orientation=="+":
                rect_height = coverage_bin[i]/self.scale_max * (box["top"]-box["bottom"])
            else:
                rect_height = coverage_bin[len(coverage_bin)-1-i]/self.scale_max * (box["top"]-box["bottom"])
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
            upperlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
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
                upperlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
                lowerlimit_string = "0" 
            else:
                lowerlimit_string = "{:.1f}".format(self.scale_max) if self.scale_max>=1 else "{:.2f}".format(self.scale_max)
                upperlimit_string = "0" 
            box["ax"].text(box["left"] - 0.5,box["top"],
                        upperlimit_string,horizontalalignment="right",verticalalignment="top",fontsize=6*self.fontscale)
            box["ax"].text(box["left"] - 0.5,box["bottom"],
                        lowerlimit_string,horizontalalignment="right",verticalalignment="bottom",fontsize=6*self.fontscale)
        

    def compute_max_regions(self,regions):
        max_cov=0
        for region in regions:
            reg = region[0]
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
        return round(max_cov)






