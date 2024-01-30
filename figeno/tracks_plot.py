import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import json
import pathlib

import importlib_resources as resources
import figeno.data

from figeno.tracks_utils import Region, Highlight, draw_highlights
from figeno.tracks_chr import chr_track
from figeno.tracks_genes import genes_track
from figeno.tracks_bed import bed_track
from figeno.tracks_bigwig import bigwig_track
from figeno.tracks_coverage import coverage_track
from figeno.tracks_alignments import alignments_track
from figeno.tracks_metfreq import metfreq_track
from figeno.tracks_HiC import hic_track
from figeno.tracks_copynumber import copynumber_track
from figeno.tracks_sv import sv_track
from figeno.tracks_ase import ase_track


matplotlib.use("svg")
matplotlib.use("pdf")
matplotlib.use("ps")
matplotlib.use("ps")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["svg.fonttype"] = "none"


class tracks_plot:
    def __init__(self,config=None,config_file=None,reference=None,genes_file=None,chrarms_file=None):
        self.total_height=0
        self.width=0
        self.tracks_list=[]
        self.regions=[]
        self.highlights = []
        self.total_width=None
        self.reference="custom"
        self.genes_file=""
        if config_file is not None:
            with open(config_file, 'r') as f:
                config = json.load(f)

        # Reference
        self.reference=reference
        self.genes_file=genes_file
        self.chrarms_file=chrarms_file
        if config is not None:
            if "reference" in config["general"]: self.reference=config["general"]["reference"]
            if "genes_file" in config["general"]: self.genes_file=config["general"]["genes_file"]
            if "chrarms_file" in config["general"]: self.chrarms_file = config["general"]["chrarms_file"]
        
        self.chr_lengths={}
        self.centromeres={}
        if self.chrarms_file is not None:
            with open(self.chrarms_file,"r") as infile2:
                for line in infile2:
                    linesplit = line.rstrip("\n").split(" ")
                    chr = linesplit[1].lstrip("chr")
                    if linesplit[0].endswith("q"):
                        self.chr_lengths[chr] = int(linesplit[3])
                    else:
                        self.centromeres[chr] = int(linesplit[4])
        elif self.reference in ["hg19","hg38"]:
            with resources.as_file(resources.files(figeno.data) / (self.reference+"_chrarms.txt")) as infile:
                with open(infile,"r") as infile2:
                    for line in infile2:
                        linesplit = line.rstrip("\n").split(" ")
                        chr = linesplit[1].lstrip("chr")
                        if linesplit[0].endswith("q"):
                            self.chr_lengths[chr] = int(linesplit[3])
                        else:
                            self.centromeres[chr] = int(linesplit[4])
        
        if config is not None:
            # Regions
            if "regions" in config:
                for s in config["regions"]:
                    if "chr" in s:
                        chr = s["chr"].lstrip("chr")
                        orientation = s["orientation"] if "orientation" in s else "+"
                        if "start" in s: start = s["start"]
                        else: start = 0
                        if "end" in s: end = s["end"]
                        else: 
                            if not chr in self.chr_lengths: 
                                raise Exception("If the end of a region is not provided, a chrarms_file must be provided (or a non-custom reference).")
                            end = self.chr_lengths[chr]
                        if end < start: 
                            start,end = end,start
                            orientation="-"
                        color=s["color"] if "color" in s else "#000000"
                        self.add_region(chr=s["chr"].lstrip("chr"),start=start,end=end,orientation=orientation,color=color) 

            # Highlights
            if "highlights" in config:
                for hl in config["highlights"]:
                    start,end = hl["start"],hl["end"]
                    if start>end: start,end = end,start
                    self.highlights.append(Highlight(hl["chr"],start,end,hl["color"],hl["alpha"]))

            # Tracks
            if "tracks" in config:
                for t in config["tracks"]:
                    if not "type" in t: raise Exception("Must specify the type for track: " +str(t))
                    track_type = t.pop("type")
                    if track_type == "alignments":
                        self.add_alignments_track(**t)
                    elif track_type=="coverage":
                        self.add_coverage_track(**t)
                    elif track_type=="bigwig":
                        self.add_bigwig_track(**t)
                    elif track_type=="bed":
                        self.add_bed_track(**t)
                    elif track_type=="Met freq" or track_type=="methylation_freq":
                        self.add_metfreq_track(**t)
                    elif track_type=="HiC":
                        self.add_hic_track(**t)
                    elif track_type=="genes":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        self.add_genes_track(**t)
                    elif track_type=="chr axis":
                        self.add_chraxis_track(**t)
                    elif track_type=="copynumber":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        t["chr_lengths"] = self.chr_lengths
                        self.add_copynumber_track(**t)
                    elif track_type=="SV":
                        self.add_sv_track(**t)
                    elif track_type=="ase":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        self.add_ase_track(**t)
                    else:
                        pass
                        #raise Exception("Unrecognized track type: " + str(track_type)+". Please use one of: alignments, coverage, bigwig, bed, methylation_freq, genes or chr_axis.")


            self.figure_type = config["general"]["figure_type"]
            self.output = None
            if "output" in config: self.output=config["output"]



    def add_alignments_track(self,**args):
        self.tracks_list.append(alignments_track(**args))
    def add_coverage_track(self,**args):
        self.tracks_list.append(coverage_track(**args))
    def add_bigwig_track(self,**args):
        self.tracks_list.append(bigwig_track(**args))
    def add_bed_track(self,**args):
        self.tracks_list.append(bed_track(**args))
    def add_metfreq_track(self,**args):
        self.tracks_list.append(metfreq_track(**args))
    def add_hic_track(self,**args):
        self.tracks_list.append(hic_track(**args))
    def add_genes_track(self,**args):
        self.tracks_list.append(genes_track(**args))
    def add_chraxis_track(self,**args):
        self.tracks_list.append(chr_track(**args))

    def add_copynumber_track(self,**args):
        self.tracks_list.append(copynumber_track(**args))
    def add_sv_track(self,**args):
        self.tracks_list.append(sv_track(**args))
    def add_ase_track(self,**args):
        self.tracks_list.append(ase_track(**args))

    def add_region(self,chr,start,end,orientation="+",color="#000000",width=20):
        self.regions.append((Region(chr,start,end,orientation,color),width))
    def update_total_width(self,total_width):
        self.total_width = total_width
    def update_regions_width(self):
        if self.total_width is not None:
            total_length= np.sum([(reg.end-reg.start) for reg,_ in self.regions])
            new_regions=[]
            for reg,_ in self.regions:
                new_regions.append((reg,self.total_width * (reg.end-reg.start) / total_length))
            self.regions = new_regions
    def update_regions_width_multilines(self):
        if self.total_width is not None:
            max_length= np.max([(reg.end-reg.start) for reg,_ in self.regions])
            new_regions=[]
            for reg,_ in self.regions:
                new_regions.append((reg,self.total_width * (reg.end-reg.start) / max_length))
            self.regions = new_regions
    def update_regions_width_tworows(self):
        if self.total_width is not None:
            total_length= np.sum([(reg.end-reg.start) for reg,_ in self.regions])
            # Assign regions to rows
            i=0
            current_length=0
            while current_length<total_length/2 and i<len(self.regions):
                current_length+=(self.regions[i][0].end - self.regions[i][0].start )
                i+=1
            length2 = current_length - (self.regions[i-1][0].end - self.regions[i-1][0].start )
            if abs(length2-total_length/2) < abs(current_length-total_length/2):
                last_row1= i-2
            else:
                last_row1 = i-1
            if last_row1<0: last_row1=0
            
            longest_length = max( np.sum([(reg.end-reg.start) for reg,_ in self.regions[:last_row1+1]]) , 
                                  np.sum([(reg.end-reg.start) for reg,_ in self.regions[last_row1+1:]]))
            
            new_regions = []
            for i in range(len(self.regions)):
                width = self.total_width * (self.regions[i][0].end-self.regions[i][0].start) / longest_length
                row = 0 if i<=last_row1 else 1
                new_regions.append((self.regions[i][0],width,row))
            self.regions = new_regions

    def update_margins(self):
        for i in range(len(self.tracks_list)):
            if isinstance(self.tracks_list[i],chr_track):
                if i>0 and self.tracks_list[i].ticks_pos=="below" and self.tracks_list[i].margin_above==0:
                    self.tracks_list[i].no_margin=True
                elif i<len(self.tracks_list)-1 and self.tracks_list[i].ticks_pos=="above" and self.tracks_list[i+1].margin_above==0:
                    self.tracks_list[i].no_margin=True
            if isinstance(self.tracks_list[i],sv_track) and isinstance(self.tracks_list[i+1],copynumber_track):
                self.tracks_list[i+1].margin_above=max(self.tracks_list[i+1].margin_above,0.3)


    def draw(self,output_config=None):
        if output_config is not None: self.output = output_config
        if self.output is None or (not "file" in self.output) or (self.output["file"]==""): raise Exception("Must specify the output file")

        if self.figure_type=="1 row": 
            self.draw_onerow(**self.output)
        elif self.figure_type=="multilines":
            self.draw_multilines(**self.output)
        elif self.figure_type=="2 rows":
            self.draw_tworows(**self.output)
        elif self.figure_type=="circle":
            self.draw_circle(**self.output)
        else:
            raise Exception("Unknown figure type: "+self.figure_type)

    def draw_onerow(self,file,width=183,dpi=150):
        self.update_total_width(width)
        self.update_regions_width()
        self.update_margins()

        #for i in range(0,len(self.tracks_list)-1):
        #    if isinstance(self.tracks_list[i],sv_track) and isinstance(self.tracks_list[i+1],copynumber_track):
        #        self.tracks_list[i+1].margin_above=max(self.tracks_list[i+1].margin_above,0.3)
        total_height = 0
        for ind,t in enumerate(self.tracks_list): 
            total_height+=t.height
            if ind>0: total_height+=t.margin_above
        total_width=0
        for _,w in self.regions: total_width+=w
        #fig,ax = plt.subplots(figsize=(total_width,total_height))
        fig,ax = plt.subplots(figsize=(width/25.4,total_height/25.4)) # convert from mm to inches
        ax.set_xlim((-5,total_width+5))
        ax.set_ylim((-4,total_height+4))

        # Highlights
        highlight_bottom=0
        highlight_top = total_height
        if isinstance(self.tracks_list[0],chr_track): 
            highlight_top+=self.tracks_list[0].get_highlights_offsets()[0]
        if isinstance(self.tracks_list[-1],chr_track): 
            highlight_bottom+=self.tracks_list[-1].get_highlights_offsets()[1]
        draw_highlights(box={"ax":ax,"bottom":highlight_bottom,"top":highlight_top,"left":0,"right":total_width},
                        highlights=self.highlights,regions=self.regions,hmargin=1.5)
        
        
        # Tracks
        current_height=0
        for t in self.tracks_list[::-1]:
            current_left=0
            #vmargin=1.5
            hmargin=1.5
            box={"ax":ax,"bottom":current_height,"top":current_height+t.height,"left":current_left,"right":current_left+total_width}
            t.draw(regions=self.regions,box=box,hmargin=hmargin)
            current_height+=t.height+t.margin_above
        plt.axis('off')
        #fig.savefig(outfile,bbox_inches="tight",pad_inches=0.0,dpi=dpi)
        plt.tight_layout(pad=0)
        fig.savefig(file,dpi=dpi)
        plt.cla() 
        plt.clf() 
        plt.close('all')
    
    def draw_multilines(self,file,width=183,dpi=150):
        self.update_total_width(width)
        self.update_regions_width_multilines()
        self.update_margins()

        #for i in range(0,len(self.tracks_list)-1):
        #    if isinstance(self.tracks_list[i],sv_track) and isinstance(self.tracks_list[i+1],copynumber_track):
        #        self.tracks_list[i+1].margin_above=max(self.tracks_list[i+1].margin_above,0.3)
        total_height = 0
        for ind,t in enumerate(self.tracks_list): 
            total_height+=t.height
            total_height+=t.margin_above
        total_width=0
        for _,w in self.regions: total_width+=w
        #fig,ax = plt.subplots(figsize=(total_width,total_height))
        fig,ax = plt.subplots(figsize=(width/25.4,total_height*len(self.regions)/25.4)) # convert from mm to inches
        ax.set_xlim((-5,width+5))
        ax.set_ylim((-4-total_height*(len(self.regions)-1),total_height+4))

        # Highlights
        #highlight_bottom=0
        #highlight_top = total_height
        #if isinstance(self.tracks_list[0],chr_track): 
        #    highlight_top+=self.tracks_list[0].get_highlights_offsets()[0]
        #if isinstance(self.tracks_list[-1],chr_track): 
        #    highlight_bottom+=self.tracks_list[-1].get_highlights_offsets()[1]
        #draw_highlights(box={"ax":ax,"bottom":highlight_bottom,"top":highlight_top,"left":0,"right":total_width},
        #               highlights=self.highlights,regions=self.regions,hmargin=1.5)
        
        
        # Tracks
        current_height=0
        for t in self.tracks_list[::-1]:
            current_left=0
            #vmargin=1.5
            hmargin=1.5
            box={"ax":ax,"bottom":current_height,"top":current_height+t.height,"left":current_left,"right":current_left+width,"total_height":total_height}
            t.draw(regions=self.regions,box=box,hmargin=hmargin)
            current_height+=t.height+t.margin_above
        plt.axis('off')
        #fig.savefig(outfile,bbox_inches="tight",pad_inches=0.0,dpi=dpi)
        plt.tight_layout(pad=0)
        fig.savefig(file,dpi=dpi)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_tworows(self,file,width=183,dpi=150):
        self.update_total_width(width)
        self.update_regions_width_tworows()
        self.update_margins()

        #for i in range(0,len(self.tracks_list)-1):
        #    if isinstance(self.tracks_list[i],sv_track) and isinstance(self.tracks_list[i+1],copynumber_track):
        #        self.tracks_list[i+1].margin_above=max(self.tracks_list[i+1].margin_above,0.3)
        total_height = 0
        for t in self.tracks_list: total_height+=t.height + t.margin_above
        #total_width=0
        #for _,w,r in self.regions: total_width+=w
        #fig,ax = plt.subplots(figsize=(total_width,total_height))
        fig,ax = plt.subplots(figsize=(width/25.4,total_height/25.4*2)) # convert from mm to inches
        ax.set_xlim((-self.total_width*0.03,self.total_width*1.03))
        ax.set_ylim((-total_height*1.02,total_height*1.02))

        # Highlights
        highlight_bottom=-total_height
        if isinstance(self.tracks_list[-1],chr_track): 
            highlight_bottom+=self.tracks_list[-1].get_highlights_offsets()[1]
        highlight_top = -highlight_bottom
        middle_pos = self.tracks_list[0].margin_above if (not isinstance(self.tracks_list[0],sv_track)) else 0
        draw_highlights(box={"ax":ax,"bottom":middle_pos,"top":highlight_top,"left":0,"right":self.total_width,"top2":-middle_pos,"bottom2":highlight_bottom},
                        highlights=self.highlights,regions=self.regions,hmargin=1.5)

        # Tracks
        current_height=0
        for index,t in enumerate(self.tracks_list):
            if index>0 or (not isinstance(t,sv_track)): current_height+=t.margin_above
            current_left=0
            #vmargin=1.5
            hmargin=1.5
            box={"ax":ax,"bottom":current_height,"top":current_height+t.height,"left":current_left,"right":current_left+self.total_width,
                 "bottom2":-current_height-t.height,"top2":-current_height}
            t.draw(regions=self.regions,box=box,hmargin=hmargin)

            current_height+=t.height
        plt.axis('off')
        plt.tight_layout(pad=0)
        fig.savefig(file,dpi=dpi)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_circle(self,file,width=183,dpi=150):
        self.update_total_width(width)
        self.update_regions_width()
        self.update_margins()
        for i in range(len(self.tracks_list)):
            if isinstance(self.tracks_list[i],sv_track):
                self.tracks_list[i].height=0
                self.tracks_list[i].margin_above=0
                self.tracks_list[i].bounding_box=False
        total_height = 0
        for t in self.tracks_list: total_height+=t.height + t.margin_above
        total_width=0
        for _,w in self.regions: total_width+=w
        #fig,ax = plt.subplots(figsize=(total_width,total_height))
        fig,ax = plt.subplots(figsize=(10,10),subplot_kw={'projection': 'polar'}) # convert from mm to inches

        current_r=0
        for t in self.tracks_list: current_r+=t.height+t.margin_above
        current_r = max(current_r* 0.7,40)
        current_r=50

        # Highlights
        highlight_bottom=current_r
        highlight_top = current_r+total_height 
        #if isinstance(self.tracks_list[0],chr_track): 
        #    highlight_top+=self.tracks_list[0].height*0.4
        if isinstance(self.tracks_list[-1],chr_track): highlight_top-= 0.85*self.tracks_list[-1].height
        draw_highlights(box={"ax":ax,"bottom":highlight_bottom,"top":highlight_top,"left":5*np.pi/2,"right":np.pi/2,"projection":"polar"},
                        highlights=self.highlights,regions=self.regions,hmargin=1.5)
        

        # Tracks
        for t in self.tracks_list:
            hmargin=1.5
            box={"ax":ax,"bottom":current_r,"top":current_r+t.height,"left":5*np.pi/2,"right":np.pi/2,"projection":"polar"}
            t.draw(regions=self.regions,box=box,hmargin=hmargin)
            current_r+=t.height+t.margin_above
        plt.axis('off')
        ax.set_rmin(0)
        ax.set_rmax(current_r)
        plt.tight_layout(pad=0)
        fig.savefig(file,dpi=dpi)
        plt.cla() 
        plt.clf() 
        plt.close('all')

