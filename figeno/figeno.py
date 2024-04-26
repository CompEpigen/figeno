import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import json
import pathlib

import importlib_resources as resources
import figeno.data

from figeno.utils import Region, Highlight, draw_highlights
from figeno.track_chr import chr_track
from figeno.track_genes import genes_track
from figeno.track_bed import bed_track
from figeno.track_bigwig import bigwig_track
from figeno.track_coverage import coverage_track
from figeno.track_alignments import alignments_track
from figeno.track_basemodfreq import basemodfreq_track
from figeno.track_hic import hic_track
from figeno.track_copynumber import copynumber_track
from figeno.track_sv import sv_track
from figeno.track_ase import ase_track


matplotlib.use("svg")
matplotlib.use("pdf")
matplotlib.use("ps")
matplotlib.use("ps")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["svg.fonttype"] = "none"


class tracks_plot:
    def __init__(self,config=None,config_file=None,reference=None,genes_file=None,cytobands_file=None):
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
        self.cytobands_file=cytobands_file
        if config is not None and "general" in config:
            if "reference" in config["general"]: self.reference=config["general"]["reference"]
            if "genes_file" in config["general"]: self.genes_file=config["general"]["genes_file"]
            if "cytobands_file" in config["general"]: self.cytobands_file = config["general"]["cytobands_file"]
        
        self.chr_lengths={}
        if self.cytobands_file is not None and self.cytobands_file!="":
            with open(self.cytobands_file,"r") as infile2:
                for line in infile2:
                    if line.startswith("#"): continue
                    linesplit = line.rstrip("\n").split("\t")
                    chr = linesplit[0].lstrip("chr")
                    end=int(linesplit[2])
                    if not chr in self.chr_lengths: self.chr_lengths[chr]=end
                    else: self.chr_lengths[chr]=max(end,self.chr_lengths[chr])
        elif self.reference in ["hg19","hg38","mm10"]:
            with resources.as_file(resources.files(figeno.data) / (self.reference+"_cytobands.tsv")) as infile:
                with open(infile,"r") as infile2:
                    for line in infile2:
                        if line.startswith("#"): continue
                        linesplit = line.rstrip("\n").split("\t")
                        chr = linesplit[0].lstrip("chr")
                        end=int(linesplit[2])
                        if not chr in self.chr_lengths: self.chr_lengths[chr]=end
                        else: self.chr_lengths[chr]=max(end,self.chr_lengths[chr])
        
        if config is not None:
            # Regions
            if "regions" in config:
                for s in config["regions"]:
                    if "chr" in s:
                        chr = str(s["chr"]).lstrip("chr")
                        if chr=="": raise Exception("Must specify the chromosome of the region.")
                        orientation = s["orientation"] if "orientation" in s else "+"
                        if "start" in s and (s["start"]is not None) and (s["start"]!=""): 
                            start = s["start"]
                            if isinstance(start,str):
                                start=start.replace(",","")
                                if start.isdigit(): start=int(start)
                                else: raise Exception("Region start should be an integer: "+start)
                        else: start = 0
                        if "end" in s and (s["end"] is not None) and (s["end"]!=""): 
                            end = s["end"]
                            if isinstance(end,str):
                                end = end.replace(",","")
                                if end.isdigit(): end=int(end)
                                else: raise Exception("Region end should be an integer: "+end)
                        else: 
                            if not chr in self.chr_lengths: 
                                raise Exception("Could not find length of chromosome \""+str(chr)+"\". If the end of a region is not provided, a cytobands_file must be provided (or use a non-custom reference).")
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
                    self.highlights.append(Highlight(str(hl["chr"]).lstrip("chr"),start,end,hl["color"],hl["opacity"]))

            # Tracks
            if "tracks" in config:
                for t in config["tracks"]:
                    if not "type" in t: raise Exception("Must specify the type for track: " +str(t))
                    track_type = t.pop("type")
                    if track_type == "alignments":
                        self.tracks_list.append(alignments_track(**t))
                    elif track_type=="coverage":
                        self.tracks_list.append(coverage_track(**t))
                    elif track_type=="bigwig":
                        self.tracks_list.append(bigwig_track(**t))
                    elif track_type=="bed":
                        self.tracks_list.append(bed_track(**t))
                    elif track_type=="basemod_freq":
                        self.tracks_list.append(basemodfreq_track(**t))
                    elif track_type=="hic":
                        self.tracks_list.append(hic_track(**t))
                    elif track_type=="genes":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        self.tracks_list.append(genes_track(**t))
                    elif track_type=="chr_axis":
                        if t["style"]=="ideogram" and "general" in config:
                            if "reference" in config["general"]:
                                t["reference"]=config["general"]["reference"]
                            if "cytobands_file" in config["general"]:
                                t["cytobands_file"] = config["general"]["cytobands_file"]
                        self.tracks_list.append(chr_track(**t))
                    elif track_type=="copynumber":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        t["chr_lengths"] = self.chr_lengths
                        self.tracks_list.append(copynumber_track(**t))
                    elif track_type=="sv":
                        self.tracks_list.append(sv_track(**t))
                    elif track_type=="ase":
                        t["reference"] = self.reference
                        t["genes_file"] = self.genes_file
                        self.tracks_list.append(ase_track(**t))
                    else:
                        print("WARNING: unrecognized track type: "+str(t))
            if "general" in config and "layout" in config["general"]:
                self.figure_layout = config["general"]["layout"]
            else:
                self.figure_layout="horizontal"
            self.output = None
            if "output" in config: self.output=config["output"]


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
    def update_regions_width_stacked(self):
        if self.total_width is not None:
            max_length= np.max([(reg.end-reg.start) for reg,_ in self.regions])
            new_regions=[]
            for reg,_ in self.regions:
                new_regions.append((reg,self.total_width * (reg.end-reg.start) / max_length))
            self.regions = new_regions
    def update_regions_width_symmetrical(self):
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
                if i>0 and self.tracks_list[i].ticklabels_pos=="below" and self.tracks_list[i].margin_above==0:
                    self.tracks_list[i].no_margin=True
                elif i<len(self.tracks_list)-1 and self.tracks_list[i].ticklabels_pos=="above" and self.tracks_list[i+1].margin_above==0:
                    self.tracks_list[i].no_margin=True
            if isinstance(self.tracks_list[i],sv_track) and i<len(self.tracks_list)-1 and isinstance(self.tracks_list[i+1],copynumber_track):
                self.tracks_list[i+1].margin_above=max(self.tracks_list[i+1].margin_above,0.3)


    def draw(self,output_config=None):
        if output_config is not None: self.output = output_config
        if self.output is None or (not "file" in self.output) or (self.output["file"]==""): raise Exception("Must specify the output file")

        if len(self.regions)==0: raise Exception("Must include at least one region.")

        if self.figure_layout=="horizontal": 
            self.draw_horizontal(**self.output)
        elif self.figure_layout=="stacked":
            self.draw_stacked(**self.output)
        elif self.figure_layout=="symmetrical":
            self.draw_symmetrical(**self.output)
        elif self.figure_layout=="circular":
            self.draw_circular(**self.output)
        else:
            raise Exception("Unknown figure layout: "+self.figure_layout)

    def draw_horizontal(self,file,width=183,dpi=150,transparent=False):
        width=float(width)
        self.update_total_width(width)
        self.update_regions_width()
        self.update_margins()


        total_height = 0
        for ind,t in enumerate(self.tracks_list): 
            total_height+=t.height
            if ind>0: total_height+=t.margin_above
        total_width=0
        for _,w in self.regions: total_width+=w

        fig,ax = plt.subplots(figsize=(width/25.4,total_height/25.4)) # convert from mm to inches
        ax.set_xlim((-10,total_width+5))
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
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')
    
    def draw_stacked(self,file,width=183,dpi=150,transparent=False):
        width=float(width)
        self.update_total_width(width)
        self.update_regions_width_stacked()
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
        ax.set_xlim((-10,width+5))
        ax.set_ylim((-4-total_height*(len(self.regions)-1),total_height+4))
        
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
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_symmetrical(self,file,width=183,dpi=150,transparent=False):
        width=float(width)
        self.update_total_width(width)
        self.update_regions_width_symmetrical()
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
        ax.set_xlim((-10,self.total_width*1.03))
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
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_circular(self,file,width=183,dpi=150,transparent=False):
        width=float(width)
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
        fig,ax = plt.subplots(figsize=(width/25.4,width/25.4),subplot_kw={'projection': 'polar'}) # convert from mm to inches

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
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

def figeno_make(config=None,config_file=None):
    if config is None and config_file is None: raise Exception("ERROR: a config or a config_file is required for figeno_make.")
    tp = tracks_plot(config=config,config_file=config_file)
    tp.draw()