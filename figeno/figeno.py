import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import json

import importlib_resources as resources
import figeno.data

from figeno.utils import Region, Highlight, draw_highlights, KnownException
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
        if config is None and config_file is None: raise KnownException("Please provide a config_file (.json) or a config as a python dictionary.")
        if config_file is not None:
            old_config=config # allow config to overwrite some parameters of config_file
            with open(config_file, 'r') as f:
                config = json.load(f)
            if old_config is not None:
                for k in old_config:
                    if not k in config: config[k]=old_config[k]
                    else:
                        for k2 in old_config[k]:
                            config[k][k2]=old_config[k][k2]

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
            if not os.path.isfile(self.cytobands_file): raise KnownException("The cytobands file could not be found: "+str(self.cytobands_file))
            with open(self.cytobands_file,"r") as infile2:
                for line in infile2:
                    if line.startswith("#"): continue
                    linesplit = line.rstrip("\n").split("\t")
                    chr = linesplit[0].lstrip("chr")
                    if self.cytobands_file.endswith(".fai"): end=int(linesplit[1])
                    else: end=int(linesplit[2])
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
                    if not "chr" in s: raise KnownException("Please provide at least the chromosome for each region (and optionally the start and end coordinates, otherwise figeno will assume than the region spans the whole chromosome).")
                    if s["chr"]=="" or s["chr"]=="chr" or s["chr"] is None: raise KnownException("Please provide at least the chromosome for each region (and optionally the start and end coordinates, otherwise figeno will assume than the region spans the whole chromosome).")
                    chr = str(s["chr"]).lstrip("chr").lstrip("Chr")
                    orientation = s["orientation"] if "orientation" in s else "+"

                    if "start" in s and (s["start"] is not None) and (s["start"]!=""): 
                        start = s["start"]
                        if isinstance(start,str):
                            start=start.replace(",","")
                            if start.isdigit(): start=int(start)
                            else: raise KnownException("Region start should be an integer: "+start)
                        elif not isinstance(start,int): raise KnownException("Region start should be an integer: "+str(start))
                    else: start = 0

                    if "end" in s and (s["end"] is not None) and (s["end"]!=""): 
                        end = s["end"]
                        if isinstance(end,str):
                            end = end.replace(",","")
                            if end.isdigit(): end=int(end)
                            else: raise KnownException("Region end should be an integer: "+end)
                        elif not isinstance(end,int): raise KnownException("Region end should be an integer: "+str(start))
                    else: 
                        if not chr in self.chr_lengths: 
                            raise KnownException("Could not find length of chromosome \""+str(chr)+"\". If the end of a region is not provided, a cytobands_file must be provided (or use a non-custom reference).")
                        end = self.chr_lengths[chr]
                    if end < start: 
                        start,end = end,start
                        orientation="-"

                    if "color" in s:
                        if matplotlib.colors.is_color_like(s["color"]): color=s["color"]
                        else: raise KnownException("The following color for a region could not be interpreted: "+str(s["color"])+". The color must either be a hexadecimal value, e.g. #FF0000, or the name of the color, e.g. red.")
                    else: color= "#000000"
                    self.add_region(chr=s["chr"].lstrip("chr"),start=start,end=end,orientation=orientation,color=color) 

            if len(self.regions)==0:
                raise KnownException("Please provide at least one region (defined by chr, start and end).")

            # Highlights
            if "highlights" in config:
                for hl in config["highlights"]:
                    if not "chr" in hl: raise KnownException("You did not specify the chromosome for a highlight.")
                    chr = hl["chr"]
                    if chr=="" or chr=="chr" or chr is None: raise KnownException("You did not specify the chromosome for one highlight.")
                    chr = str(chr).lstrip("chr")

                    if "start" in hl and (hl["start"] is not None) and (hl["start"]!=""): 
                        start = hl["start"]
                        if isinstance(start,str):
                            start=start.replace(",","")
                            if start.isdigit(): start=int(start)
                            else: raise KnownException("Highlight start should be an integer: "+start)
                        elif not isinstance(start,int): raise KnownException("Highlight start should be an integer: "+str(start))
                    else: start = 0

                    if "end" in hl and (hl["end"] is not None) and (hl["end"]!=""): 
                        end = hl["end"]
                        if isinstance(end,str):
                            end = end.replace(",","")
                            if end.isdigit(): end=int(end)
                            else: raise KnownException("Highlight end should be an integer: "+end)
                        elif not isinstance(end,int): raise KnownException("Highlight end should be an integer: "+str(start))
                    else: 
                        if not chr in self.chr_lengths: 
                            raise KnownException("Could not find length of chromosome \""+str(chr)+"\". If the end of a highlight is not provided, a cytobands_file must be provided (or use a non-custom reference).")
                        end = self.chr_lengths[chr]
                    if end < start: 
                        start,end = end,start

                    if "color" in hl:
                        if matplotlib.colors.is_color_like(hl["color"]): color=hl["color"]
                        else: raise KnownException("The following color for a region could not be interpreted: "+str(s["color"])+". The color must either be a hexadecimal value, e.g. #FF0000, or the name of the color, e.g. red.")
                    else: color="#eba434"

                    if "opacity" in hl:
                        try:
                            opacity=float(hl["opacity"])
                        except: raise KnownException("The opacity for a highlight should be a float between 0 and 1, but you specified: "+str(hl["opacity"]))
                        if opacity<0 or opacity>1: raise KnownException("The opacity for a highlight should be a float between 0 and 1, but you specified: "+str(hl["opacity"]))
                    else: opacity=0.3
                    self.highlights.append(Highlight(chr,start,end,color,opacity))

            # Tracks
            if "tracks" in config:
                for t in config["tracks"]:
                    if not "type" in t: raise KnownException("Must specify the type for track: " +str(t))
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
                        if "style" in t and t["style"]=="ideogram" and "general" in config:
                            if "reference" in config["general"]:
                                t["reference"]=config["general"]["reference"]
                            if "cytobands_file" in config["general"] and (not config["general"]["cytobands_file"].endswith(".fai")):
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
                        raise KnownException("Unrecognized track type: "+str(t))
            if len(self.tracks_list)==0: raise KnownException("Please provide at least one track.")
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
            total_length=0
            for x in [(reg.end-reg.start) for reg,_ in self.regions]: total_length+=x
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
            total_length=0
            for x in [(reg.end-reg.start) for reg,_ in self.regions]: total_length+=x
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
            
            length1=0
            for x in [(reg.end-reg.start) for reg,_ in self.regions[:last_row1+1]]: length1+=x
            length2=0
            for x in [(reg.end-reg.start) for reg,_ in self.regions[last_row1+1:]]: length2+=x
            longest_length = max(length1,length2)
            
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
        if self.figure_layout=="circular" and isinstance(self.tracks_list[-1],chr_track): self.tracks_list[-1].margin_above=max(0.3,self.tracks_list[-1].margin_above)


    def draw(self,output_config=None,warnings=[]):
        if output_config is not None: self.output = output_config
        if self.output is None or (not "file" in self.output) or (self.output["file"]==""): raise KnownException("Please provide an output file.")
        
        if os.path.dirname(self.output["file"])!="" and not os.path.isdir(os.path.dirname(self.output["file"])): 
            try: os.makedirs(os.path.dirname(self.output["file"]))
            except: raise KnownException("The directory for the output file does not exist ("+os.path.dirname(self.output["file"])+") and could not be created. Please make sure that you have write permissions.")

        if len(self.regions)==0: raise KnownException("Please include at least one region.")

        if not "." in self.output["file"]:
            warnings.append("No file extension was included for the output file. This was automatically set to .png.")
            self.output["file"]+=".png"

        # Show warnings if some margin between chr_axis and copynumber, or between copynumber and sv.
        for i in range(1,len(self.tracks_list)):
            if self.figure_layout!="circular" and isinstance(self.tracks_list[i],chr_track) and isinstance(self.tracks_list[i-1],copynumber_track) and self.tracks_list[i].margin_above>0:
                warnings.append("You might want to set the 'margin_above' parameter of the chr_axis track to 0 when it is below a copynumber track.")
            if isinstance(self.tracks_list[i],copynumber_track) and isinstance(self.tracks_list[i-1],sv_track) and self.tracks_list[i].margin_above>0:
                warnings.append("You might want to set the 'margin_above' parameter of the copynumber track to 0 when it is below a sv track.")
            if self.figure_layout=="circular" and isinstance(self.tracks_list[i],sv_track):
                warnings.append("For circular layouts, the sv track is generally expected to be placed at the top, which is not the case here.")

        
        #Ensure that all regions have similar sizes, otherwise print a warning
        min_size=-1
        max_size=-1
        for r,_ in self.regions:
            size=abs(r.end-r.start)
            if min_size==-1:
                min_size=size
                max_size=size
            else:
                min_size=min(size,min_size)
                max_size=max(size,max_size)
        if max_size>100*min_size:
            warnings.append("You used regions with very different sizes, so the smaller regions may not be visible.")

        # Check if group autoscale was used, in which case compute them.
        self.compute_scales()

        if self.figure_layout=="horizontal": 
            self.draw_horizontal(**self.output,warnings=warnings)
        elif self.figure_layout=="stacked":
            self.draw_stacked(**self.output,warnings=warnings)
        elif self.figure_layout=="symmetrical":
            self.draw_symmetrical(**self.output,warnings=warnings)
        elif self.figure_layout=="circular":
            self.draw_circular(**self.output,warnings=warnings)
        else:
            raise KnownException("Unknown figure layout: "+self.figure_layout+". Please select one from horizontal, symmetrical, circular or stacked.")

    def draw_horizontal(self,file,width=183,dpi=150,transparent=False,warnings=[]):
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
        for i,t in enumerate(self.tracks_list[::-1]):
            if i==len(self.tracks_list)-1 and i!=0 and isinstance(t,chr_track) and t.ticklabels_pos=="below":
                warnings.append("In the chr_axis track, you might want to change the ticklabels_pos parameter to \"above\" instead of \"below\" for better results, since you put this axis at the top of the figure. "\
                                 "The \"below\" option is intended for when the chr_axis is at the bottom of the figure.")
            current_left=0
            #vmargin=1.5
            hmargin=1.5
            box={"ax":ax,"bottom":current_height,"top":current_height+t.height,"left":current_left,"right":current_left+total_width}
            t.draw(regions=self.regions,box=box,hmargin=hmargin,warnings=warnings)
            current_height+=t.height+t.margin_above
        plt.axis('off')
        #fig.savefig(outfile,bbox_inches="tight",pad_inches=0.0,dpi=dpi)
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')
    
    def draw_stacked(self,file,width=183,dpi=150,transparent=False,warnings=[]):
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
            t.draw(regions=self.regions,box=box,hmargin=hmargin,warnings=warnings)
            current_height+=t.height+t.margin_above
        plt.axis('off')
        #fig.savefig(outfile,bbox_inches="tight",pad_inches=0.0,dpi=dpi)
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_symmetrical(self,file,width=183,dpi=150,transparent=False,warnings=[]):
        width=float(width)
        self.update_total_width(width)
        self.update_regions_width_symmetrical()
        self.update_margins()
        if len(self.regions)==1: warnings.append("The symmetrical layout is intended for displaying more than one region.")

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
            t.draw(regions=self.regions,box=box,hmargin=hmargin,warnings=warnings)

            current_height+=t.height
        plt.axis('off')
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def draw_circular(self,file,width=183,dpi=150,transparent=False,warnings=[]):
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
            current_r+=t.margin_above
            box={"ax":ax,"bottom":current_r,"top":current_r+t.height,"left":5*np.pi/2,"right":np.pi/2,"projection":"polar"}
            t.draw(regions=self.regions,box=box,hmargin=hmargin,warnings=warnings)
            current_r+=t.height 
        plt.axis('off')
        ax.set_rmin(0)
        ax.set_rmax(current_r)
        plt.tight_layout(pad=0.1)
        fig.savefig(file,dpi=float(dpi),transparent=transparent)
        plt.cla() 
        plt.clf() 
        plt.close('all')

    def compute_scales(self):
        self.compute_scales_instance(bigwig_track,"bigwig")
        self.compute_scales_instance(coverage_track,"coverage")

    def compute_scales_instance(self,instance=bigwig_track,instance_name="bigwig"):
        group2scales={}
        group2scaletype={}
        for t in self.tracks_list:
            if isinstance(t,instance):
                if t.scale=="custom":
                    if t.scale_max=="": t.scale_max=None
                    if t.scale_max is None: raise KnownException("Please provide the scale_max parameter if you use a custom scale.")
                    if isinstance(t.scale_max,str) and "," in t.scale_max:
                        try:t.scale_max=[float(x) for x in t.scale_max.split(",")]
                        except: raise KnownException("The scale_max parameter in a "+instance_name+" track should be a number (or a list of numbers separated by commas): "+str(self.scale_max))
                    else:
                        try: t.scale_max=[float(t.scale_max)]
                        except: raise KnownException("The scale_max parameter in a "+instance_name+" track should be a number: "+str(self.scale_max))
                elif t.scale=="auto":
                    t.scale_max=t.compute_max(self.regions,per_region=False)
                elif t.scale=="group auto":
                    if not hasattr(t,"group"): raise KnownException("Please provide the group parameter, if you use the group auto scale.")
                    maximum=t.compute_max(self.regions,per_region=False)
                    if t.group in group2scaletype:
                        if group2scaletype[t.group]!="group auto": raise KnownException("Please use the same scale (group auto or group auto per region) for all bigwig tracks of the same group.")
                    else: group2scaletype[t.group]="group auto"
                    if t.group in group2scales: group2scales[t.group] = max(maximum,group2scales[t.group])
                    else: group2scales[t.group] = maximum
                elif t.scale=="auto per region":
                    t.scale_max= t.compute_max(self.regions,per_region=True)
                elif t.scale=="group auto per region":
                    if not hasattr(t,"group"): raise KnownException("Please provide the group parameter, if you use the group auto per region scale.")
                    if t.group in group2scaletype:
                        if group2scaletype[t.group]!="group auto per region": raise KnownException("Please use the same scale (group auto or group auto per region) for all bigwig tracks of the same group.")
                    else: group2scaletype[t.group]="group auto per region"
                    maxima=t.compute_max(self.regions,per_region=True)
                    if t.group in group2scales: group2scales[t.group] = [max(group2scales[t.group][i],maxima[i]) for i in range(len(maxima))]
                    else: group2scales[t.group] = maxima
                else:
                    raise KnownException("Invalid scale for "+instance_name+" track: "+str(t.scale)+". Must be auto, auto per region, group auto, group auto per region, or custom.")

        for t in self.tracks_list:
            if isinstance(t,instance):
                if t.scale=="group auto" or t.scale=="group auto per region":
                    t.scale_max=group2scales[t.group]
                if len(t.scale_max)>1 and t.scale_pos!="none": t.scale_pos="corner all"
                


def figeno_make(config=None,config_file=None,warnings=[]):
    if config is None and config_file is None: raise Exception("ERROR: a config or a config_file is required for figeno_make.")
    tp = tracks_plot(config=config,config_file=config_file)
    tp.draw(warnings=warnings)