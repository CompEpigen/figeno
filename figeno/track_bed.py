import os
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

from figeno.utils import split_box, draw_bounding_box, KnownException
from figeno.track_genes import reorder_balanced

from collections import namedtuple
BedRect = namedtuple('BedRect', 'chr start end name strand color')


class bed_track:
    def __init__(self,file,color="#051650",label="",label_rotate=False,show_names=False,show_strand=False,use_file_colors=False,collapsed=False,fontscale=1,bounding_box=False,height=3,margin_above=1.5,**kwargs):
        self.file = file
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.show_names = show_names
        self.show_strand=show_strand
        self.use_file_colors=use_file_colors
        self.collapsed=collapsed
        self.fontscale=float(fontscale)
        self.bounding_box=bounding_box
        self.height = float(height)
        self.margin_above= float(margin_above)
        self.kwargs=kwargs

    def draw(self, regions, box ,hmargin,warnings):
        # Check that the bed file exists
        if self.file=="" or self.file is None:
            raise KnownException("Please provide a file for the bed track.")
        self.file=str(self.file)
        if not os.path.isfile(self.file):
            raise KnownException("The following bed file could not be found: "+self.file)
        
        # Make sure that the file can be read
        try:
            if self.file.endswith(".gz"):
                with gzip.open(self.file,"rt") as infile:
                    line=infile.readline()
            else:
                with open(self.file,"r") as infile:
                    line=infile.readline()
        except:
            raise KnownException("The following bed file could not be read: "+self.file)
        
        boxes = split_box(box,regions,hmargin)
        lines_regions=self.read_records_regions(regions)
        if (box["top"]-box["bottom"])/len(lines_regions[0])<3.5 and self.show_names:
            warnings.append("The bed track is very dense ({} lines for a height of {} mm). You might want to increase the height for this track or to not show names.".format(len(lines_regions[0]),self.height))
        for i in range(len(regions)):
            self.draw_region(regions[i][0],boxes[i],lines_regions[i],warnings)
        self.draw_title(box)

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the bed track because it is not one of the accepted parameters.")

    def draw_region(self,region,box,lines,warnings):
        if self.bounding_box: draw_bounding_box(box)
        
        if self.show_names: 
            rect_height = (box["top"]-box["bottom"]) * 0.4 / len(lines)
            #rect_bottom = box["bottom"] +  (box["top"]-box["bottom"]) * 0.08
        else:
            rect_height = (box["top"]-box["bottom"]) * 0.7 / len(lines)
            #rect_bottom = box["bottom"] +  (box["top"]-box["bottom"]) * 0.15

        # Params for strand
        self.arrow_gapwidth=3 # gap between arrows
        self.arrow_height=min(0.8*rect_height,1.4)
        self.arrow_width=min(0.5,1.5*self.arrow_height)

        # Params for text
        fontsize = 8*self.fontscale if len(lines)<=1 else 6*self.fontscale

        for line_index,line in enumerate(lines):
            line_bottom = box["bottom"] + (box["top"]-box["bottom"])/len(lines) * line_index + (box["top"]-box["bottom"]) * 0.15 / len(lines)
            for bedrect in line:
                if region.orientation=="+":
                    converted_start = box["left"] + (bedrect.start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                    converted_end = box["left"] + (bedrect.end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                else:
                    converted_start = box["right"] - (bedrect.end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                    converted_end = box["right"] - (bedrect.start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                converted_start = min(box["right"],max(box["left"],converted_start))
                converted_end = min(box["right"],max(box["left"],converted_end))
                rect_width=converted_end-converted_start

                rect = patches.Rectangle((converted_start,line_bottom),rect_width,rect_height,color=bedrect.color,lw=0.1)
                box["ax"].add_patch(rect)

                if self.show_names and bedrect.name!="":
                    box["ax"].text((converted_start+converted_end)/2,(line_bottom+min(rect_height*1.1,rect_height+1)),bedrect.name,
                               horizontalalignment="center",verticalalignment="bottom",fontsize=fontsize)
                if self.show_strand and bedrect.strand in ["+","-"]:
                    arrow_center=line_bottom+rect_height/2
                    arrow_bottom=arrow_center-self.arrow_height/2
                    arrow_top=arrow_center+self.arrow_height/2
                    if rect_width>=1.2*self.arrow_width: 
                        # Number of arrows that can fit in the rectangle
                        n_arrows= 1+ int((rect_width-3*self.arrow_width)/(self.arrow_gapwidth+self.arrow_width))

                        if rect_width<3*self.arrow_width:
                            current_left_pos=(converted_start+converted_end)/2 - self.arrow_width/2
                        else:
                            arrows_width= 3*self.arrow_width + (n_arrows-1)*(self.arrow_gapwidth+self.arrow_width)
                            current_left_pos=(converted_start+converted_end)/2 -arrows_width/2 + self.arrow_width
                        for i in range(n_arrows):
                            if bedrect.strand=="+":
                                box["ax"].plot([current_left_pos,current_left_pos+self.arrow_width,current_left_pos],[arrow_top,arrow_center,arrow_bottom],color="white",lw=0.4,zorder=3)
                            else:
                                box["ax"].plot([current_left_pos+self.arrow_width,current_left_pos,current_left_pos+self.arrow_width],[arrow_top,arrow_center,arrow_bottom],color="white",lw=0.4,zorder=3)
                            current_left_pos+=self.arrow_gapwidth+self.arrow_width


    def draw_title(self,box):
        if len(self.label)>0:
            self.label = self.label.replace("\\n","\n")
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 1,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
            
    def get_record_name(self,linesplit):
        if self.show_names:
            if len(linesplit)>3:
                return linesplit[3]
            else: 
                raise KnownException("For the bed file "+self.file+", you ticked the show_names option, but no names were found for line "+"\t".join(linesplit) +" (should be the 4th column).")
        else: return ""


    def get_record_strand(self,linesplit):
        if self.show_strand:
            if len(linesplit)>5:
                strand=linesplit[5]
                if not strand in ["+","-","."]: raise KnownException("Unknown strand information for line "+"\t".join(linesplit) +" (should be the '+', '-', or '.' in the 6th column).")
                return strand
            else: raise KnownException("For the bed file "+self.file+", you ticked the show_strand option, but no strand information was found for line "+"\t".join(linesplit) +" (should be the 6th column).")
        else: return "."
            
    def get_record_color(self,linesplit):
        if self.use_file_colors:
            if len(linesplit)>=9:
                return "#{:02x}{:02x}{:02x}".format(*[int(x) for x in linesplit[8].split(",")])
            else: raise KnownException("For the bed file "+self.file+", you ticked the use_file_colors option, but no color information was found for line "+"\t".join(linesplit) +" (should be the 9th column).")
        else: return self.color

    def read_records_regions(self,regions):
        regions = [reg[0] for reg in regions]
        lines_regions=[]
        for region in regions:
            lines_region=[]
            if self.collapsed: lines_region.append([])
            infile = gzip.open(self.file,"rt") if self.file.endswith(".gz") else open(self.file,"r")
            for line in infile:
                if line.startswith("#") or line.startswith("track") or line.startswith("browser"): continue
                linesplit = line.rstrip("\n").split("\t")
                if len(linesplit)<3: raise KnownException("The following bed file has the wrong format: "+self.file+".\nThe file should be tab-separated, without header, with at least 3 columns: chr, start, end.\nProblematic line:\n"+line[:70])
                if linesplit[0].lstrip("chr")==region.chr.lstrip("chr"):
                    try:
                        start,end = int(linesplit[1]),int(linesplit[2])
                    except:
                        raise KnownException("The following bed file has the wrong format: "+self.file+".\nThe file should be tab-separated, without header, with at least 3 columns: chr, start, end (where start and end should be integers). \nProblematic line:\n"+line[:70])
                    if start <=region.end and end >= region.start:
                        name=self.get_record_name(linesplit)
                        strand=self.get_record_strand(linesplit)
                        color=self.get_record_color(linesplit)
                        bedrect=BedRect(region.chr,start,end,name,strand,color)
                        if self.collapsed: lines_region[0].append(bedrect)
                        else: add_BedRect_pile(bedrect,lines_region,margin=0)
            infile.close()
            lines_regions.append(lines_region)
        max_nlines = max([len(l) for l in lines_regions])        

        new_line_order = reorder_balanced(0,max_nlines)

        for i in range(len(regions)):
            lines_regions[i]+=[[] for j in range(0,max_nlines-len(lines_regions[i]))]
            lines_regions[i] = [lines_regions[i][j] for j in new_line_order]
        return lines_regions

        



def add_BedRect_pile(bedrect,l,margin=0):
    # Assign rectangles to lines (to avoid overlap)
    i=0
    while i<len(l):
        last_rect = l[i][-1]
        if last_rect.end+margin<=bedrect.start:
            l[i].append(bedrect)
            return
        i+=1
    l.append([bedrect])