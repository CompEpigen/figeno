import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

from figeno.utils import split_box, draw_bounding_box, KnownException


class bed_track:
    def __init__(self,file,color="#051650",label="",label_rotate=False,show_names=False,show_strand=False,fontscale=1,bounding_box=False,height=3,margin_above=1.5,**kwargs):
        self.file = file
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.show_names = show_names
        self.show_strand=show_strand
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
            with open(self.file,"r") as infile:
                line=infile.readline()
        except:
            raise KnownException("The following bed file could not be read: "+self.file)
        
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            self.draw_region(regions[i][0],boxes[i],warnings)
        self.draw_title(box)

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the bed track because it is not one of the accepted parameters.")

    def draw_region(self,region,box,warnings):
        if self.bounding_box: draw_bounding_box(box)
        
        if self.show_names: 
            rect_height = (box["top"]-box["bottom"]) * 0.4
            rect_bottom = box["bottom"] +  (box["top"]-box["bottom"]) * 0.08
        else:
            rect_height = (box["top"]-box["bottom"]) * 0.7
            rect_bottom = box["bottom"] +  (box["top"]-box["bottom"]) * 0.15

        # Params for strand
        self.arrow_gapwidth=3 # gap between arrows
        self.arrow_height=min(0.8*rect_height,1.4)
        self.arrow_width=min(0.5,1.5*self.arrow_height)


        with open(self.file,"r") as infile:
            for line in infile:
                if line.startswith("#"): continue
                linesplit = line.rstrip("\n").split("\t")
                if len(linesplit)<3: raise KnownException("The following bed file has the wrong format: "+self.file+".\nThe file should be tab-separated, without header, with at least 3 columns: chr, start, end.\nProblematic line:\n"+line[:70])
                if linesplit[0].lstrip("chr")==region.chr.lstrip("chr"):
                    try:
                        start,end = int(linesplit[1]),int(linesplit[2])
                    except:
                        raise KnownException("The following bed file has the wrong format: "+self.file+".\nThe file should be tab-separated, without header, with at least 3 columns: chr, start, end (where start and end should be integers). \nProblematic line:\n"+line[:70])
                    if start <=region.end and end >= region.start:
                        if region.orientation=="+":
                            converted_start = box["left"] + (start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                            converted_end = box["left"] + (end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                        else:
                            converted_start = box["right"] - (end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                            converted_end = box["right"] - (start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                        converted_start = min(box["right"],max(box["left"],converted_start))
                        converted_end = min(box["right"],max(box["left"],converted_end))
                        rect_width=converted_end-converted_start

                        rect = patches.Rectangle((converted_start,rect_bottom),rect_width,rect_height,color=self.color)
                        box["ax"].add_patch(rect)


                        if self.show_names:
                            if len(linesplit)>3:
                                box["ax"].text((converted_start+converted_end)/2,(rect_bottom+rect_height+1),linesplit[3],
                                           horizontalalignment="center",verticalalignment="bottom",fontsize=8*self.fontscale)
                            else: 
                                warnings.append("For the bed file "+self.file+", you ticked the show_names option, but no names were found for line "+line +" (should be the 4th column).")
                        

                        if self.show_strand:
                            arrow_center=rect_bottom+rect_height/2
                            arrow_bottom=arrow_center-self.arrow_height/2
                            arrow_top=arrow_center+self.arrow_height/2
                            if len(linesplit)>5:
                                strand=linesplit[5]
                                if not strand in ["+","-"]: continue
                                if rect_width<1.2*self.arrow_width: continue
                                # Number of arrows that can fit in the rectangle
                                n_arrows= 1+ int((rect_width-3*self.arrow_width)/(self.arrow_gapwidth+self.arrow_width))

                                if rect_width<3*self.arrow_width:
                                    current_left_pos=(converted_start+converted_end)/2 - self.arrow_width/2
                                else:
                                    arrows_width= 3*self.arrow_width + (n_arrows-1)*(self.arrow_gapwidth+self.arrow_width)
                                    current_left_pos=(converted_start+converted_end)/2 -arrows_width/2 + self.arrow_width
                                for i in range(n_arrows):
                                    if strand=="+":
                                        box["ax"].plot([current_left_pos,current_left_pos+self.arrow_width,current_left_pos],[arrow_top,arrow_center,arrow_bottom],color="white",lw=0.4,zorder=3)
                                    else:
                                        box["ax"].plot([current_left_pos+self.arrow_width,current_left_pos,current_left_pos+self.arrow_width],[arrow_top,arrow_center,arrow_bottom],color="white",lw=0.4,zorder=3)
                                    current_left_pos+=self.arrow_gapwidth+self.arrow_width
                            else:
                                warnings.append("For the bed file "+self.file+", you ticked the show_strand option, but no strand information was found for line "+line +" (should be the 6th column).")


    def draw_title(self,box):
        if len(self.label)>0:
            self.label = self.label.replace("\\n","\n")
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 1,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)

