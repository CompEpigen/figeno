import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

from figeno.utils import split_box, draw_bounding_box


class bed_track:
    def __init__(self,file,color="#051650",label="",label_rotate=False,show_names=False,fontscale=1,bounding_box=False,height=6,margin_above=1.5):
        self.file = file
        self.color=color
        self.label=label
        self.label_rotate=label_rotate
        self.show_names = show_names
        self.fontscale=fontscale
        self.bounding_box=bounding_box
        self.height = height
        self.margin_above=margin_above

    def draw(self, regions, box ,hmargin):
        boxes = split_box(box,regions,hmargin)
        for i in range(len(regions)):
            self.draw_region(regions[i][0],boxes[i])
        self.draw_title(box)

    def draw_region(self,region,box):
        if self.bounding_box: draw_bounding_box(box)

        rect_height = (box["top"]-box["bottom"]) * 0.7
        rect_bottom = box["bottom"] +  (box["top"]-box["bottom"]) * 0.15

        with open(self.file,"r") as infile:
            for line in infile:
                linesplit = line.rstrip("\n").split("\t")
                if linesplit[0].lstrip("chr")==region.chr.lstrip("chr"):
                    start,end = int(linesplit[1]),int(linesplit[2])
                    if start <=region.end and end >= region.start:
                        if region.orientation=="+":
                            converted_start = box["left"] + (start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                            converted_end = box["left"] + (end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                        else:
                            converted_start = box["right"] - (end-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                            converted_end = box["right"] - (start-region.start) / (region.end-region.start) * (box["right"] - box["left"])
                        converted_start = min(box["right"],max(box["left"],converted_start))
                        converted_end = min(box["right"],max(box["left"],converted_end))

                        rect = patches.Rectangle((converted_start,rect_bottom),converted_end-converted_start,rect_height,color=self.color)
                        box["ax"].add_patch(rect)

                        if self.show_names and len(linesplit)>2:
                            box["ax"].text((converted_start+converted_end)/2,(rect_bottom+rect_height/2),linesplit[2],
                                           horizontalalignment="center",verticalalignment="center",fontsize=8*self.fontscale)
    def draw_title(self,box):
        if len(self.label)>0:
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 1,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)

