import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import importlib_resources as resources
import figeno.data
from figeno.genes import read_transcripts
from figeno.utils import KnownException, correct_region_chr, split_box, draw_bounding_box, interpolate_polar_vertices, compute_rotation_text, polar2cartesian, cartesian2polar


class genes_track:
    def __init__(self,reference="custom",genes_file="",style="default",collapsed=True,only_protein_coding=True,genes="auto",show_gene_names=True,exon_color="#4a69bd",fontscale=1,bounding_box=False,height=12,margin_above=1.5,
                 label="",label_rotate=False,**kwargs):
        self.reference=reference
        self.genes_file = genes_file
        self.style=style
        self.collapsed=collapsed
        self.only_protein_coding=only_protein_coding
        self.genes=genes
        self.show_gene_names=show_gene_names
        self.exon_color=exon_color
        self.fontscale=float(fontscale)
        self.bounding_box=bounding_box
        self.height = float(height)
        self.margin_above= float(margin_above)
        self.label=label
        self.label_rotate=label_rotate
        self.kwargs=kwargs

    def draw(self, regions, box ,hmargin=0,warnings=[]):
        boxes = split_box(box,regions,hmargin)
        self.margin_between_genes = 1.5*np.sum([abs(r[0].end-r[0].start) for r in regions]) / abs(box["right"]-box["left"]) # in bp
        lines_regions = self.read_transcripts_lines_regions(regions,warnings=warnings)
        max_region_size=0
        for i in range(len(regions)):
            max_region_size=max(max_region_size,abs(regions[i][0].end-regions[i][0].start))
            self.draw_region(regions[i][0],boxes[i],lines_regions[i])
        self.draw_title(box)

        if max_region_size>1e7:
            warnings.append("The genes track is intended for rather small regions (<2Mb), but here you used very large regions.")

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the genes track because it is not one of the accepted parameters.")

    def draw_region(self,region,box,lines,xy_scale=1):
        if self.bounding_box: draw_bounding_box(box)
        def transform_coord(pos):
            if region.orientation=="+":
                transformed_pos = box["left"] + (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])
            else:
                transformed_pos = box["right"] - (pos-region.start) / (region.end-region.start) *(box["right"]-box["left"])
            if box["left"] < box["right"]:
                if transformed_pos <=box["left"]: return box["left"]
                elif transformed_pos >=box["right"]: return box["right"]
            else:
                if transformed_pos >=box["left"]: return box["left"]
                elif transformed_pos <=box["right"]: return box["right"]
            return transformed_pos

        n_lines = max(1,len(lines))
        height_transcript=(box["top"]-box["bottom"])/n_lines * 0.01
        height_exon=(box["top"]-box["bottom"])/n_lines * 0.42
        arrow_width = 1.5
        if "projection" in box and box["projection"]=="polar": arrow_width = -0.02

        for i in range(len(lines)):
            y=box["bottom"] + (box["top"]-box["bottom"])/n_lines * i + height_exon*0.7
            for transcript in lines[i]:
                if transcript.chr!=region.chr or (transcript.start > region.end) or transcript.end <=region.start: continue
                # Transcript
                #if region.orientation=="-": transcript = (transcript[1],transcript[0],transcript[2])
                start_coord = transform_coord(transcript.start) if transcript.strand=="+" else transform_coord(transcript.end)
                end_coord = transform_coord(transcript.end) if transcript.strand=="+" else transform_coord(transcript.start)
                vertices = [(start_coord,y-height_transcript/2) , (start_coord,y+height_transcript/2) , (end_coord,y+height_transcript/2) , (end_coord,y-height_transcript/2)]
                if "projection" in box and box["projection"]=="polar":
                    vertices = interpolate_polar_vertices(vertices)
                polygon = patches.Polygon(vertices,color="gray")
                box["ax"].add_patch(polygon)
   
                fontsize = 8*self.fontscale if n_lines<=1 else 6*self.fontscale
                if self.show_gene_names: 
                    if "projection" in box and box["projection"]=="polar":
                        angle_text,halign,valign = compute_rotation_text((start_coord+end_coord)/2)
                        box["ax"].text((start_coord+end_coord)/2,y+height_exon*1.4,transcript.name,rotation=angle_text,horizontalalignment="center",verticalalignment="center",fontsize=fontsize)
                        #box["ax"].text((start_coord+end_coord)/2,y+height_exon*0.6,transcript.name,rotation=angle_text,horizontalalignment=halign,verticalalignment=valign,fontsize=fontsize)
                    else:
                        box["ax"].text((start_coord+end_coord)/2,y+height_exon*0.6,transcript.name,horizontalalignment="center",verticalalignment="bottom",fontsize=fontsize)
                # Exons
                exons = sorted(transcript.exons)
                #exons_tmp = []
                #for exon in exons:
                #    if exon[1]>=region.start and exon[0]<=region.end: exons_tmp.append(exon)
                #exons = exons_tmp
                exon_count=0
                #transcript_orientation_region = "+" if (transcript[2]=="+" and region.orientation=="+") or (transcript[2]=="-" and region.orientation=="-") else "-"
                for exon in exons:
                    if exon[1]>=region.start and exon[0]<=region.end:

                        start_coord = transform_coord(exon[0])
                        end_coord = transform_coord(exon[1])
                        if region.orientation=="-": start_coord, end_coord = end_coord,start_coord
                        if self.style=="default" and \
                            ((exon_count==len(exons)-1 and transcript.strand=="+" and region.orientation=="+") or (exon_count==0 and transcript.strand=="-" and region.orientation =="-")) :
                            vertices=[(start_coord,y+height_exon/2) , (end_coord,y+height_exon/2), (end_coord+arrow_width,y) ,
                                                    (end_coord,y-height_exon/2) , (start_coord,y-height_exon/2)]
                            if "projection" in box and box["projection"]=="polar":
                                vertices = interpolate_polar_vertices(vertices)
                            polygon = patches.Polygon(vertices,color=self.exon_color)
                            box["ax"].add_patch(polygon)
                        elif self.style=="default" and \
                            ((exon_count==len(exons)-1 and transcript.strand=="+" and region.orientation=="-") or (exon_count==0 and transcript.strand=="-" and region.orientation =="+")) :
                            vertices=[(start_coord-arrow_width,y) , (start_coord,y+height_exon/2) , (end_coord,y+height_exon/2),
                                                    (end_coord,y-height_exon/2) , (start_coord,y-height_exon/2)]
                            if "projection" in box and box["projection"]=="polar":
                                vertices = interpolate_polar_vertices(vertices)
                            polygon = patches.Polygon(vertices,color=self.exon_color)
                            box["ax"].add_patch(polygon)
                        else:
                            vertices = [(start_coord,y-height_exon/2) , (start_coord,y+height_exon/2) , (end_coord,y+height_exon/2) , (end_coord,y-height_exon/2)]
                            if "projection" in box and box["projection"]=="polar":
                                vertices = interpolate_polar_vertices(vertices)
                            polygon = patches.Polygon(vertices,color=self.exon_color)
                            box["ax"].add_patch(polygon)
                        exon_count+=1

                # Arrow at TSS
                if self.style=="TSS_arrow":
                    if (transcript.strand=="+" and region.orientation=="+") or (transcript.strand=="-" and region.orientation=="-"):
                        left_coord = transform_coord(transcript.start) if transcript.strand=="+" else transform_coord(transcript.end) 
                        rect_width = (box["top"]-box["bottom"])*0.1
                        right_coord = left_coord + rect_width
                        height1 = height_exon*0.5
                        right_coord2 = left_coord+ height1 / xy_scale
                        bottom_coord = y + height_exon*0.7
                        top_coord = bottom_coord+height1
                        rect_height = rect_width * xy_scale
                        bottom_coord2 = top_coord - rect_height
                        rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color=self.exon_color)
                        box["ax"].add_patch(rect)
                        rect = patches.Rectangle((left_coord,bottom_coord2),right_coord2-left_coord,rect_height,color=self.exon_color)
                        box["ax"].add_patch(rect)
                        triangle = plt.Polygon([[right_coord2+rect_width*1.3,(bottom_coord2+top_coord)/2],[right_coord2,bottom_coord2-rect_height*0.8],[right_coord2,top_coord+rect_height*0.8]], color=self.exon_color)
                        box["ax"].add_patch(triangle)
                    else:
                        right_coord =  transform_coord(transcript.start) if transcript.strand=="+" else transform_coord(transcript.end) 
                        rect_width =(box["top"]-box["bottom"])*0.1
                        left_coord = right_coord - rect_width
                        height1 = height_exon*0.5
                        left_coord2 = right_coord - height1 / xy_scale
                        bottom_coord = y + height_exon*0.7
                        top_coord = bottom_coord+height1
                        rect_height = rect_width * xy_scale
                        bottom_coord2 = top_coord - rect_height
                        rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color=self.exon_color)
                        box["ax"].add_patch(rect)
                        rect = patches.Rectangle((left_coord2,bottom_coord2),right_coord-left_coord2,rect_height,color=self.exon_color)
                        box["ax"].add_patch(rect)
                        triangle = plt.Polygon([[left_coord2-rect_width*1.3,(bottom_coord2+top_coord)/2],[left_coord2,bottom_coord2-rect_height*0.8],[left_coord2,top_coord+rect_height*0.8]], color=self.exon_color)
                        box["ax"].add_patch(triangle)

    def draw_title(self,box):
        if "projection" in box and box["projection"]=="polar": 
            if len(self.label)>0:
                self.label = self.label.replace("\\n","\n")
                rotation = 90 if self.label_rotate else 0
                x,y= polar2cartesian((box["left"],(box["top"]+box["bottom"])/2))
                theta,r = cartesian2polar((x-1,y))

                box["ax"].text(theta,r,
                            self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=10*self.fontscale)
        else:
            if len(self.label)>0:
                self.label = self.label.replace("\\n","\n")
                rotation = 90 if self.label_rotate else 0
                box["ax"].text(box["left"] - 1.0,(box["top"]+box["bottom"])/2,
                            self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
                
    def read_transcripts_lines_regions(self,regions,warnings=[]):
        regions = [reg[0] for reg in regions]
        lines_regions=[]
        max_nlines=0
        for region in regions:
            if self.genes_file is None or self.genes_file=="":
                if self.reference in ["hg19","hg38","mm10"]:
                    with resources.as_file(resources.files(figeno.data) / (self.reference+"_genes.txt.gz")) as infile:
                        transcripts = read_transcripts(infile,region.chr,region.start,region.end,self.genes,collapsed=self.collapsed,only_protein_coding=self.only_protein_coding)
                else:
                    raise KnownException("When using a custom reference genome, you have to provide a genes file if you want to display a genes track. See https://figeno.readthedocs.io/en/latest/content/describe_figure.html#general for the format of this file.")
            else:
                transcripts = read_transcripts(self.genes_file,region.chr,region.start,region.end,self.genes,
                                               collapsed=self.collapsed,only_protein_coding=self.only_protein_coding,warnings=warnings)
            
            # Assign transcripts to lines, so that they do not overlap
            transcripts= sorted(transcripts,key=lambda x:x.start)
            lines = []
            if self.style=="default":
                for transcript in transcripts:
                    add_transcript_pile(transcript,lines,margin=self.margin_between_genes)
            else:
                lines.append([])
                for transcript in transcripts:
                    lines[0].append(transcript)
            lines_regions.append(lines)
            max_nlines = max(max_nlines,len(lines))

        # Reorder lines to distribute transcripts more evenly
        new_line_order = reorder_balanced(0,max_nlines)

        for i in range(len(regions)):
            lines_regions[i]+=[[] for j in range(0,max_nlines-len(lines_regions[i]))]
            lines_regions[i] = [lines_regions[i][j] for j in new_line_order]
        return lines_regions



def add_transcript_pile(transcript,l,margin=10000):
    # Assign transcripts to lines (to avoid overlap)
    i=0
    while i<len(l):
        last_transcript = l[i][-1]
        if last_transcript.end+margin<transcript.start:
            l[i].append(transcript)
            return
        i+=1
    l.append([transcript])


def merge_interleave(l1,l2):
    l=[]
    for i in range(max(len(l1),len(l2))):
        if i<len(l1): l.append(l1[i])
        if i<len(l2): l.append(l2[i])
    return l

def reorder_balanced(start,end):
    if start>end: return []
    if end-start==0: return []
    if end-start==1: return [start]
    mid=(end-start)//2
    if (end-start)%2==1:
        return [start+mid] + merge_interleave(reorder_balanced(start,start+mid),reorder_balanced(start+mid+1,end))
    else:
        return  merge_interleave(reorder_balanced(start,start+mid),reorder_balanced(start+mid,end))
