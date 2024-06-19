import numpy as np
from matplotlib import patches
from matplotlib.collections import PatchCollection
from collections import namedtuple
Region = namedtuple('Region', 'chr start end orientation color')
Highlight = namedtuple('Highlight', 'chr start end color alpha')


class KnownException(Exception):
    pass

def correct_region_chr(region,chromosomes,file=""):
    if region.chr in chromosomes:
        return region
    elif "chr"+region.chr in chromosomes:
        return Region("chr"+region.chr,region.start,region.end,region.orientation,region.color)
    elif "Chr"+region.chr in chromosomes:
        return Region("Chr"+region.chr,region.start,region.end,region.orientation,region.color)
    elif region.chr.lstrip("chr") in chromosomes:
        return Region(region.chr.lstrip("chr"),region.start,region.end,region.orientation,region.color)
    else: 
        error_message="Could not find chromosome "+region.chr
        if file!="": error_message+=" in file "+file+"."
        else: error_message+="."
        error_message+=" Only the following chromosome names were found: "+", ".join(chromosomes)+" (the chr prefix can be omitted)."
        raise KnownException(error_message)

def split_box(box,regions,hmargin):
    boxes=[]
    current_left=box["left"]
    additional_width = hmargin/len(regions) # compensate the last margin
    if "projection" in box and box["projection"]=="polar":
        total_width = np.sum([x[1] for x in regions]) + hmargin * (len(regions)-1)
        for region in regions:
            width = region[1]
            width_angle = (box["right"] +0.25 - box["left"]) * width / total_width 
            margin_angle = (box["right"] +0.25  - box["left"]) * hmargin / total_width 
            box1= {"ax":box["ax"],"bottom":box["bottom"],"top":box["top"],"left":current_left,"right":current_left+width_angle,"projection":"polar"}
            current_left+=width_angle + margin_angle
            boxes.append(box1)
    elif "bottom2" in box: # 2 rows
        # Count number of regions in each row
        count0,count1=0,0
        length0,length1=0,0
        for reg,width,row in regions:
            if row==0: 
                count0+=1
                length0+=width
            else:
                count1+=1
                length1+=width
        if length0>length1: additional_width = hmargin/count0
        else: additional_width = hmargin/count1

        current_row=0
        for reg,width,row in regions:
            if row==0:
                box1= {"ax":box["ax"],"bottom":box["bottom"],"top":box["top"],"left":current_left,"right":current_left+width+additional_width-hmargin,"upside_down":True}
            else:
                if current_row==0:
                    current_row=1
                    current_left=box["left"]
                box1= {"ax":box["ax"],"bottom":box["bottom2"],"top":box["top2"],"left":current_left,"right":current_left+width+additional_width-hmargin}
            current_left+=width+additional_width
            boxes.append(box1)
    elif "total_height" in box:
        for i, region in enumerate(regions):
            box1= {"ax":box["ax"],"bottom":box["bottom"]-i*box["total_height"],"top":box["top"]-i*box["total_height"],"left":0,"right":region[1]}
            boxes.append(box1)
    else:  
        for region in regions:
            width = region[1]
            box1= {"ax":box["ax"],"bottom":box["bottom"],"top":box["top"],"left":current_left,"right":current_left+width+additional_width-hmargin}
            current_left+=width+additional_width
            boxes.append(box1)
    return boxes

def draw_bounding_box(box):
    #rect = patches.Rectangle((box["left"],box["bottom"]),(box["right"]-box["left"])*1.00,box["top"]-box["bottom"],color="black",fc='none',lw=0.5,zorder=2)
    #box["ax"].add_patch(rect)
    w=0.3
    if "projection" in box and box["projection"]=="polar":
        w=0.2
        # Left
        theta,r1,r2 = box["left"] , box["bottom"] - w , box["top"] + w
        x1,y1 = polar2cartesian((theta,r1))
        x2,y2 = x1 - w*np.sin(theta) , y1 + w*np.cos(theta)
        x4,y4 = polar2cartesian((theta,r2))
        x3,y3 =x4- w*np.sin(theta) , y4 + w*np.cos(theta)

        left_vertices = [cartesian2polar((x1,y1)) , cartesian2polar((x2,y2)) , cartesian2polar((x3,y3)) ,cartesian2polar((x4,y4))]
        #left_vertices = interpolate_polar_vertices(left_vertices)
        left = patches.Polygon(left_vertices,color="black",lw=0,zorder=2)

        # Right
        theta,r1,r2 = box["right"] , box["bottom"] - w , box["top"] + w
        x1,y1 = polar2cartesian((theta,r1))
        x2,y2 = x1 + w*np.sin(theta) , y1 - w*np.cos(theta)
        x4,y4 = polar2cartesian((theta,r2))
        x3,y3 =x4+ w*np.sin(theta) , y4 - w*np.cos(theta)
        right_vertices = [cartesian2polar((x1,y1)) , cartesian2polar((x2,y2)) , cartesian2polar((x3,y3)) ,cartesian2polar((x4,y4))]
        #right_vertices = interpolate_polar_vertices(right_vertices)
        right = patches.Polygon(right_vertices,color="black",lw=0,zorder=2)

        # Top
        top_vertices = [(box["left"]+0.001,box["top"]) , (box["left"]+0.001 , box["top"]+w) , (box["right"]-0.001 , box["top"] + w) , (box["right"]-0.001,box["top"])]
        top_vertices = interpolate_polar_vertices(top_vertices)
        top = patches.Polygon(top_vertices,color="black",lw=0,zorder=2)

        # Bottom
        bottom_vertices = [(box["left"]+0.001,box["bottom"]) , (box["left"]+0.001 , box["bottom"]-w) , (box["right"]-0.001 , box["bottom"] - w) , (box["right"]-0.001,box["bottom"])]
        bottom_vertices = interpolate_polar_vertices(bottom_vertices)
        bottom = patches.Polygon(bottom_vertices,color="black",lw=0,zorder=2)

        coll = PatchCollection([left,right,top,bottom],match_original=True,zorder=5)
        box["ax"].add_collection(coll)
    
    else:
        rect_left = patches.Rectangle((box["left"]-w,box["bottom"]-w),width=w,height=box["top"]-box["bottom"]+2*w,color="black",zorder=2,lw=0)
        rect_right = patches.Rectangle((box["right"],box["bottom"]-w),width=w,height=box["top"]-box["bottom"]+2*w,color="black",zorder=2,lw=0)
        rect_top = patches.Rectangle((box["left"]-w,box["top"]),width=box["right"]-box["left"]+2*w,height=w,color="black",zorder=2,lw=0)
        rect_bottom = patches.Rectangle((box["left"]-w,box["bottom"]-w),width=box["right"]-box["left"]+2*w,height=w,color="black",zorder=2,lw=0)
        coll = PatchCollection([rect_left,rect_top,rect_right,rect_bottom],match_original=True,zorder=5)
        box["ax"].add_collection(coll)
        #box["ax"].add_patch(rect_left)

def draw_highlights(box,regions,highlights,hmargin):
    boxes = split_box(box,regions,hmargin)
    for highlight in highlights:
        for i in range(len(regions)):
            box = boxes[i]
            region = regions[i][0]
            if highlight.end>=region.start and highlight.start<=region.end:
                if region.orientation=="+":
                    left = box["left"] + (highlight.start-region.start) / (region.end-region.start) * (box["right"]-box["left"])
                    right = box["left"] + (highlight.end-region.start) / (region.end-region.start) * (box["right"]-box["left"])
                else:
                    left = box["right"] - (highlight.start-region.start) / (region.end-region.start) * (box["right"]-box["left"])
                    right = box["right"] - (highlight.end-region.start) / (region.end-region.start) * (box["right"]-box["left"])
                vertices = [(left,box["top"]),(right,box["top"]) , (right,box["bottom"]) , (left,box["bottom"])]
                if box["left"] < box["right"]:
                    vertices = [(max(box["left"],min(box["right"],x)),y) for (x,y) in vertices]
                else:
                    vertices = [(min(box["left"],max(box["right"],x)),y) for (x,y) in vertices]
                if "projection" in box and box["projection"]=="polar":
                    vertices = interpolate_polar_vertices(vertices)
                box["ax"].add_patch(patches.Polygon(vertices,color=highlight.color,alpha=highlight.alpha,lw=0,zorder=100))
            


def interpolate_polar_vertices(vertices):
    new_vertices = [vertices[0]]
    for i in range(len(vertices)):
        theta1,r1 = vertices[i]
        if i+1<len(vertices):
            theta2,r2 = vertices[i+1]
        else:
            theta2,r2 = vertices[0]
        n_points = round(abs(theta1-theta2) * 100)
        if n_points>0:
            x_list = np.linspace(theta1,theta2,n_points+2)
            y_list = np.linspace(r1,r2,n_points+2)
            new_vertices = new_vertices + list(zip(x_list,y_list))[1:-1]
        if i+1<len(vertices):
            new_vertices.append(vertices[i+1])
    return new_vertices

def compute_rotation_text(angle_text):
    while angle_text>2*np.pi: angle_text-=2*np.pi
    while angle_text<0: angle_text+=2*np.pi
    angle_text_deg = angle_text *180 / np.pi # in [0 , +360]
    if angle_text_deg>180:
        rotation_text = angle_text_deg+90
        valign="top"
    else:
        rotation_text = angle_text_deg-90
        valign="bottom"

    if angle_text_deg<45 or angle_text_deg>315 or (angle_text_deg>135 and angle_text_deg < 225):
        valign="center"
    if angle_text_deg <90 or angle_text_deg>270:
        halign="left"
    else:
        halign="right"
    return rotation_text,halign,valign

def polar2cartesian(vertex):
    theta,r = vertex
    return ( r*np.cos(theta) , r*np.sin(theta) )

def cartesian2polar(vertex):
    x,y = vertex
    r = np.sqrt(x*x + y*y)
    theta = np.arctan(y/x)
    if x<0: theta+=np.pi
    if theta<np.pi/2: theta+=2*np.pi
    return (theta,r)
            
def chr_to_int(chr):
    chr=str(chr)
    if chr.isdigit():
        return int(chr)
    elif chr=="X":
        return 23
    else:
        return 24
    
def set_box_defaults(box):
    if not "bottom" in box: box["bottom"]=0
    if not "top" in box: box["top"]=1
    if not "left" in box: box["left"]=0
    if not "right" in box: box["right"]=1
    return box