import cooler
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm 
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PatchCollection

from figeno.utils import KnownException, correct_region_chr, split_box, draw_bounding_box, interpolate_polar_vertices

class hic_track:
    def __init__(self,file,max_dist=700,color_map="red",scale_max_percentile=95,interactions_across_regions=True,double_interactions_across_regions=True,
                 extend=True,upside_down=False,pixel_border=False,show_colorbar=False,scale="auto",scale_min=1.0,scale_max=1.12,
                 rasterize=True,label="",label_rotate=False,fontscale=1,bounding_box=True,height=50,margin_above=1.5,**kwargs):
        self.file = file # must be in cool format
        self.max_dist = int(float(max_dist) *1000)
        self.color_map=color_map
        self.scale_max_percentile = float(scale_max_percentile)
        self.interactions_across_regions=interactions_across_regions
        if isinstance(interactions_across_regions,str): self.interactions_across_regions= (interactions_across_regions.lower()=="true")
        self.double_interactions_across_regions=double_interactions_across_regions
        if isinstance(double_interactions_across_regions,str): self.double_interactions_across_regions= (double_interactions_across_regions.lower()=="true")
        self.extend=extend
        self.upside_down=upside_down
        self.pixel_border=pixel_border
        self.show_colorbar = show_colorbar
        self.scale=scale
        self.scale_min=float(scale_min)
        self.scale_max= float(scale_max)
        
        self.rasterize=rasterize

        self.label=label
        self.label_rotate=label_rotate
        self.fontscale = float(fontscale)
        self.bounding_box=bounding_box
        self.height=float(height)
        self.margin_above=float(margin_above)
        self.kwargs=kwargs
        

    def draw(self, regions, box ,hmargin,warnings=[]):
        if "projection" in box and box["projection"]=="polar": raise KnownException("The circular layout does not support hic tracks.")
        self.draw_title(box)
        boxes = split_box(box,regions,hmargin)
        if self.bounding_box:
            if self.interactions_across_regions: draw_bounding_box(box)
            else:
                for box1 in boxes:  draw_bounding_box(box1)

        regions= [reg[0] for reg in regions]

        if self.file=="" or self.file is None:
            raise KnownException("Please provide a file for the hic track. This file must be in .cool format.")
        try:
            if self.file.endswith(".mcool"):
                f=h5py.File(self.file,'r')
                resolutions=list(f["resolutions"].keys())
                selected_resolution=resolutions[0]
                warnings.append("You used a .mcool file without specifying the resolution ("+self.file+"), so the resolution was automatically set to "+str(selected_resolution)+". This .mcool file contains the following resolutions: "+\
                                ", ".join([str(x) for x in resolutions])+". You can select the resolution by adding ::resolutions//"+str(selected_resolution)+" to the filename.")
                self.file+="::resolutions//"+str(selected_resolution)
               
            c=cooler.Cooler(self.file)
        except: 
            raise KnownException("Failed to open file for the hic track: "+self.file+".\nThis file must be in .cool format.\nIf you have data in .hic format,"\
                                 "please convert it to .cool format first. See https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html.")
        resolution = c.binsize
        total_dist = np.sum([(region.end-region.start) for region in regions])
        self.max_dist=min(self.max_dist,total_dist)
        max_bindist=self.max_dist//resolution /2
        n = len(regions)

        self.balanced= "weight" in c.bins()[:].columns
        if not self.balanced:
            warnings.append("No balancing weights were found for file "+self.file+". Figeno will show the raw, unnormalized data. If you want to show normalized data, please run: \"cooler balance "+self.file+"\" and re-generate the figure."\
                            " See https://cooler.readthedocs.io/en/latest/cli.html#cooler-balance for more information and options.")

        angle = self.find_angle(regions,boxes,max_bindist)
        

        if self.scale=="auto":
            vmin,vmax = self.get_min_max_values(regions,0,self.scale_max_percentile)
        else: vmin,vmax = np.exp(self.scale_min),np.exp(self.scale_max)

        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        if self.color_map=="red" or self.color_map=="Red":
             self.color_map = LinearSegmentedColormap.from_list('interaction',
                    ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
        else:
            self.color_map = cm.RdYlBu_r
       
        colormap = cm.ScalarMappable(norm=norm, cmap=self.color_map)

        # Draw colorbar 
        if self.show_colorbar:
            colorbar_matrix = np.array([x for x in np.linspace(0,1,100)]).reshape((100,1))
            colorbar_top = box["top"] -(box["top"]-box["bottom"])*0.05
            colorbar_bottom=box["top"] -(box["top"]-box["bottom"])*0.35
            vertices=[(-8,colorbar_top) , (-2,colorbar_top) ,(-2,colorbar_bottom) , (-8,colorbar_bottom) ]
            polygon = patches.Polygon(vertices,lw=0.6,edgecolor="black",facecolor="none",zorder=3)
            #box["ax"].add_patch(polygon)
            box["ax"].text(-8.4,colorbar_top,"{:10.4f}".format(np.log(vmax)),ha="right",va="top",fontsize=6*self.fontscale)
            box["ax"].text(-8.4,colorbar_bottom,"{:10.4f}".format(np.log(vmin)),ha="right",va="bottom",fontsize=6*self.fontscale)
            box["ax"].imshow(colorbar_matrix,cmap=self.color_map,extent=(-8,-2,colorbar_top,colorbar_bottom))

        patches_list = []

        # For each pair of region
        for a in range(n):
            for b in range(a,n):
                if (not self.interactions_across_regions) and a!=b: continue
                region1,region2 = correct_region_chr(regions[a],c.chromnames), correct_region_chr(regions[b],c.chromnames,file=self.file)
                box1,box2 = boxes[a],boxes[b]
                start1,end1 = region1.start, region1.end
                start2,end2 = region2.start, region2.end
                left_offset,right_offset=0,0

                if self.extend: # Extend regions so that we can show a rectangle instead of triangles.
                    if a==0 or (not self.interactions_across_regions):
                        if region1.orientation=="+":
                            new_start = max(0,region1.start-self.max_dist)
                            left_offset = (region1.start-new_start) // resolution
                            start1 = region1.start - left_offset *resolution
                        else:
                            new_end = min(c.chromsizes[region1.chr],region1.end+self.max_dist)
                            left_offset = (new_end-region1.end) // resolution
                            end1 = region1.end + left_offset * resolution
                    if b==len(regions)-1 or (not self.interactions_across_regions):
                        if region2.orientation=="+":
                            new_end = min(c.chromsizes[region2.chr],region2.end+self.max_dist)
                            right_offset = (new_end-region2.end) // resolution
                            end2 = region2.end + right_offset * resolution
                        else:
                            new_start = max(0,region2.start-self.max_dist)
                            right_offset = (region2.start-new_start) // resolution
                            start2 = region2.start - right_offset *resolution

                #start1,end1,left_offset = region1.start, region1.end,0
                mat = c.matrix(balance=self.balanced).fetch(region1.chr+":"+str(start1)+"-"+str(end1),
                                                region2.chr+":"+str(start2)+"-"+str(end2))
                mat[np.isnan(mat)]=0
                if region1.orientation=="-": mat = mat[::-1,:]
                if region2.orientation=="-": mat = mat[:,::-1]


                # Translocations have a copy number 1/2 of intra-chromosomal interactions, so double them.
                if a!=b and self.double_interactions_across_regions: mat = 1+2*mat
                else: mat = 1+mat

                # Pixel sizes
                width1 = (box1["right"]-box1["left"]) / (mat.shape[0]-left_offset)
                height1 = width1 * angle
                width2 = (box2["right"]-box2["left"]) / (mat.shape[1]-right_offset)
                height2 = width2 * angle

                # For each pair of bin from each region
                for i in range(0,mat.shape[0]):
                    for j in range(mat.shape[1]):
                        if a==b and j<i-left_offset: continue

                        x1 = box1["left"] + (i-left_offset)*width1
                        x2 = box2["left"] + (j)*width2
                        x=(x1+x2)/2
                        if box["left"]<box["right"]:
                            if x+width1/2+width2/2<boxes[0]["left"] or x>boxes[-1]["right"]: continue
                        else:
                            if x+width1/2+width2/2>boxes[0]["left"] or x<boxes[-1]["right"]: continue
                        y = box1["bottom"] + abs(x-x1) * abs(angle) 
                        if self.upside_down: y = box1["top"] - abs(x-x1) * abs(angle)
                        if ((not self.upside_down) and y-height2>box1["top"]) or (self.upside_down and (y+height2<box["bottom"])): continue

                        vertices = [(x,y) , (x+width2/2, y+height2/2) , (x+width2/2 +width1/2, y+height2/2-height1/2) , (x+width1/2, y-height1/2)]

                        if box["left"]<box["right"]:
                            if vertices[2][0]<box1["left"]+1e-6: continue
                            if vertices[0][0]>box2["right"]-1e-6:continue
                        else:
                            if vertices[2][0]>box1["left"]-1e-6: continue
                            if vertices[0][0]<box2["right"]+1e-6:continue
                        if vertices[-1][1]>box["top"]-1e-6:continue
                        if vertices[1][1]<box["bottom"]+1e-6: continue
                        
                        if box["left"]<box["right"]:
                            vertices = [(min(box2["right"],max(box1["left"],z2)),max(box1["bottom"]-0.1,min(box1["top"]+0.1,z))) for (z2,z) in vertices] # TODO -0.05
                        else:
                            vertices = [(max(box2["right"],min(box1["left"],z2)),max(box1["bottom"]-0.1,min(box1["top"]+0.1,z))) for (z2,z) in vertices]
                        
                        if vertices[0][0]==vertices[1][0] and vertices[2][0]==vertices[3][0]: continue
                        #if abs(vertices[1][0]-vertices[0][0])<1e-6: vertices = vertices[1:]
                        
                        
                        # Adjust the coordinates for the extremities of the triangle.
                        if a!=b:
                            if j==0:
                                if abs(vertices[0][1]-box["top"])<1e-6 and (not self.upside_down):
                                    y_left= abs(vertices[-1][1]-vertices[0][1])
                                    vertices[0]= (max(vertices[0][0],vertices[-1][0]-y_left/angle),vertices[0][1])
                                if abs(vertices[0][1]-box["bottom"])<1e-6 and  self.upside_down:
                                    y_left= abs(vertices[1][1]-vertices[0][1])
                                    vertices[0]= (max(vertices[0][0],vertices[-1][0]-y_left/angle),vertices[0][1])
                            elif i==mat.shape[0]-1:
                                if abs(vertices[0][1]-box["top"])<1e-6 and (not self.upside_down):
                                    y_right= abs(vertices[-1][1]-vertices[2][1])
                                    vertices[2]= (min(vertices[2][0],vertices[-1][0]+y_right/angle),vertices[2][1])
                                if abs(vertices[0][1]-box["bottom"])<1e-6 and  self.upside_down:
                                    y_right= abs(vertices[1][1]-vertices[2][1])
                                    vertices[2]= (min(vertices[2][0],vertices[-1][0]+y_right/angle),vertices[2][1])


                        # Left side: triangle instead
                        if abs(vertices[1][0]-vertices[0][0])<1e-6: vertices = vertices[1:]

                        if "projection" in box and box["projection"]=="polar":
                            vertices = interpolate_polar_vertices(vertices)
                        
                        
                        if self.pixel_border:
                            polygon = patches.Polygon(vertices,lw=0.2,facecolor=colormap.to_rgba(mat[i,j]),edgecolor="black")
                        else:
                            polygon = patches.Polygon(vertices,lw=0.2,color=colormap.to_rgba(mat[i,j]))
                        #box1["ax"].add_patch(polygon)
                        patches_list.append(polygon)
        patches_coll = PatchCollection(patches_list,rasterized=self.rasterize,match_original=True)
        box1["ax"].add_collection(patches_coll)

        for x in self.kwargs:
            warnings.append(x+" parameter was ignored in the hic track because it is not one of the accepted parameters.")

    def get_min_max_values(self,regions,low_percentile,high_percentile):
        c=cooler.Cooler(self.file)
        l=[]
        for a in range(len(regions)):
            for b in range(a,len(regions)):
                if (not self.interactions_across_regions) and a!=b: continue
                region1,region2 = correct_region_chr(regions[a],c.chromnames), correct_region_chr(regions[b],c.chromnames,file=self.file)
                mat = c.matrix(balance=self.balanced).fetch(region1.chr+":"+str(region1.start)+"-"+str(region1.end),
                                                region2.chr+":"+str(region2.start)+"-"+str(region2.end))
                mat[np.isnan(mat)]=0
                if a!=b and self.double_interactions_across_regions: mat = 1+2*mat
                else: mat = 1+mat
                l.append(mat.flatten())
        l = np.concatenate(l)
        return np.percentile(l,low_percentile),np.percentile(l,high_percentile)
    
    def find_angle(self,regions,boxes,max_bindist):
        c=cooler.Cooler(self.file)
        min_angle=100
        if boxes[0]["left"]>boxes[0]["right"]: min_angle=-100
        for a in range(len(regions)):
            region1 = correct_region_chr(regions[a],c.chromnames,file=self.file)
            if not region1.chr in c.chromnames: raise KnownException("Could not find chromosome "+region1.chr+" in the .cool file. Please make sure that you specified the correct chromosome name (the chr prefix can be omitted). "\
                                                                     "Only the following chromosome names were found in the .cool file: "+", ".join(c.chromnames))
            box1 = boxes[a]
            try:
                mat = 1+c.matrix(balance=self.balanced).fetch(region1.chr+":"+str(region1.start)+"-"+str(region1.end))
            except ValueError:
                raise KnownException("Could not retrieve region "+region1.chr+":"+str(region1.start)+"-"+str(region1.end)+" in file "+self.file+"."\
                                     " Make sure that the region that you specified does not extend beyond the chromosome length, and that you did not subset your .cool file.")
            except Exception as e:
                raise e
            width = (box1["right"]-box1["left"]) / mat.shape[0]
            height = (box1["top"]-box1["bottom"]) / max_bindist
            angle = height / width
            if box1["left"]<box1["right"]:
                min_angle=min(min_angle,angle)
            else:
                min_angle=max(min_angle,angle)
        return min_angle
    
    def draw_title(self,box):
        if len(self.label)>0:
            self.label = self.label.replace("\\n","\n")
            rotation = 90 if self.label_rotate else 0
            box["ax"].text(box["left"] - 1.0,(box["top"]+box["bottom"])/2,
                        self.label,rotation=rotation,horizontalalignment="right",verticalalignment="center",fontsize=7*self.fontscale)
        





def compute_dist_reg(regions,a,b,x1,x2):
    if a==b:
        return abs(x1-x2)
    else:
        s=0
        if regions[a].orientation=="+": s+= abs(regions[a].end-x1)
        else:s+= abs(regions[a].start-x1)
        if regions[b].orientation=="+": s+= abs(regions[b].start-x2)
        else:s+= abs(regions[b].end-x2)
        for i in range(a+1,b):
            s+=regions[i].end-regions[i].start
        return s