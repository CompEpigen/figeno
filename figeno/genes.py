import gzip

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pathlib

from collections import namedtuple
Transcript = namedtuple('Transcript', 'name chr start end strand exons')
Gene = namedtuple('Gene', 'gene_id gene_name chr start end strand biotype')





def read_transcript_name_gtf(gtf_file,gene_name,transcript_name=None):
    """
    Used for plotting.
    """
    exons=[]
    if transcript_name is None:
        if gene_name=="MECOM":
            transcript_name = "MECOM-201"
        else:
            transcript_name = gene_name+"-001"
    transcript_found=False
    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if not "transcript_name \""+transcript_name+"\"" in linesplit[8]: continue
        if linesplit[2]=="transcript":
            chr = linesplit[0]
            transcript_end=int(linesplit[4]) 
            transcript_start = int(linesplit[3])
            transcript_orientation = linesplit[6]
            transcript_found=True
        if linesplit[2]=="exon":
            exons.append((int(linesplit[3]),int(linesplit[4])))
    infile.close()
    if not transcript_found: 
        if transcript_name==gene_name+"-001":
            return read_transcript_name_gtf(gtf_file,gene_name,gene_name+"-002")
        return False
    return Transcript(gene_name,chr,transcript_start,transcript_end,transcript_orientation,exons)


def read_transcripts_region_gtf(gtf_file,chr,start,end,transcript_name=None):
    gene_names=set()
    gene_name2transcript={}
    transcript2exons={}
    transcript2infos={}

    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        
        
        if linesplit[0].lstrip("chr") != chr.lstrip("chr"): continue
        if linesplit[2]=="gene":
            gene_name = ""
            filter_gene=False
            for x in linesplit[8].split(";"):
                if x.startswith(" gene_name"):
                    x = x[x.find("\"")+1:]
                    gene_name = x[:x.find("\"")]
                if x.startswith(" gene_biotype"):
                    x = x[x.find("\"")+1:]
                    biotype = x[:x.find("\"")]
                    if biotype!="protein_coding": filter_gene=True
            if not filter_gene:
                gene_names.add(gene_name)
        elif linesplit[2]=="transcript":
            transcript_end=int(linesplit[4]) 
            transcript_start = int(linesplit[3])
            if transcript_start>end or transcript_end < start:continue
            transcript_orientation = linesplit[6]
            # Find transcript and gene names 
            gene_name,transcript_name = "",""
            for x in linesplit[8].split(";"):
                if x.startswith(" gene_name"):
                    x = x[x.find("\"")+1:]
                    gene_name = x[:x.find("\"")]
                elif x.startswith(" transcript_name"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
            transcript2infos[transcript_name] = Transcript(gene_name,chr,transcript_start,transcript_end,transcript_orientation,[])
            if not transcript_name in transcript2exons: transcript2exons[transcript_name] = []
            if not gene_name in gene_name2transcript:
                gene_name2transcript[gene_name] = transcript_name
            elif transcript_name.endswith("001"): gene_name2transcript[gene_name] = transcript_name
            elif transcript_name.endswith("002") and (not gene_name2transcript[gene_name].endswith("001")): gene_name2transcript[gene_name] = transcript_name
            
        elif linesplit[2]=="exon":
            exon_end=int(linesplit[4]) 
            exon_start = int(linesplit[3])
            if exon_start>end or exon_end < start:continue
            transcript_name = ""
            for x in linesplit[8].split(";"):
                if x.startswith(" transcript_name"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
            if not transcript_name in transcript2exons: transcript2exons[transcript_name] = []
            transcript2exons[transcript_name].append((exon_start,exon_end))
    infile.close()
    if "MECOM" in gene_name2transcript: gene_name2transcript["MECOM"] = "MECOM-201"
    result=[]
    for gene_name in gene_names:
        if not gene_name in gene_name2transcript: continue
        transcript_name = gene_name2transcript[gene_name]
        tx= transcript2infos[transcript_name]
        exons = transcript2exons[transcript_name]
        transcript = Transcript(tx.name,tx.chr,tx.start,tx.end,tx.strand,exons)
        result.append(transcript)

    return result


def read_transcripts_region(file,chr,start,end,transcript_name=None):
    transcripts=[]

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath):
        if file.name.endswith(".gz") : infile = gzip.open(file,"rt")
        else: infile =open(file,"r")
    elif file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        
        if linesplit[2].lstrip("chr") != chr.lstrip("chr"): continue
        txStart = int(linesplit[4])
        txEnd = int(linesplit[5])
        if start<=txEnd and end>=txStart:
            name=linesplit[12]
            strand = linesplit[3]
            exons=[]
            exon_starts = linesplit[9].split(",")
            exon_ends = linesplit[10].split(",")
            exons = []
            for i in range(len(exon_starts)-1):
                exons.append((int(exon_starts[i]),int(exon_ends[i])))
            transcripts.append(Transcript(name,chr,txStart,txEnd,strand,exons))
    infile.close()
    return transcripts

def read_transcripts_names(file,gene_names):
    transcripts=[]
    genes_found=set()

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath):
        if file.name.endswith(".gz") : infile = gzip.open(file,"rt")
        else: infile =open(file,"r")
    elif file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[12] in gene_names and not linesplit[12] in genes_found:
            genes_found.add(linesplit[12])
            txStart = int(linesplit[4])
            txEnd = int(linesplit[5])
            name=linesplit[12]
            strand = linesplit[3]
            exons=[]
            exon_starts = linesplit[9].split(",")
            exon_ends = linesplit[10].split(",")
            transcripts.append(Transcript(name,linesplit[2].lstrip("chr"),txStart,txEnd,strand,exons))
    infile.close()
    return transcripts
      

def read_transcripts(file,chr=None,start=None,end=None,gene_names="auto"):
    """Wrapper depending on the file type and whether a region or gene_names are provided"""
    if gene_names=="auto":
        if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath) or file.endswith("txt.gz") or file.endswith(".txt"):
            return read_transcripts_region(file,chr,start,end)
        else:
            return read_transcripts_region_gtf(file,chr,start,end)
    else:
        if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath) or file.endswith("txt.gz") or file.endswith(".txt"):
            return read_transcripts_names(file,gene_names)
        else:
            l=[read_transcript_name_gtf(file,gene_name) for gene_name in gene_names]
            l2=[]
            for x in l:
                if x!=False: l2.append(x)
            return l2



def find_gene_id(gene_name,gtf_file):
    gene_id=None
    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2]=="gene" and "gene_name \""+gene_name+"\"" in linesplit[8]: 
            for x in linesplit[8].split(";"):
                if "gene_id" in x:
                    gene_id = x[x.find("\"")+1:-1]
    infile.close()
    return gene_id

def find_gene_name(gene_id,gtf_file):
    gene_id=None
    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2]=="gene" and "gene_id\""+gene_id+"\"" in linesplit[8]: 
            for x in linesplit[8].split(";"):
                if "gene_name" in x:
                    gene_name = x[x.find("\"")+1:-1]
    infile.close()
    return gene_name

def draw_genes(box,segment,gtf,gene_names,xy_scale=1,fontscale=1):
    """
    box is a dictionary with keys: ax, left, right, top, bottom . 
    Tested with top-bottom = 2* right-left
    """
    def transform_coord(pos):
        if segment.orientation=="+":
            transformed_pos = box["left"] + (pos-segment.start) / (segment.end-segment.start) *(box["right"]-box["left"])
        else:
            transformed_pos = box["right"] - (pos-segment.start) / (segment.end-segment.start) *(box["right"]-box["left"])
        if transformed_pos <=box["left"]: return box["left"]
        elif transformed_pos >=box["right"]: return box["right"]
        return transformed_pos

    height_transcript=(box["top"]-box["bottom"]) * 0.01
    height_exon=(box["top"]-box["bottom"]) * 0.42
    y= box["bottom"] + height_exon*0.7

    if gene_names == "auto":
        genes_transcripts_exons = read_transcripts_region_gtf(gtf_file=gtf,chr=segment.chr,start=segment.start,end=segment.end)
    else:
        genes_transcripts_exons = [read_transcript_name_gtf(gtf,gene_name) for gene_name in gene_names]

    for gene,transcript,exons in genes_transcripts_exons:
        if gene.chr!=segment.chr or gene.start > segment.end or gene.end <=segment.start: continue
        # Transcript
        if segment.orientation=="-": transcript = (transcript[1],transcript[0],transcript[2])
        start_coord = transform_coord(transcript[0])
        end_coord = transform_coord(transcript[1])
        rect = patches.Rectangle((start_coord,y-height_transcript/2),end_coord-start_coord,height_transcript,color="gray")
        box["ax"].add_patch(rect)
        box["ax"].text((start_coord+end_coord)/2,y+height_exon*0.6,gene.gene_name,horizontalalignment="center",verticalalignment="bottom",fontsize=12*fontscale)
        # Exons
        for exon in exons:
            if exon[1]>=segment.start and exon[0]<=segment.end:
                start_coord = transform_coord(exon[0])
                end_coord = transform_coord(exon[1])
                rect = patches.Rectangle((start_coord,y-height_exon/2),end_coord-start_coord,height_exon,color="#4a69bd")
                box["ax"].add_patch(rect)

        # Arrow at TSS
        if (transcript[2]=="+" and segment.orientation=="+") or (transcript[2]=="-" and segment.orientation=="-"):
            left_coord = transform_coord(transcript[0])
            rect_width = (box["right"]-box["left"])*0.005
            right_coord = left_coord + rect_width
            height1 = height_exon*0.5
            right_coord2 = left_coord+ height1 / xy_scale
            bottom_coord = y + height_exon*0.7
            top_coord = bottom_coord+height1
            rect_height = rect_width * xy_scale
            bottom_coord2 = top_coord - rect_height
            rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color="#4a69bd")
            box["ax"].add_patch(rect)
            rect = patches.Rectangle((left_coord,bottom_coord2),right_coord2-left_coord,rect_height,color="#4a69bd")
            box["ax"].add_patch(rect)
            triangle = plt.Polygon([[right_coord2+rect_width*1.3,(bottom_coord2+top_coord)/2],[right_coord2,bottom_coord2-rect_height*0.8],[right_coord2,top_coord+rect_height*0.8]], color="#4a69bd")
            box["ax"].add_patch(triangle)
        else:
            right_coord = transform_coord(transcript[1])
            rect_width = (box["right"]-box["left"])*0.005
            left_coord = right_coord - rect_width
            height1 = height_exon*0.5
            left_coord2 = right_coord - height1 / xy_scale
            bottom_coord = y + height_exon*0.7
            top_coord = bottom_coord+height1
            rect_height = rect_width * xy_scale
            bottom_coord2 = top_coord - rect_height
            rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color="#4a69bd")
            box["ax"].add_patch(rect)
            rect = patches.Rectangle((left_coord2,bottom_coord2),right_coord-left_coord2,rect_height,color="#4a69bd")
            box["ax"].add_patch(rect)
            triangle = plt.Polygon([[left_coord2-rect_width*1.3,(bottom_coord2+top_coord)/2],[left_coord2,bottom_coord2-rect_height*0.8],[left_coord2,top_coord+rect_height*0.8]], color="#4a69bd")
            box["ax"].add_patch(triangle)

def draw_genes_transcript(box,gtf,chr,start,end,gene_name,transcript_name,xy_scale=20,fontscale=1):
    """
    box is a dictionary with keys: ax, left, right, top, bottom . 
    Tested with top-bottom = 2* right-left
    """
    def transform_coord(pos):
        transformed_pos = box["left"] + (pos-start) / (end-start) *(box["right"]-box["left"])
        if transformed_pos <=box["left"]: return box["left"]
        elif transformed_pos >=box["right"]: return box["right"]
        return transformed_pos

    height_transcript=(box["top"]-box["bottom"]) * 0.1
    height_exon=(box["top"]-box["bottom"]) * 0.3
    y= box["bottom"] + height_exon

    transcript = read_transcript_name_gtf(gtf,gene_name,transcript_name=transcript_name)

    # Transcript
    start_coord = transform_coord(transcript[0])
    end_coord = transform_coord(transcript[1])
    rect = patches.Rectangle((start_coord,y-height_transcript/2),end_coord-start_coord,height_transcript,color="#4a69bd")
    box["ax"].add_patch(rect)
    box["ax"].text((start_coord+end_coord)/2,y+height_exon*0.56,gene_name,horizontalalignment="center",verticalalignment="bottom",fontsize=14*fontscale)
    # Exons
    for exon in transcript.exons:
        if exon[1]>=start and exon[0]<=end:
            start_coord = transform_coord(exon[0])
            end_coord = transform_coord(exon[1])
            rect = patches.Rectangle((start_coord,y-height_exon/2),end_coord-start_coord,height_exon,color="#4a69bd")
            box["ax"].add_patch(rect)

    # Arrow at TSS
    if transcript[2]=="+":
        left_coord = transform_coord(transcript[0])
        rect_width = (box["right"]-box["left"])*0.005
        right_coord = left_coord + rect_width
        height1 = height_exon*1.0
        right_coord2 = left_coord+ height1 / xy_scale
        bottom_coord = y + height_exon*0.8
        top_coord = bottom_coord+height1
        rect_height = rect_width * xy_scale
        bottom_coord2 = top_coord - rect_height
        rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color="#4a69bd")
        box["ax"].add_patch(rect)
        rect = patches.Rectangle((left_coord,bottom_coord2),right_coord2-left_coord,rect_height,color="#4a69bd")
        box["ax"].add_patch(rect)
        triangle = plt.Polygon([[right_coord2+rect_width*1.3,(bottom_coord2+top_coord)/2],[right_coord2,bottom_coord2-rect_height*0.8],[right_coord2,top_coord+rect_height*0.8]], color="#4a69bd")
        box["ax"].add_patch(triangle)
    else:
        right_coord = transform_coord(transcript[1])
        rect_width = (box["right"]-box["left"])*0.005
        left_coord = right_coord - rect_width
        height1 = height_exon*1.0
        left_coord2 = right_coord - height1 / xy_scale
        bottom_coord = y + height_exon*0.8
        top_coord = bottom_coord+height1
        rect_height = rect_width * xy_scale
        bottom_coord2 = top_coord - rect_height
        rect = patches.Rectangle((left_coord,bottom_coord),right_coord-left_coord,height1,color="#4a69bd")
        box["ax"].add_patch(rect)
        rect = patches.Rectangle((left_coord2,bottom_coord2),right_coord-left_coord2,rect_height,color="#4a69bd")
        box["ax"].add_patch(rect)
        triangle = plt.Polygon([[left_coord2-rect_width*1.3,(bottom_coord2+top_coord)/2],[left_coord2,bottom_coord2-rect_height*0.8],[left_coord2,top_coord+rect_height*0.8]], color="#4a69bd")
        box["ax"].add_patch(triangle)
        



#fig,ax = plt.subplots(figsize=(10,5))
#ax.set_xlim(0,1)
#ax.set_ylim(0,10)

#gtf = "/home/e840r/Documents/Data/GRCh37/RNAseq/Homo_sapiens.GRCh37.87.gtf"
#draw_ax({"ax":ax,"left":0,"right":1,"bottom":0,"top":1},"7",156790000,156806000)
#draw_genes({"ax":ax,"left":0,"right":1,"bottom":1,"top":3},gtf,"7",156790000,156806000,["MNX1"])
#plt.axis('off')
#plt.show()
