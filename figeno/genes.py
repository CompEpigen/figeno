import gzip
import importlib_resources as resources
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pathlib

import figeno.data

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



      

def read_transcripts(file,chr=None,start=None,end=None,gene_names="auto",collapsed=True,only_protein_coding=True):
    """Wrapper depending on the file type and whether a region or gene_names are provided"""
    if gene_names =="auto": gene_names = None
    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath) or file.endswith("txt.gz") or file.endswith(".txt"):
        return read_genes_refseq(file,chr,start,end,gene_names,collapsed=collapsed,only_protein_coding=only_protein_coding)
    else:
        return read_genes_gtf(file,chr,start,end,gene_names,collapsed=collapsed,only_protein_coding=only_protein_coding) 
        

def add_exon(exons,exon):
    start2,end2=exon
    for i in range(len(exons)):
        start1,end1 = exons[i]
        if end2<start1: # insert here
            exons.insert(i,(start2,end2))
            return
        elif start2<end1: # merge here
            exons[i] = (min(start1,start2),max(end1,end2))
            return
    # Add at the end if not added before.
    exons.append((start2,end2))


def merge_transcripts(transcript1,transcript2):
    exons=[x for x in transcript1.exons]
    for exon in transcript2.exons:
        add_exon(exons,exon)
    transcript=Transcript(transcript1.name,transcript1.chr,min(transcript1.start,transcript2.start),max(transcript1.end,transcript2.end),transcript1.strand,exons)
    return transcript


# Consider all exons of all transcripts for a given gene.
def read_genes_refseq(file,chr=None,start=None,end=None,gene_names=None,collapsed=True,only_protein_coding=True):
    transcripts={}

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath):
        if file.name.endswith(".gz") : infile = gzip.open(file,"rt")
        else: infile =open(file,"r")
    elif file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        
        # Filter based on coordinates
        if (chr is not None) and linesplit[2].lstrip("chr") != chr.lstrip("chr"): continue
        txStart = int(linesplit[4])
        txEnd = int(linesplit[5])
        if (start is not None) and (end is not None) and (start>txEnd or end<txStart): continue

        # Filter based on gene name
        gene_name=linesplit[12]
        if (gene_names is not None) and (not gene_name in gene_names):continue

        # Filter based on protein coding
        if only_protein_coding and linesplit[13] in ["none","unk"]: continue

        strand = linesplit[3]
        exons=[]
        exon_starts = linesplit[9].split(",")
        exon_ends = linesplit[10].split(",")
        exons = []
        for i in range(len(exon_starts)-1):
            exons.append((int(exon_starts[i]),int(exon_ends[i])))
        transcript = Transcript(gene_name,linesplit[2].lstrip("chr"),txStart,txEnd,strand,exons)

        # Collapse all transcripts corresponding to the same gene or not.
        if collapsed: name=gene_name
        else: name = linesplit[1]
        if name in transcripts:
            transcripts[name] = merge_transcripts(transcripts[name],transcript)
        else:
            transcripts[name] = transcript
    infile.close()
    transcripts = [transcripts[name] for name in transcripts]
    return transcripts



def read_genes_gtf(gtf_file,chr=None,start=None,end=None,gene_names=None,collapsed=True,only_protein_coding=True):
    transcripts={}
    name2exons={}

    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        
        if chr is not None and linesplit[0].lstrip("chr") != chr.lstrip("chr"): continue
        elif linesplit[2]=="transcript":
            transcript_end=int(linesplit[4]) 
            transcript_start = int(linesplit[3])
            if start is not None and end is not None and (transcript_start>end or transcript_end < start):continue
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
            if (gene_names is not None) and (not gene_name in gene_names): continue
            if only_protein_coding and (not "protein_coding" in linesplit[8]): continue
            if collapsed: name = gene_name
            else: name=transcript_name
            transcript = Transcript(gene_name,linesplit[0].lstrip("chr"),transcript_start,transcript_end,transcript_orientation,[])
            if name in transcripts: transcripts[name] = merge_transcripts(transcripts[name],transcript)
            else: transcripts[name] = transcript
        elif linesplit[2]=="exon":
            exon_end=int(linesplit[4]) 
            exon_start = int(linesplit[3])
            if (start is not None) and (end is not None) and (exon_start>end or exon_end < start):continue
            gene_name,transcript_name = "",""
            for x in linesplit[8].split(";"):
                if x.startswith(" transcript_name"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
                if x.startswith(" gene_name"):
                    x = x[x.find("\"")+1:]
                    gene_name = x[:x.find("\"")]
            if (gene_names is not None) and (not gene_name in gene_names): continue
            if collapsed: name = gene_name
            else: name=transcript_name
            if not name in name2exons: name2exons[name] = []
            add_exon(name2exons[name],(exon_start,exon_end))
    infile.close()

    result=[]
    for x in transcripts:
        if x in name2exons: exons=name2exons[x]
        else: exons = []
        result.append(Transcript(transcripts[x].name,transcripts[x].chr,transcripts[x].start,transcripts[x].end,
                                 transcripts[x].strand,exons))
    return result


def find_genecoord_refseq_wrapper(gene_name,reference,genes_file=None):
    
    if genes_file is None or genes_file=="":
        if reference in ["hg19","hg38","mm10"]:
            with resources.as_file(resources.files(figeno.data) / (reference+"_genes.txt.gz")) as infile:
                return find_genecoord_refseq(gene_name,infile)
        else:
            raise Exception("Using a custom reference genome, but no gene file was provided.")
    else:
        return find_genecoord_refseq(gene_name,genes_file)

def find_genecoord_refseq(gene_name,file=None):
    chr=""
    min_coord=1e10
    max_coord=0

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath):
        if file.name.endswith(".gz") : infile = gzip.open(file,"rt")
        else: infile =open(file,"r")
    elif file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")

    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[12]==gene_name:
            chr=linesplit[2].lstrip("chr")
            min_coord=min(min_coord,int(linesplit[4]))
            max_coord=max(max_coord,int(linesplit[5]))

    if chr!="":
        length=max_coord-min_coord
        min_coord-= max(10,int(0.05*length))
        max_coord+= max(10,int(0.05*length))

        
    return (chr,min_coord,max_coord)







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