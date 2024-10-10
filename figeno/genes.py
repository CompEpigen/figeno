import os
import gzip
import importlib_resources as resources
import pathlib

import figeno.data
from figeno.utils import KnownException

from collections import namedtuple
Transcript = namedtuple('Transcript', 'name chr start end strand exons')
Gene = namedtuple('Gene', 'gene_id gene_name chr start end strand biotype')
      

def read_transcripts(file,chr=None,start=None,end=None,gene_names="auto",collapsed=True,only_protein_coding=True,warnings=[]):
    """Wrapper depending on the file type"""
    if gene_names =="auto" or gene_names=="": gene_names = None
    if gene_names is not None: 
        if not isinstance(gene_names,str): raise KnownException("gene_names should be a string of comma-separated gene names.")

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath) or file.endswith("txt.gz") or file.endswith(".txt"):
        return read_genes_refseq(file,chr,start,end,gene_names,collapsed=collapsed,only_protein_coding=only_protein_coding)
    elif file.endswith("gff3") or file.endswith("gff3.gz"):
        return read_genes_gff3(file,chr,start,end,gene_names,collapsed=collapsed,only_protein_coding=only_protein_coding) 
    elif file.endswith("gtf") or file.endswith("gtf.gz"):
        return read_genes_gtf(file,chr,start,end,gene_names,collapsed=collapsed,only_protein_coding=only_protein_coding,warnings=warnings) 
    else:
        raise KnownException("The extension for genes file "+str(file)+" was not recognized. The filename must end with .txt(.gz) for RefSeq format, .gtf(.gz) for gtf or .gff3(.gz) for gff3.")
        

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
    if gene_names is not None: gene_names=gene_names.upper()
    transcripts={}

    if isinstance(file,pathlib.PosixPath) or isinstance(file,pathlib.WindowsPath):
        if file.name.endswith(".gz") : infile = gzip.open(file,"rt")
        else: infile =open(file,"r")
    elif file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if len(linesplit)<14: raise KnownException("Wrong format for the genes file "+str(file)+" (there should be at least 14 columns per line).")
        
        # Filter based on coordinates
        if (chr is not None) and linesplit[2].lstrip("chr") != chr.lstrip("chr"): continue
        txStart = int(linesplit[4])
        txEnd = int(linesplit[5])
        if (start is not None) and (end is not None) and (start>txEnd or end<txStart): continue

        # Filter based on gene name
        gene_name=linesplit[12]
        if (gene_names is not None) and (not gene_name.upper() in gene_names):continue

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



def read_genes_gtf(gtf_file,chr=None,start=None,end=None,gene_names=None,collapsed=True,only_protein_coding=True,warnings=[]):
    if gene_names is not None: gene_names=gene_names.upper()
    transcripts={}
    name2exons={}

    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if len(linesplit)<9: raise KnownException("Wrong format for the genes file (there should be at least 9 columns per line).")
        
        if chr is not None and linesplit[0].lstrip("chr") != chr.lstrip("chr"): continue
        elif linesplit[2]=="transcript":
            transcript_end=int(linesplit[4]) 
            transcript_start = int(linesplit[3])
            if start is not None and end is not None and (transcript_start>end or transcript_end < start):continue
            transcript_orientation = linesplit[6]
            # Find transcript and gene names 
            gene_name,transcript_name = "",""
            for x in linesplit[8].split(";"):
                if x.lstrip(" ").startswith("gene_name"):
                    x = x[x.find("\"")+1:]
                    gene_name = x[:x.find("\"")]
                elif x.lstrip(" ").startswith("transcript_name"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
                elif transcript_name=="" and x.lstrip(" ").startswith("transcript_id"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
            if (gene_names is not None) and (not gene_name.upper() in gene_names): continue
            if only_protein_coding and (not "protein_coding" in linesplit[8]):
                if ("gene_biotype" in linesplit[8]): continue
                else: 
                    warn_message="Could not filter for protein coding genes because this information was not provided in the gtf file."
                    if not warn_message in warnings: warnings.append(warn_message) 
                
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
                if x.lstrip(" ").startswith("transcript_name"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
                elif transcript_name=="" and x.lstrip(" ").startswith("transcript_id"):
                    x = x[x.find("\"")+1:]
                    transcript_name = x[:x.find("\"")]
                elif x.lstrip(" ").startswith("gene_name"):
                    x = x[x.find("\"")+1:]
                    gene_name = x[:x.find("\"")]
            
            if (gene_names is not None) and (not gene_name.upper() in gene_names): continue
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

def read_genes_gff3(gff3_file,chr=None,start=None,end=None,gene_names=None,collapsed=True,only_protein_coding=True):
    if gene_names is not None: gene_names=gene_names.upper()
    transcriptID2transcript={}
    transcriptID2exons={}
    geneID2transcriptIDs={}
    geneID2gene_name={}

    if gff3_file.endswith(".gz") : infile = gzip.open(gff3_file,"rt")
    else: infile =open(gff3_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        if line=="\n": continue
        linesplit = line.rstrip("\n").split("\t")
        if len(linesplit)<9: raise KnownException("Wrong format for the genes file (there should be at least 9 columns per line).")
        
        if chr is not None and linesplit[0].lstrip("chr") != chr.lstrip("chr"): continue
        elif linesplit[2] in ["transcript","mRNA","rRNA","miRNA"]:
            transcript_end=int(linesplit[4]) 
            transcript_start = int(linesplit[3])
            if start is not None and end is not None and (transcript_start>end or transcript_end < start):continue
            transcript_orientation = linesplit[6]
            # Find transcript and gene names 
            gene_id,transcript_id, transcript_name = "","",""
            for x in linesplit[8].split(";"):
                if x.startswith("Parent"):
                    gene_id = x[x.find("=")+1:]
                elif x.startswith("Name") or x.startswith("transcript_name"):
                    transcript_name = x[x.find("=")+1:]
                elif x.startswith("ID"):
                    transcript_id = x[x.find("=")+1:]
            transcript={"id":transcript_id,"name":transcript_name,"chr":linesplit[0].lstrip("chr"),"start":transcript_start,"end":transcript_end,"orientation":transcript_orientation}
            transcriptID2transcript[transcript_id]=transcript
            if not gene_id in geneID2transcriptIDs: geneID2transcriptIDs[gene_id]=[]
            geneID2transcriptIDs[gene_id].append(transcript_id)

        elif linesplit[2]=="exon":
            exon_end=int(linesplit[4]) 
            exon_start = int(linesplit[3])
            if (start is not None) and (end is not None) and (exon_start>end or exon_end < start):continue
            transcript_id = ""
            for x in linesplit[8].split(";"):
                if x.startswith("Parent"):
                    transcript_id = x[x.find("=")+1:]
            if not transcript_id in transcriptID2exons: transcriptID2exons[transcript_id]=[]
            add_exon(transcriptID2exons[transcript_id],(exon_start,exon_end))

        elif linesplit[2] in ["gene", "miRNA_gene","rRNA_gene"]:
            if chr is not None and linesplit[0].lstrip("chr") != chr.lstrip("chr"): continue
            gene_end=int(linesplit[4]) 
            gene_start = int(linesplit[3])
            if start is not None and end is not None and (gene_start>end or gene_end < start):continue
            if only_protein_coding:
                if "biotype" in linesplit[8]:
                    if not "biotype=protein_coding" in linesplit[8]: continue
                elif "gene_type" in linesplit[8]:
                    if not "gene_type=protein_coding" in linesplit[8]: continue
                else:
                    if linesplit[2]!="gene": continue
            gene_id,gene_name="",""
            for x in linesplit[8].split(";"):
                if x.startswith("ID"):
                    gene_id = x[x.find("=")+1:]
                elif x.startswith("Name") or x.startswith("gene_name"):
                    gene_name = x[x.find("=")+1:]
            geneID2gene_name[gene_id]=gene_name
    infile.close()

    result=[]
    for gene_id in geneID2gene_name:
        if gene_id in geneID2transcriptIDs:
            gene_name=geneID2gene_name[gene_id]
            if gene_name=="": gene_name = gene_id
            if (gene_names is not None) and (not gene_name.upper() in gene_names): continue
            if collapsed:
                transcripts_chr=""
                transcripts_start=10000000000
                transcripts_end=0
                transcripts_orientation=""
                exons=[]
                for transcript_id in geneID2transcriptIDs[gene_id]:
                    transcript=transcriptID2transcript[transcript_id]
                    transcripts_start= min(transcripts_start,transcript["start"])
                    transcripts_end=max(transcripts_end,transcript["end"])
                    transcripts_orientation=transcript["orientation"]
                    transcripts_chr=transcript["chr"]
                    for exon in transcriptID2exons[transcript_id]:
                        add_exon(exons,exon)
                result.append(Transcript(gene_name,transcripts_chr,transcripts_start,transcripts_end,transcripts_orientation,exons))
            else:
                for transcript_id in geneID2transcriptIDs[gene_id]:
                    transcript=transcriptID2transcript[transcript_id]
                    result.append(Transcript(transcript["name"],transcript["chr"],transcript["start"],transcript["end"],transcript["orientation"],transcriptID2exons[transcript_id]))
    return result
                



def find_genecoord_wrapper(gene_name,reference,genes_file=None):
    
    if genes_file is None or genes_file=="":
        if reference in ["hg19","hg38","mm10"]:
            
            with resources.as_file(resources.files(figeno.data) / (reference+"_genes.txt.gz")) as infile:
                (chr,min_coord,max_coord) = find_genecoord_refseq(gene_name,infile,custom_ref=False)
                with resources.as_file(resources.files(figeno.data) / (reference+"_cytobands.tsv")) as infile2: #Ensure the end is smaller than the chr size.
                    chr_length=-1
                    with open(infile2,"r") as f:
                        for line in f:
                            if line.startswith("#"): continue
                            linesplit=line.split("\t")
                            if linesplit[0].lstrip("chr")==chr:
                                chr_length=max(chr_length,int(linesplit[2]))
                    if chr_length>=0:
                        min_coord=min(min_coord,chr_length)
                        max_coord=min(max_coord,chr_length)
                return (chr,min_coord,max_coord) 
        else:
            raise KnownException("You are using a custom reference genome, but did not provide a genes file. See https://figeno.readthedocs.io/en/latest/content/describe_figure.html#general for the format of a genes file.")
    else:
        if not os.path.isfile(genes_file): raise KnownException("The provided genes file does not exist: "+str(genes_file)+".")
        
        if genes_file.endswith(".gtf.gz") or genes_file.endswith(".gtf"):
            return find_genecoord_gtf(gene_name,genes_file)
        elif genes_file.endswith(".gff3.gz") or genes_file.endswith(".gff3"):
            return find_genecoord_gff3(gene_name,genes_file)
        else:
            return find_genecoord_refseq(gene_name,genes_file)

def find_genecoord_refseq(gene_name,file=None,custom_ref=True):
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
        if linesplit[12].upper()==gene_name.upper():
            chr_tmp=linesplit[2].lstrip("chr")
            if (custom_ref) or (not "_" in chr_tmp): # avoid alternative contigs...
                chr=linesplit[2].lstrip("chr")
                min_coord=min(min_coord,int(linesplit[4]))
                max_coord=max(max_coord,int(linesplit[5]))

    if chr!="":
        length=max_coord-min_coord
        min_coord-= max(10,int(0.05*length))
        max_coord+= max(10,int(0.05*length))

    if min_coord<=0: min_coord=0
    return (chr,min_coord,max_coord)

def find_genecoord_gtf(gene_name,file):
    if file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")
    current_chr=""
    current_start=-1
    current_end=-1

    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2] in ["gene","miRNA_gene","rRNA_gene"]:
            for x in linesplit[8].split(";"):
                if x.lstrip(" ").startswith("gene_name"):
                    x = x[x.find("\"")+1:]
                    name = x[:x.find("\"")]
                    if gene_name.upper()==name.upper():
                        chr=linesplit[0].lstrip("chr")
                        start=int(linesplit[3])
                        end=int(linesplit[4])
                        length=end-start
                        start-=max(10,int(0.05*length))
                        end+=max(10,int(0.05*length))
                        if start<=0: start=0
                        return (chr,start,end)
                    else: break
        if linesplit[2] in ["transcript"]:
            for x in linesplit[8].split(";"):
                if x.lstrip(" ").startswith("gene_name"):
                    x = x[x.find("\"")+1:]
                    name = x[:x.find("\"")]
                    if gene_name.upper()==name.upper():
                        chr=linesplit[0].lstrip("chr")
                        start=int(linesplit[3])
                        end=int(linesplit[4])
                        length=end-start
                        start-=max(10,int(0.05*length))
                        end+=max(10,int(0.05*length))
                        if start<=0: start=0
                        if current_chr=="":
                            current_chr,current_start,current_end=chr,start,end
                        else:
                            current_start=min(start,current_start)
                            current_end=max(end,current_end)
                    else: break
    return (current_chr,current_start,current_end)

def find_genecoord_gff3(gene_name,file):
    if file.endswith(".gz") : infile = gzip.open(file,"rt")
    else: infile =open(file,"r")

    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2] in ["gene","miRNA_gene","rRNA_gene"]:
            for x in linesplit[8].split(";"):
                if x.startswith("Name") or x.startswith("ID") or x.startswith("gene_name"):
                    name=x[x.find("=")+1:]
                    if gene_name.upper()==name.upper():
                        chr=linesplit[0].lstrip("chr")
                        start=int(linesplit[3])
                        end=int(linesplit[4])
                        length=end-start
                        start-=max(10,int(0.05*length))
                        end+=max(10,int(0.05*length))
                        if start<=0: start=0
                        return (chr,start,end)
                    else: continue
    return ("",0,1)








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