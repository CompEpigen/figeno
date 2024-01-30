import pysam

from figeno.tracks_utils import correct_region_chr


def add_read_pile(read,l,margin=0):
    i=0
    while i<len(l):
        last_read = l[i][-1]
        if last_read.reference_end + margin<read.reference_start:
            l[i].append(read)
            return
        i+=1
    l.append([read])

def read_read_groups(samfile,region):
    region = correct_region_chr(region,samfile.references)
    reads_HP1=[]
    reads_HP2=[]
    reads_unphased=[]
    reads_all=[]
    margin = (region.end-region.start) / 500
    for read in samfile.fetch(region.chr,region.start , region.end):
        if "H" in read.cigarstring or read.is_secondary or read.is_duplicate: continue
        add_read_pile(read,reads_all,margin)
        if read.has_tag("HP"):
            if read.get_tag("HP")==1: add_read_pile(read,reads_HP1,margin)
            else: add_read_pile(read,reads_HP2,margin)
        else: add_read_pile(read,reads_unphased,margin)
    return reads_HP1, reads_HP2, reads_unphased, reads_all

def decode_read_basemod(read,base="C",mod="m"):
    if "H" in read.cigarstring: return [],0
    index_seq2ref_coord = {}
    for x in read.get_aligned_pairs():
        if x[1] is not None:
            index_seq2ref_coord[x[0]] = x[1]
    if (base, 1, mod) in read.modified_bases: 
        methyl = read.modified_bases[(base, 1, mod)]
        strand=1
    elif (base, 0, mod) in read.modified_bases: 
        methyl = read.modified_bases[(base, 0, mod)]
        strand=0
    else:
        return [],0
    l=[]
    for pos,lik in methyl:
        if pos in index_seq2ref_coord:
            if lik<110: l.append((index_seq2ref_coord[pos],index_seq2ref_coord[pos]+1,0))
            else: l.append((index_seq2ref_coord[pos],index_seq2ref_coord[pos]+1,1))
    return l,strand

def decode_read_basemods2(read,basemods):
    """basemods is a list of 1 or 2 tuples: (base,mod,color)"""
    if "H" in read.cigarstring: return [],0
    index_seq2ref_coord = {}
    for x in read.get_aligned_pairs():
        if x[1] is not None:
            index_seq2ref_coord[x[0]] = x[1]

    if (basemods[0][0], 1, basemods[0][1]) in read.modified_bases: 
        strand=1
    elif (basemods[0][0], 0, basemods[0][1]) in read.modified_bases: 
        strand=0
    else:
        return []
    
    basemods1 = read.modified_bases[(basemods[0][0],strand,basemods[0][1])]
    d={}
    for pos,lik in basemods1:
        if pos in index_seq2ref_coord:
            if lik<110: d[index_seq2ref_coord[pos]] = 0
            else: d[index_seq2ref_coord[pos]] = basemods[0][2]

    if len(basemods)>1:
        basemods2 = read.modified_bases[(basemods[1][0],strand,basemods[1][1])]
        for pos,lik in basemods2:
                if pos in index_seq2ref_coord:
                    if lik<110 and (not index_seq2ref_coord[pos] in d): d[index_seq2ref_coord[pos]] = 0
                    elif lik >=110: d[index_seq2ref_coord[pos]] = basemods[1][2]

    l=[]
    for x in sorted(d.keys()):
        l.append((x,x+1,d[x]))
    return l


def decode_read_basemods(read,basemods):
    """basemods is a list of 1 or 2 tuples: (base,mod,color)"""
    if "H" in read.cigarstring: return [],0
    index_seq2ref_coord = {}
    read_length = read.infer_read_length()
    for x in read.get_aligned_pairs():
        if (x[1] is not None) and (x[0] is not None):
            if (read.flag&16)==0:
                index_seq2ref_coord[x[0]] = x[1]
            else:
                index_seq2ref_coord[read_length-x[0]-1] = x[1]

    MMs = read.get_tag("MM").rstrip(";").split(";")
    ML=read.get_tag("ML")
    ML_index=0
    seq=read.get_forward_sequence()
    d={} # reference position to color of base modification
    for MM in MMs:
        MM = MM.split(",")
        base,mod = (MM[0][0],MM[0][2]) # Assume MM[0] is: base+mod?
        color=None
        for b,m,c in basemods:
            if b==base and m==mod: color=c
        if color is None: # Ignore this basemod.
            ML_index+= len(MM)-1
            continue
        seq_index=0
        for toskip in MM[1:]:
            toskip = int(toskip)
            current_base = seq[seq_index]
            while toskip>0 or current_base!=base:
                if current_base==base:
                    toskip-=1
                seq_index+=1
                current_base=seq[seq_index]
            if seq_index in index_seq2ref_coord:
                ref_index= index_seq2ref_coord[seq_index]
                if ML[ML_index]>120: d[ref_index] = color
                elif not ref_index in d: d[ref_index] = 0
            ML_index+=1
            seq_index+=1
    l=[]
    for x in sorted(d.keys()):
        l.append((x,x+1,d[x]))
    return l

def decode_read_m_hm(read):
    index_seq2ref_coord = {}
    for x in read.get_aligned_pairs():
        if x[1] is not None:
            index_seq2ref_coord[x[0]] = x[1]
    if ("C", 1, "m") in read.modified_bases: 
        strand=1
    elif ("C", 0, "m") in read.modified_bases: 
        strand=0
    else:
        return [],0
    l=[]
    for i in range(len(read.modified_bases[("C",strand,"m")])):
        pos = read.modified_bases[("C",strand,"m")][i][0]
        lik_m = read.modified_bases[("C",strand,"m")][i][1]
        lik_h = read.modified_bases[("C",strand,"h")][i][1]
        if pos in index_seq2ref_coord:
            if lik_m>120: l.append((index_seq2ref_coord[pos],index_seq2ref_coord[pos]+1,1))
            elif lik_h>120: l.append((index_seq2ref_coord[pos],index_seq2ref_coord[pos]+1,2))
            else: l.append((index_seq2ref_coord[pos],index_seq2ref_coord[pos]+1,0))
    return l,strand