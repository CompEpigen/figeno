from figeno.utils import KnownException, correct_region_chr, chr_to_int
import numpy as np
import pysam
import copy

from collections import namedtuple
Read = namedtuple('Read', 'read qstart qend region_index')
Breakpoint = namedtuple('Breakpoint', 'chr1 pos1 strand1 chr2 pos2 strand2')
# qstart and qend are the start and end position in the query, if the alignment is in forward orientation and if we do not ignore hard alignments.
# This is a way to have the relative positioning of different split reads in the original fragment.


# Read placement
# For each group (a group is either all reads, or the reads corresponding to a haplotype):
#   - assign each read to the first row where the read would not overlap with a read already assigned to this row
#     (optionally with some minimum margin between reads)
#     Each group is implemented by a list of lists called a pile: where pile[row][i] is the i-th read in the row-th row (the reads are sorted within a row)
# In case we group splitreads:
# We first place the reads with supplementary alignments, and place the other reads after (so splitreads will be towards the top of the figure).
# When a read has supplementary alignments:
# Compute the coordinates of this alignment in the original query (in forward orientation): qstart and qend.
# Order the reads based on qstart
# Group together alignments which are adjacent (ie qend of the first is close to qstart of the next)
# Try to place all grouped alignments in the same row. This is not always possible, eg in case of foldback inversion two alignments 
# of the same read will align to some overlapping positions. When this is not possible, try to minimize the row distance between adjacent alignments.
# For now, this is implemented with a rather simple greedy algorithm. This seems to work well for simple cases, 
# but has not been extensively tested for reads with a large number of chimeric alignments.


def add_read_pile(read,l,margin=0):
    i=0
    while i<len(l):
        last_read = l[i][-1]
        if last_read.reference_end + margin<read.reference_start:
            l[i].append(read)
            return
        i+=1
    l.append([read])

def add_read_pile_greedy(read,pile,ideal_row,margin=10):
    # Try to add the read to ideal_row, or to the closest row if that's not possible.
    # Return the row where the read was inserted.
    distance_row=0 # distance from the current row to the ideal row
    while distance_row<1000:
        # row = ideal_row - distance_row
        row = ideal_row-distance_row
        if row>=0:
            if row>len(pile)-1:
                for x in range(len(pile),row+1): pile.append([])
            row_index=0
            # The rows are assumed to be sorted. Try to insert the read in this row.
            while row_index<len(pile[row]) and pile[row][row_index].reference_start<=read.reference_end+margin: row_index+=1
            if row_index==0:
                pile[row] = [read]+pile[row]
                return row
            else:
                if read.reference_start-margin>pile[row][row_index-1].reference_end:
                    pile[row].insert(row_index,read)
                    return row
        # row = ideal_row + distance_row
        if distance_row>0:
            row=ideal_row+distance_row
            if row>len(pile)-1:
                for x in range(len(pile),row+1): pile.append([])
            row_index=0
            # The rows are assumed to be sorted. Try to insert the read in this row.
            while row_index<len(pile[row]) and pile[row][row_index].reference_start<=read.reference_end+margin: row_index+=1
            if row_index==0:
                pile[row] = [read]+pile[row]
                return row
            else:
                if read.reference_start-margin>pile[row][row_index-1].reference_end:
                    pile[row].insert(row_index,read)
                    return row
        distance_row+=1
    return -1

def add_read_pile_nopreference(read,pile,margin=10,row_offset=0):
    row = row_offset
    while row<len(pile):
        row_index=0
        # The rows are assumed to be sorted. Try to insert the read in this row.
        while row_index<len(pile[row]) and pile[row][row_index].reference_start<=read.reference_end+margin: row_index+=1
        if row_index==0:
            pile[row] = [read]+pile[row]
            return
        else:
            if read.reference_start-margin>pile[row][row_index-1].reference_end:
                pile[row].insert(row_index,read)
                return
        row+=1
    pile.append([read])

def add_splitreads_piles(reads,piles_list,margin=10,max_depth=None):
    # reads is a list of read with(read,qstart,qend,region_index)
    #piles: for each region, one pile per group (either all, or HP1, HP2, unphased...)
    # if len(pile)==1: all, if 2: HP1 and HP2. If 3: HP1, HP2 and unphased.
    best_score = 100000
    best_piles=[]
    max_len = np.max([np.max([len(p) for p in piles]) for piles in piles_list])
    n_piles = len(piles_list[0]) # 1 if ungrouped, 2 if HP1+HP2, 3 if HP1+HP2+unphased

    # Greedy approach: try all possible values for the row of the first read, and then select the rows which work best for the next reads.
    i=0 # i: row of the first read
    while i<max_len+2:
        current_piles_list = copy.deepcopy(piles_list)
        current_ideal_row=i
        score=0
        for read_index in range(len(reads)):
            region_index = reads[read_index].region_index
            group_index = 0
            if n_piles>1:
                if reads[read_index].read.has_tag("HP"):
                    if reads[read_index].read.get_tag("HP")==2: group_index=1
                else: group_index=2

            row = add_read_pile_greedy(reads[read_index].read,current_piles_list[region_index][group_index],current_ideal_row,margin=margin)
            score+=abs(row-current_ideal_row)
            current_ideal_row=row
        if score<best_score:
            best_score = score 
            best_piles = current_piles_list
        i+=1
    return best_piles


def read2query_start_end(read):
    cigar = read.cigarstring
    length=0
    current_start=0
    current_pos=1
    qstart=0
    while current_pos<len(cigar):
        while cigar[current_pos].isdigit(): current_pos+=1
        if current_start==0 and cigar[current_pos] in ["S","H"]:
            qstart = int(cigar[:current_pos])
        if cigar[current_pos] in ["I","M"]:
            length+=int(cigar[current_start:current_pos])
        current_start = current_pos+1
        current_pos = current_start+1

    current_pos = len(cigar)-2
    while cigar[current_pos].isdigit(): current_pos-=1
    if cigar[-1] in ["S","H"]:qend = int(cigar[current_pos+1:-1])
    else: qend=0

    if (read.flag&16)!=0:
        qstart= qend
    qend = qstart+length
    return qstart,qend



def find_splitreads(samfile,regions,keep_unphased=True,min_splitreads_breakpoints=2):
    d={} # First find reads which have an SA tag
    for region_index,region in enumerate(regions):
        for read in samfile.fetch(region.chr,region.start,region.end):
            if (not keep_unphased) and (not read.has_tag("HP")): continue
            if read.is_secondary: continue
            if read.has_tag("SA"): 
                qstart,qend = read2query_start_end(read)
                read2 = Read(read,qstart,qend,region_index)
                if read.query_name in d: d[read.query_name].append(read2)
                else: d[read.query_name] = [read2]
    
    d_filtered={}  # Only keep reads for which the supplementary alignments are in the regions.
    for x in d:
        if len(d[x])>1: 
            # Sort the reads by qstart
            l = sorted(d[x],key=lambda r:r.qstart)
            # Create groups where the alignments were adjacent in the original fragment.
            max_grouplength=0
            adjacent_groups=[]
            current_group=[l[0]]
            current_index=1
            while current_index<len(l):
                if abs(l[current_index-1].qend-l[current_index].qstart)<150: current_group.append(l[current_index])
                else:
                    adjacent_groups.append(current_group)
                    max_grouplength= max(max_grouplength,len(current_group))
                    current_group = [l[current_index]]
                current_index+=1
            adjacent_groups.append(current_group)
            max_grouplength= max(max_grouplength,len(current_group))
            if max_grouplength>1:
                d_filtered[x] = adjacent_groups

    # Only keep splitreads corresponding to recurrent breakpoints
    query_qpos_breakpoint, bp_counts= splitreads2breakpoints(d_filtered)
    if min_splitreads_breakpoints>0:
        d_filtered2={}
        for qname in d_filtered:
            keep_qname=False
            for g in d_filtered[qname]:
                for r in g[1:]:
                    bp = query_qpos_breakpoint[r.read.query_name][r.qstart]
                    if bp_counts[bp]>min_splitreads_breakpoints: keep_qname=True
            if keep_qname: d_filtered2[qname] = d_filtered[qname]
        d_filtered=d_filtered2



    
    return d_filtered, query_qpos_breakpoint, bp_counts

def round_bp(bp,bp_counts):
    # if a similar breakpoint is already similar, return it instead
    for bp2 in bp_counts:
        if bp.chr1==bp2.chr1 and bp.chr2==bp2.chr2 and bp.strand1==bp2.strand1 and bp.strand2==bp2.strand2 and \
            abs(bp.pos1-bp2.pos1)<30 and abs(bp.pos2-bp2.pos2)<30:
            return bp2
    return bp

def splitreads2breakpoints(splitreads):
    query_qpos_breakpoint={}
    bp_counts={}
    for query_name in splitreads:
        for group in splitreads[query_name]:
            for i in range(1,len(group)):
                r1 = group[i-1]
                r2 = group[i]
                qpos = r1.qend
                qpos2 = r2.qstart
                chr1 = r1.read.reference_name
                pos1 = r1.read.reference_end if (r1.read.flag&16)==0 else r1.read.reference_start
                strand1 = "-" if (r1.read.flag&16)==0 else "+"
                chr2 = r2.read.reference_name
                pos2 = r2.read.reference_start if (r2.read.flag&16)==0 else r2.read.reference_end
                strand2 = "+" if (r2.read.flag&16)==0 else "-"
                if chr_to_int(chr1)>chr_to_int(chr2) or (chr_to_int(chr1)==chr_to_int(chr2) and pos2<pos1):
                    chr1,chr2 = chr2,chr1
                    pos1,pos2 = pos2,pos1
                    strand1,strand2=strand2,strand1
                if not query_name in query_qpos_breakpoint: query_qpos_breakpoint[query_name] = {}
                bp = Breakpoint(chr1,pos1,strand1,chr2,pos2,strand2)
                # Make sure that a similar breakpoint is not already present
                bp = round_bp(bp,bp_counts)
                query_qpos_breakpoint[query_name][qpos] = bp
                if qpos2!=qpos:query_qpos_breakpoint[query_name][qpos2] = bp
                if bp in bp_counts: bp_counts[bp]+=1
                else: bp_counts[bp]=1
    
    return query_qpos_breakpoint, bp_counts

    





    
def add_reads_to_piles(samfile,piles_list,regions,splitreads={},margin=10,only_show_splitreads=False,only_one_splitread_per_row=True):
    # First add the split reads
    for queryname in splitreads:
        for group in splitreads[queryname]:
            piles_list = add_splitreads_piles(group,piles_list,margin=margin)

    piles_lengths=[]
    for reg_piles in piles_list:
        lengths=[]
        for pile in reg_piles:
            lengths.append(len(pile))
        piles_lengths.append(lengths)
    if not only_show_splitreads:
        # Then add the other reads
        for region_index,region in enumerate(regions):
            for read in samfile.fetch(region.chr,region.start,region.end):
                if read.query_name in splitreads: continue
                if read.is_secondary: continue
                group_index=0
                if len(piles_list[region_index])>1:
                    if read.has_tag("HP"):
                        if read.get_tag("HP")==2: group_index=1
                    else:
                        group_index=2
                if group_index < len(piles_list[region_index]):
                    row_offset=piles_lengths[region_index][group_index] if only_one_splitread_per_row else 0
                    add_read_pile_nopreference(read,piles_list[region_index][group_index],margin=margin,row_offset=row_offset)
    return piles_list
    
    


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


def reverse_complement(seq):
    l=""
    for x in seq[::-1]:
        if x=="A": l+="T"
        elif x=="C": l+="G"
        elif x=="G": l+="C"
        elif x=="T": l+="A"
        else: l+=x
    return l

def fix_hardclipped_read(read,samfile):
    # It is not possible to get methylation information (MM/ML tags) for hardclipped alignments, because part of the sequence information is missing.
    # This function will try to find the full sequence (in the primary alignment), and update the read accordingly
    if not "H" in read.cigarstring: return read
    if not read.has_tag("SA"): return read
    SA_list = read.get_tag("SA").split(";")
    for SA in SA_list:
        SA = SA.split(",")
        if len(SA)<=1: continue
        for record in samfile.fetch(SA[0],int(SA[1])-1,int(SA[1])+1):
            if record.query_name==read.query_name:
                if not "H" in record.cigarstring:
                    seq=record.query_sequence
                    if (read.flag&16)!=(record.flag&16):
                        seq = reverse_complement(seq)
                    read.cigartuples = [(op,length) if op!=5 else (4,length) for (op,length) in read.cigartuples] # Replace H with S in the cigar
                    read.query_sequence=seq
                    return read
    return read

                    


def decode_read_basemods(read,basemods,samfile,fix_hardclip_basemod=True,warnings=[]):
    """basemods is a list of 1 or 2 tuples: (base,mod,color)"""
    if fix_hardclip_basemod: read = fix_hardclipped_read(read,samfile)
    if "H" in read.cigarstring: return []
    index_seq2ref_coord = {}
    read_length = read.infer_read_length()
    for x in read.get_aligned_pairs():
        if (x[1] is not None) and (x[0] is not None):
            if (read.flag&16)==0:
                index_seq2ref_coord[x[0]] = x[1]
            else:
                index_seq2ref_coord[read_length-x[0]-1] = x[1]

    if (not read.has_tag("MM")) or (not read.has_tag("ML")):
        raise KnownException("MM and ML tags are missing from the bam file, but are required in order to visualize base modifications.\n\n"\
                            "For ONT data, these tags should automatically be added if you run dorado with a modified bases model "\
                            "(see https://github.com/nanoporetech/dorado#modified-basecalling).")
    MMs = read.get_tag("MM").rstrip(";").split(";")
    ML=read.get_tag("ML")
    ML_index=0
    seq=read.get_forward_sequence()
    if seq is None: return []
    basemods_in_bam=[]
    basemods_missing=set()
    for b,m,c in basemods: basemods_missing.add((b,m))
    d={} # reference position to color of base modification
    for MM in MMs:
        MM = MM.split(",")
        base,mod = (MM[0][0],MM[0][2]) # Assume MM[0] is: base+mod?
        basemods_in_bam.append("("+base+","+mod+")")
        color=None
        for b,m,c in basemods:
            if b==base and m==mod:
                color=c
                basemods_missing.remove((b,m))
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

    if len(basemods_missing)>0:
        basemods_missing=["("+b+","+m+")" for b,m in basemods_missing]
        warning_message="The following base modifications were not found in the bam file: "+", ".join(basemods_missing)+".\n"
        warning_message+="Only the following base modifications were found in the bam file: "+", ".join(basemods_in_bam)+"."
        if not warning_message in warnings:
            warnings.append(warning_message)
    l=[]
    for x in sorted(d.keys()):
        l.append((x,x+1,d[x]))
    return l
