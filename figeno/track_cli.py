import json
from figeno.figeno import tracks_plot


#/home/e840r/Documents/Scripts/pyGenViz/figeno/config_16PB3075_TEKT1_enhancer2.json
with open('/home/e840r/Documents/tests/config_GDM1.json', 'r') as f:
    data = json.load(f)
    tp = tracks_plot()
    for s in data["regions"]:
        if "chr" in s:
            orientation = s["orientation"] if "orientation" in s else "+"
            if s["end"] < s["start"]: 
                s["start"],s["end"] = s["end"],s["start"]
                orientation="-"
            tp.add_region(chr=s["chr"],start=s["start"],end=s["end"],orientation=orientation)
        else: # Region can be a string formatted as  chr:start-end
            chr,coord=s.split(":")
            start,end=coord.split("-")
            start,end=int(start),int(end)
            if start>end:
                start,end=end,start
                orientation="-"
            else: orientation="+"
            tp.add_region(chr=chr,start=start,end=end,orientation=orientation)
    for t in data["tracks"]:
        if not "type" in t: raise Exception("Must specify the type for track: " +str(t))
        track_type = t.pop("type")
        if track_type == "alignments":
            tp.add_alignments_track(**t)
        elif track_type=="coverage":
            tp.add_coverage_track(**t)
        elif track_type=="bigwig":
            tp.add_bigwig_track(**t)
        elif track_type=="bed":
            tp.add_bed_track(**t)
        elif track_type=="Met freq" or track_type=="methylation_freq":
            tp.add_metfreq_track(**t)
        elif track_type=="HiC":
            tp.add_hic_track(**t)
        elif track_type=="genes":
            tp.add_genes_track(**t)
        elif track_type=="chr axis":
            tp.add_chraxis_track(**t)
        elif track_type=="copynumber":
            tp.add_copynumber_track(**t)
        elif track_type=="sv":
            tp.add_sv_track(**t)
        else:
            pass
            #raise Exception("Unrecognized track type: " + str(track_type)+". Please use one of: alignments, coverage, bigwig, bed, methylation_freq, genes or chr_axis.")
    
    if (not "output" in data) or (not "file") in data["output"]: raise Exception("Must specify the outfile.")
    tp.draw(**data["output"])

