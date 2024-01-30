import vcfpy

def read_phased_vcf(vcf,chr,start,end,selected_SNPs=None):
    """Returns a list of (pos, allele in HP1, allele in HP2)"""
    l=[] 
    reader = vcfpy.Reader.from_path(vcf)
    for record in reader.fetch(chr,start,end):
        if (selected_SNPs is None) or record.POS in selected_SNPs:
            if record.calls[0].data.get("GT")=="0|1":
                l.append((record.POS-1,record.REF,record.ALT[0].value))
            elif record.calls[0].data.get("GT")=="1|0":
                l.append((record.POS-1,record.ALT[0].value,record.REF))
    return l
