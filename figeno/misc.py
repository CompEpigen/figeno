

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