def f1(filepath):
    fileHandle=open(filepath,'r+')
    aminoAcidDict={'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',       
                 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                 'TTA': 'L', 'TCA': 'S', 'TAA': 'STOP', 'TGA': 'STOP',
                 'TTG': 'L', 'TCG': 'S', 'TAG': 'STOP', 'TGG': 'W',
                 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                 'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                 'ATG': 'START', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G','ATG':'M'}
               
    for aminoacid in fileHandle:
        aminoAcidDetail=aminoacid.split('')
        if aminoAcidDetail[0] in ['Start', 'Stop']:
            for code in aminoAcidDetail[1:]:
                aminoAcidDict[code] = [aminoAcidDetail[0]]
        else:
            for code in aminoAcidDetail[3:]:
                code = code.replace('\n', '')
                aminoAcidDict[code] = [aminoAcidDetail[0], aminoAcidDetail[1], aminoAcidDetail[2]]
    return f1
