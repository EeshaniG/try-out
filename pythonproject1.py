import operator
#################################
"""DNA TEMPLATE TO RNA STRAND"""
#################################

#function which has the dictionary to convert the DNA template to a complementary strand
def complement(dna):
    basecomplement={'A':'T','C':'G','T':'A','G':'C'}
    template=list(dna)
    template1=[basecomplement[base] for base in template]
    return ''.join(template1)
    
tmp='GTGTACCCCAAACCCATCGGGATTACTAAA'

cpm=complement(tmp)
print 'the complementary dna stand for the input template is:', cpm

#to get the final RNA strand corresponding to the DNA strand
rna=cpm.replace('T', 'U')
print 'the corresponding translated rna strand for the input dna template is:', rna

#function to split the DNA complementary strand into codons
def codons(cpm):
    end=len(cpm)-(len(cpm)%3)-1
    codons=[cpm[i:i+3] for i in range (0,end,3)]
    return codons
    
cd=codons(cpm)
print 'the dna strand stripped into codons looks like:', cd

#dictionary which gives amino acids in fasta format for the codons

amino_acids={'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',       
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
               

startTags = ['ATG'] #amino acid assembly begins once he start codon is encountered
stopTags = ['TAA', 'TGA', 'TAG'] #amino acid assembly stops once the stop codons are encountered
bStart = False
bStop = False
completeStr = ''
concatanteStr = ''
for v in cd:
    if v in startTags:
        bStart = operator.xor(1, bStart)
        if not bStart:
            concatanteStr = ''#it is required to empty the string if a relevant stop tag is not found after the start tag
            bStart = 1
        #bStart = True
    elif v in stopTags:
        if bStart:
            bStop = True
            completeStr = completeStr + ' ' + concatanteStr
            concatanteStr = ''
            bStart = False
            bStop = False
    elif bStart:
        concatanteStr = concatanteStr + v
if concatanteStr and not bStop:#it means we have not been able to find a valid 'Stop' for a 'Start'
    concatanteStr = ''
print "the codons which occur within the start and stop tags are:", completeStr

#the required codons have been isolatedd from within the start and the stop codons

def finalcodons(sampleCode):#this function splits the entire string into groups of 3 and creates the codes
    sampleCode = sampleCode.strip() + ' '#this extra space being added is very important as else, after it has read the last character
    #since there are no more char within the string, it comes out of the for loop and you will lose the last 3 characters. in this case: TAA
    aminoAcidCodesWithinInput = []
    counter = 0
    aminoAcidCode = ''
    for charVal in sampleCode:
        if counter < 3:#keep adding the char unitl the count is less than 3
            aminoAcidCode = aminoAcidCode + charVal
            counter = counter + 1
        else:#the moment the counter goes above 3, add the 3 char string formed above into a list, reset counter to 1 because you have already read the next char
            aminoAcidCodesWithinInput.append(aminoAcidCode)
            aminoAcidCode = charVal
            counter = 1
    return aminoAcidCodesWithinInput

c=finalcodons(completeStr)
print c   
t1=[amino_acids[base] for base in c]
d=''.join(t1)
print 'the final amino acid string is:',d   




##############################################
"""AMINO ACIDS TO DNA"""
##############################################

def protein(cd):
    AA={'F':'TTT', 'S':'TCT', 'Y':'TAT', 'C':'TGT',
                 'F':'TTC', 'S':'TCC', 'Y':'TAC', 'C':'TGC',
                 'L':'TTA', 'S':'TCA', 'STOP':'TAA','STOP':'TGA',
                 'L':'TTG', 'S':'TCG', 'STOP':'TAG', 'W':'TGG',
                 'L':'CTT', 'P':'CCT', 'H':'CAT', 'R':'CGT',
                 'L':'CTC', 'P':'CCC', 'H':'CAC', 'R':'CGC',
                 'L':'CTA', 'P':'CCA', 'Q':'CAA', 'R':'CGA', 
                 'L':'CTG', 'P':'CCG', 'Q':'CAG', 'R':'CGG',
                 'I':'ATT', 'T':'ACT', 'N':'AAT', 'S':'AGT',
                 'I':'ATC', 'T':'ACC', 'N':'AAC', 'S':'AGC',
                 'M':'ATA', 'T':'ACA', 'K':'AAA', 'R':'AGA',
                 'START':'ATG', 'T':'ACG', 'K':'AAG', 'R':'AGG',
                 'V':'GTT', 'A':'GCT', 'D':'GAT', 'G':'GGT', 
                 'V':'GTC', 'A':'GCC', 'D':'GAC', 'G':'GGC',
                 'V':'GTA','A':'GCA', 'E':'GAA', 'G':'GGA',
                 'V':'GTG', 'A':'GCG', 'EE':'GAG', 'G':'GGG', 'M':'ATG'}
                 
                                                    
    template=list(cd)
    template1=[AA[base] for base in template]
    return ''.join(template1)
    
aasequence='MMMYG'
dnacmp=protein(aasequence)
print 'The DNA strand corresponding to the amino acid sequence is:',dnacmp
                                