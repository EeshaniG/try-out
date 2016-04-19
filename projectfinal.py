import operator
import string
import random

#references: 1. http://www.proteinstructures.com/Structure/Structure/amino-acids.html


#################################
"""DNA TEMPLATE TO RNA STRAND"""
#################################

#function which has the dictionary to convert the DNA template to a complementary strand
def complement(dna):
    basecomplement={'A':'T','C':'G','T':'A','G':'C'}
    template=list(dna)
    template1=[basecomplement[base] for base in template]
    return ''.join(template1)
    
tmp='GTGTACCCCAAACCCTTTATCGGGATTACTAAA'
l=len(tmp)


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

print '###############################################################'

#################################################################################
'Predefined mutation in the string and the resultant amino acid sequence'
#################################################################################

d=completeStr.replace('A', 'G')
e=completeStr.replace('G', 'C')

print 'randomly mutated string is:', e


f=finalcodons(e)
print f



g=[amino_acids[base] for base in f]
h=''.join(g)

print 'the new amino acid sequence for the mutated string is:', h


print '###############################################################'

#################################################################################
'Random mutation in the string and the resultant amino acid sequence'
#################################################################################


def id_generator(size=6, chars=string.ascii_uppercase):
    return ''.join(random.choice(chars) for _ in range(size))

s=id_generator(6, "ATCG")

print 'the initial template is', s

s1=finalcodons(s)


p0=[amino_acids[base] for base in s1]
p1=''.join(p0)


n=0
p2=id_generator(6, "ATCG")


while p2!=p1:
   
    w1=id_generator(6, "ATCG")

 
    w2=finalcodons(w1)
    

    w3=[amino_acids[base] for base in w2]
    p2=''.join(w3)
    
    n=n+1
   
print 'the number of iterations after which we will get the same protein string is:', n


print '###################################################################'
#############################################################################
'identifying the hydrophobic and hydrophillic parts of the amino acid string'
#############################################################################

hydrophobic_hydrophillic={'R':'hp', 'K':'hp', 'D':'hp', 'E':'hp', 'Q':'hp', 'N':'hp', 'H':'hp', 'S':'hp', 'T':'hp', 'Y':'hp', 'C':'hp', 'M':'hp', 'W':'hp',
                          'A':'ph', 'I':'ph', 'L':'ph', 'F':'ph', 'V':'ph', 'P':'ph', 'G':'ph', 'STOP':'-'}
                         
#hp=hydrophillic: includes charged amino acids and polar amino acids
#ph=hydrophobic: includes non polar amino acids

solubility=[hydrophobic_hydrophillic[base] for base in t1]
solubility1=''.join(solubility)
print 'the nature of the amino acid string corresponding to the initial template with regards to solubility is:', solubility1 

phobic=solubility.count('ph')
phillic=solubility.count('hp')

print 'the number of hyrdophobic amino acids is:', phobic
print 'the number of hyrdophillic amino acids is:', phillic

if phobic>phillic:
    print 'the protein string has more hydrophobic nature'
    
elif phillic>phobic:
    print 'the protein string has more hydrophillic nature'
    
elif phillic==phobic:
    print 'the protein string has equal hydrophobic and hydrophillic nature'
print '#############################################################'    
##############################################################################
'sexual reproduction using DNA strand'
##############################################################################

parent1=id_generator(12, "ATCG")
parent2=id_generator(12, "ATCG")  

print 'one parent is:',parent1
print 'other parent is:',parent2

pp1=finalcodons(parent1)
pp2=[amino_acids[base] for base in pp1]
pp3=''.join(pp2)

print 'the amino acid sequence for the first parent is:', pp3

solp12=[hydrophobic_hydrophillic[base] for base in pp2]
solp13=''.join(solp12)
print 'the nature of the amino acid string corresponding to the first parent with regards to solubility is:', solp13 

phobicp12=solp12.count('ph')
phillicp12=solp12.count('hp')

print 'the number of hyrdophobic amino acids in the first parent is:', phobicp12
print 'the number of hyrdophillic amino acids in the first parent is:', phillicp12

if phobicp12>phillicp12:
    print 'the protein string of the first parent has more hydrophobic nature'
    
elif phillicp12>phobicp12:
    print 'the protein string of the first parent has more hydrophillic nature'

elif phillicp12==phobicp12:
    print 'the protein string has equal hydrophobic and hydrophillic nature'
print '#############################################################' 
############################################################################


mm1=finalcodons(parent2)
mm2=[amino_acids[base] for base in mm1]
mm3=''.join(mm2)


print 'the amino acid sequence for the second parent is:', mm3


solp22=[hydrophobic_hydrophillic[base] for base in mm2]
solp23=''.join(solp22)
print 'the nature of the amino acid string corresponding to the second parent with regards to solubility is:', solp23 

phobicp22=solp22.count('ph')
phillicp22=solp22.count('hp')

print 'the number of hyrdophobic amino acids in the second parent is:', phobicp22
print 'the number of hyrdophillic amino acids in the second parent is:', phillicp22

if phobicp22>phillicp22:
    print 'the protein string of the first parent has more hydrophobic nature'
    
elif phillicp22>phobicp22:
    print 'the protein string of the first parent has more hydrophillic nature'

elif phillicp22==phobicp22:
    print 'the protein string has equal hydrophobic and hydrophillic nature'
print '#############################################################' 
############################################################################


def chunk(in_string,num_chunks):
    chunk_size = len(in_string)//num_chunks
    if len(in_string) % num_chunks: chunk_size += 1
    iterator = iter(in_string)
    for _ in range(num_chunks):
        accumulator = list()
        for _ in range(chunk_size):
            try: accumulator.append(next(iterator))
            except StopIteration: break
        yield ''.join(accumulator)
        
parent1parts=list(chunk(parent1,2))
parent2parts=list(chunk(parent2,2))

print 'the first parent string split into two equal parts:', parent1parts
print 'the second parent string split into two equal parts:', parent2parts

p1p1=parent1parts[0]
p1p2=parent1parts[1]
p2p1=parent2parts[0]
p2p2=parent2parts[1]



child1=p1p1+p2p2

child2=p2p1+p1p2


print 'the first offspring is:',child1
print 'the second offspring is:', child2

print '#############################################################' 
o1=finalcodons(child1)
o2=[amino_acids[base] for base in o1]
o3=''.join(o2)

print 'the amino acid sequence for the first offspring is:', o3


c1=finalcodons(child2)
c2=[amino_acids[base] for base in c1]
c3=''.join(c2)

print 'the amino acid sequence for the second offspring is:', c3
print '#############################################################' 
########################################################################
sol=[hydrophobic_hydrophillic[base] for base in o2]
sol1=''.join(sol)
print 'the nature of the amino acid string corresponding to the first offspring with regards to solubility is:', sol1 

phobic1=sol.count('ph')
phillic1=sol.count('hp')

print 'the number of hyrdophobic amino acids in the first offspring is:', phobic1
print 'the number of hyrdophillic amino acids in the first offspring is:', phillic1

if phobic1>phillic1:
    print 'the protein string of the first offspring has more hydrophobic nature'
    
elif phillic1>phobic1:
    print 'the protein string of the first offspring has more hydrophillic nature'

elif phillic1==phobic1:
    print 'the protein string has equal hydrophobic and hydrophillic nature'
print '#############################################################' 
########################################################################
sol2=[hydrophobic_hydrophillic[base] for base in c2]
sol3=''.join(sol2)
print 'the nature of the amino acid string corresponding to the second offspring with regards to solubility is:', sol3 

phobic2=sol2.count('ph')
phillic2=sol2.count('hp')

print 'the number of hyrdophobic amino acids in the second offspring is:', phobic2
print 'the number of hyrdophillic amino acids in the second offspring is:', phillic2

if phobic2>phillic2:
    print 'the protein string of the second offspring has more hydrophobic nature'
    
elif phillic2>phobic2:
    print 'the protein string of the second offspring has more hydrophillic nature'

elif phillic2==phobic2:
    print 'the protein string has equal hydrophobic and hydrophillic nature'
print '#############################################################' 
############################################################################
                     