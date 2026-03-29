from pathlib import Path
import json
import os

def infer_aa(result, cdn, notEnoughInfoThreshold=3, highConfThreshold=0.9999):
    """
    Inputs:
        result: a dictionary with decoding probabilities
        cdn: three letter codon
        notEnoughInfoThreshold: a threshold for not enough info
        highConfThreshold: a threshold for high confidence
    Outputs:
        - upper case letter for high confidence aa
        - '?' for insufficient info
        - '?' for uncertain (below high confidence threshold)
    """
    if result[cdn]["n_"]<=notEnoughInfoThreshold:
        inferredAA = '?'
    else:
        probs = [result[cdn][aa] for aa in aas]
        bestAA = aas[probs.index(max(probs))]
        p = probs[aas.index(bestAA)]
        if p>highConfThreshold:
            if bestAA=='?':
                inferredAA = '?'
            else:
                inferredAA = bestAA
        else:
            inferredAA = '?'
    return inferredAA

def count_all(stTT,infTT):
    out = {}
    names = ['# amino acid inferred stop codons', 
             '# uninferred stop codons', 
             '# expected amino acid inferred sense codons',
             '# unexpected amino acid infered sense codons',
             '# uninferred sense codons',
             'is interesting']
    out['# amino acid inferred stop codons'] = 0
    out['# uninferred stop codons'] = 0
    out['# expected amino acid inferred sense codons'] = 0
    out['# unexpected amino acid infered sense codons'] = 0
    out['# uninferred sense codons'] = 0
    out['is interesting'] = "FALSE"
    for stAA, infAA, codon in zip(stTT, infTT, codons):
        if stAA=='?':
            if infAA!='?':
                out['# amino acid inferred stop codons'] += 1
                out['is interesting'] = "TRUE"
            else:
                out['# uninferred stop codons'] += 1
        else:
            if infAA=='?':
                out['# uninferred sense codons'] += 1
            else:
                if infAA==stAA:
                    out['# expected amino acid inferred sense codons'] += 1
                else:
                    out['# unexpected amino acid infered sense codons'] += 1
                    if not ((codon=='AGA' or codon=='AGG') and infAA=='K'):
                         out['is interesting'] = "TRUE"
    return names,out
    

codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 
          'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
          'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 
          'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 
          'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 
          'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 
          'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 
          'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
standardTT = "FFLLSSSSYY??CC?WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
aas = "ACDEFGHIKLMNPQRSTVWY?"

dirsToProcess = ["results/"]
fileOutBase = "summary_"
categories = ["genomic"]

header = "assembly,timestamp,inferred genetic code,"
names,_ = count_all("","")
for name in names:
    header += name+','

for category in categories:
    outDicts = []
    fileOut = fileOutBase + category + ".csv"
    fOut = open(fileOut,'w')
    fOut.write(header[:-1]+'\n')
    for d in dirsToProcess:
            
        print(f"Processing directory {d} for {category}")

        fileList = os.listdir(d)
        fileListProcessed = []
        for fname in fileList:
            if not fname.endswith(category+".json"): continue
            f = open(d+'/'+fname,'r')
            try:
                result = json.load(f)
            except:
                print(f"Unable to load {fname}")
                f.close()
                continue
            f.close()
            inferredTT = ""
            for codon in codons:
                inferredTT += infer_aa(result,codon)
            countNames,counts = count_all(standardTT,inferredTT)
            line = fname.split('_')[0] + '_' + fname.split('_')[1] + ','
            line += str(result['timestamp']) + ','
            line += inferredTT + ','
            for name in countNames:
                line += str(counts[name])+','
            fOut.write(line[:-1]+'\n')

    fOut.close()

