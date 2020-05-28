import os



def searchFunction(ListOfString,target):
    hitsAll = []
    for stringIdx in range(len(ListOfString)):
        A_string = ListOfString[stringIdx]
        i = len(target)-1
        hitsReplace = []
        hits = []
        for k in range(len(target)-1,len(A_string)):
                if A_string[k] in target[-1]:
                    hitsReplace.append(k)
       
        if len(hitsReplace) ==  len(A_string) - (len(target)-1):
            for k in target:
                if not k in target[-1]:
                    hits  = []
                    continueTF =True
            if continueTF:
                continue
            hits = []
            for hit in hitsReplace:
                hits.append(hit - len(target))
            hitsAll.append([stringIdx,hits])
            continue
            
        while i >= 0:
            i -= 1
            
            hits = hitsReplace
            if not hitsReplace:
                break
            hitsReplace = []
           # print(hits)
            for hit in hits:
                if A_string[hit-1] in target[i]:
                        hitsReplace.append(hit-1)
            #print(hitsReplace,i)
                
                    
        if hits:
            hitsAll.append([stringIdx,hits])
   
    return hitsAll


def findMatchInProtien(names,sequences,key,sizeMotif):
    hits = searchFunction(sequences,key)
    matches = []
    #check position between 2 names then assign the last name to the motif
    #this should be replaced with a binary search or radix sort
    
    for hit in hits:
        for i in hit[1]:
          
            matches.append([names[hit[0]],sequences[hit[0]][i:i+sizeMotif]])

    return matches
def generateMotifTSVFiles(fileName,proteinsMotif):

    writeFile = "motifTSV/" + fileName.split(".")[0]
    f = open(writeFile.split(".")[0] + ".tsv", 'w')
    f.write("Sequence:\tIdentifier:\n")
    for i in range(0, len(proteinsMotif)):
      
        f.write(proteinsMotif[i][0] + "\t" + proteinsMotif[i][1]+'\n')
    f.close()
def openFastaInFolder(path,key,numResidues):
    for filename in os.listdir(path):
        with open(path +"/"+filename) as file:
            print(filename)
            lines = "\n"+ file.read()
            sequences = []
            names = []
            nameSeq = lines.split("\n>")
            
            for seq in nameSeq:
                if len(seq.split("\n")) > 1:
                    sequences.append("".join(seq.split("\n")[1:]))
                   
                    names.append(seq.split("\n")[0])
            
            motifMatchingProteins = findMatchInProtien(names,sequences,key,numResidues)
           
            generateMotifTSVFiles(filename,motifMatchingProteins)
        
    
            
    
