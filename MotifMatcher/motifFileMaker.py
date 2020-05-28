import os
import xml.etree.cElementTree as ET
import lxml.etree as etree
#written by Jordan Valgardson

def motif_Finder(inTSVFolder,motif_type,motifLen,outFile,refOrganism):
    # takes an input folder and a reference organism computes all of the matches between
    # the reference organism motifs and the motifs of the reference organisms
    allEndings = {}
    motifEndings = {}
    writeFile = open(outFile,"wb")
    refEndings = set()
    # gets all of the unique protein endings for each organism in the tsv folder saves them to a dictionary
    for file in os.listdir(inTSVFolder):
        with open(inTSVFolder + file, 'r') as readFile:
            allEndings[file] = set()
            
            for line in readFile:
                tempSequence = line.split('\t')[1].replace('\n','')
                if len(tempSequence) == motifLen:
                    allEndings[file].add(tempSequence)
    # checks if the protein endings fufill the motif criteria
    for organism in allEndings:  
        motifEndings[organism] = set()
        for motif in allEndings[organism]:
           
            if all(motif[key] in motif_type[key] for key in motif_type):
                 motifEndings[organism].add(motif)

            
    # removes the reference organism from the rest of the organisms
    refEndings = refEndings | motifEndings[refOrganism]
    del motifEndings[refOrganism]
    #initializes the xml tree output
    
    root = etree.Element("root")
    counter = 0
    summary =  etree.SubElement(root, "Summary",name = refOrganism+","+str(len(refEndings))+";" +";".join([",".join([key,str(len(value))]) for key, value in motifEndings.items()]))
    node = etree.SubElement(summary,"node")
    tree = etree.ElementTree(root)
    beginingAndEnd = etree.tostring(tree, pretty_print=True)
    beginingAndEnd = beginingAndEnd.split(b"<node/>")
    refNum = 0
    writeFile.write(beginingAndEnd[0])
    for Ref_motif in refEndings:
        refEnding  = etree.SubElement(node, "RefSequence",name = Ref_motif)
        refNum += 1
        # calculates the score for every unique motif in the organism list against the reference organism
        for organism in motifEndings:
            counter += 1
          
            tempDict = {n:set() for n in range(1,motifLen+1)}
            for motif in motifEndings[organism]:
                score = sum([1 if motif[i] == Ref_motif[i] else 0 for i in  range(motifLen)])
                if score != 0:
                    tempDict[score].add(motif)
            NonRefOrganism = etree.SubElement(refEnding, "NonRefOrganism",name = organism)
            # add to tree
            for score in range(motifLen,0,-1):
                matchScore = etree.SubElement(NonRefOrganism, "match"+str(score)).text = ",".join(tempDict[score])
        if counter > 10000:
            tree = etree.ElementTree(node)
            writeFile.write(etree.tostring(tree, pretty_print=True).replace(b"<node>",b"").replace(b"</node>",b""))
            node.clear()
            counter = 0
    tree = etree.ElementTree(node)
    writeFile.write(etree.tostring(tree, pretty_print=True).replace(b"<node>",b"").replace(b"</node>",b""))
    #print tree to file
    writeFile.write(beginingAndEnd[1])
    writeFile.close()
    



        
        
    
    
