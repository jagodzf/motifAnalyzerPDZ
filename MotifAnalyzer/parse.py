import os
import requests
from fileOutput import writeRawTSV
from dl import getFasta

def readFile(struct, silent=True):
    struct.organism = getFasta(struct.organism, silent)
    
    if not struct.organism:
        return False

    filename = 'fastas/{}'.format(struct.organism)

    with open(filename, 'r') as f:
        #generate a list of all proteins and headers,structure: [[header,protein],etc]
        lines = f.read()
        proteins, proteinsLastN = [], []

        # open fasta file and append sequences to a proteins and last N residues to last proteins last N
        for i in lines.split('\n>'):
            x = i.split('\n',1)
            if len(x) > 1 and len(x[1]) > struct.numResidues:
                sequence = x[1].replace('\n','')
                proteins.append([x[0],sequence])
                proteinsLastN.append([x[0],sequence[-struct.numResidues:]])

        # check if raw tsv organism files already exist if not make them
        tsv_file = 'rawTSV/{}.tsv'.format(struct.organism.split('.')[0])
        if not os.path.exists(tsv_file):
            writeRawTSV(struct.organism, proteinsLastN)

        # remove the proteins that dont match the given motif
        motifMatchingProteins = [[i[0],i[1]]  for i in proteinsLastN if all([i[1][-(1+struct.pList[k])] in struct.importantPositions[k] for k in range(len(struct.pList))])]

        # only runs if motifMatchingProteins is not empty to prevent error in processFile().
        if motifMatchingProteins:
            processFile(proteinsLastN, motifMatchingProteins, struct)
            struct = findEnrichment(struct)
        return struct


def processFile(proteinsLastN,motifMatchingProteins,currentStructure):
    # record the freq of each amino acid for all the proteins
    for protein in proteinsLastN:
        char =  protein[1][-(1+currentStructure.searchPosition)]
        currentStructure.freqAllPositional[char] += 1
    for protein in motifMatchingProteins:
        char =  protein [1][-(1+currentStructure.searchPosition)]
        currentStructure.freqMotifPositional[char] += 1

    # record the freq of each amino acid for all the non redundant proteins
    for protein in set(list(zip(*proteinsLastN))[1]):
        char = protein[-(1+currentStructure.searchPosition)]
        currentStructure.freqNonRedAllPositional[char] += 1

    LabelAndSequences = list(zip(*motifMatchingProteins)) # transpose the motifMatchingProteins list
    # record the freq of each amino acid for all the motif matching non redundant proteins

    for protein in set(LabelAndSequences[1]):

        char = protein[-(1+currentStructure.searchPosition)]
        currentStructure.freqNonRedMotifPositional[char] += 1

    # save the redundant motif information to the current Structure
    currentStructure.sequences = LabelAndSequences[1]
    currentStructure.sequenceLabels = LabelAndSequences[0]
    return currentStructure

def findEnrichment(currentStructure):
    for k in currentStructure.freqNonRedAllPositional:
        acidFreqAll = currentStructure.freqNonRedAllPositional[k]
        acidFreqPos = currentStructure.freqNonRedMotifPositional[k]
        sumFreqAll = sum(currentStructure.freqNonRedAllPositional.values())
        sumFreqPos = sum(currentStructure.freqNonRedMotifPositional.values())

        if acidFreqAll != 0 and sumFreqAll != 0 and sumFreqPos != 0:
            currentStructure.enrichment[k] = round((acidFreqPos / sumFreqPos) / (acidFreqAll / sumFreqAll),3)

    return currentStructure

