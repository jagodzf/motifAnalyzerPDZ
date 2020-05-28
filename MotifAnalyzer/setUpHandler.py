import argparse
import os
import time

from parse import readFile
from fileOutput import writeCSV, writeSequenceLists


def getImportantPositions(motifs):
    '''
    Extract the positions and amino acids from the motifs entered at run time.
    :param motifs: (list) the motifs supplied (e.g. P0:ILVF)
    :return: a list of lists containing the position followed by the amino acid letters.
    '''
    importantPositions = []
    for m in motifs:
        position, aminoacids = m.split(':')
        currentPositions = [int(position[1:])] + list(aminoacids)
        importantPositions.append(currentPositions)
    return importantPositions


class structure(object):
    def __init__(self):
        aminoacids = ["G", "A", "V", "L", "I", "M", "S", "C", "T", "P", "F", "Y", "W", "H", "K", "R", "D", "E", "N", "Q", "X", "U", "Z", "B", "O"]
        self.orgFile = ""
        self.organism = "" # holds fasta name
        self.numResidues = 0 
        self.searchPosition = 0 
        self.pList = []
        self.importantPositions = []
        self.sequences = [] # holds the *redundant* sequences that also match the motifs.
        self.sequenceLabels = [] # labels for the above sequences from the fasta
        self.freqAllPsnl = {a:0 for a in aminoacids} 
        self.freqLastN = {a:0 for a in aminoacids}
        self.freqMotifPositional  = {a:0 for a in aminoacids} 
        # freqPositional is the same as freqAllPsnl but with redundancy removed
        self.freqAllPositional = {a:0 for a in aminoacids}
        self.freqNonRedMotifPositional  = {a:0 for a in aminoacids}
        self.freqNonRedAllPositional = {a:0 for a in aminoacids}
        self.enrichment = {a:0 for a in aminoacids}

    # prints are for debugging purposes
    def printStructure(self):
        print("Organism: " + self.organism)
        print("numResidues: " + self.numResidues)
        print("searchPostion: " + self.numResidues)
        print("pList: " + self.pList)
        print("importantPositions: " + self.importantPositions)
    
    def printFrequencies(self):
        print("freqAllPsnl: " + self.freqAllPsnl)
        print("freqLastN: " + self.freqLastN)
        print("freqPositional: " + self.freqNonRedMotifPositional)
        print("enrichment: " + self.enrichment)



def setUpStructure(filenames, numResidues, searchPosition, pList, importantPositions, o):
    '''
    Create the structure containing the organism, number of residues, given position,
    and the motif information. Calculate and write output to the required directories.
    :param filenames: (list) organisms provided by the user
    :param numResidues: (int) the number of residues to provide statistics on.
    :param searchPosition: (int) the current position of interest.
    :param pList: (list) list of positions provided by the user.
    :param importantPositions: (list) list of lists containing the position followed by the amino acid letters.
    :o: (bool) flag to save tsv files containing the sequences.
    '''
    for curr_file in filenames:
        #set up class
        struct = structure()
        struct.organism = curr_file
        struct.numResidues = numResidues
        struct.searchPosition = searchPosition
        struct.pList = pList
        struct.importantPositions = importantPositions
        
        org = curr_file.split('.')[0]
        
         # parse file (if file cannot be found, add to unfound organisms)
        ### readFile() is how parse.py gets called ###
        struct = readFile(struct)

        if not (struct):
            with open('unfoundOrganisms.txt', 'a') as f:
                f.write('{}\n'.format(org))
            continue

        # file writing
        csv_file = 'csv/{}{}.csv'.format(struct.organism.split('.')[0], str(struct.searchPosition))
        writeCSV(csv_file, struct)

        if o:
            seq_file = 'sequenceLists/{}.tsv'.format(struct.organism.split('.')[0])
            writeSequenceLists(seq_file, struct)
        


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Process C-terminal decameric matchingSequences')
    parser.add_argument('file', type=str,
        help='The file to be parsed')
    parser.add_argument('position', type=int,
        help='The position to search in')
    parser.add_argument('numResidues', type=int,
        help='The number of residues to provide statistics on.')
    parser.add_argument('sort', type=int,
        help='How to sort final output')
    parser.add_argument('-o', action='store_true',
        help='Add this flag to enable file output.')
    parser.add_argument('-q', action='store_true',
        help='Add this flag to silence printout of ratios')
    parser.add_argument('motifs', nargs='+',
        help='The matching positions with desired amino acids. Usage: P0:ILVF P2:ST ... PX:X')
    return parser.parse_args()


def main():
    args = parseArgs()

    importantPositions = getImportantPositions(args.motifs)
    pList = []
    for i in importantPositions:
        pList.append(i[0])

    start = time.time()

    setUpStructure(args.file.split(','), args.numResidues, args.position, 
                   pList, importantPositions, args.o)

    if not args.q:
        organism = '_'.join([args.file.split('_')[0], args.file.split('_')[1]])
        print('Position {} for organism {} completed in {} seconds'.format(args.position, organism, round(time.time() - start, 2)))


if __name__ == '__main__':
    main()
