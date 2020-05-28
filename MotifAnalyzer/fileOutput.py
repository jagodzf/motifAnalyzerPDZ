from datetime import datetime
import os

def parseFileNames(organismfile):
    ''' takes the organisms file and generates the names to use for the downloaded fastas.
    
        Pre conditions:
            organismfile: string with the directory for the organism file
        
        Post conditions:
            returns files: a list containing the names to be used for the fastas that will be downloaded. '''
    files = []
    with open(organismfile, 'r') as f:
        for line in f.readlines():
            line = line.strip().lower()
            organism, taxon, reviewed = line.split(' ')
            files.append('{}_TaxID_{}_R_{}.fasta'.format(organism, taxon, reviewed))
    files = list(set(files))
    return files

def makeFolders(path='./'):
    ''' creates the necessary directories for the results to be stored in.

        Pre conditions:
            path is a string with the location for the directories to be made
        
        Post conditions:
            the directories have been created '''
    folders = ['csv', 'fastas', 'sequenceLists', 'rawTSV']
    for folder in folders:
        filepath = os.path.join(path, folder)
        if not os.path.exists(filepath):
            os.makedirs(filepath)

def writeSummaryFile(rawFolder, filteredFolder, outfile, q=False):
    ''' writeSummaryFile generates the summary file containing organism name, total motif matching proteins, and total proteins.

        Pre conditions:
            rawFolder and filteredFolder are the locations of the raw tsv's and filtered tsv's
            outfile is the location to write the summary file
        
        Post conditions:
            Summary file exists in the location specified by outfile '''
    with open(outfile, 'w') as f:
        rawLengths, filteredLengths = {}, {}

        for curr_file in os.listdir(rawFolder):
            with open(os.path.join(rawFolder, curr_file), 'r') as r:
                rawLengths[curr_file] = len(r.readlines()) - 1

        for curr_file in os.listdir(filteredFolder):
            with open(os.path.join(filteredFolder, curr_file), 'r') as r:
                filteredLengths[curr_file] = len(r.readlines()) - 1

        f.write(',motif Proteins, all proteins\n')
        for organism in filteredLengths:
            output = '{},{},{}\n'.format(organism, str(filteredLengths[organism]),
                                         str(rawLengths[organism]))
            f.write(output)


def writeCSV(filename, struct, q=False):
    headers = ['Frequencies of all amino acids at search position:',
    'Frequencies of all amino acids in the last {} positions (non-redundant):'.format(str(struct.numResidues)),
    'Frequencies of all amino acids at the given position (non-redundant)']

    dicts = [struct.freqMotifPositional,
             struct.freqNonRedMotifPositional,
             struct.freqAllPositional]

    with open(filename, 'w') as f:
        f.write('File name: {},,\n'.format(filename.split('/')[1]))

        # TODO delete this code
        lineToWrite = 'Requirements: '
        for item in struct.importantPositions:
            lineToWrite += 'P{}:{} '.format(str(item[0]), ''.join(str(item[x]) for x in range(1, len(item))))
        lineToWrite += ',,\n'
        f.write(lineToWrite)

        f.write('Comparison Position: {},,\n'.format(str(struct.searchPosition)))
        f.write('Timestamp: {},,\n\n'.format(datetime.now().strftime("%m/%d/%Y %H:%M")))
        f.write('Amino Acid,# Occurrences,Percentage\n')

        for header, aminoacidDict in zip(headers, dicts):
            f.write('\n{},,\n'.format(header))

            totalChars = sum(aminoacidDict.values())
            aminoacidDict = sorted(aminoacidDict.items(), key=lambda x:x[1], reverse=True)
            for a in aminoacidDict:
                if totalChars > 0:
                    percentString = str(round(a[1] / totalChars * 100, 4)) + '%'
                    f.write('{},{},{},\n'.format(a[0], str(a[1]), percentString))

        header = 'Enrichment ratios for each amino acid at given position'
        aminoacidDict = struct.enrichment

        f.write('\n{},,\n'.format(header))
        
        totalChars = sum(aminoacidDict.values())
        aminoacidDict = sorted(aminoacidDict.items(), key=lambda x:x[0], reverse=False)        
        for a in aminoacidDict:
            if totalChars > 0:
                f.write('{},{},,\n'.format(a[0], str(a[1])))


def writeSequenceLists(filename, struct):
    filename = 'sequenceLists/{}.tsv'.format(struct.organism.split('.')[0])
    with open(filename, 'w') as f:
        f.write('Sequence:\tIdentifier:\n')
        for i in range(0, len(struct.sequences)):
            f.write('{}\t{}\n'.format(struct.sequences[i], struct.sequenceLabels[i]))

def writeRawTSV(filename, proteinsLastN):
    writefile = 'rawTSV/{}.tsv'.format(filename.split('.')[0])
    with open(writefile, 'w') as f:
        f.write('Sequence:\tIdentifier:\n')
        for i in range(0, len(proteinsLastN)):
            f.write('{}\t{}\n'.format(proteinsLastN[i][0], proteinsLastN[i][1]))
