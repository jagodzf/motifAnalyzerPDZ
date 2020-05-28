import argparse
import time
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from fileOutput import parseFileNames, makeFolders, writeSummaryFile

def distributeWork(positions, fileNames, motifs, numResidues, q=False):
    ''' 
    Calls setUpHandler.py for each position and organism supplied with the motifs and numResidues.
    The threads are run in parallel to reduce time requirements.
    '''
    motifs = ' '.join(motifs).upper()

    # array of tasks to be completed
    futures = []

    # open thread pool
    with ThreadPoolExecutor(max_workers=8) as executor:
        for x in positions:
            for i in fileNames:
                #submit tasks to be completed by threads
                silence = '-q' if q else ''
                setup = ' setUpHandler.py {} {} {} 2 -o {} {}'.format(
                        i, str(x), str(numResidues), silence, motifs)
                futures.append(executor.submit(os.system, (sys.executable + setup)))

def createHeatmaps(organisms, motifID, positions):
    '''
    Grabs the generated csv files and runs both extractCSV.py and createHeatmap.py.
    Saves the generated heat map pngs.
    '''
    for p in positions:
        infiles = '/csv/*{}.csv'.format(p)
        
        pickle = '{}/position{}.p'.format(motifID, p)
        title = '{}_position{}'.format(motifID, p)
        outfile = '{}/{}.png'.format(motifID, title)

        extract = ' extractCSV.py --files {} -outfile {}'.format(infiles, pickle)
        create = ' createHeatMap.py --enrichment {} --out {} -title {} -organisms {}'.format(
            pickle, outfile, title, organisms)

        os.system(sys.executable + extract)
        os.system(sys.executable + create)


def writeData(read, write, index):
    with open(read, 'r') as fr, open(write, 'a') as fw: 
        lines = fr.readlines()
        if index < len(lines):
            line = lines[index]
            fw.write(line.split('\n')[0])
            if len(line.strip('\n')) > 0:
                fw.write(',,')
            else:
                fw.write(',,,,')
        else:
            fw.write(',,,,')

def combineData(filenames, positions):
    '''
    Creates combinedCSVs directory and reads in the data from the separate csv files. 
    '''
    folderPos = os.getcwd()
    directory = '{}/combinedCSVs'.format(folderPos)
    if not os.path.exists(directory):
        os.makedirs(directory)

    for pos in positions.split(','):
        write_file = '{}/combined{}.csv'.format(directory, pos)
        open(write_file, 'w').close()
        for index in range(115):
            for org_file in filenames:
                read_file = '{}/csv/{}{}.csv'.format(folderPos, org_file.split('.')[0], pos)
                if os.path.exists(read_file):
                    writeData(read_file, write_file, index)
            with open(write_file, 'a') as f:
                f.write('\n')


def checkInput(args):
    '''
    Ensure the arguments entered by the user are valid, exit the program otherwise
    '''
    assert os.path.exists(args.organisms), 'Could not find file: {}'.format(args.organisms)

    pos = args.positions.split(',')
    for p in pos:
        assert p.isdigit(), 'Invalid position entered: {}'.format(p)

    for mot in args.motifs:
        m = mot.split(':')
        assert len(m) == 2, 'Invalid motif entered: {}'.format(mot)
        assert m[1] != '', 'Invalid motif entered: {}'.format(mot)
        assert m[0][0].lower() == 'p', 'Invalid motif entered: {}'.format(mot)
        assert m[0][1:].isdigit(), 'Invalid motif entered: {}'.format(mot)


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Run C-terminal decameric sequence processing on many files simultaneously')
    parser.add_argument('-organisms', required=True, type=str,
        help='The text file containing latin names of organisms. (e.g. organisms.txt)')
    parser.add_argument('-positions', required=True, type=str,
        help='The positions to search over, delimited with commas. (e.g. 1,3,4,5)')
    parser.add_argument('-numResidues', required=True, type=int,
        help='The number of residues to provide statistics on. (e.g. 6)')
    parser.add_argument('-motifs', nargs='+',
        help='The matching positions with desired amino acids. (e.g. P0:ILVF P2:ST ... PX:X)')
    parser.add_argument('-c', action='store_true',
        help='Add this flag to create combined CSVs.')
    parser.add_argument('-q', action='store_true',
        help='Add this flag to quiet the command line arguments')
    parser.add_argument('-heatmaps', action='store_true',
        help='Add this flag to make enrichment heat maps for csv data.')
    parser.add_argument('-motifID', type=str, default='motif',
        help='If heat maps are created, use this naming convention. [default=motif]')
    return parser.parse_args()


def main():
    args = parseArgs()
    checkInput(args)
    start = time.time()
    positions = args.positions.split(',')

    makeFolders()
    with open('unfoundOrganisms.txt', 'w') as f:
        f.write('Organisms not found in uniprot:\n\n')
    
    organism_filenames = parseFileNames(args.organisms)
    distributeWork(positions, organism_filenames, args.motifs, args.numResidues, q=args.q)

    writeSummaryFile('rawTSV/', 'sequenceLists/', 'summaryFile.csv', q=args.q)

    # if the user supplies -heatmap flag, create heatmaps for each position provided
    if args.heatmaps:
        print('Create heat maps...')
        if not os.path.exists(args.motifID):
            os.makedirs(args.motifID)
        createHeatmaps(args.organisms, args.motifID, positions)

    # if the user supplies -c flag, combine CSVs with combineData.py
    if args.c:
        print('Combining CSVs...')
        combineData(organism_filenames, args.positions)
    print('All tasks completed in {} seconds'.format(str(round(time.time() - start, 2))))

if __name__ == '__main__':
    main()

