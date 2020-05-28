import argparse
import sys
import motifFileMaker
import motifCounterXMLParse
import os
import time
import makeXMLWrapper
import PDZGUI
import motifIsolation
def autoCompleteFileName(fileName,folder):
    for file in os.listdir(folder):
        if fileName in file:
            return file
    raise FileNotFoundError

def parseArgs():
    parser = argparse.ArgumentParser(
        description='TODO')
    parser.add_argument('-refOrganism', required=True, type=str,
                         help = "Provide the file name of the organism that the other organisms will \
                            be compared with. Usage: /path/folder/organism.tsv")
    parser.add_argument('-sequenceFolder',required=True, type=str,
                        help="The path of the folder that contains tsv files that \
                        contain the last N many residues of all the proteins in the proteomes in question.\
                        Usage: /path/folder")
    parser.add_argument('-numResidues', required=True, type=int, 
                        help='The number of residues to provide statistics on. Usage: 6')
    parser.add_argument('-motifs', required=True,nargs='+',
                        help='The matching positions with desired amino acids. Usage: P0:ILVF P2:ST ... PX:X')
    parser.add_argument('-out', required=True, type=str,
                        help='out put file directory')
    parser.add_argument('-significance', required=True, type=float,
                        help='significance cut off value example 0.95 corresponds to P<0.05 or 5 percent chance the item was flagged in error')
    parser.add_argument("-decouple", action='store_true',required = False, default = False,
                        help = "If a xml file has all ready been generated use  this flag to skip xml generation")
    parser.add_argument("-pyth", action='store_false',required = False, default = True,
                        help = "force the use of motifFileMaker over makeXML. Using this flag will slowdown exicution.")
    parser.add_argument("-GUI", action='store_true',required = False, default = False,
                        help = "Launch GUI automatically after program exicution.")
    parser.add_argument("-full", action='store_true',required = False, default = False,
                        help = "Allow the program to identify motifs in the internals of the protein.")

    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs()
    start = time.clock()
    folderName = os.path.normpath(args.sequenceFolder)
    motifName ='-'.join(args.motifs).replace(':','_')
    motif = {(args.numResidues-1-int(k[0])):set(list(k[1])) for k in[i.replace('p','').split(':') for i in args.motifs]}
    if args.full:
        file_path = str(os.getcwd().replace("\\", "/")) + ("/%s" % "motifTSV")
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        motif1 = [set(list("QWERTYIPASDFGHKLXCVNM")) for i in range(args.numResidues)]
        for item in motif:
            motif1[item] = motif[item]
        motifIsolation.openFastaInFolder(folderName,motif1,args.numResidues)
        folderName = "motifTSV/"
    xmlFileName = args.out +".xml"
    tsvFileName = args.out +".tsv"
    refOrgFile = autoCompleteFileName(args.refOrganism,folderName)
    if args.decouple:
        motifCounterXMLParse.generateOutput(xmlFileName,tsvFileName,motif,args.numResidues,args.significance)
    else:
        # get the operating system
        OpSys = sys.platform
        # decide whether to using the makeXML or the python motifFileMaker
        if OpSys in ["win32","linux2","linux"] and args.pyth:
            print('make')
            makeXMLWrapper.callMakeXML(folderName+'/',motif,args.numResidues,xmlFileName,refOrgFile)
        else:
            print("pyth")
            motifFileMaker.motif_Finder(folderName+'/',motif,args.numResidues,xmlFileName,refOrgFile)
        motifCounterXMLParse.generateOutput(xmlFileName,tsvFileName,motif,args.numResidues,args.significance)
    print("completed in: " +str(time.clock()-start))
    if args.GUI:
        PDZGUI.PDZGUI_wrapper(folderName + '/',xmlFileName)

