"""
extractCSV.py

With one or more provided csv files, read each line of each file and produce nested 
dictionaries containing all information. Data structure produced:
{all_f:  {organism: {aminoacid:[# Occurrences, Percentage]}}    // all positions freqs
 last_f: {organism: {aminoacid:[# Occurrences, Percentage]}}    // last 6 positions freqs
 pos_f:  {organism: {aminoacid:[# Occurrences, Percentage]}}    // current position freqs
 pos_e:  {organism: {aminoacid:[EnrichmentValue]}}              // current position enrichment
}

python3 extractCSV.py --files "bioParser/csv/*.csv" -outfile model.p
arguments:
--files FILEPATTERN     Uses glob to list all filenames matching a specified pattern
                        The pattern provided must have quotes around it
-outfile FILENAME       The filename to write condensed CSV information to in pickle format
-headerlines NUM        The number of header lines in the csv file if different from default
"""

import argparse
import glob
import os
import re
from os import path
import csv
import pickle

# Mapping the long CSV category names to shorter names for dictionary references
filelists = {"Frequencies of all amino acids in the last 6 positions (non-redundant):":'last_f', "Enrichment ratios for each amino acid at given position":'pos_e', "Frequencies of all amino acids at the given position (non-redundant)":'pos_f', "Frequencies of all amino acids at search position:":'all_f'}

# Mapping the value type to the category names (pos_e does not contain the Percentage category)
enrichment = ["# Occurrences"]
frequencies = ["# Occurrences", "Percentage"]
freqlists = {"pos_e":enrichment, "pos_f":frequencies, "last_f":frequencies, "all_f":frequencies}

"""
readDirectory
input: list of csv file names and the number of header lines to skip over
output: dictionary of organisms containing the four data cells (all_f, last_f, pos_e, and pos_e)
each data cell contains the Percentage/AminoAcidName/#Occurrences for each amino acid present
"""
def readDirectory(files, headerlines):
    data = {}
    for f in files:
        filesize = os.stat(f).st_size
       
        if os.path.isfile(f) and filesize > 50:   # Only read files that are 1kb or greater
            name, info = readFile(f, headerlines)
            
            if name != "":
                data[name.lower()] = info
    return data   

"""
readFile
input: csv file name and the number of header lines present in the file to skip over
output: organism name and the information for each of the four lists within the csv file
"""
def readFile(inputfile, headerlines):
    index = ""
    info = {}
    with open(inputfile) as csvfile:
        organism = csvfile.readline().split()[2][:-2][:-4]
        organism = '_'.join([organism.split('_')[0], organism.split('_')[1]])
        for i in range(headerlines-1): 
            csvfile.readline()
        reader = csv.DictReader(csvfile, delimiter=',')
        fields = reader.fieldnames
        for row in reader:
            if not row[fields[2]] and not row[fields[1]]:
                index = filelists[row[fields[0]]]
                
            else:
                info[index] = info[index] if index in info.keys() else []
                if index == "pos_e":
                    row = {k: v for k,v in row.items() if k != "Percentage"}
                info[index].append(row)
    return organism, info

"""
createDictionary
input: the dictionary produced with readDirectory
output: a new, rearranged dictionary in the following order: {cell:{organism:{aminoacid:[data]}}}
example: {'all_f':{'homo_sapiens':{'A':[<Occ>, <Freq>]}}, 'pos_e':{'homo_sapiens':{'A':[<Enri>]}}}
"""
def createDictionary(data):
    dicts = {}
    for identifier, entries in freqlists.items():
        all_amino = {aminoacid:{} for aminoacid in data.keys()}
        for organism, values in data.items():
            keys = list(values.keys())
            if len(keys) == 4: # Ensure all values are present
                for pos in values[identifier]:
                    all_amino[organism][pos["Amino Acid"]] = [pos[name] for name in entries]
        dicts[identifier] = all_amino
    return dicts

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', required=True, type=str,
        help="Path to csv files to read in. (e.g. 'csv/*.csv')")
    parser.add_argument('-outfile', type=str, default="model.p",
        help="File to write condensed enrichment CSV information to [default=model.p].")
    parser.add_argument('-headerlines', type=int, default=5,
        help="Number of header lines to skip over [default=5].")
    return parser.parse_args()

def findFiles(searchPara):
    fileList = []
    pathAndSearchPara = os.path.split(searchPara)
 
    path  = os.getcwd().replace("\\", "/")
    return [x.path for x in os.scandir(path + pathAndSearchPara[0]+"/") if x.name[-len(pathAndSearchPara[1][1:]):] == pathAndSearchPara[1][1:]]
    

if __name__ == "__main__":
    # Get all arguments from commandline
    args = parseArguments()
    files = findFiles(args.files)

    # Read each csv file provided and extract all information
    data = readDirectory(files, args.headerlines)
  
    # Store information in nested dictionary structure described above
    dicts = createDictionary(data)
    with open(args.outfile, 'wb') as f: 
        pickle.dump(dicts, f)
