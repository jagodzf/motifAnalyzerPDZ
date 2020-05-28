import ctypes
import os
import time
import sys
def callMakeXML(inTSVFolder,motif_type,motifLen,outFile,refOrganism):
    cwd = os.getcwd()
    OpSys= sys.platform
    if OpSys == "win32":
        findFunc = ctypes.CDLL(cwd + "/lib/makeXML.dll")
    elif OpSys == "linux2":
        findFunc = ctypes.CDLL(cwd + "/lib/makeXML.so")
    elif OpSys == "linux":
        findFunc = ctypes.CDLL(cwd + "/lib/makeXML.so")
    else:
        raise AttributeError("unrecognized operating system")
    function = findFunc.makeXML
    function.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_char_p)),ctypes.POINTER(ctypes.c_int),ctypes.c_int, ctypes.POINTER(ctypes.c_char_p),ctypes.c_int,ctypes.c_int,ctypes.c_char_p,ctypes.POINTER(ctypes.c_char_p),ctypes.c_char_p]
    # takes an input folder and a reference organism computes all of the matches between
    # the reference organism motifs and the motifs of the reference organisms
    # takes an input folder and a reference organism computes all of the matches between
    # the reference organism motifs and the motifs of the reference organisms
    allEndings = {}
    refEndings = set()
    # gets all of the unique protein endings for each organism in the tsv folder saves them to a dictionary
    for file in os.listdir(inTSVFolder):
        with open(inTSVFolder + file, 'r') as readFile:
            allEndings[file.encode("utf-8")] = set()
            
            for line in readFile:
                motif = line.split('\t')[1].replace('\n','')
                if (len(motif) == motifLen) and all(motif[key] in motif_type[key] for key in motif_type):
                    allEndings[file.encode("utf-8")].add(motif.encode("utf-8"))

            
    # removes the reference organism from the rest of the organisms
    refEndings = refEndings | allEndings[refOrganism.encode("utf-8")]
    del allEndings[refOrganism.encode("utf-8")]
    #convert python data structures to c data structures 
    orgList = []
    EntireArray = []
    EALen = []
    for org in allEndings:
        orgList.append(org)
        
        EntireArray.append(list(allEndings[org]))
        EALen.append(len(allEndings[org]))
    
    P_ARRAY_ARRAY_CHAR = ctypes.POINTER(ctypes.c_char_p)*len(orgList)
    ptrEA = P_ARRAY_ARRAY_CHAR()
    for org in range(len(orgList)):
        ptrEA[org] = (ctypes.c_char_p*len(EntireArray[org]))()
        for j in range(len(EntireArray[org])):
            ptrEA[org][j] = EntireArray[org][j]
    
    
    EntireArrayLen = ((ctypes.c_int)*len(EALen))()
    for j in range(len(EALen)):
        EntireArrayLen[j] = EALen[j]
    
    refArray = (ctypes.c_char_p*len(refEndings))()
    refArray[:] = list(refEndings)
    organisms = (ctypes.c_char_p*len(orgList))()
    organisms[:] = orgList
    summary = (refOrganism.encode("utf-8")+b","+str(len(refEndings)).encode("utf+8")+b";" +b";".join([b",".join([key,str(len(value)).encode("utf-8")]) for key, value in allEndings.items()]))
    #pass arguments to to makeXML
    function(ptrEA,EntireArrayLen,len(orgList),refArray,len(refEndings),motifLen,outFile.encode("utf-8"),organisms,summary)
    
