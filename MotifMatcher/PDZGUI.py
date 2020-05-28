# written by Michael Lee

# 12/28/2018
# updated 1/1/2019 (bug fixes, comments)

# This program takes the tsv folder and xml file as arguments and provides a GUI to search for
# proteins by match num, ref sequence, other organism sequence, and Latin name.


'''
PDZGUI
python PDZGUI.py -tsvFolder /path/folder -xmlFile /path/folder/filename.xml
'''

import tkinter as tk
import os
import xml.etree.ElementTree as ET
import argparse

class Backend: # the fun stuff happens here...
    def __init__(self,all_seqs,organisms,root,matchValues,tsvs_folder,refOrg,\
                 org_tsv=None,ref_seq=None,org_seq=None,lat_name=None,\
                 console_print=None,file_write=None):
        
        # if a refseq was specified, find its location in the list of all refOrg sequences
        if ref_seq:
            hseq_location = all_seqs.index(ref_seq)

        # if a lat name was specified, find its location in the list of all orgs
        # and specify tsv path
        if lat_name:
            org_location = organisms.index(lat_name)
            tsv_address = tsvs_folder + "/" + org_tsv
        
        if ref_seq and lat_name:
            # find relevant org sequences from xml
            sequences = self.sequences_from_xml(hseq_location,org_location,root,matchValues)

            # find relevant ref proteins from tsv
            ref_matches = self.find_ref_proteins_tsv(tsvs_folder,[ref_seq],refOrg)

            # find relevant org proteins from tsv (based on xml data)
            org_matches = {i:[] for i in matchValues}
            for i in matchValues:
                org_matches[i] = self.find_org_protein(sequences[i],tsvs_folder,tsv_address)
                
            self.output(ref_matches,org_matches,console_print,file_write,ref_seq,org_seq,lat_name)
            
        elif org_seq and lat_name:
            # find relevant org sequences from tsv
            org_matches = self.find_org_protein([org_seq],tsvs_folder,tsv_address)

            # find relevant ref proteins from xml
            ref_proteins = self.find_ref_proteins_xml(matchValues,root,org_location,org_seq)

            # find relevant ref proteins from tsv (based on xml data)
            ref_matches = {i:[] for i in matchValues}
            for i in matchValues:
                ref_matches[i] = self.find_ref_proteins_tsv(tsvs_folder,ref_proteins[i],refOrg)
                
            self.output(ref_matches,org_matches,console_print,file_write,ref_seq,org_seq,lat_name)
            
        elif ref_seq and org_seq:
            # find relevant ref sequences from tsv
            ref_matches = self.find_ref_proteins_tsv(tsvs_folder,ref_seq,refOrg)

            # find relevant organisms from xml
            org_list = self.orgs_from_xml(hseq_location,org_seq,root)
            
            # make a list of all tsvs of all relevant orgs (based on xml)
            orgs = []
            for item in os.listdir(tsvs_folder):
                for org in org_list:
                    if org in item.lower():
                        orgs.append(item)

            # parse all relevant tsvs and put relevant lines into list (org_matches)
            org_matches = []
            for tsv in orgs:
                tsv_address = tsvs_folder + "/" + tsv
                matches = self.find_org_protein([org_seq],tsvs_folder,tsv_address)
                org_matches.extend(matches)
                
            self.output(ref_matches,org_matches,console_print,file_write,ref_seq,org_seq,lat_name)

    def none_to_string(self,a_string):
        # when ET creates a tree from an xml file, the contents of empty nodes are type None.
        # Empty strings are more tractable. If this function is passed a string, it returns the string.
        # If this function is passed None, it returns an empty string.
        return a_string or ''

    def sequences_from_xml(self,hseq_location,org_location,root,matchValues):
        # this function returns all org sequences for a given ref sequence and a given organism
        location = root[0][hseq_location][org_location]

        # the keys are the matchvalues of interest. The values are lists of sequences
        sequences = {i:self.none_to_string(location[6-i].text).split(",") for i in matchValues}
        return sequences

    def orgs_from_xml(self,hseq_location,org_seq,root):
        # this function parses the xml file and returns all orgs with a given sequence,
        # given a specific ref sequence
        ref_seq = root[0][hseq_location]
        org_list = []
        for nonRefOrganism in ref_seq:
            # make a list of all sequences for a given ref seq and given organism name
            for i in range(6):
                seq_list = []
                seq_string = nonRefOrganism[i].text
                seq_list.extend(self.none_to_string(seq_string).split(","))

            # if org_seq is in seq_list, add org name to org_list (in all lowercase,
            # without the '.tsv' extension on the end)
                if org_seq in seq_list:
                    org_list.append(list(nonRefOrganism.attrib.values())[0][:-4].lower())
        return org_list

    def find_ref_proteins_tsv(self,tsvs_folder,ref_seq_list,refOrg):
        # this function returns all proteins from the ref tsv whose sequence is in
        # ref_seq_list
        protein_matches = []
        filename = tsvs_folder + "/" + refOrg
        with open(filename,'r') as tsv_file:
            
            for line in tsv_file:
                if line.split("\t")[1].replace("\n","") in ref_seq_list:
                    protein_matches.append(line)
        return protein_matches

    def find_ref_proteins_xml(self,matchValues,root,org_location,org_seq):
        # given organism sequence and match values of interest, find ref sequences
        match_dict = {i:[] for i in matchValues}
        for i in matchValues:
            
            # in the xml file, matches are presented in descending order (i.e. 6-matches are presented first).
            match_num = -i+6
            
            for hum_seq in root[0]:
                protein_string = hum_seq[org_location][match_num].text # a string of all proteins for a given hum_seq, organism, and matchValue
                protein_list = self.none_to_string(protein_string).split(",") # same thing, in list form
                for protein in protein_list:
                    if protein == org_seq:
                        match_dict[i].append(list(hum_seq.attrib.values())[0])  # if the specified organism sequence is in protein_list, then
                                                                                # the ref protein (hum_seq) is interesting and should
                                                                                # eventually be outputted, so it's saved to a dictionary.
                        continue
        return match_dict

    def find_org_protein(self,org_seq_list,tsvs_folder,tsv_address):
        # for a given list of organism sequences, return lines from org tsv file that
        # have that organism sequence in them
        count=0
        protein_matches = []
        with open(tsv_address,'r') as tsv_file:
            org_seq_set = set(org_seq_list)
            for line in tsv_file:
                if line.replace("\n","").split("\t")[1] in org_seq_set: # sequence is second item in line
                    protein_matches.append(line)
        return protein_matches

    def output(self,ref_matches,org_matches,console_print,file_write,ref_seq,org_seq,lat_name):
        to_output = "Reference proteins matching the input criteria:\n"
        
        if type(ref_matches) is list:
            # if a ref sequence was provided in the initial input, then ref_matches is list.
            for match in ref_matches:
                ref_matches_pretty = [(match.split("\t")[1].replace("\n",""),match.split("\t")[0]) for match in ref_matches]
            for match in ref_matches_pretty:
                to_output+=match[0]+"\t"+match[1]+"\n\n"
                
        elif type(ref_matches) is dict:
            # if a ref sequence was NOT provided in the initial input, then ref_matches is dict.
            for match_num in ref_matches:
                to_output+="{0} match:\n".format(match_num)
                ref_matches_pretty = [(match.split("\t")[1].replace("\n",""),match.split("\t")[0]) for match in ref_matches[match_num]]
                for match in ref_matches_pretty:
                    to_output+=match[0]+"\t"+match[1]+"\n\n"
                if len(ref_matches_pretty) == 0:
                    to_output+="[no matches]\n\n"
                    
        to_output+="\n\nother proteins matching the input criteria:\n\n"
        if type(org_matches) is dict:
            # if an org sequence was NOT provided in the initial input, then org_matches is dict.
            for match_num in org_matches:
                to_output+="{0} matches:\n".format(match_num)
                org_matches_pretty = [(o.split("\t")[1].replace("\n",""),o.split("\t")[0]) for o in org_matches[match_num]]
                for org in org_matches_pretty:
                    to_output+=org[0]+"\t"+org[1]+"\n\n"
                if len(org_matches_pretty) == 0:
                    to_output+="[no matches]\n\n"
                    
        elif type(org_matches) is list:
            # if an org sequence was provided in the initial input, then org_matches is list.
            org_matches_pretty = [(o.split("\t")[1].replace("\n",""),o.split("\t")[0]) for o in org_matches]
            for match in org_matches_pretty:
                to_output+=match[0]+"\t"+match[1]+"\n\n"
                
        to_output+="\n\n--end of data printout--"
        
        if console_print:
            print(to_output)
        if file_write:
            file_name = "{0} {1} {2}.txt".format(ref_seq,org_seq,lat_name[:-4])
            with open(file_name,'w') as file_object:
                file_object.write(to_output)
                print("--task completed--")


class GUI:
    def __init__(self,master,refOrg,tsvs_folder,xml_path):
    
        # drawing window
        self.master = master
        master.title("protein retrieval tool")

        self.instructions = tk.Label(master, text="""Fill data into two of the following\
 three boxes and select
 "write to file" or "print to console".

 If you provide a Latin name, it's fine if you just put in the first
 few letters (but be unambiguous!)""")
        self.instructions.grid(columnspan=8,sticky=tk.W)

        # Labels
        self.ref_seq_label = tk.Label(master, text="Reference sequence:")
        self.ref_seq_label.grid(columnspan=5,sticky = tk.W)

        self.org_seq_label = tk.Label(master, text="other org sequence:   ")
        self.org_seq_label.grid(columnspan=5,sticky = tk.W)

        self.lat_name_label = tk.Label(master,text="Lat. name:")
        self.lat_name_label.grid(columnspan=3,sticky=tk.W)

        # Entry boxes (Human sequence, org sequence, Lat name)
        self.ref_entry = tk.Entry(master,width=15)
        self.ref_entry.grid(columnspan=4,row=1,column=4,sticky=tk.W)

        self.org_entry = tk.Entry(master,width=15)
        self.org_entry.grid(columnspan=4,row=2,column=4,sticky=tk.W)

        self.lat_entry = tk.Entry(master,width=30)
        self.lat_entry.grid(columnspan=7,row=3,column=4)

        # Vars for checkboxes (1-6)
        self.var6 = tk.IntVar()
        self.var5 = tk.IntVar()
        self.var4 = tk.IntVar()
        self.var3 = tk.IntVar()
        self.var2 = tk.IntVar()
        self.var1 = tk.IntVar()

        # Checkboxes (1-6)
        self.button6 = tk.Checkbutton(master,text="6",variable=self.var6)
        self.button6.grid(row=4,column=0)
        self.button6.select()

        self.button5 = tk.Checkbutton(master,text="5",variable=self.var5)
        self.button5.grid(row=4,column=1)
        self.button5.select()

        self.button4 = tk.Checkbutton(master,text="4",variable=self.var4)
        self.button4.grid(row=4,column=2)
        self.button4.select()

        self.button3 = tk.Checkbutton(master,text="3",variable=self.var3)
        self.button3.grid(row=4,column=3)

        self.button2 = tk.Checkbutton(master,text="2",variable=self.var2)
        self.button2.grid(row=4,column=4)

        self.button1 = tk.Checkbutton(master,text="1",variable=self.var1)
        self.button1.grid(row=4,column=5)

        self.match_label = tk.Label(master, text="matches")
        self.match_label.grid(row=4,column=6)

        # vars for checkboxes ("write to file", "print to console")
        self.file_write = tk.IntVar()
        self.console_print = tk.IntVar()

        # Checkboxes ("write to file", "print to console")
        self.file_write_button = tk.Checkbutton(master,text="write to file",variable=self.file_write)
        self.file_write_button.grid(columnspan=3,row=5)

        self.console_print_button = tk.Checkbutton(master,text="print to console",variable=self.console_print)
        self.console_print_button.grid(columnspan=3,column=3,row=5)
        self.console_print_button.select()

        # find tsv and xml files
        self.refOrg =refOrg
        self.tsvs_folder = tsvs_folder 
        self.xml_path  = xml_path
    
            
        # parse xml file
        root,all_seqs,organisms = self.parse_xml(xml_path)

        # submit button
        self.submit = tk.Button(master,text="submit",command=lambda:self.get_entries_and_run(organisms,tsvs_folder,refOrg,all_seqs,root))
        self.submit.grid(columnspan=6,row=7)

    def parse_xml(self,xml_path):
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # The refsequences are the children of root[0].
        # All of the refsequences in the xml are gathered into a list.
        all_seqs = [list(child.attrib.values())[0] for child in root[0]]
        
        # root[0] contains a comma- and semicolon- separated string of all nonref orgs.
        # A list is created with the names of all nonref orgs.
        organisms = list(root[0].attrib.values())[0].split(";")
        organisms = [item.split(",") for item in organisms]
        organisms = [organisms[x][0] for x in range(1,len(organisms))] # start with 1 because homo_sapiens is in position 0
        organisms = [org.lower() for org in organisms] # lowercase is better
        return root,all_seqs,organisms

    def get_entries_and_run(self,organisms,tsvs_folder,refOrg,all_seqs,root):
        # gets data from input boxes as strings. Empty boxes become blank strings "".
        self.ref_seq = self.ref_entry.get()
        self.org_seq = self.org_entry.get()
        self.lat_name = self.lat_entry.get()
        self.lat_name = self.lat_name.lower().strip().replace(" ","_") # lowercase is better. Spaces are bad. Underscores are friendly.
        
        if len(self.lat_name)!=0: # if there's something in the box...
            org_tsv = [file for file in os.listdir(tsvs_folder) if self.lat_name.lower() in file.lower()]
            if len(org_tsv)>1:
                raise RuntimeError("More than one organism tsv file with this Latin name!")
            elif len(org_tsv)==0:
                raise RuntimeError("No org tsv file with this Latin name!")
            else:
                org_tsv=org_tsv[0]
                self.lat_name = [org for org in organisms if self.lat_name in org]
                if len(self.lat_name)>1:
                    raise RuntimeError("More than one org in xml file with this Latin name!")
                elif len(self.lat_name)==0:
                    raise RuntimeError("This org not in xml file!")
                else:
                    self.lat_name=self.lat_name[0]
                
        matchValues = []
        print_to_console,write_to_file = [bool(val) for val in [self.console_print.get(),self.file_write.get()]]
        
        if not print_to_console and not write_to_file:
            raise RuntimeError("No output method specified!")

        # populates matchValues with the user-specified match values of interest
        button_data = [self.var6.get(),self.var5.get(),self.var4.get(),self.var3.get(),self.var2.get(),self.var1.get()]
        for i in range(6):
            if button_data[i] == 1:
                matchValues.append(-i+6)
                
        if self.check_for_valid(): # if user input is ok...
            if self.ref_seq != "" and self.lat_name != "":
                computation = Backend(all_seqs,organisms,root,matchValues,tsvs_folder,refOrg,\
                                      org_tsv,self.ref_seq,lat_name=self.lat_name,\
                                      console_print=print_to_console,file_write=write_to_file)
            elif self.org_seq != "" and self.lat_name != "":
                computation = Backend(all_seqs,organisms,root,matchValues,tsvs_folder,refOrg,\
                                      org_tsv,org_seq=self.org_seq,lat_name=self.lat_name,console_print=print_to_console,\
                                      file_write=write_to_file)
            elif self.ref_seq != "" and self.org_seq != "":
                computation = Backend(all_seqs,organisms,root,matchValues,tsvs_folder,refOrg,\
                                      ref_seq=self.ref_seq,org_seq=self.org_seq,\
                                      console_print=print_to_console,file_write=write_to_file)
            else:
                pass # the program shouldn't ever get here, because check_for_valid would have raised an error
        else:
            pass # the program shouldn't ever get here because check_for_valid raises an error if False

    def check_for_valid(self):
        amino_acids = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
        refSeqOk = -1
        orgSeqOk = -1
        latNameOk = -1
        error_message = "" # stays blank unless user input is invalid
        
        if len(self.ref_seq) == 6 and all(x.upper() in amino_acids for x in self.ref_seq):
            refSeqOk = True
            self.ref_seq = self.ref_seq.upper()
        elif self.ref_seq == "":
            pass
        else:
            refSeqOk = False
            error_message+="Human sequence not OK. "
        if len(self.org_seq) == 6 and all(x.upper() in amino_acids for x in self.org_seq):
            orgSeqOk = True
            self.org_seq = self.org_seq.upper()
        elif self.org_seq == "":
            pass
        else:
            orgSeqOk = False
            error_message+="Org sequence not ok."

        # lat_name already checked in get_entries_and_run, so no checks needed here
        if self.lat_name != "":
            latNameOk = True
        else:
            pass
        
        if all([refSeqOk,orgSeqOk,latNameOk]): # i.e. if none of them are False...
            if sum([refSeqOk,orgSeqOk,latNameOk])==1:# i.e. if exactly two boxes are filled (1+1-1 == 1)
                return True
            else:
                raise RuntimeError("Make sure that exactly two boxes are filled!")
        else:
            raise ValueError(error_message)
def PDZGUI_wrapper(tsvs_folder, xml_path):
    with open(xml_path,'rb') as f:
        tree = ET.iterparse(f,events=("start", "end"))
        for event,element in tree:
            if element.tag == "Summary":
                motifRuningSumAndProteinTotals = element.attrib['name']
                refOrg = [protein.split(',') for protein in motifRuningSumAndProteinTotals.split(";")][0][0]
                break
    
    source = tk.Tk()
    source.geometry("350x240")
    source.grid_rowconfigure(7,minsize=40)
    protein_gui = GUI(source, refOrg, tsvs_folder, xml_path)
    source.mainloop()
    
def parseArgs():
    parser = argparse.ArgumentParser(description="protein retrieval")
    parser.add_argument('-tsvFolder',required=True,type=str,
                        help="Provide address of folder with tsv files.\
Usage: /path/folder")
    parser.add_argument('-xmlFile',required=True,type=str,
                        help="Provide address of xml file.\
Usage: /path/folder/file.xml")
    return parser.parse_args()

if __name__ == "__main__":
    args=parseArgs()
    tsvs_folder=os.path.normpath(args.tsvFolder)
    xml_path=os.path.normpath(args.xmlFile)
    if not os.path.exists(tsvs_folder):
        raise FileNotFoundError("tsv folder location invalid.")
    if not os.path.exists(xml_path):
        raise FileNotFoundError("xml file location invalid.")
    PDZGUI_wrapper(tsvs_folder, xml_path)
    
