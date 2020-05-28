# PDZBioParser
# BioParser

BioParser is a tool used to process C-terminal decameric sequences,
with the evaluation of PDZ binding domains as is primary
goal. BioParser will filter and analyze the proteomes of as many
organisms as are supplied to it, and will filter based off as many
motifs as are supplied to it.

Motifs refers to specifying specific amino acids at certain
positions. For example, if you wanted to analyze the enrichment of the
amino acids S and T at position 2, you could supply the -motifs
argument with ``` P2:ST ```

BioParser also creates detailed data on the proteomes that is stored
in the CSV folder that is created during the programs execution. The
-c argument combines all CSV's for one organism, which can be usefull
when large numbers of positions are being processed.

# MotifMatcher
# Reqirements
Python3

# Usage
```
main.py [-h] -refOrganism REFORGANISM -sequenceFolder SEQUENCEFOLDER
               -numResidues NUMRESIDUES -inputValues INPUTVALUES
               [INPUTVALUES ...] -out OUT -significance SIGNIFICANCE
Arguments:
               
   -refOrganism REFORGANISM
                        Provide the file name of the organism that the other
                        organisms will be compared with. Usage:
                        /path/folder/organism.tsv
  -sequenceFolder SEQUENCEFOLDER
                        The path of the folder that contains tsv files that
                        contain the last N many residues of all the proteins
                        in the proteomes in question. Usage: /path/folder
  -numResidues NUMRESIDUES
                        The number of residues to provide statistics on.
                        Usage: 6
  -inputValues INPUTVALUES [INPUTVALUES ...]
                        The matching positions with desired amino acids.
                        Usage: P0:ILVF P2:ST ... PX:X
  -out OUT              Output file directory and name. Do not include file
                        extensions. Usage: /path/folder/fileName_motif_example
  -significance SIGNIFICANCE
                        significance cut off value example 0.9999 corresponds
                        to P<0.0001 or 0.1 percent chance the item was flagged
                        in error
  -decouple             If a xml file has all ready been generated use this
                        flag to skip xml generation
  -pyth                 force the use of motifFileMaker over makeXML. Using
                        this flag will slowdown exicution.
  -GUI                  Launch GUI automatically after program exicution.

     
Example Usage:
Linux
python3 main.py -refOrganism homo_sapiens_TaxID_9606_R_reviewed.tsv -sequenceFolder rawTSV/ -numResidues 6 -inputValues P0:WYIAFLVM P2:ST -out motif_2_human -significance 0.9999
Windows
main.py -refOrganism homo_sapiens_TaxID_9606_R_reviewed.tsv -sequenceFolder rawTSV/ -numResidues 6 -inputValues P0:WYIAFLVM P2:ST -out motif_2_human -significance 0.9999

```
