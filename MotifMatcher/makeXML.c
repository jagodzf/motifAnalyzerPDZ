#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void makeXML(const char *** EntireArray,int * EntireArrayLen, int EntireLen, const char ** refArray,int refLen, int motifLen,const char * outFile,const char ** organisms, const char * summary){
    FILE * writeFile = fopen(outFile,"w");
    fprintf(writeFile,"<root>\n");
    fprintf(writeFile,"\t<Summary name=\"%s\">\n",summary);
    int maxSoFar = 0;
    for (int i = 0; i < EntireLen;i++){
        if (EntireArrayLen[i]>maxSoFar ){
            maxSoFar = EntireArrayLen[i];
        }
    }
    
    int score;
    int l;
    char scoreArray[motifLen+1][maxSoFar*(motifLen+1)+1];
    
    int indexArray[motifLen+1];
    /* initialize scoreArray and indexArray*/
    for (l =0; l<motifLen;l++){
            
                for (int s =0 ;s<motifLen+1; s++){
                    indexArray[s] = 0;
                    memset(scoreArray[s],'\0',sizeof(char)*(maxSoFar*(motifLen+1)+1));
                    
                }
            }
    /*begin main loog*/
    for (int g =0;g<refLen;g++){
        
        fprintf(writeFile,"\t<RefSequence name=\"%s\">\n",refArray[g]);
        
        
        for (int i = 0; i < EntireLen;i++){
            
            for (int k = 0; k < EntireArrayLen[i];k++){
                /* calculate score*/
                score = 0;
                for (int j =0;j<motifLen;j++){
                    if (EntireArray[i][k][j] == refArray[g][j]){
                        score++;
                    }
                }
                /* sort the array into buckets based on the scores*/
                if (score != 0){
                    for (l =0; l<motifLen;l++){
                        
                        scoreArray[score][indexArray[score]] = EntireArray[i][k][l];
                        indexArray[score]++;
                        
                    }
                    scoreArray[score][indexArray[score]] = ',';
                    indexArray[score]++;
            }
            }
            /*write the sequences to the file */
            
            fprintf(writeFile,"\t\t\t<NonRefOrganism name=\"%s\">\n",organisms[i]);
            for (int s = motifLen;s >0;s-- ){
                if (scoreArray[s][0] != '0'){
                    scoreArray[s][indexArray[s]-1]  = '\0';
                    fprintf(writeFile,"\t\t\t\t<match%d>%s</match%d>\n",s,scoreArray[s],s);
                }
                
                
            }
            fprintf(writeFile,"\t\t\t</NonRefOrganism>\n");
            /* clear memory for next cycle*/
            for (l =0; l<motifLen;l++){
            
                for (int s =0 ;s<motifLen+1; s++){
                    
                    memset(scoreArray[s],'\0',sizeof(char)*(indexArray[s]+1));
                    indexArray[s] = 0;
                }
            }
        }
        fprintf(writeFile,"\t\t</RefSequence>\n");
    }
    fprintf(writeFile,"\t</Summary>\n");
    fprintf(writeFile,"</root>\n");
    fclose(writeFile);
}
