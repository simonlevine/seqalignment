
import sys
import numpy as np

#python NW.py example.fasta to run

def ReadFASTA(filename):
    fp=open(filename, 'r')
    Sequences={}
    tmpname=""
    tmpseq=""
    for line in fp:
        if line[0]==">":
            if len(tmpseq)!=0:
                Sequences[tmpname]=tmpseq
            tmpname=line.strip().split()[0][1:]
            tmpseq=""
        else:
            tmpseq+=line.strip()
    Sequences[tmpname]=tmpseq
    fp.close()
    return Sequences

# helper functions for Needleman-Wunsch algorithm here



def score_chars(char1,char2): 
    if char1 == char2:
        return 0  # match
    elif char1 == '-' or char2 == '-':
        return -1 #gap w/ one character
    else:
        return -2 #mismatch


def F_matrix(seq1,seq2):
    #mismatch = -2, gap = -1, match = 0
    F=np.empty([len(seq1)+1,len(seq2)+1])
    
    d = -1 #gap cost
    for i in range(0,len(seq1)+1):
        F[i][0] = d*i
    for j in range(0,len(seq2)+1):
        F[0][j] = d*j


    for i in range(1, len(seq1)+1):
        for j in range(1,len(seq2)+1):
            match = F[i-1][j-1] + score_chars(seq1[i-1],seq2[j-1])
             #score: cost of matching characters: 0 if same, -2 if mismatch
            delete = F[i-1][j] + d #gap cost
            insert = F[i][j-1] + d

            F[i][j] = max(match,insert,delete)
    print(F)
    return F


def traceback(F, seq1, seq2):
    
    d = -1
    
    i = len(seq1)
    j = len(seq2)
    Alignment1 = ""
    Alignment2 = ""


    score = 0
    # score = score_chars(seq1[len(seq1)-1],seq2[len(seq2)-1])

    while (i>0 or j>0): #stay in matrix bounds

        

        if F[i][j] == F[i-1][j-1] + score_chars(seq1[i-1], seq2[j-1]) :
            #match or mismatch case
            Alignment1 = seq1[i-1] + Alignment1
            Alignment2 = seq2[j-1] + Alignment2

            score += score_chars(seq1[i-1], seq2[j-1])
            # print(i,j)
            i += -1
            j += -1


        elif ((i > 0) and (F[i][j] == F[i-1][j] + d)):
            #gap case, seq2 
            Alignment1 = seq1[i-1] + Alignment1
            Alignment2 = "-" + Alignment2 #gap character
            score += -1
            # print(i,j)
            i += -1

        else: #gap for seq2
            Alignment1 = "-" + Alignment1
            Alignment2 = seq2[j-1] + Alignment2
            score += -1
            # print(i,j)
            j += -1
    
    #finally, score the alignment
    # score = 99
    return Alignment1, Alignment2, score



# main function signature
def needleman_wunsch(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: a string of three lines like so:
    '<alignment score>\n<alignment in seq1>\n<alignment in seq2>'
    """
    F = F_matrix(seq1,seq2)
    A1,A2,score = (traceback(F, seq1, seq2))

    # out_str = "strscore}\n{A1}\n{A2}"
    out_str = str(score)+"\n"+str(A1)+"\n"+str(A2)
    return out_str

if __name__=="__main__":
    Sequences=ReadFASTA(sys.argv[1])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq2=Sequences[list(Sequences.keys())[0]]
    seq1=Sequences[list(Sequences.keys())[1]]
    # print(needleman_wunsch('AGCTA', 'AGGTCA'))
    print(needleman_wunsch(seq1,seq2))

