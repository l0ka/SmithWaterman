#!/usr/bin/env python
# alessio.locallo@studenti.unitn.it

###############################################################################
## Building blocks
# seqA, seqB                --> 1. input sequences
# match, mismatch, gap      --> 2. scoring function
# rows, columns             --> 3. two matrices: scores and paths 
# highest score             --> 4. starting position for traceback
# traceback                 --> 5. alignment(s)

###############################################################################
## Get input sequences
import sys
if len(sys.argv) == 3:
    seq1, seq2 = sys.argv[1], sys.argv[2]
else:
    print 'Please provide exactly 2 DNA sequences as input.\n', 'Exit!\n'
    sys.exit()

##  Check if the input strings are valid DNA sequences
valid = 'ATCG'
if any(i not in valid for i in seq1 + seq2):
    print 'Please provide 2 valid DNA sequences as input.\n', 'Exit!\n'
    sys.exit()

## Define scoring function 
match    =  3
mismatch = -3
gap      = -2

###############################################################################
## Implement Smith-Waterman algorithm
def s_w(seqA, seqB):
    ## Initialize variables and matrices
    cols      = len(seqA)
    rows      = len(seqB)
    matrix    = [[0 for row in range(rows+1)] for col in range(cols+1)]
    paths     = [[0 for row in range(rows+1)] for col in range(cols+1)]
    max_score = 0
    s1, s2 = [], []
    ## Fill the scoring matrix
    for i in range(cols):
        for j in range(rows):
            if seqA[i] == seqB[j]:
                diag = matrix[i][j] + match
            else:
                diag = matrix[i][j] + mismatch
            up    = matrix[i + 1][j] + gap
            left  = matrix[i][j + 1] + gap
            score = max(0,diag, up, left)
            matrix[i+1][j+1] = score
            if score > max_score:
                max_score = score
                start_pos = [i+1, j+1]
    ## Fill the paths matrix
            if matrix[i+1][j+1]   == diag and matrix[i+1][j+1] != 0:
                paths[i+1][j+1] = 'diag'
            elif matrix[i+1][j+1] == up   and matrix[i+1][j+1] != 0:
                paths[i+1][j+1] = 'up'
            elif matrix[i+1][j+1] == left and matrix[i+1][j+1] != 0:
                paths[i+1][j+1] = 'left'
    ## Traceback
    i, j = start_pos
    start_path = paths[i][j]
    while start_path != 0:
        if start_path == 'diag':
            s1.append(seqA[i-1])
            s2.append(seqB[j-1])
            i, j = i-1, j-1
        elif start_path == 'up':
            s1.append('-')
            s2.append(seqB[j-1])
            j = j-1
        else:
            s1.append(seqA[i-1])
            s2.append('-')
            i = i-1
        start_path = paths[i][j]
    ## Return the optimal local alignment
    return ''.join(list(reversed(s1))), ''.join(list(reversed(s2)))

###############################################################################    
## Run the function and print the result
aln1, aln2 = s_w(seq1, seq2)
print 'Optimal local alignment\n', aln1, '\n', aln2, '\n'
