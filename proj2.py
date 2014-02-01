#!/bin/python

from math import *
import sys, re

GAP='-'
GAP_SCORE = -5

# Assumes only two fasta entries
def read_fasta(fasta):
    a=''
    a_name=''
    b=''
    b_name=''
    state=0
    for line in open(fasta, 'r'):
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        if 0 == state:
            if None == re.search('^>', line):
                print('error: not valid fasta format')
            a_name = line[1:]
            state = state + 1
        elif 1 == state:
            if None == re.search('^>', line):
                a = a + line
            else:
                state = state + 1
                b_name = line[1:]
        else:
            b = b + line
    return (a,a_name,b,b_name)

def read_scoring(scoring):
    state = 0
    header = []
    sm = {}
    for line in open(scoring, 'r'):
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        if None == re.search('^#', line):
            a = line.split(' ')
            if 0 == state:
                header = a
                state = state + 1
            else:
                if len(a) != 1 + len(header):
                    print('error: the number of entries did not match the header')
                    sys.exit(2)
                for i in range(len(header)):
                    sm[a[0] + header[i]] = int(a[i+1])
   
    return sm

def align(a, b, sm):
    n=len(a) # rows
    m=len(b) # cols

     # Initialize the score and traceback matrices
    globalScore = []
    globalTraceback = []
    localScore = []
    localTraceback = []
    for r in range(n+1):
        globalScore.append([])
        globalTraceback.append([])
        localScore.append([])
        localTraceback.append([])
        for c in range(m+1):
            globalScore[r].append(0)
            globalTraceback[r].append(0)
            localScore[r].append(0)
            localTraceback[r].append(0)

     # Set up score row gap penalties
    for r in range(n+1):
        globalScore[r][0] = r * GAP_SCORE
        localScore[r][0] = 0

    # Set up score column gap penalties
    for c in range(m+1):
        globalScore[0][c] = c * GAP_SCORE
        localScore[0][c] = 0

    # Set up traceback matrices
    globalTraceback[0][0] = 'done'
    localTraceback[0][0] = 'done'

    for r in range(1, n+1):
        globalTraceback[r][0] = 'up'
        localTraceback[r][0] = 'done'

    for c in range(1, m+1):
        globalTraceback[0][c] = 'left'
        localTraceback[0][c] = 'done'

    # Set variable to keep track of the current highest local score.
    # First entry is the value, second is row position, third is column position
    currHighestLocalScore = [0, 0, 0]

    for i in range(1, n+1):
        for j in range(1, m+1):
            globalDiag = globalScore[i-1][j-1] + sm[a[i-1]+b[j-1]]
            globalUp = globalScore[i-1][j] + GAP_SCORE
            globalLeft = globalScore[i][j-1] + GAP_SCORE

            localDiag = localScore[i-1][j-1] + sm[a[i-1]+b[j-1]]
            localUp = localScore[i-1][j] + GAP_SCORE
            localLeft = localScore[i][j-1] + GAP_SCORE
            
            globalScore[i][j] = max(globalDiag, globalUp, globalLeft)
            localScore[i][j] = max(0, localDiag, localUp, localLeft)
            
            # Global score comparisons
            if globalScore[i][j] == globalDiag:
                globalTraceback[i][j] = 'diag'
            elif globalScore[i][j] == globalUp:
                globalTraceback[i][j] = 'up'
            else:
                globalTraceback[i][j] = 'left'

            # Local score comparisons
            if localScore[i][j] > currHighestLocalScore[0]:
                currHighestLocalScore[0] = localScore[i][j]
                currHighestLocalScore[1] = i
                currHighestLocalScore[2] = j

            if localScore[i][j] == localDiag:
                localTraceback[i][j] = 'diag'
            elif localScore[i][j] == localUp:
                localTraceback[i][j] = 'up'
            elif localScore[i][j] == localLeft:
                localTraceback[i][j] = 'left'
            else:
                localTraceback[i][j] = 'done'

    #Global alignment calculation
    ag = ''
    bg = ''
    globalTraceRow = n
    globalTraceCol = m

    # Set the tracback to start at the bottom right corner
    currentGlobalTrace = globalTraceback[globalTraceRow][globalTraceCol]
    while currentGlobalTrace != 'done':
        if currentGlobalTrace == 'diag':
            globalTraceRow = globalTraceRow - 1
            globalTraceCol = globalTraceCol - 1
            ag = a[globalTraceRow] + ag
            bg = b[globalTraceCol] + bg
            currentGlobalTrace = globalTraceback[globalTraceRow][globalTraceCol]
        elif currentGlobalTrace == 'left':
            ag = GAP + ag
            globalTraceCol = globalTraceCol - 1
            bg = b[globalTraceCol] + bg
            currentGlobalTrace = globalTraceback[globalTraceRow][globalTraceCol]
        else:
            bg = GAP + bg
            globalTraceRow = globalTraceRow - 1
            ag = a[globalTraceRow] + ag
            currentGlobalTrace = globalTraceback[globalTraceRow][globalTraceCol]

    # Get the global score
    gs = globalScore[n][m]

    # Local alignment calculation
    al = ''
    bl = ''
    localTraceRow = currHighestLocalScore[1]
    localTraceCol = currHighestLocalScore[2]

    # Set traceback to start at the location of the highest local score in the matrix
    currentLocalTrace = localTraceback[localTraceRow][localTraceCol]
    while currentLocalTrace != 'done':
        if currentLocalTrace == 'diag':
            localTraceRow = localTraceRow - 1
            localTraceCol = localTraceCol - 1
            al = a[localTraceRow] + al
            bl = b[localTraceCol] + bl
            currentLocalTrace = localTraceback[localTraceRow][localTraceCol]
        elif currentLocalTrace == 'left':
            al = GAP + al
            localTraceCol = localTraceCol - 1
            bl = b[localTraceCol] + bl
            currentLocalTrace = localTraceback[localTraceRow][localTraceCol]
        else:
            bl = GAP + bl
            localTraceRow = localTraceRow - 1
            al = a[localTraceRow] + al
            currentLocalTrace = localTraceback[localTraceRow][localTraceCol]
    
    # Get the local score
    ls = currHighestLocalScore[0]


    return (ag, bg, gs, al, bl, ls)


def main():
    # get arguments
    if 3 != len(sys.argv):
        print('error: two input files required')
        sys.exit(2)
    fasta = sys.argv[1]
    scoring = sys.argv[2]

    # read in input files
    (a,a_name,b,b_name)=read_fasta(fasta)
    (sm)=read_scoring(scoring)

    (ag, bg, gs, al, bl, ls) = align(a, b, sm)

    print '#GLOBAL SCORE=' + str(gs)
    print a_name + '\t' + ag
    print b_name + '\t' + bg
    print '#LOCAL SCORE=' + str(ls)
    print a_name + '\t' + al
    print b_name + '\t' + bl

if __name__ == '__main__':
    main()