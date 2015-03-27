#!/usr/bin/env python3
# Christopher Lee (chmalee) Rolando Perez (rcperez)
# Group Members : Rolando Perez (rcperez)

import sys
class FastAreader :
    
    """
    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):
 
    object attributes:
    fname: the initial file name
 
    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick, modifications by Christopher Lee(May 2014)
    Date: April 19, 2013

    readHMMER(): returns specified information from hmmer3 file
    Author: Christopher Lee
    Date: May 31, 2014
    """
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence. If filename was not provided,
        stdin is used. Whitespace is removed, and sequence is upcased.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            # initialize return containers
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
 
            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence
    
    def readHMMER(self, posSeqs, maxScore = 0.01):

        '''
        using filename given in init, compute list of headers, and
        dictionary of scores associated with headers. 
        If filename was not provided, stdin is used. 
        '''
        dnaScores = {} #evalue:header
        aaScores = {}  #header:evalue
        with self.doOpen() as fileH:
            #go through each line individually
            for line in fileH:
                #remove whitespace
                trimmed = line.strip()
                
                #make array to get each indiv value of table
                tmplist = trimmed.split()
                
                """
                examine list, check if first elem is an E-Value.
                If it is, check if there is a fasta header within the line. 
                If so, check if the fasta header is in the list of possible 
                headers. When a corresponding header is found, add the header
                as the key and the score as the value in the score dict.
               
                This creates a dictionary of headers that the user wants.
                """
                if tmplist:
                    try:
                        eval = float(tmplist[0])
                        #meets cutoff value, if it doesn't we do nothing
                        if (eval < maxScore):
                            #determine if nuc or aa
                            if len(tmplist) > 7: #protein
                                try:
                                    header = tmplist[8] #can only refer to data we want
                                    if header in posSeqs:
                                        aaScores[header] = eval
                                    else:
                                        for i in range(len(posSeqs)):
                                            cmp = header.find(posSeqs[i])
                                            if cmp != -1:
                                                aaScores[posSeqs[i]] = eval
                                except IndexError: #floats outside of table
                                    continue
                            else: #nucleotide
                                try:
                                    fasta = tmplist[3]
                                    start = tmplist[4]
                                    end = tmplist[5]
                                    dnaList = []
                                    dnaList.extend((fasta,start,end))
                                    if fasta in posSeqs:
                                        dnaScores[eval] = dnaList
                                    else:
                                        for i in range(len(posSeqs)):
                                            cmp = fasta.find(posSeqs[i])
                                            if cmp != -1:
                                                dnaScores[eval] = dnaList
                                except IndexError: #must be in non-table part of file
                                    continue
                    except ValueError: #catches instances of blank lines
                        continue
        if dnaScores:
            yield dnaScores
        else:
            yield aaScores
