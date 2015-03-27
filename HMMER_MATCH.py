#!/usr/bin/env python3
#Christopher Lee (chmalee) Rolando Perez (rcperez)
#Group Members: Rolando Perez (rcperez)

"""
Program Overview:
    Parses files for FASTA header and returns corresponding sequence.

    Optional args are:

        -s --score     Only sequences better than a certain score are returned

    Bugs and Limitations:

Class CommandLine:
    Parses commandline for optional arguments. Only an input HMMER3.out
    file and corresponding FASTA file are required.

Class TargetMatcher:
    Makes use of fastaReader module to parse both FASTA files and HMMER3 files.
    After storing all the relevant information from HMMER3 file, read in FASTA
    file. Then check for matches of FASTA headers in both files, outputting the
    matches, thus allowing one visualization of homologous genes of interest.
"""
class CommandLine():
    def __init__(self, inOpts=None) :
        import argparse
        self.parser = argparse.ArgumentParser(description = "",
                                             epilog = "",
                                             add_help = True,
                                             prefix_chars = "-",
                                             usage = ""
                                            )
        self.parser.add_argument('hmmer3', action = 'store', help='input HMMER3 file') #query HMMER3 file
        self.parser.add_argument('fasta', action = 'store', nargs = '+', help = 'input FASTA file') #allows for multiple FASTA files

        #following line specify optional params
        self.parser.add_argument('-s', '--score', action='store', nargs='?', default=0.01,help='E-value')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


import fastaReader
"""Reads in files, does matching, returns appropriate info"""
class TargetMatcher():

    """class constructor"""
    def __init__(self, commandLine):
        self.hmmerFile = commandLine.args.hmmer3
        self.hmmerReader = fastaReader.FastAreader(self.hmmerFile)
        self.fastaFile = commandLine.args.fasta
        try:
            self.score = float(commandLine.args.score)
        except TypeError:
            self.score = 0.01
        self.dnaFlag = False
        self.seqDict = {}
        self.Error = False

    """reads in the appropriate lines from the HMMER3 file"""
    def findHMMERSeqs(self):
        #open hmmer3.out file specified by myCommandLine.args.hmmer3
        #read in lines of sequences that meet cutoff score if specified
        #and create a list of header sequences
        for scoreDict in self.hmmerReader.readHMMER(self.posSeqs, self.score):
            self.scores = scoreDict
        keylist = list(self.scores.keys())
        print(keylist)
        if keylist:
            try:
                """when target FASTA sequences are DNA,
                   a dictionary of scores:sequence,startPos,endPos
                   is returned. When the target FASTA sequences
                   are aa, a dictionary of FASTA:sequence is returned."""
                if float(keylist[0]):
                    self.dnaFlag = True
            except ValueError:
                self.dnaFlag = False
        else:
            self.Error = True

    """searches for appropriate sequences in FASTA file"""
    def findFastas(self):
        #given a sorted list of possible sequence headers, iterate through sorted FASTA headers
        #once a match is found
        #Create dictionary with header as key and sequence as value
        self.posSeqs = []
        for file in range(len(self.fastaFile)):
            myfastaReader = fastaReader.FastAreader(self.fastaFile[file])
            for header,sequence in myfastaReader.readFasta():
                self.posSeqs.append(header)
                self.seqDict[header] = sequence
        self.posSeqs.sort()
        self.findHMMERSeqs()

    """Finding and printing appropriate sequences"""
    def getSeqs(self):
        #check if HMMER3 seqs actually in possible fasta headers
        if self.Error:
            print("No Matches Found")
        else:

            #check if dict returned was DNA or AA
            if self.dnaFlag:

                #self.scores = {eval: [header,start,stop]}
                for eval in self.scores:

                    #get slice positions
                    startPos = int(self.scores[eval][1])
                    endPos = int(self.scores[eval][2])

                    #check if on reverse strand, in which case look for backwards string
                    #and do reverse complement on it
                    if startPos > endPos:
                        reverse = True
                        #these flags hold what we want to reverse complement
                        revStart = endPos
                        revEnd = startPos
                    else:
                        reverse = False

                    header = self.scores[eval][0]
                    if reverse:
                        reverseComp = self.reverseComplement(self.seqDict[header][revStart:revEnd])
                    if reverse:
                        #print the reverse start positions, which are the original positions
                        print(">%s %d-%d" % (header,startPos, endPos))
                    else:
                        print(">%s %d-%d" % (header, startPos, endPos))

                    #slice string
                    seq = self.seqDict[header]
                    if reverse:
                        print("\t%s" % reverseComp)
                    else:
                        print("\t%s" % seq[startPos:endPos])
            else: #Amino Acid, no special specifics needed
                for key in self.seqDict:
                    print(">%s" % key)
                    print("\t%s" % self.seqDict[key])

    """Method reverse complements the string if HMMER searched the reverse complement"""
    def reverseComplement(self,seqString):
        reversed = seqString[::-1]
        reverseDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        return ''.join(reverseDict[key] for key in reversed)


def main(myCommandLine = None):
    #parse commandline
    myCommandLine = CommandLine(myCommandLine)

    #read in files and print out homologous seqs
    #if nothing found the self.error flag prints: "No Matches"
    #and program exits
    myTargetMatcher = TargetMatcher(myCommandLine)
    myTargetMatcher.findFastas()
    myTargetMatcher.getSeqs()


if __name__ == "__main__":
    main()
