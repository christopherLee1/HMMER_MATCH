#!/usr/bin/env python3
#Christopher Lee(chmalee)
#Group Members: Rolando Perez (rcperez)

"""
Parses commandline for optional arguments. Only an input HMMER3.out 
file and corresponding FASTA file are required.
"""


class CommandLine():
    def __init__(self, inOpts=None) :
        import argparse
        self.parse = argparse.ArgumentParser(description = "",
                                             epilog = "",
                                             add_help = True,
                                             prefix_chars = "-",
                                             usage = ""
                                            )
        self.parser.add_argument('hmmer3', action = 'store', help='input HMMER3 file')
        self.parser.add_argument('fasta', action = 'store', help = 'input FASTA file')
        #self.parser.add_argument('outFile', action = 'store', default=TargetMatch.out, help='output file name')
        self.parser.add_argument('-n', '--noDups', action = 'store', nargs='?', help='no duplicate sequences')
        self.parser.add_argument('-m', '--mRNA', action = 'store', nargs='?', help='Only mRNA returned')
        self.parser.add_argument('-d', '--dna', nargs='?', action = 'store', help='Only DNA returned')
        self.parser.add_argument('-a', '--aminoAcid', action = 'append', nargs='?', help='Only amino acid sequence returned')
        self.parser.add_argument('-s', '--score', action='store', nargs='?',help='E-value')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

