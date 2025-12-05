#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the FIRST Advanced Question of the Hidden Markov Models
    assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via
    corresponding Canvas assignment. Note this script will be graded manually,
    if and only if your "hmm.py" script succesfully implements Baum-Welch
    training! Continuous Feedback will not be available for this script.

AUTHOR:
    Finn van Vlaanderen
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from hmm_utility import load_tsv
from numpy.random import choice



def parse_args():
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline. Take a look at the
    # argparse documentation, the parser in hmm_utility.py or align.py
    # (from the Dynamic Programming exercise) for hints on how to do this.

    parser = ArgumentParser(prog='python3 sequence_generator.py',
                            formatter_class=RawTextHelpFormatter,
                            description=
                            '   Generate random sequences with given transition and\n'
                            '   emission matrices\n\n'
                            '   Example syntax:\n'
                            '       python3 sequence_generator.py A.tsv E.tsv\n'
                            '       python3 sequence_generator.py -n 20 A.tsv E.tsv')

    # Positionals
    parser.add_argument('transition', help='path to TSV formatted transition matrix')
    parser.add_argument('emission', help='path to TSV formatted emission matrix')

    # Optionals
    parser.add_argument('-n', dest='number', type=int, default=1, help=
                        'number of sequences to generate')
    parser.add_argument('-o', dest='out_dir', help=
                        'path to a directory where output files are saved.')
    
    #####################
    #  END CODING HERE  #
    #####################

    return parser.parse_args()


def generate_sequence(A,E):
    #####################
    # START CODING HERE #
    #####################
    # Implement a function that generates a random sequence using the choice()
    # function, given a Transition and Emission matrix.
    
    # Look up its documentation online:
    # https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.random.choice.html
            
    sequence = ''
    allStates = list(A.keys())
    emittingStates = list(E.keys())
    allSymbols = list(list(E.values())[0].keys())

    state = 'B'
    while state != 'E':
        transition_p = list(A[state].values())
        state = choice(allStates, p=transition_p)
        if state in emittingStates:
            emission_p = list(E[state].values())
            symbol = choice(allSymbols, p=emission_p)
            sequence = sequence + symbol
    
    #####################
    #  END CODING HERE  #
    #####################
    
    return sequence



def main():
    args = parse_args()
    #####################
    # START CODING HERE #
    #####################
    # Uncomment and complete (i.e. replace '?' in) the lines below:
            
    N = args.number                # The number of sequences to generate
    if args.out_dir:
        out_file = args.out_dir    # The file path to which to save the sequences
    else:
        out_file = 'output.fasta'
    A = load_tsv(args.transition)  # Transition matrix
    E = load_tsv(args.emission)    # Emission matrix
    with open(out_file,'w') as f:
        for i in range(N):
            seq = generate_sequence(A, E)
            f.write('>random_sequence_%i\n%s\n' % (i,seq))
        
    #####################
    #  END CODING HERE  #
    #####################
    


if __name__ == "__main__":
    main()
