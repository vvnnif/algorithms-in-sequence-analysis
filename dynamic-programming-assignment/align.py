#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Name and student ID here!>
"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2



def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix = []
    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)
    
    if strategy == 'global':
        for i in range(M):
            score_matrix[i][0] -= i * gap_penalty
        for i in range(N):
            score_matrix[0][i] -= i * gap_penalty
    else:
        for i in range(M):
            score_matrix[i][0] = 0
        for i in range(N):
            score_matrix[0][i] = 0


    ### 2: Fill in Score Matrix
 
    #####################
    # START CODING HERE #
    #####################

    def dp_function(i, j):
        if strategy == 'global':
            if i == j == 0:
                return 0, ['stop'] # For global alignment, only stop at (0,0)
        else:
            if i == 0 or j == 0:
                return 0, ['stop'] # For local and semiglobal alignment, stop at first row or column
        options = {}
        if i > 0:
            options.update({'up': score_matrix[i-1][j] - gap_penalty})
        if i > 0 and j > 0:
            options.update({'diagonal': score_matrix[i-1][j-1] + substitution_matrix[seq1[i-1]][seq2[j-1]]})
        if j > 0:
            options.update({'left': score_matrix[i][j-1] - gap_penalty})
        if strategy == 'local':
            options.update({'stop': 0}) # For local alignment, also stop at a "fat zero"
        max_value = max(options.values())
        return max_value, [k for k, v in options.items() if v == max_value]


    for i in range(1,M):
        for j in range(1,N):
            max_value, _ = dp_function(i, j)
            score_matrix[i][j] = max_value
    

    #####################
    #  END CODING HERE  #
    #####################   
    
    ''' Returns the position in the score matrix
        from whence to start the traceback,
        which depends on the algorithm used.
        The output will consist of two integers,
        the first being the row index and the second
        being the column index.
    '''
    def get_traceback_start():
        if strategy == 'global':
            return M-1, N-1 # Bottom-right corner
        if strategy == 'local':
            # Return the position with maximum score. If
            # there are multiple such positions,
            # choose the upper-rightmost one.
            maxscore = float("-inf")
            best_i, best_j = N, 0
            for i in range(M):
                for j in range(N):
                    val = score_matrix[i][j]
                    if (val > maxscore
                        or (val == maxscore and (j > best_j or (j == best_j and i < best_i)))):                        
                        maxscore = val
                        best_i, best_j = i, j
            return best_i, best_j
        
        if strategy == 'semiglobal':
            max_score = float("-inf")
            start_positions = []
    
            # Check last row (seq1 fully aligned, seq2 may have suffix gaps)
            i = M - 1
            for j in range(N):
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    start_positions = [(i, j)]
                elif score_matrix[i][j] == max_score:
                    start_positions.append((i, j))

            # Check last column (seq2 fully aligned, seq1 may have suffix gaps)  
            j = N - 1
            for i in range(M):
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    start_positions = [(i, j)]
                elif score_matrix[i][j] == max_score:
                    start_positions.append((i, j))

            # Apply "high road": prefer higher row (lower i index)
            if start_positions:
                start_positions.sort(key=lambda x: (x[0], x[1]))  # Sort by i, then j
                return start_positions[0]  # Lowest i index
        

    ''' Dynamic programming traceback starting at matrix index
        (start_i, start_j). 
        Immediately returns the resulting aligned sequences.
    '''
    def traceback(start_i, start_j):
        aligned_seq1 = aligned_seq2 = ''
        i, j = start_i, start_j
        while True:
            _, options = dp_function(i, j)
            if 'up' in options: # Corresponds to a gap in sequence 2
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                i -= 1
            elif 'diagonal' in options: # Corresponds to no gap
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                i -= 1
                j -= 1
            elif 'left' in options: # Corresponds to a gap in sequence 1
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
            elif 'stop' in options: # For global alignment, this happens when i=j=0.
                                    # For local alignment, this happens when i=0, j=0,
                                    # or traversing a "fat zero" is the only remaining option.
                                    # For semiglobal alignment, this happens when i=0 or j=0.
                break

        if strategy == 'semiglobal': # For semiglobal alignment, we still have to add trailing
                                     # parts of the sequences after alignment.
            unaligned_seq1_left = i # Number of still unaligned letters right of the alignment
            unaligned_seq1_right = M-1-start_i # Number of still unaligned letters left of the alignment
            unaligned_seq2_left = j # Same for seq2
            unaligned_seq2_right = N-1-start_j  

            # Should do this because s[-0:] = s, not ''
            last_seq1_letters = '' if unaligned_seq1_right == 0 else seq1[-unaligned_seq1_right:]
            last_seq2_letters = '' if unaligned_seq2_right == 0 else seq1[-unaligned_seq2_right:]

            # Fill up both aligned sequences with unaligned letters and gaps
            aligned_seq1 =  (
                                ('-' * unaligned_seq2_left)   + 
                                seq1[:unaligned_seq1_left]    + 
                                aligned_seq1                  + 
                                last_seq1_letters             +
                                ('-' * unaligned_seq2_right)  
                            )
            aligned_seq2 =  (
                                ('-' * unaligned_seq1_left)  + 
                                seq2[:unaligned_seq2_left]   + 
                                aligned_seq2                 + 
                                last_seq2_letters            +
                                ('-' * unaligned_seq1_right)
                            )

        return aligned_seq1, aligned_seq2
    

    start_i, start_j = get_traceback_start()
    align_score = score_matrix[start_i][start_j]
    aligned_seq1, aligned_seq2 = traceback(start_i, start_j)
    alignment = (aligned_seq1, aligned_seq2, align_score)
    return (alignment, score_matrix)



def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()