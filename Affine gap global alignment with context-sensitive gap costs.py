
# coding: utf-8

# ## PROBLEM 1: Affine gap global alignment with context-sensitive gap costs (60 points)
# 
# Thus far, the algorithms that we have seen for alignment use the same gap score parameters regardless of where the gap is located.  It is plausible that insertions and deletions may be more likely to occur within certain contexts in a sequence.  In fact, in the human-chimp alignnment that we have been analyzing in class, gaps are much more frequently flanked by `A` and `T` and than by `G` or `C`, even after controlling for the fact that `A` and `T` are generally more frequent in those sequences.
# 
# In this problem you will implement the dynamic programming algorithm for global alignment with affine gap penalties and with gap scores that vary depending on the context of the gap.  Specifically, we will allow the gap score, $g$, to vary depending on the character immediately preceding and the character immediately following the gap.  In other words, interpreting each gap as a deletion, the gap score will be dependent on the single character before and the single character after the deletion, in the sequence in which the deletion occurred.  This variable gap score will be represented by a matrix, $g(a, b)$ with rows corresponding to the preceding character, $a$, and columns corresponding to the following character, $b$.  With this formulation, the primary dynamic programming equations are as follows:
# 
# $M(i, j) = \max\left\{
# \begin{array}{l}
# M(i - 1, j - 1) + S(x_i, y_j) \\
# I_x(i - 1, j - 1) + S(x_i, y_j) \\
# I_y(i - 1, j - 1) + S(x_i, y_j) \\
# \end{array}
# \right.$
# 
# $I_x(i, j) = \max\left\{
# \begin{array}{l}
# M(i - 1, j) + g(y_j, y_{j+1}) + s \\
# I_x(i - 1, j) + s \\
# \end{array}
# \right.$
# 
# $I_y(i, j) = \max\left\{
# \begin{array}{l}
# M(i, j - 1) + g(x_i, x_{i+1}) + s \\
# I_y(i, j - 1) + s \\
# \end{array}
# \right.$
# 
# When a gap occurs at the beginning of a sequence, we will consider the preceding character to be the empty string, and when a gap occurs at the end of a sequence, we will consider the following character to be the empty string. Correspondingly, in the above equations, we can consider $y_j = empty\_string$ if $j < 1$ or $j > length(y)$, and similarly for $x_i$.
# 
# Implement the algorithm as a function `align_global_affine_cs_gaps` below, that takes as input two sequences, `x` and `y`, a substitution matrix, a gap score matrix, and a space score.  Both matrices will be reprsented as dictionaries with two-element tuples, `(a, b)`, as keys and scores as values.  Your function should output a tuple of two elements, the first being the score of an optimal alignment, and the second being a single alignment that obtains that score. Your alignment should be represented as a list of two strings.  See the "Tests for PROBLEM 1" section at the bottom of this notebook for examples of the inputs and outputs.  
# 
# In the case that there are multiple optimal alignments, you should output the **highroad** alignment.  Note that in the case of affine gap alignment, the highroad alignment corresponds to preferring to traceback to the $I_x$ matrix first, the $M$ matrix second, and the $I_y$ matrix third, in the case that there are ties.  

# In[14]:


NEGATIVE_INFINITY = float("-inf")

def initiate(score_M, score_X, score_Y, trace_M, trace_X, trace_Y, gap_score, space_score, row, clm):
    """Initiate the score and trace matrix of X, Y and Z with corresponding values"""
    
    # Initiate score matrix of M, X, and Y
    for i in range(row):
        score_M.append([NEGATIVE_INFINITY] * clm)
        score_X.append([NEGATIVE_INFINITY] * clm)
        score_Y.append([NEGATIVE_INFINITY] * clm)
    score_M[0][0] = 0
    
    # Initiate base conditions in score_X and trace_X 
    for i in range(row):
        score_X[i][0] = gap_score + space_score * i
        trace_X[(i+1, 0)] = ("Ix", (i, 0))
    
    # Initiate base conditions in score_Y and trace_Y
    for j in range(clm):
        score_Y[0][j] = gap_score + space_score * j
        trace_Y[(0, j+1)] = ("Iy", (0, j))
        
def find_index(current_MXY):
    """Return the index to next cell following high-road criteria"""
    
    index_list = []
    max_value = max(current_MXY.values())
    
    # Find index with tie max_value
    for index, value in current_MXY.items():
        if value == max_value:
            index_list.append(index)
    
    # Return index following high-road criteria
    if "Ix" in index_list:
        return "Ix"
    elif "M" in index_list:
        return "M"
    else:
        return "Iy"

def find_path(score_M, score_X, score_Y, trace_M, trace_X, trace_Y, submatrix, gap_score_matrix, space_score, x, y):
    """Find a optimal path to the end of alignment by filling out the score and trace matrix"""
    
    # Initiate a temporary dict to compute the alignment score for the current cell
    current_M = {}
    current_X = {}
    current_Y = {}
    
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
    
            # Fill out M score and trace matrix
            current_M['M']  = score_M[i-1][j-1] + submatrix[(x[i-1], y[j-1])]
            current_M['Ix'] = score_X[i-1][j-1] + submatrix[(x[i-1], y[j-1])]
            current_M['Iy'] = score_Y[i-1][j-1] + submatrix[(x[i-1], y[j-1])]           
            score_M[i][j] = max(current_M.values())
            trace_M[(i,j)] = (find_index(current_M), (i-1, j-1))
            
            # Fill out Ix score and trace matrix          
            if j >= len(y):
                current_X['M'] = score_M[i-1][j] + gap_score_matrix[(y[j-1], "")] + space_score
            else:
                current_X['M'] = score_M[i-1][j] + gap_score_matrix[(y[j-1], y[j])] + space_score
            current_X['Ix'] = score_X[i-1][j] + space_score
            score_X[i][j] = max(current_X.values())
            trace_X[(i,j)] = (find_index(current_X), (i-1, j))
            
            ## Fill out Iy score and trace matrix
            if i >= len(x):
                current_Y['M'] = score_M[i][j-1] + gap_score_matrix[(x[i-1], "")] + space_score                
            else:
                current_Y['M'] = score_M[i][j-1] + gap_score_matrix[(x[i-1], x[i])] + space_score
            current_Y['Iy'] = score_Y[i][j-1] + space_score
            score_Y[i][j] = max(current_Y.values())
            trace_Y[(i,j)] = (find_index(current_Y), (i, j-1))

def traceback(pre_cell, trace_M, trace_X, trace_Y, x, y):
    """Trace back following the trace matrix to get the best global alignment"""

    # Initiate alignments for x and y sequences
    align_x = ""
    align_y = ""

    # Trace back to the top left cell
    while pre_cell[1] != (0, 0):
        
        trace_index = pre_cell[0]
        trace_coor = pre_cell[1]
        
        if trace_index == "M":
            align_x = x[trace_coor[0] - 1] + align_x
            align_y = y[trace_coor[1] - 1] + align_y
        if trace_index == "Ix":
            align_x = x[trace_coor[0] - 1] + align_x
            align_y = "-" + align_y
        if trace_index == "Iy":
            align_x = "-" + align_x
            align_y = y[trace_coor[1] - 1] + align_y
        pre_cell = get_pre_cell(pre_cell, trace_M, trace_X, trace_Y)
        
    return [align_x, align_y]

def get_pre_cell(pre_cell, trace_M, trace_X, trace_Y):
    """Get to the previous cell that points to the current cell"""
    
    trace_index = pre_cell[0]
    trace_coor = pre_cell[1]
    if trace_index == "M":
        return trace_M[trace_coor]
    if trace_index == "Ix":
        return trace_X[trace_coor]
    if trace_index == "Iy":
        return trace_Y[trace_coor]
            
def align_global_affine_cs_gaps(x, y, submatrix, gap_score_matrix, space_score):
    """Returns a global alignment of X with Y using the given subsitution matrix, SUBMATRIX, and
    an affine gap scoring, with the given context-sensitive gap score (a function of the character before
    and after a deletion), and the given space score.
    """

    # Initiate M, X, Y score_matrix and trace_matrix
    # Score matrix stores optimal alignment score
    # Trace matrix stores where the algnment score comes from
    row = len(x) + 1
    clm = len(y) + 1
    score_M = []
    score_X = []
    score_Y = []
    trace_M = {}
    trace_X = {}
    trace_Y = {}
    initiate(score_M, score_X, score_Y, trace_M, trace_X, trace_Y, gap_score_matrix[('','')], space_score, row, clm)
    
    # Fill out M, X, Y score_matrix and trace_matrix
    find_path(score_M, score_X, score_Y, trace_M, trace_X, trace_Y, submatrix, gap_score_matrix, space_score, x, y)
    
    # Get the global alignment score and starting point from the bottom right cell
    bottom_right = {"M": score_M[row-1][clm-1], "Ix": score_X[row-1][clm-1], "Iy": score_Y[row-1][clm-1]}
    alignment_score = max(bottom_right.values())
    trace_start = [find_index(bottom_right), (len(x), len(y))]

    # Get global alignment by tracing back from the starting point
    optimal_alignment = traceback(trace_start, trace_M, trace_X, trace_Y, x, y)

    return (alignment_score, optimal_alignment)


# ## Tests for PROBLEM 1

# ### Substitution and gap score matrices for testing

# In[18]:


DNA = list("ACGT")
DNA_AND_EMPTY = DNA + [""]

# A simple match=+1 and mismatch=-1 matrix
basic_submatrix = {(a, b): 1 if a == b else -1 for a in DNA for b in DNA}

# A gap score matrix with the same penalty for all contexts
uniform_gap_score_matrix = {(a, b): -2 for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}

# A gap score matrix with smaller penalties for gaps flanked by A and/or T
biased_gap_score_matrix = {(a, b): -1 if 'T' in (a, b) or 'A' in (a, b) else -2 
                           for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}

# A gap score matrix with all zeros.  Equivalent to a linear gap penalty function
zero_gap_score_matrix = {(a, b): 0 for a in DNA_AND_EMPTY for b in DNA_AND_EMPTY}

# A utility function for displaying a substitution or context-sensitive gap score matrix
def print_dict_matrix(m, width=5):
    """Prints the dict-based matrix m with the specified width for each column."""
    def print_row(fields):
        print("".join(["{:>{width}}".format(field, width=width) for field in fields]))
    labels = sorted({c for pair in m for c in pair})
    print_row([""] + labels)
    for a in labels:
        print_row([a] + [m.get((a, b), "") for b in labels])


# ### Functions for scoring an a global alignment with affine and context-sensitive gap scores

# In[19]:


def score_alignment(alignment, submatrix, gap_score_matrix, space_score):
    """Returns the score of the alignment given a substitution matrix, gap score matrix, and space score."""
    aligned_pair_scores = [submatrix[pair] for pair in zip(*alignment) if '-' not in pair]
    gap_x_scores, gap_y_scores = [gap_scores(s, gap_score_matrix, space_score) for s in alignment]
    return sum(aligned_pair_scores) + sum(gap_x_scores) + sum(gap_y_scores)
    
def gap_scores(s, gap_score_matrix, space_score):
    """Returns a list of the scores of the gaps within the aligned sequence s."""
    return [gap_score_matrix[context] + space_score * length for length, context in gaps(s)]

import re
def gaps(s):
    """Returns an iterator over the gaps of the string s.
    Each element of the iterator has the form (length, context) 
    where length is the length of the gap and context is a tuple of the two characters flanking the gap."""
    gap_pattern = re.compile('-+')
    for match in gap_pattern.finditer(s):
        yield (match.end() - match.start(), 
               (s[match.start() - 1: match.start()], s[match.end(): match.end() + 1]))


# ###  Test cases and testing functions

# In[20]:


test_case_inputs = {
    'small_1':  ("AGTA", "AGA",  basic_submatrix, uniform_gap_score_matrix, -1),
    'small_2':  ("AGT",  "AGT",  basic_submatrix, uniform_gap_score_matrix, -1),
    'small_3':  ("AGT",  "G",    basic_submatrix, uniform_gap_score_matrix, -1),
    'small_4':  ("G",    "AGT",  basic_submatrix, uniform_gap_score_matrix, -1),
    'small_5':  ("AGT",  "",     basic_submatrix, uniform_gap_score_matrix, -1),
    'small_6':  ("",     "AGT",  basic_submatrix, uniform_gap_score_matrix, -1),
    'small_7':  ("A",    "",     basic_submatrix, uniform_gap_score_matrix, -1),
    'small_8':  ("",     "",     basic_submatrix, uniform_gap_score_matrix, -1),
    'small_9':  ("CTC",  "TGCT", basic_submatrix, zero_gap_score_matrix,    -1),
    'small_10': ("TATA", "TTAC", basic_submatrix, biased_gap_score_matrix,  -1),
    'large_1':  ("TCATTCTGTTTATACTATCTTACTGGTTACCTTAATAATACAATCAGAATCGTAATTCGTCCTGTTCGT",
                 "TGTATAACTATGTCATCTAACCCCAAGCTTATCACTGCTTACGGAGGACAG",
                                 basic_submatrix, biased_gap_score_matrix,  -1)
}

test_case_correct_outputs = {
    'small_1': (0, ['AGTA', 
                    'AG-A']),
    'small_2': (3, ['AGT', 
                    'AGT']),
    'small_3': (-5, ['AGT', 
                     'G--']),
    'small_4': (-5, ['--G', 
                     'AGT']),
    'small_5': (-5, ['AGT', 
                     '---']),
    'small_6': (-5, ['---', 
                     'AGT']),
    'small_7': (-3, ['A', 
                     '-']),
    'small_8': (0, ['', 
                    '']),
    'small_9': (-1, ['--CTC', 
                     'TGCT-']),
    'small_10': (-1, ['TATA-', 
                      'T-TAC']),
    'large_1': (-15, ['TCATTCTGTTTATACTATCTTACTGGTTACCTTAATAATACAATCAGAATCGTAATTCGTCCTGTTCGT',
                      'T------GTATA-ACTATGTCA-TC-TAACCCCAAGCTTA---TCACT-GCTTA---CGGA--GGACAG'])
}

import numbers
def check_valid_alignment_result(result, x, y):
    """Checks that the alignment result is valid for sequences x and y."""
    assert isinstance(result, tuple), "Output is not a tuple"
    assert len(result) == 2, "Output does not have exactly two elements"
    score, alignment = result
    assert isinstance(alignment, list), "Alignment is not a list"
    assert isinstance(score, numbers.Number), "Score is not a number"
    assert len(alignment) == 2, "Alignment does not have exactly two elements"
    assert all(isinstance(element, str) for element in alignment), "Alignment elements are not strings"
    assert len(alignment[0]) == len(alignment[1]), "Alignment strings do not have the same length"
    assert alignment[0].replace('-', '') == x, "First string of alignment is not x"
    assert alignment[1].replace('-', '') == y, "Second string of alignment is not y"

def check_valid_alignment_score(result, submatrix, gap_score_matrix, space_score):
    """Checks that the computed score of the alignment is equal to the score given in the result."""
    score, alignment = result
    computed_score = score_alignment(alignment, submatrix, gap_score_matrix, space_score)
    assert computed_score == score, "Computed score ({}) does not equal the returned score ({})".format(computed_score,
                                                                                                        score)
def check_test_case(case_name, test_name=None, valid_result=True, valid_score=True, correct_alignment=True):
    inputs = test_case_inputs[case_name]
    correct_output = test_case_correct_outputs[case_name]
    result = align_global_affine_cs_gaps(*inputs)
    if valid_result:
        check_valid_alignment_result(result, *inputs[:2])
    if valid_score:
        check_valid_alignment_score(result, *inputs[2:])
    if correct_alignment:
        assert result == correct_output
    print("SUCCESS:", test_name if test_name else case_name, "passed!")


# ### Visible tests

# In[21]:


# TEST: small_1 valid output
check_test_case("small_1", test_name="small_1 valid output", valid_score=False, correct_alignment=False)


# In[22]:


# TEST: small_1 valid score
check_test_case("small_1", test_name="small_1 valid score", correct_alignment=False)


# In[23]:


# TEST: small_2
check_test_case("small_2")


# In[24]:


# TEST: small_3
check_test_case("small_3")


# In[25]:


# TEST: small_4
check_test_case("small_4")


# In[26]:


# TEST: small_5
check_test_case("small_5")


# In[27]:


# TEST: small_6
check_test_case("small_6")


# In[28]:


# TEST: small_7
check_test_case("small_7")


# In[29]:


# TEST: small_8
check_test_case("small_8")


# In[30]:


# TEST: small_9
check_test_case("small_9")


# In[31]:


# TEST: small_10
check_test_case("small_10")


# In[32]:


# TEST: large_1
check_test_case("large_1")

