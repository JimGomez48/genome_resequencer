
#=============================================================
# Alignment Parameters
#=============================================================

class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, match, mismatch, gap_extend, gap_start=0):
        self.gap_start = gap_start
        self.gap_extend = gap_extend
        self.match = match
        self.mismatch = mismatch

    def match_char(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap_extend
        )


#=============================================================
# Printing and utility functions
#=============================================================

INFINITY = float('inf')


def __make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]


def __print_matrix(x, y, A):
    """
    Print the matrix with the (0,0) entry in the top left corner. Will label
    the rows by the sequence and add in the 0-row if appropriate.
    """
    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print "%5s" % (" "),
    else:
        print "%5s %5s" % (" ","*"),
        y = "*" + y
    # print the top row
    for c in x:
        print "%5s" % (c),
    print
    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print


#=============================================================
# Sequence Alignment
#=============================================================

def global_align(s1, s2, score=ScoreParam(4, -2, -1, -3)):
    """
    Uses the needleman-wunsch algorithm to globally align two sequences
    :param s1: the first sequence
    :param s2: the second sequence
    :return: the aligned sequences as a tuple
    """
    ROWS = len(s1) + 1
    COLS = len(s2) + 1
    UP = 'u'
    LEFT = 'l'
    DIAG = 'd'
    align1 = ''
    align2 = ''
    c_matrix = __make_matrix(ROWS, COLS)
    t_matrix = __make_matrix(ROWS, COLS)
    last_was_indel = False
    #initialize first column
    for row in xrange(ROWS):
      c_matrix[row][0] = row * score.gap_extend
    #initalize first row
    for col in xrange(COLS):
      c_matrix[0][col] = col * score.gap_extend
    # compute cost matrix
    for row in xrange(1, ROWS):
        for col in xrange(1, COLS):
            if s1[row-1] == s2[col-1]:
                top_left = c_matrix[row-1][col-1] + score.match
            else:
                top_left = c_matrix[row-1][col-1] + score.mismatch
            if t_matrix[row][col-1] == LEFT:
                left = c_matrix[row][col-1] + score.gap_extend
            else:
                left = c_matrix[row][col-1] + score.gap_start
            if t_matrix[row-1][col] == UP:
                top = c_matrix[row-1][col] + score.gap_extend
            else:
                top = c_matrix[row-1][col] + score.gap_start
            scores = (top_left, left, top)
            max_pos = scores.index(max(scores))
            c_matrix[row][col] = scores[max_pos]
            if max_pos == 0:
                t_matrix[row][col] = DIAG
            elif max_pos == 1:
                t_matrix[row][col] = LEFT
            elif max_pos == 2:
                t_matrix[row][col] = UP
    # print c_matrix
    # print t_matrix

    # traceback to get alignments
    row = len(t_matrix) - 1
    col = len(t_matrix[0]) - 1
    while row > 0 and col > 0:
        current = t_matrix[row][col]
        if current == DIAG:
            row -= 1
            col -= 1
            align1 += s1[row]
            align2 += s2[col]
        elif current == LEFT:
            col -= 1
            align1 += '-'
            align2 += s2[col]
        elif current == UP:
            row -= 1
            align1 += s1[row]
            align2 += '-'
    if row > 0:
        align1 += s1[row]
        align2 += '-'
    elif col > 0:
        align1 += '-'
        align2 += s2[col]
    return align1[::-1], align2[::-1]


def local_align(s1, s2, score=ScoreParam(4, -2, -1, -3)):
    """
    Uses the smith-waterman algorithm to locally align two sequences
    :param s1: the first sequence
    :param s2: the second sequence
    :return: the aligned sequences as a tuple
    """
    ROWS = len(s1) + 1
    COLS = len(s2) + 1
    UP = 'u'
    LEFT = 'l'
    DIAG = 'd'
    align1 = ''
    align2 = ''
    max_row = 0
    max_col = 0

    c_matrix = __make_matrix(ROWS, COLS)
    t_matrix = __make_matrix(ROWS, COLS)
    # compute cost matrix
    for row in xrange(1, ROWS):
        for col in xrange(1, COLS):
            if s1[row-1] == s2[col-1]:
                top_left = c_matrix[row-1][col-1] + score.match
            else:
                top_left = c_matrix[row-1][col-1] + score.mismatch
            if t_matrix[row][col-1] == LEFT:
                left = c_matrix[row][col-1] + score.gap_extend
            else:
                left = c_matrix[row][col-1] + score.gap_start
            if t_matrix[row-1][col] == UP:
                top = c_matrix[row-1][col] + score.gap_extend
            else:
                top = c_matrix[row-1][col] + score.gap_start
            scores = (top_left, left, top)
            max_pos = scores.index(max(scores))
            c_matrix[row][col] = max(0, scores[max_pos])
            if max_pos == 0:
                t_matrix[row][col] = DIAG
            elif max_pos == 1:
                t_matrix[row][col] = LEFT
            elif max_pos == 2:
                t_matrix[row][col] = UP
            if c_matrix[row][col] >= c_matrix[max_row][max_col]:
                max_row = row
                max_col = col
    # print c_matrix
    # print t_matrix

    # traceback pre-alignment
    while 1:
        if row == max_row and col == max_col:
            break
        if row > max_row and col > max_col:
            row -= 1
            col -= 1
            align1 += s1[row]
            align2 += s2[row]
        else:
            if row > max_row:
                row -= 1
                align1 += s1[row]
                align2 += '-'
            elif col > max_col:
                col -= 1
                align1 += '-'
                align2 += s2[col]
    # traceback local alignment
    while 1:
        if t_matrix[row][col] == 0:
            break
        if t_matrix[row][col] == DIAG:
            row -= 1
            col -= 1
            align1 += s1[row]
            align2 += s2[col]
        elif t_matrix[row][col] == LEFT:
            col -= 1
            align1 += '-'
            align2 += s2[col]
        elif t_matrix[row][col] == UP:
            row -= 1
            align1 += s1[row]
            align2 += '-'
    # traceback post-alignment
    while 1:
        if row == 0 and col == 0:
            break
        if row > 0 and col > 0:
            row -= 1
            col -= 1
            align1 += s1[row]
            align2 += s2[row]
        else:
            if row > 0:
                row -= 1
                align1 += s1[row]
                align2 += '-'
            elif col > 0:
                col -= 1
                align1 += '-'
                align2 += s2[col]
    return align1[::-1], align2[::-1]


def affine_align(s1, s2, score=ScoreParam(4, -2, 0, -3)):
    """Global alignment with affine penalties. We assume we are maximizing."""
    M = __make_matrix(len(s1) + 1, len(s2) + 1)
    X = __make_matrix(len(s1) + 1, len(s2) + 1)
    Y = __make_matrix(len(s1) + 1, len(s2) + 1)

    # initialize the columns
    for i in xrange(1, len(s1)+1):
        M[i][0] = -INFINITY
        X[i][0] = -INFINITY
        Y[i][0] = score.gap_start + i * score.gap_extend
    # initialize the rows
    for i in xrange(1, len(s2)+1):
        M[0][i] = -INFINITY
        X[0][i] = score.gap_start + i * score.gap_extend
        Y[0][i] = -INFINITY
    # fill scoring matrices
    for i in xrange(1, len(s1)+1):  # for every column
        for j in xrange(1, len(s2)+1):  # for every row
            # update the top level score
            M[i][j] = score.match_char(s1[i-1], s2[j-1]) + max(
                    M[i-1][j-1],
                    X[i-1][j-1],
                    Y[i-1][j-1]
            )
            # update the horizontal level score
            X[i][j] = max(
                    score.gap_start + score.gap_extend + M[i][j-1],
                    score.gap_extend + X[i][j-1],
                    score.gap_start + score.gap_extend + Y[i][j-1]
            )
            # update the vertical level score
            Y[i][j] = max(
                    score.gap_start + score.gap_extend + M[i-1][j],
                    score.gap_start + score.gap_extend + X[i-1][j],
                    score.gap_extend + Y[i-1][j]
            )

    return __affine_traceback(s1, s2, M, X, Y)
    # The best alignment comes from the matrix with the highest score
    # opt = max(M[len(x)][len(y)], X[len(x)][len(y)], Y[len(x)][len(y)])
    # print "x = %s & y = %s" % (x,y)
    # print "Scoring:", str(score)
    # print "M matrix ="
    # print_matrix(x,y,M)
    # print "X matrix ="
    # print_matrix(x,y,X)
    # print "Y matrix ="
    # print_matrix(x,y,Y)
    # print "Optimal =", opt


def __affine_traceback(s1, s2, M, X, Y):
    align1 = ''
    align2 = ''
    matrix = {
        'M': M,
        'X': X,
        'Y': Y,
    }
    i = len(M) - 1
    j = len(M[0]) - 1
    curr_matrix = 'M'
    last_matrix = curr_matrix
    # move through the matrices
    while i > 0 or j > 0:
        curr_matrix = __determine_next_matrix(last_matrix, M[i][j], X[i][j], Y[i][j])
        # if we're in the same matrix as the last iteration
        if curr_matrix == last_matrix:
            if curr_matrix == 'M':
                i -= 1
                j -= 1
                align1 += s1[i]
                align2 += s2[j]
            elif curr_matrix == 'X':
                j -= 1
                align1 += '-'
                align2 += s2[j]
            elif curr_matrix == 'Y':
                i -= 1
                align1 += s1[i]
                align2 += '-'
        else:  # if we hopped to a new matrix
            last_matrix = curr_matrix
    return align1[::-1], align2[::-1]


def __determine_next_matrix(last_matrix, M, X, Y):
    maximum = max(M, X, Y)
    if last_matrix == 'M'and maximum == M:
        return 'M'
    elif last_matrix == 'X'and maximum == X:
        return 'X'
    elif last_matrix == 'Y'and maximum == Y:
        return 'Y'
    if maximum == M:
        return 'M'
    elif maximum == X:
        return 'X'
    elif maximum == Y:
        return 'Y'


def main():
    # u = 'ACGTA'
    # v = 'ATA'
    u = 'TTAGAAATTGC'
    v = 'TAGATGCC'
    # u = 'GTGACCCTGGGGGAGTTCCCATCTTAGATTGTATTGTCTCGTCATGCATCGAGGCTAGGTGCCCGAGATCCCATGGCAAATAAACTTGTTCAACGCGGGA'
    # v = 'GTGACCCTGGGGGAGTTCCCATCTTAGATTGTATTGTCTCGTCATGCATCGAGGCTAGGTGCCCGAGAGATCATGGCAAATAAACTTGTTCAACGCGGGA'
    score_scheme = ScoreParam(
        match=4,
        mismatch=-2,
        gap_extend=-1,
        gap_start=-3,
    )
    # affine_align(u, v)
    # a, b = affine_align(u, v, score_scheme)
    # a, b = global_align(u, v, score=score_scheme)
    a, b = local_align(u, v, score=score_scheme)
    print a
    print b


if __name__ == '__main__':
    main()