#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    pass
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))
import math


class GeneSequencing:
    def __init__( self ):
        self.d = 3

    def align_all( self, sequences, banded, align_length ):
        #print(banded)
        #print(align_length)
        results = []
        if not banded:
            for i in range(len(sequences)):
                jresults = []
                for j in range(len(sequences)):
                    if j < i:
                        s = {'align_cost':i+j,
                            'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                                 len(sequences[i]), align_length, ',BANDED' if banded else ''),
                            'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                                 len(sequences[j]), align_length, ',BANDED' if banded else '')}
                    else:
                        value, top, left = self.compare_sequence_unbanded(sequences[i],sequences[j],align_length)
                        s = {'align_cost':value,'seqi_first100':top[:100],'seqj_first100':left[:100]}
                    jresults.append(s)
                results.append(jresults)
            return results
        else:
            for i in range(len(sequences)):
                jresults = []
                for j in range(len(sequences)):
                    if j < i:
                        s = {'align_cost':i+j,
                            'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                                 len(sequences[i]), align_length, ',BANDED' if banded else ''),
                            'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                                 len(sequences[j]), align_length, ',BANDED' if banded else '')}
                    else:
                        value, top, left = self.compare_banded(sequences[i],sequences[j],align_length)
                        s = {'align_cost':value,'seqi_first100':top[:100],'seqj_first100':left[:100]}
                    jresults.append(s)
                results.append(jresults)
            return results

    #space: O(1) time O(n+m)
    def finalize_unbanded(self,value_array, prev_array, width, height):
        for x in range(2, width):
            value_array[1][x] = value_array[1][x - 1] + 5
            prev_array[1][x] = 'left'
        for y in range(2, height):
            value_array[y][1] = value_array[y-1][1] + 5
            prev_array[y][1] = 'up'

    #space: O(1) time:(1)
    def finalize_banded(self,value_array, prev_array, width, height):
        value_array[1][1] = 0
        for i in range(1,self.d + 1):
            value_array[1][1+i] = value_array[1][i] + 5
            value_array[1+i][1] = value_array[i][1] + 5


    #space: O(mn) time: O(nm)
    def initialize(self, seq1, seq2, width, height, banded):
        # for both the banded and unbanded this takes space O(mn) time O(mn) but since this is just set up and doesn't
        # involve any actual calculations I am not going to count it
        if not banded:
            value_array = [[0 for x in range(width)] for y in range(height)]
        else:
            value_array = [[math.inf for x in range(width)] for y in range(height)]
        # space O(mn) time O(mn) but like the previous entry this is merely set up and involves on calcs so not counted
        prev_array = [['' for x in range(width)] for y in range(height)]

        #set first entry in the top string to be a hyphen, used for alignment purposes
        value_array[0][1] = '-'
        #places the string across the top of the value and backtrace array
        for x in range(width - 2):
            value_array[0][x + 2] = seq1[x]
            prev_array[0][x + 2] = seq1[x]

        #set first entry in left string to be a hyphen used for alignment purposes
        value_array[1][0] = '-'
        #places the string across the left side of both the value and the backtrace array
        for y in range(height - 2):
            value_array[y + 2][0] = seq2[y]
            prev_array[y + 2][0] = seq2[y]

        # if the algorithm is unbanded place values increasing in increments of 5 across
        # the top and left side of the value array, set prev values accordingly
        # time O(n+m)
        if not banded:
            self.finalize_unbanded(value_array,prev_array,width,height)
        # if banded only initialize the first 6 cells in the array, the first 3 to the right of 1,1 and the first 3
        # below 1,1 time O(1)
        else:
            self.finalize_banded(value_array, prev_array, width, height)
        # return the initialized arrays
        return value_array, prev_array

    # time: O(n+m)
    def get_alignment(self, prev_array, height, width):
        y = height - 1
        x = width - 1
        top_align = ''
        left_align = ''
        # iterates through the previous array following the backtraces until it arrives at 1,1.
        # While doing this it constructs the strings according to the rules
        # will take at most time: O(n+m) if the trace takes it across the bottom of the entire alignment
        # then across the side, should be impossible but you never know
        while x > 1 and y > 1:
            # gets the value of the back trace from the prev array
            string = prev_array[y][x]

            # if this space was acquired from the space to its diagonal adjust the x and y values according and add the
            # corresponding characters to the alignment string
            if string == 'diag':
                top_align = prev_array[0][x] + top_align
                left_align = prev_array[y][0] + left_align
                x = x-1
                y = y-1

            # if the cell was derived from the cell above it place a hyphen in the alignment for the top string and
            # the char in the left string
            elif string == 'up':
                top_align = '-' + top_align
                left_align = prev_array[y][0] + left_align
                y = y-1

            # if the cell was derived from the cell to its left add a hyphen to left alignment
            # and char in the top alignment
            else:
                top_align = prev_array[0][x] + top_align
                left_align = '-' + left_align

            x = x - 1

        return top_align, left_align

    # compare the sequences using the unbounded algorithm
    # time O(mn) space O(mn)
    def compare_sequence_unbanded(self, seq1, seq2, align_length):
        width = len(seq1[:align_length]) + 2
        height = len(seq2[:align_length]) + 2

        # space: O(nm) time: O(mn)
        # I am not really considering the time here, insignificant compared to that of the algorithm
        value_array,prev_array = self.initialize(seq1,seq2,width,height, False)

        # iterates over the entire array, will always be O(mn)
        # time: O(mn) space: (1)
        for y in range(2,height):
            for x in range(2,width):

                # time: O(1) space: O(1)
                self.adjust_cell(y,x,value_array,prev_array)

        # gets the alignment for the given sequences
        # time: O(n+m) space: (
        top_align, left_align = self.get_alignment(prev_array, height, width)
        return value_array[height -1][width -1], top_align, left_align

    # compares the values surrounding the given cell, updates its value and places the backtrace values
    # time: O(1) space: O(1)
    def adjust_cell(self,y, x, value_array, prev_array):
        # get the values to be compared
        left = value_array[y][x - 1] + 5
        up = value_array[y - 1][x] + 5
        diag = value_array[y - 1][x - 1]
        # if the letters in the sequence match score of -3 is subtraced from the diagonal squares value
        # otherwise a 1 is added
        if value_array[0][x] == value_array[y][0]:
            diag = diag - 3
        else:
            diag = diag + 1

        # get the min value from the 3
        new_val = min(left, up, diag)
        # set the min value to the newly acquired value
        value_array[y][x] = new_val
        # check which square the value came from, prioritize diagonals if possible
        # unimportant merely the decision I made
        if new_val == diag:
            prev_array[y][x] = 'diag'
        elif new_val == up:
            prev_array[y][x] = 'up'
        else:
            prev_array[y][x] = 'left'

    # implementation of the banded algorithm, runs significantly faster than the unbanded for the same length strings
    # time: O(n+m) space: O(nm)
    def compare_banded(self,seq1, seq2, align_length):
        width = len(seq1[:align_length]) + 2
        height = len(seq2[:align_length]) + 2
        # time: O(nm) but I am not really including this since it is initialization overhead
        # space: O(nm)
        value_array, prev_array = self.initialize(seq1,seq2,width,height,True)
        # time(n) for sequences long than 4 character i.e., all of them the inner 4 loop is overtaken
        # by the outer because it stays constant
        for y in range(2,height):
            # time O(1)
            for i in range(4):
                new_x = y + i
                if new_x < width:
                    self.adjust_cell(y,new_x,value_array,prev_array)
                new_y = y + i
                if new_y < height and y < width:
                    self.adjust_cell(new_y,y,value_array,prev_array)
        # checks to see if the bottom right-most cell has been set, if not returns no alignment
        # in the absolute worst case the alignment will be retrieved in time: O(n+m)
        if not value_array[height - 1][width - 1] == math.inf:
            top_align, left_align = self.get_alignment(prev_array,height,width)
            return value_array[height - 1][width - 1], top_align, left_align
        else:
            return math.inf, 'No Alignment Possible', 'No Alignment Possible'
