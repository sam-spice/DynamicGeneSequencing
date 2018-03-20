#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np



class GeneSequencing:
    def __init__( self ):
        pass

    def align_all( self, sequences, banded, align_length ):
        #print(banded)
        #print(align_length)
        results = []
        print(len(sequences[2]))
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                s = {'align_cost':i+j,
                     'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                         len(sequences[i]), align_length, ',BANDED' if banded else ''),
                     'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                         len(sequences[j]), align_length, ',BANDED' if banded else '')}
                jresults.append(s)
            results.append(jresults)
        return results

    def initialize_unbanded(self, seq1, seq2, width, height):
        value_array = [[0 for x in range(width)] for y in range(height)]
        prev_array = [['' for x in range(width)] for y in range(height)]
        value_array[0][0] = ' '
        value_array[0][1] = '-'
        for x in range(len(seq1)):
            value_array[0][x + 2] = seq1[x]
            prev_array[0][x + 2] = seq1[x]
        value_array[1][0] = '-'
        for y in range(len(seq2)):
            value_array[y + 2][0] = seq2[y]
            prev_array[y + 2][0] = seq2[y]
        for x in range(2, width):
            value_array[1][x] = value_array[1][x - 1] + 5
            prev_array[1][x] = 'left'
        for y in range(2, height):
            value_array[y][1] = value_array[y-1][1] + 5
            prev_array[y][1] = 'up'
        return value_array, prev_array

    def get_alignment(self,seq1,seq2,prev_array,height,width):
        y = height
        x = width
        top_align = ''
        left_align = ''
        while(x > 2 and y > 2):
            string = prev_array[y][x]
            if string == 'diag':
                pass
            elif string == 'up':
                pass
            else:
                pass


    def compare_sequence_unbanded(self, seq1, seq2, align_length):
        width = len(seq1[:align_length]) + 2
        height = len(seq2[:align_length]) + 2
        value_array,prev_array = self.initialize_unbanded(seq1,seq2,width,height)
        for y in range(2,height):
            for x in range(2,width):
                left = value_array[y- 1][x] + 5
                up = value_array[y-1][x] + 5
                diag = value_array[y - 1][x - 1]
                if value_array[0][x] == value_array[y][0]:
                    diag = diag - 3
                else:
                    diag = diag + 1
                new_val = min(left, up, diag)
                value_array[y][x] = new_val
                if new_val == diag:
                    prev_array[y][x] = 'diag'
                elif new_val == up:
                    prev_array[y][x] = 'up'
                else:
                    prev_array[y][x] = 'left'
        for row in range(height):
            print(value_array[row])
        for row in range(height):
            print(prev_array[row])

poly = 'polynomial'
exp = 'exponential'
temp = GeneSequencing()
temp.compare_sequence_unbanded(poly,exp,1000)