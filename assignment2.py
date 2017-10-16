'''
Author: Jason Li
Assignment: Scientific Computing HW1
Professor: DuPont
Date: 10/15/2017
'''

import matplotlib.pyplot as plt
import numpy as np
import math

def a(i):
    return 0.03*i

def b(i):
    return 0.01*i*i

def d(i):
    return 2.1 + 0.1*i

row1 = [d(0), a(0), 0, 0, 0]
row2 = [b(0), d(1), a(1), 0, 0]
row3 = [0, b(1), d(2), a(2), 0]
row4 = [0, 0, b(2), d(3), a(3)]
row5 = [0, 0, 0, b(3), d(4)]

x = [1, 2, 3, 4, 5]
arrayA = [row1, row2, row3, row4, row5]

def multiply(arrayA, vectorV):
    arrayA= np.array(arrayA)
    vectorV= np.array(vectorV)
    return arrayA*vectorV

testArray = [[2, 1, 0], [3, 6, 2], [0, 2, 5]]

def factor(arrayA):

    rowCount = 0

    for index in range(len(arrayA)-1):
        # factor row
        factor = 1/arrayA[index][rowCount]
        conceptRow1 = list(arrayA[rowCount])
        conceptRow2 = list(arrayA[rowCount+1])

        for idx, coeff in enumerate(conceptRow1):
            conceptRow1[idx] = coeff*factor

        variable = conceptRow2[rowCount]
        # changes conceptRow2 to the difference
        for i, coeff in enumerate(conceptRow2):
            conceptRow2[i] = conceptRow2[i]-conceptRow1[i]*variable

        #changes real row1
        for idx, val in enumerate(arrayA[rowCount]):

            print(conceptRow1)
            print(arrayA[0])

            if (idx > rowCount):
                arrayA[rowCount][idx] = conceptRow1[idx]

        #changes real row2
        for idx, val in enumerate(arrayA[rowCount+1]):
            if idx > rowCount:
                arrayA[rowCount+1][idx] = conceptRow2[idx]

        rowCount += 1

    return arrayA

# print(np.array(arrayA))
print("final", factor(testArray))
