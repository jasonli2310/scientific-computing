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

def factor(arrayA):

    rowCount = 0

    for row in arrayA:

        # factor row
        factor = 1/row[rowCount]

        print(factor)

    # for element in row1:
    #     element = element*factor
    #
    #
    #
    # row1[0]
