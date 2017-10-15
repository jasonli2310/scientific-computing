'''
Author: Jason Li
Assignment: Scientific Computing HW1
Professor: DuPont
Date: 10/2/2017
'''

import matplotlib.pyplot as plt

import numpy
import math

def function1(x):
	return x**2

def function2(x):
	return math.sqrt(math.sin(math.pi*(x/2)))

def function3(x):
	return x**4

def function4(x):
	if(x < 1.0/3.0):
		return 1
	else:
		return 0

maxLevel = 0
numFunctCalls = 0
goodInterval = []
goodMidpoint = []

def recursion(mathFunct, a, b, tolerance, level, coarseArea=None):
	global numFunctCalls
	global maxLevel
	global goodInterval
	global goodMidpoint
	midpoint = (a + b) / 2.0
	if coarseArea is None:
		coarseArea = mathFunct(midpoint) * (b - a)
		numFunctCalls += 1

	fineLeft = mathFunct((a + midpoint) / 2.0) * (midpoint - a)
	fineRight = mathFunct((b + midpoint) / 2.0) * (b - midpoint)
	fineArea = fineLeft + fineRight
	interval = b-a
	goodInterval.append(interval)
	goodMidpoint.append(midpoint)
	numFunctCalls += 2
	level += 1
	if level >= 30:
		maxLevel = 30
		return fineArea
	if level < 4:
		return recursion(mathFunct, a, midpoint, tolerance / 2.0, level, coarseArea=fineLeft) + recursion(mathFunct, midpoint, b, tolerance / 2.0, level, coarseArea=fineRight)
	if abs(float(coarseArea) - float(fineArea)) <= tolerance:
		if level >= maxLevel:
			maxLevel = level
		return fineArea
	else:
		return recursion(mathFunct, a, midpoint, tolerance / 2.0, level, coarseArea=fineLeft) + recursion(mathFunct, midpoint, b, tolerance / 2.0, level, coarseArea=fineRight)

def plotEvals(mathFunct, a, b):
	global numFunctCalls
	global maxLevel
	plotData = []
	tols = [0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003, 0.00001, 0.000003, 0.000001]
	for tol in tols:
		plotData.append([recursion(mathFunct, a, b, tol, 0), numFunctCalls])
		numFunctCalls = 0
		maxLevel = 0

	logNumFunctCalls = [math.log(x[1]) for x in plotData]
	logTol = [math.log(1/x) for x in tols]
	plt.plot(logTol, logNumFunctCalls)
	plt.title(mathFunct)
	plt.xlabel("log of inverse tolerance")
	plt.ylabel("log of inverse of num function evals")
	plt.show()

plotEvals(function1, 0, 2)
plotEvals(function2, 0, 1)
plotEvals(function3, -1, 1)
plotEvals(function4, 0, 1)

def plotError(mathFunct, a, b):
	global numFunctCalls
	global maxLevel
	logErrors = []
	plotData = []
	lowTol = 0.0000000001
	tols = [0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003, 0.00001, 0.000003, 0.000001]
	actualValue = recursion(mathFunct, a, b, lowTol, 0)
	print(actualValue)
	for tol in tols:
		actError = abs (recursion(mathFunct, a, b, tol, 0) - actualValue)
		plotData.append(actError)
		print(actError)
		numFunctCalls = 0
		maxLevel = 0

	for x in plotData:
		try:
			logErrors.append(math.log(1.0/x))
		except:
			logErrors.append(numpy.inf)

	logTol = [math.log(1.0/x) for x in tols]


	plt.plot(logTol, logErrors)
	plt.title(mathFunct)
	plt.xlabel("log of inverse tolerance")
	plt.ylabel("log of inverse actual error")
	plt.show()

plotError(function1, 0, 2)
plotError(function2, 0, 1)
plotError(function3, -1, 1)
plotError(function4, 0, 1)

def plotGoodInterval(mathFunct, a, b):
	global goodInterval
	global goodMidpoint
	hardTol = 0.0003

	goodInterval = []
	goodMidpoint = []
	logInterval = []

	recursion(mathFunct, a, b, hardTol, 0)
	for x in goodInterval:
		logInterval.append(-math.log(x))

	plt.plot(goodMidpoint, goodInterval)
	plt.title(mathFunct)
	plt.xlabel("x position")
	plt.ylabel("negative log of good interval")
	plt.show()

plotGoodInterval(function1, 0, 2)
plotGoodInterval(function2, 0, 1)
plotGoodInterval(function3, -1, 1)
plotGoodInterval(function4, 0, 1)
