'''
Author: Jason Li
Assignment: CMSC 285 Project 3 - Runge Kutta 4
Professor: DuPont
Date: 10/31/2017
'''

import matplotlib.pyplot as plt
import numpy as np
import math


def f(t, y):
    return (-1 * np.array(y)) / (euclideanNorm(y))**3

def hermiteCubicSolve(aX, bX, aY, bY, aPrime, bPrime, x):

    term1 = (1+ 2*((x-aX)/(bX - aX))) * ((bX-x)/(bX-aX))**2 * np.array(aY)
    term2 = (x-aX) * ((bX-x)/(bX-aX))**2 * np.array(aPrime)
    term3 = (1+ 2*((bX-x)/(bX - aX))) * ((aX-x)/(aX-bX))**2 * np.array(bY)
    term4 = (x-bX) * ((bX-x)/(aX-bX))**2 * np.array(bPrime)

    return (term1 + term2 + term3 + term4)


def euclideanNorm(matrix):
    result = 0
    for term in matrix:
        result += term**2.0
    return math.sqrt(result)



class ComputeYnew:
    """docstring for step."""
    def __init__(self, t_now, Ynow, dt, fnow):

        #self.k1 = f(t_now, Ynow)
        self.k1 = fnow

        self.k2 = f(t_now + dt/2, np.sum([Ynow, np.array(self.k1)*dt/2], axis=0))

        self.k3 = f(t_now + dt/2, np.sum([Ynow, np.array(self.k2)*dt/2], axis=0))
        self.k4 = f(t_now + dt, np.sum([Ynow, np.array(self.k3)*dt], axis=0))
        self.dY = (np.array(self.k1) + 2*np.array(self.k2)+ 2.0*np.array(self.k2) + np.array(self.k4)) * dt/2
        self.Ynew = np.sum([Ynow + self.dY], axis=0)


class ComputeError:
    """docstring for Error."""

    def __init__(self, t_now, dt, Ynow, Ynew, k1):

        self.t_new = t_now+dt
        self.f_new = f(self.t_new,Ynew)
        self.c1 = 0.4-math.sqrt(0.06)
        self.c2 = 0.4+math.sqrt(0.06)

        self.t_c1 = t_now + dt*self.c1
        self.t_c2 = t_now + dt*self.c2

        self.Yc1 = hermiteCubicSolve(t_now, self.t_new, Ynow, Ynew, k1, self.f_new, self.t_c1)
        self.Yc2 = hermiteCubicSolve(t_now, self.t_new, Ynow, Ynew, k1, self.f_new, self.t_c2)

        self.Fc1 = f(self.t_c1, self.Yc1)
        self.Fc2 = f(self.t_c2, self.Yc2)

        self.w1 = (3*self.c2-1)/(6*(self.c1-self.c2)*(self.c1-1))
        self.w2 = (3*self.c1-1)/(6*(self.c2-self.c1)*(self.c2-1))
        self.wnew = 1 - self.w1 - self.w2
        self.dY_alt = dt * (self.w1 * np.array(self.Fc1) + self.w2 * np.array(self.Fc2) + self.wnew * np.array(self.f_new))



def plot(listF, tol):
    x = []
    y = []
    for value in listF:
        x.append(listF[0])
        y.append(listF[1])

    plt.plot(x, y)
    plt.title('trajectory with tolerance', tol)
    plt.show()


xVals = []
yVals = []
tVals = []
energys = []
lengthTime = 3 * 2 * math.pi / (2-math.sqrt(0.3))**1.5
timeFrame = [0,lengthTime]
time = timeFrame[0]
endTime = timeFrame[1]
Ynow = [1, 0]
dt = 0.05
dtmin = 0.01
dtmax = 0.1
agrow = 1.25
ashrink = 0.8
tol = 0.01
fnow = [0, 0.3]
listF = []

def takeStep(tol, timeFrame, time, endTime, Ynow, dt, dtmin, dtmax, agrow, ashrink, fnow, listF, xVals, yVals, tVals, energys):

    while time < endTime:
        tryStep = ComputeYnew(time, Ynow, dt, fnow)
        errorStep = ComputeError(time, dt, Ynow, tryStep.Ynew, tryStep.k1)
        ei = euclideanNorm(tryStep.dY - errorStep.dY_alt)/ dt

        # print("ei is ", ei)
        # print('y new is', tryStep.Ynew)

        if (ei < tol or dt == dtmin):

            time += dt
            Ynow = tryStep.Ynew
            print("yup")

            print('Energy is', 0.5 * (euclideanNorm(f(time, Ynow)))**2 - 1/euclideanNorm(Ynow))
            print('Y is', Ynow)
            print('t is', time)

            '''compile Ynow for trajectory plot'''
            xVals.append(Ynow[0])
            yVals.append(Ynow[1])
            tVals.append(time)
            energys.append(0.5 * (euclideanNorm(f(time, Ynow)))**2 - 1/euclideanNorm(Ynow))


            '''altering dt for better estimations'''
            if ei < tol/4 and dt*agrow < dtmax:
                dt *= agrow
            elif ei > 0.75* tol:
                dt *= ashrink

            '''making sure dt fits at the end of tolerance'''
            if time+dt > endTime:
                dt = endTime-time;
            elif time + 2*dt > endTime:
                dt = (endTime-time)/2

            fnow = f(time+dt, tryStep.Ynew)

        else:
            print("nope")
            if dt/2 < dtmin:
                dt = dtmin
            else:
                dt = dt/2
            print(dt)


    print(xVals)
    print(yVals)

    #plt.plot(tVals, xVals)
    #plt.title('trajectory with tolerance', tol)

    plt.plot(tVals, energys)
    plt.show()


takeStep(tol, timeFrame, time, endTime, Ynow, dt, dtmin, dtmax, agrow, ashrink, fnow, listF, xVals, yVals, tVals, energys)
