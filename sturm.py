
import numpy as np
import math

def signChange(a,b):
    if a*b < 0:
        return 1
    return 0

def newtons(p, pder, guess):
    np.polyval(p,guess)
    nextterm = np.polyval(p,guess)/np.polyval(pder,guess)
    guess = guess - nextterm
    if(abs(nextterm) < .001):
        return guess
    else:
        return newtons(p,pder,guess)

def sturm(x):
    p = []
    p.append(np.poly1d(x))
    p.append(np.polyder(p[0]))
    while True:
        remainder = np.polydiv(p[len(p)-2],p[len(p)-1])[1]
        if(remainder == np.poly1d(0)):
            break
        p.append(-remainder)
    return p

def nRoots(p):
    changes = 0
    bound = 9999999999
    for x in range(len(p)-1):
        changes += signChange(np.polyval(p[x],-bound),np.polyval(p[x+1],-bound))
    for x in range(len(p)-1):
        changes -= signChange(np.polyval(p[x],bound),np.polyval(p[x+1],bound))
    return changes

def intervals(p):
    numRoot = nRoots(p)
    num = numRoot
    interval = []
    upperbound = 1000
    for x in range(num-1):
        changes = numRoot - x

        if(x == 0):
            lowerbound = -1000
        else:
            lowerbound = interval[x-1]

        while(changes == numRoot-x):
            changes = 0
            lowerbound += .5
            for y in range(len(p)-1):
                changes += signChange(np.polyval(p[y],lowerbound),np.polyval(p[y+1],lowerbound))
            for y in range(len(p)-1):
                changes -= signChange(np.polyval(p[y],upperbound),np.polyval(p[y+1],upperbound))
        interval.append(lowerbound-.5)

    numRoot = 0
    for y in range(len(p)-1):
        numRoot += signChange(np.polyval(p[y],upperbound),np.polyval(p[y+1],upperbound))
    changes = numRoot
    while(changes == numRoot):
        changes = 0
        upperbound -= 1
        for y in range(len(p)-1):
            changes += signChange(np.polyval(p[y],upperbound),np.polyval(p[y+1],upperbound))
    interval.append(upperbound)
    interval.append(upperbound+1)
    return(interval)

def roots(poly, interval):
    p = np.poly1d(poly)
    pder = np.polyder(p)
    root = []
    for x in range(len(interval)-1):
        mid = (interval[x]+interval[x+1])/2
        root.append(newtons(p,pder,mid))
    return root

def complexroots(poly,roots):
    p = np.poly1d(poly)
    for x in roots:
        y = np.poly1d([1, -x])
        p = np.polydiv(p, y)[0]
    a = p[0]
    b = p[1]
    c = p[2]
    disc = np.lib.scimath.sqrt(b*b-4*a*c)
    root1 = (-b+(disc))/(2*a)
    root2 = (-b-(disc))/(2*a)
    return root1, root2

f0 = [1.0,-5.0,3.0,5.0,9.0,1.0,-2.0]
print("there are ", nRoots(sturm(f0)), "real roots")
print("the intervals are ", intervals(sturm(f0)))
print("the real roots are ", roots(f0,intervals(sturm(f0))))
print("the complex roots are", complexroots(f0, roots(f0,intervals(sturm(f0)))))
