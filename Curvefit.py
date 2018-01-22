#!/usr/bin/env python
#<examples/doc_model1.py>
from numpy import sqrt, pi, exp, linspace, loadtxt, log, diff, zeros, float, clip

from lmfit import  Model, Parameters

import lmfit
import numpy
import matplotlib.pyplot as plt
#plt.rcParams["figure.figsize"] = (15,8)
import constants

data = loadtxt('C:\curvefit\stmfull.dat',delimiter=',',skiprows=1)
x = data[:, 0]
y = data[:, 1]
l = data[:, 2]
j = data[:, 3]

newl = []
newlindice = []
t=0
for k in y:
	if k > -2.7E-6 and k < 0:
		newlindice.append(t)
		newl.append(k)
	t += 1
	
dy = zeros(y.shape,float)
dy[0:-1] = diff(x)/diff(j)
dyclip = clip(dy, 1E-7,1.86E-7)

truncdy = []
t = 0
for i in dy:
	if t >= newlindice[0] and t <= newlindice[-1]:
		truncdy.append(i)
	t += 1

def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))

def schottky(V, phiB, n):
	"Ideal Diode equation"
	return (constants.A*constants.Ar*constants.T**2 * exp((-constants.q*phiB)/(constants.k*constants.T)) * exp((constants.q*V)/(n*constants.k*constants.T)) * (1-exp((-constants.q*V)/(constants.k*constants.T))))
	
def linear(P,m,b):
	"Linear Fitting"
	#y = pars['y']
	#m = pars['n']
	return (m*P+b)
	
def logschottky(V,I,n):
	"Log fitting linear"
	return (log(I)+((constants.q)/(n*constants.k*constants.T))*V)

pars = Parameters()
pars.add('y', min=9.32E-8, max=1.88E-7)
pars.add('dy', min=0.162, max=1.157)
pars.add('n', value=1 ,min=1)

#j = pars['y']
#dy= pars['dy']
vlin = linspace(-2,2,200)
#gmodel = Model(gaussian)
#result = gmodel.fit(y, x=x, amp=5, cen=5, wid=1)

#gmodel = Model(schottky)
#result = gmodel.fit(y, V=x, phiB=constants.phiB, n=constants.n)

#gmodel = Model(linear)
#result = gmodel.fit(truncdy, P=newl, m=constants.n, b=0)

slope=(truncdy[-1]-truncdy[0])/(newl[-1]-newl[0])

#gmodel = Model(logschottky)
#result = gmodel.fit(y, V=x, I=constants.phiB, n=constants.n)
#print(result.fit_report())

plt.plot(y, dy,         'bo')
#plt.plot(x, result.init_fit, 'k--')
#plt.plot(newl, result.best_fit, 'r-')
#print (truncdy)
#print (newl)
#print (diff(x))
#print (diff(j))
#axes = plt.gca()
#axes.set_xlim([1E-7,1.75E-7])
#axes.set_ylim([0,2])
#plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')

# Four figure plot of schottky data
#plt.rcParams["figure.figsize"] = (15,8)
#fig, ax = plt.subplots(nrows=2, ncols=2)
#plt.figure(200)
#plt.subplot(2,2,1)
#plt.plot(y, dy, 'r-')
#plt.title('dV/dln|I|')
#plt.subplot(2,2,2)
#plt.plot(x, j, "r-")
#plt.title('ln|I|') 
#plt.subplot(2,2,3)
#plt.plot(x, y, "r-")
#plt.title('V vs I') 
#plt.subplot(2,2,4)
#plt.plot(x, l, "r-")
#plt.title('V vs |I|') 
plt.show()

# Using below to display 2 figure windows.
#plt.figure(300)
#plt.plot(x, j, "r-")
#plt.show(block=False)
#plt.show()
#<end examples/doc_model1.py>