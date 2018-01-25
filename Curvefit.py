#!/usr/bin/env python
#<examples/doc_model1.py>
from numpy import sqrt, pi, exp, linspace, loadtxt, log, diff, zeros, float, clip

from lmfit import  Model, Parameters

import lmfit
import numpy
import matplotlib.pyplot as plt
#plt.rcParams["figure.figsize"] = (15,8)
import constants

# ------------------------------------------------------------
# Setting defaults for plots
# ------------------------------------------------------------

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# ------------------------------------------------------------
# Importing data from file
# ------------------------------------------------------------

data = loadtxt('diode.dat',delimiter=',',skiprows=1)
x = data[:, 0]
y = data[:, 1]
l = data[:, 2]
j = data[:, 3]

# ------------------------------------------------------------
# Creating arrays for differentiation
# ------------------------------------------------------------
	
dy = zeros(y.shape,float)
dy[0:-1] = diff(x)/diff(j)
H = zeros(y.shape,float)
H[0:] = (x-(((constants.k*constants.T)/constants.q)*log((l)/(constants.A*constants.Ar*(constants.T)*(constants.T)))))

lower_lim=-1.5E-5
upper_lim=-5.99E-3

# ------------------------------------------------------------
# Finding the indices of the arrays to the linear approximations
# ------------------------------------------------------------

newl = []
newlindice = []
truncdy = []
truncH = []
newHl = []
newHlindice = []
t = 0
for k in y:
	if k <=lower_lim and k >=upper_lim:
	#if k <=lower_lim:
		newlindice.append(t)
		newl.append(k)
	t += 1

t = 0
for i in dy:
	if t >= newlindice[0] and t <= newlindice[-1]:
		truncdy.append(i)
	t += 1

t = 0
for k in y:
	if k <=lower_lim and k >=upper_lim:
	#if k <=lower_lim:
		newHlindice.append(t)
		newHl.append(k)
	t += 1
t = 0
for i in H:
	if t >= newHlindice[0] and t <= newHlindice[-1]:
		truncH.append(i)
	t += 1

# ------------------------------------------------------------
# Definitions for model
# ------------------------------------------------------------
	
def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))

def schottky(V, phiB, n):
	"Ideal Diode equation"
	return (constants.A*constants.Ar*constants.T**2 * exp((-constants.q*phiB)/(constants.k*constants.T)) * exp((constants.q*V)/(n*constants.k*constants.T)) * (1-exp((-constants.q*V)/(constants.k*constants.T))))
	
def linear(P,Rs,n):
	"Linear Fitting"
	#y = pars['y']
	#n = pars['n']
	return (Rs*P+n*constants.k*constants.T/(constants.q))
	
def Hlinear(P,phiB, Rs):
	"Linear Fitting"
	#y = pars['y']
	#m = pars['n']
	return ((Rs)*P+n_fit*phiB)
	
def logschottky(V,I,n):
	"Log fitting linear"
	return (log(I)+((constants.q)/(n*constants.k*constants.T))*V)

# ------------------------------------------------------------
# Evaluating the models
# ------------------------------------------------------------	

#gmodel = Model(schottky)
#result = gmodel.fit(y, V=x, phiB=constants.phiB, n=constants.n)

# ------------------------------------------------------------	
# --Model of the derivative of voltage and Ln I---------------	
# ------------------------------------------------------------	
gmodel = Model(linear)
Gresult = gmodel.fit(truncdy, P=newl, Rs=0, n=0)
Rs_fit = Gresult.best_values["Rs"]
n_fit = abs(Gresult.best_values["n"])

# ------------------------------------------------------------	
# --Model of the Cheung Function to extract Schottky Barrier--	
# ------------------------------------------------------------	
HI = zeros(y.shape,float)
HI[0:] = n_fit*constants.phiB*constants.q + y*Rs_fit
Hmodel = Model(Hlinear)
Hresult = Hmodel.fit(truncH,Rs=Rs_fit, P=newHl, phiB=0.85)

phiB_fit = abs(Hresult.best_values["phiB"])
phiB_final = phiB_fit
#n_calc=n_fit*constants.q/(constants.k*constants.T)
slope=(truncdy[-1]-truncdy[0])/(newl[-1]-newl[0])

# ------------------------------------------------------------
# Saving data to files
# ------------------------------------------------------------

f=open("stmfull_out.dat","a+")
t=0
f.write("!Current(A),|Current(A)|,dV/dlnI(V),H(I)\n")
for i in y:

	f.write('%s' % y[t]+","+'%s' % l[t]+","+'%s' % dy[t]+","+'%s' % H[t]+"\n")
	t += 1
	
f.close()

f=open("stmfull_linear_out.dat","a+")
t=0
f.write("|Current(A)|,dx/dlnI(V),H(I)\n")
for i in newl:

	f.write('%s' % newl[t]+","+'%s' % truncdy[t]+","+'%s' % truncH[t]+"\n")
	t += 1
	
f.close()

f=open("stmfull_linear_fits.dat","a+")
f.write('%s' % Gresult.fit_report())
f.write(Hresult.fit_report())	
f.close()

#gmodel = Model(logschottky)
#result = gmodel.fit(y, V=x, I=constants.phiB, n=constants.n)

# ------------------------------------------------------------
# Plotting the fit results
# ------------------------------------------------------------

print(Gresult.fit_report())
print("----")
print(Gresult.best_values)
print("----")
print(Hresult.fit_report())
print("----")
print(phiB_final)
#print("----")
#print(Hresult.params)


##plt.plot(newl, truncdy,         'bo')
plt.plot(x,y,'go')
#plt.plot(newl, truncdy,         'bo')
#plt.plot(newHl,truncH,'ro')
#plt.plot(newHl,Hresult.best_fit, 'r-' )
#plt.plot(newl, Gresult.best_fit, 'b-')
#plt.plot(newl, Gresult.init_fit, 'k--')

#plt.plot(y, HI, 'ro')
##plt.plot(y, H, 'go')
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