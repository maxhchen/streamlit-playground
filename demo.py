import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import altair as alt

# HW 5 Demo

from math import log10
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

st.title("Damping in RLC Circuits")

## CHANGELOG ##
# 2019/02/13: initial version (code based on slider_demo.py in the matplotlib documentation). JR <jr@berkeley.edu>
# 2019/09: Adapted for Fall 2019
# 2021/05: Streamlit

j=1j

class structtype: # MATLAB-like struct, useful for organization
    pass

## SETUP
parms = structtype()


cond = st.sidebar.radio("Damping Condition:",
('Playground', 'Undamped', 'Underdamped', 'Critically Damped', 'Overdamped'))

## DEFAUT
parms.R = 1
parms.C = 10e-9 # 10nF
parms.L = 25e-6 # 25uH
parms.f0 = 1/np.sqrt(parms.L*parms.C)/2/np.pi
parms.tauC = parms.R*parms.C
parms.tauL = parms.L/parms.R
parms.alpha = parms.R/(2*parms.L)
Cmin = parms.C/3 * 10 ** 8
Cmax = 3*parms.C * 10 ** 8

## Undamped
if cond == 'Undamped':
    parms.R = 0
    parms.C = 10e-9 # 10nF
    parms.L = 25e-6 # 25uH
    parms.f0 = 1/np.sqrt(parms.L*parms.C)/2/np.pi
    parms.tauC = parms.R*parms.C
    parms.tauL = np.inf # parms.L/parms.R
    parms.alpha = parms.R/(2*parms.L)
    Cmin = parms.C/3 * 10 ** 8
    Cmax = 3*parms.C * 10 ** 8

## Overdamped
if cond == 'Overdamped':
    parms.R = 1000
    parms.C = 10e-9 # 10nF
    parms.L = 25e-6 # 25uH
    parms.f0 = 1/np.sqrt(parms.L*parms.C)/2/np.pi
    parms.tauC = parms.R*parms.C
    parms.tauL = parms.L/parms.R
    parms.alpha = parms.R/(2*parms.L)
    Cmin = parms.C/3 * 10 ** 8
    Cmax = 3*parms.C * 10 ** 8

## Underdamped
if cond == 'Underdamped':
    parms.R = 1
    parms.C = 10e-9 # 10nF
    parms.L = 25e-6 # 25uH
    parms.f0 = 1/np.sqrt(parms.L*parms.C)/2/np.pi
    parms.tauC = parms.R*parms.C
    parms.tauL = parms.L/parms.R
    parms.alpha = parms.R/(2*parms.L)
    Cmin = parms.C/3 * 10 ** 8
    Cmax = 3*parms.C * 10 ** 8

## Critically Damped
if cond == 'Critically Damped':
    parms.C = 10e-9 # 10nF
    parms.L = 25e-6 # 25uH
    parms.R = 2 * np.sqrt(parms.L / parms.C)
    parms.f0 = 1/np.sqrt(parms.L*parms.C)/2/np.pi
    parms.tauC = parms.R*parms.C
    parms.tauL = parms.L/parms.R
    parms.alpha = parms.R/(2*parms.L)
    Cmin = parms.C/3 * 10 ** 8
    Cmax = 3*parms.C * 10 ** 8

# st.write(parms.R)
# st.write(parms.C)
# st.write(parms.L)

def vc_of_t(t, initcond): # solution for vc of t given an initial condition
    """
        vc(t) = K1 e^{lambda_1 t} + K2 e^{lambda_1 t}, where:
            lambda1/2 = -alpha +- sqrt(alpha^2 - 4 pi^2 f0^2)
            K1 = vc(0)/(1-lambda1/lambda2)
            K2 = vc(0) - K1
    """
    sqrtterm = np.sqrt(parms.alpha**2 - 4 * np.pi**2 * parms.f0**2 + 0j)
    lambda1 = -parms.alpha + sqrtterm
    lambda2 = -parms.alpha - sqrtterm
    K1 = initcond/(1-lambda1/lambda2)
    K2 = initcond - K1
    retval = np.real(K1*np.exp(lambda1*t) + K2*np.exp(lambda2*t))
    return(retval)

def vcenvelope_of_t(t, initcond): # approximate envelope for vc of t given an initial condition
    """
        vc(t) = K1 e^{lambda_1 t} + K2 e^{lambda_1 t}, where:
            lambda1/2 = -alpha +- sqrt(alpha^2 - 4 pi^2 f0^2)
            K1 = vc(0)/(1-lambda1/lambda2)
            K2 = vc(0) - K1
    """
    if (4 * np.pi**2 * parms.f0**2 < parms.alpha**2):
        sqrtterm = np.sqrt(parms.alpha**2 - 4 * np.pi**2 * parms.f0**2 + 0j)
        lambda1 = -parms.alpha + sqrtterm
        lambda2 = -parms.alpha - sqrtterm
        K1 = initcond/(1-lambda1/lambda2)
        K2 = initcond - K1
        retval = np.real(0 + K1*np.exp(lambda1*t) + K2*np.exp(lambda2*t))
    else:
        sqrtterm = 0
        lambda1 = -parms.alpha
        K1 = initcond
        retval = np.real(K1*np.exp(lambda1*t)); # - lambda1*t*np.exp(lambda1*t)); 
    return(retval)

n_cycles = 20
pts_per_cycle = 50
N = pts_per_cycle*n_cycles
ts = 1/parms.f0*n_cycles*np.array(range(0,N+1))/N

vcs = vc_of_t(ts, initcond=1)
vcenvelope = vcenvelope_of_t(ts, initcond=1)

## PLOTTING

fig, ax = plt.subplots(); # like MATLAB's figure()
plt.subplots_adjust(top=0.75); # changes the margins of the plot

l, = plt.plot(ts, vcs, linewidth=1, color='red', linestyle='-')
lenv, = plt.plot(ts, vcenvelope, color='green', linewidth=3, linestyle=':')

plt.grid(color='k', linestyle='-.', linewidth=0.5)
plt.xlabel('time (seconds)'); # 
plt.ylabel('vc(t) (Volts)');
plt.title('vc(t), series RLC circuit');

N_C = 100
sC = st.sidebar.slider('C (10^8)', Cmin, Cmax); # streamlit slider
# if cond == 'Playground':
    # sC = st.sidebar.slider('C (10^8)', Cmin, Cmax); # streamlit slider
# else:
    # sC = st.sidebar.slider('C (10^8)', parms.C, parms.C); # streamlit slider
    # sC = parms.C

Rmin = 1e-2
Rmax = 1e4
N_R = 100
# sR = st.sidebar.slider('log(R)', log10(Rmin), log10(Rmax)); # streamlit slider
if cond == 'Playground':
    sR = st.sidebar.slider('log(R)', log10(Rmin), log10(Rmax)); # streamlit slider
elif cond == 'Undamped':
    sR = -np.inf
else:
    sR = log10(parms.R)

textstr = '\n'.join((
    r'C=%g' % (parms.C, ),
    r'R=%g' % (parms.R, ) ))

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textH = ax.text(0.0, 1.2, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

## UPDATE PLOT

C = sC / 10 ** 8
R = 10 ** sR

parms.R = R
parms.C = C
parms.tauC = parms.R*parms.C
if cond == 'Undamped':
    parms.tauL = np.inf
else:
    parms.tauL = parms.L/parms.R
parms.alpha = parms.R/(2*parms.L)
parms.f0 = 1/(2*np.pi*np.sqrt(parms.L*parms.C))

vcs = vc_of_t(ts, 1)
vcenvelope = vcenvelope_of_t(ts,1)

l.set_ydata(vcs)
lenv.set_ydata(vcenvelope)
ax.relim(True)
ax.autoscale_view()

textstr = '\n'.join((
    r'C=%g' % (parms.C, ),
    r'R=%g' % (parms.R, ) ))
textH.set_text(textstr)

st.pyplot(fig)