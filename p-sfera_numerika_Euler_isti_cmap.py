# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 11:54:04 2019

@author: Admin
"""

import numpy as np
from numpy import sin, cos, reshape, transpose
from scipy.integrate import quadrature
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
import time
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import scipy
plt.rcParams.update({'font.size': 14})

pratim_te = 0 # brojac
z_da_ne = 0 # 0 ili 1
primjer = 2
parametri = 1

p = 4
n = 8
T = 10

if(parametri == 1):
    a = 1
    b = 2
    L = (b**3 - a**3) / 3
    R = 1.
    c_d = 1.
    c_0 = 1.
    c_v = 1.
    j_I = 1.
    #c_a = 0
    lam = -2
    mi_r = 1.
    mi = 3
    kappa = 0.024
    d = z_da_ne*1 # ?
    delta = z_da_ne*1 # ?
    
    def fja_r(v1, v2, v3):
        eps = 0.2 # 10
        m = 2
        if(z_da_ne):
            return eps * (v1**(m-1)) * (v3**m) * np.exp((v2-1.)/(eps * v2))
        else:
            return 0.
elif(parametri == 2):
    a = 1.
    b = 2.
    L = (b**3 - a**3) / 3.
    R = 10.
    c_d = 3.
    c_0 = -2.
    c_v = 20.
    j_I = 8.
    #c_a = 0
    lam = -2
    mi_r = 5.
    mi = 3.
    kappa = 5.
    d = z_da_ne*1 # ?
    delta = z_da_ne*1 # ?
    
    def fja_r(v1, v2, v3):
        eps = 0.2 # 10
        m = 2
        if(z_da_ne):
            return eps * (v1**(m-1)) * (v3**m) * np.exp((v2-1.)/(eps * v2))
        else:
            return 0.
elif(parametri == 3):
    a = 1.
    b = 2.
    L = (b**3 - a**3) / 3.
    R = 10.
    c_d = 3.
    c_0 = -2.
    c_v = 20.
    j_I = 8.
    #c_a = 0
    lam = -2
    mi_r = 15.
    mi = 5.
    kappa = 5.
    d = z_da_ne*1 # ?
    delta = z_da_ne*1 # ?
    
    def fja_r(v1, v2, v3):
        eps = 0.2 # 10
        m = 2
        if(z_da_ne):
            return eps * (v1**(m-1)) * (v3**m) * np.exp((v2-1.)/(eps * v2))
        else:
            return 0.

indeksi = np.array(list(range(n+1)))
indeksi2 = np.array(list(range(n+1))).reshape(n+1, -1)


# INTEGRACIJA
deg = 20
cvorovi, tezine = leggauss(deg)
cvorovi = 0.5 * (cvorovi + 1)
tezine = 0.5 * tezine

#------------------------------------------------------------------------------
# POCETNI UVJETI

## Primjer 1 -----------------------------------
if(primjer == 1):
    def ro_0(x):
        return 1.
    def int1krp_0(x):
        return x

    def r_0(x):
        return np.cbrt(a**3 + 3 * L * int1krp_0(x))
    def dr_0(x):
        return L / ( ro_0(x) * (r_0(x))**2)

    v_pocetni = np.zeros(n)
    v_pocetni[0] = 1.

    omega_pocetni = np.zeros(n)
    omega_pocetni[1] = 1.

    teta_pocetni = np.zeros(n+1)
    teta_pocetni[0] = 2.
    teta_pocetni[1] = 1.

    z_pocetni = np.zeros(n+1)
    z_pocetni[0] = z_da_ne*1.

    q_pocetni = np.zeros(n)
    q2_pocetni = np.zeros(n**2)
    q3_pocetni = np.zeros(n**3)    

## Primjer 2 -------------------------------------
#Pr1 iz rad27
elif(primjer == 2):
    def ro_0(x):
        return 1.
    def int1krp_0(x):
        return x
    
    def r_0(x):
        return np.cbrt(a**3 + 3 * L * int1krp_0(x))
    def dr_0(x):
        return 3 * L / (3 * ro_0(x) * (r_0(x))**2)
    
    v_pocetni = np.zeros(n)
    def v_poc_f(x, i):
        return (x-x*x) * (2. * sin(np.pi * i * x))
    for i in range(n):
        v_pocetni[i] = quadrature(v_poc_f, 0, 1, i+1)[0]
    
    omega_pocetni = np.zeros(n)
    def om_poc_f(x, i):
        return (x*x*(1-x*x)) * (2. * sin(np.pi * i * x))
    for i in range(n):
        omega_pocetni[i] = quadrature(om_poc_f, 0, 1, i+1)[0]
    
    teta_pocetni = np.zeros(n+1)
    def teta_poc_f(x, i):
        return (2.*x*x*x - 3.*x*x + 2.) * (2. * cos(np.pi * i * x))
    for i in range(n+1):
        teta_pocetni[i] = quadrature(teta_poc_f, 0, 1, i)[0]
    teta_pocetni[0] = teta_pocetni[0]/2.
    
    z_pocetni = np.zeros(n+1)
    # def z_poc_f(x, i):
    #     return (2.*x*x*x - 3.*x*x + 1.) * (2. * cos(np.pi * i * x))
    # for i in range(n+1):
    #     z_pocetni[i] = quadrature(z_poc_f, 0, 1, i)[0]
    # z_pocetni[0] = z_pocetni[0]/2.
    #z_pocetni[0] = 1.
    
    q_pocetni = np.zeros(n)
    q2_pocetni = np.zeros(n**2)
    q3_pocetni = np.zeros(n**3)

# ## Primjer 3 -----------------------------------
# elif(primjer == 3):
#     def ro_0(x):
#         return np.abs(x*x-0.25)+ 1.
#     def jkro_0(x):
#         return 1./(np.abs(x*x-0.25)+ 1.)
#     def int1krp_0(x):
#         return quadrature(jkro_0, 0, x)[0]
    
#     def r_0(x):
#         return np.cbrt(a**3 + 3 * L * int1krp_0(x))
#     def dr_0(x):
#         return 3 * L / (3 * ro_0(x) * (r_0(x))**2)
    
#     v_pocetni = np.zeros(n)
#     def v_poc_f(x, i):
#         return 0. * (2. * sin(np.pi * i * x))
#     for i in range(n):
#         v_pocetni[i] = quadrature(v_poc_f, 0, 1, i+1)[0]
            
#     omega_pocetni = np.zeros(n)
#     def om_poc_f(x, i):
#         return 4.*(x*x*(1-x*x)) * (2. * sin(np.pi * i * x))
#     for i in range(n):
#         omega_pocetni[i] = quadrature(om_poc_f, 0, 1, i+1)[0]
            
#     teta_pocetni = np.zeros(n+1)
#     def teta_poc_f(x, i):
#         return 0.1 * (2. * cos(np.pi * i * x))
#     for i in range(n+1):
#         teta_pocetni[i] = quadrature(teta_poc_f, 0, 1, i)[0]
#     teta_pocetni[0] = teta_pocetni[0]/2.
        
#     z_pocetni = np.zeros(n+1)
#     def z_poc_f(x, i):
#         vrati = np.zeros(np.shape(x))
#         jj = 0
#         for xx in x:
#             if(xx < 0.2):
#                 vrati[jj] = 1. * (2. * cos(np.pi * i * xx))
#             elif(xx > 0.7):
#                 vrati[jj] = 0. * (2. * cos(np.pi * i * xx))
#             else:
#                 vrati[jj] = (-2*xx+1.4) * (2. * cos(np.pi * i * xx))
#             jj = jj + 1
#         return vrati
#     for i in range(n+1):
#         z_pocetni[i] = quadrature(z_poc_f, 0, 1, i, maxiter = 1000)[0]
#     z_pocetni[0] = z_pocetni[0]/2.
            
#     q_pocetni = np.zeros(n)
#     q2_pocetni = np.zeros(n**2)
#     q3_pocetni = np.zeros(n**3)
# Fromiraj pocetni----------------------------------------------------

pocetni = np.concatenate((v_pocetni, omega_pocetni, teta_pocetni, z_pocetni, q_pocetni, q2_pocetni, q3_pocetni))

#------------------------------------------------------------------------------
# JEDNADZBE
def r(x, t, q):
	return r_0(x) + np.sum( q * np.sin(np.pi * indeksi[1:] * x))
def dr(x, t, q):
	return dr_0(x) + np.pi * np.sum( q * indeksi[1:] * np.cos(np.pi * indeksi[1:] * x))
def ro(x, t, q, q2, q3):
    n = np.size(q)
    p1 = 2 * r_0(x) * dr_0(x) * np.sum(q * np.sin(np.pi * indeksi[1:] * x)) + r_0(x)**2  * np.pi * np.sum(q * indeksi[1:] * np.cos(np.pi * indeksi[1:] * x))
    p2 = 2 * dr_0(x) * np.sum( reshape(q2, (n, -1)) * np.sin(np.pi * indeksi2[1:] * x) * transpose(np.sin(np.pi * indeksi2[1:] * x))) + 2 * r_0(x) * np.pi * np.sum( reshape(q2, (n, -1)) * (indeksi2[1:] * np.cos(np.pi * indeksi2[1:] * x)) * transpose(np.sin(np.pi * indeksi2[1:] * x)) ) + 2 * r_0(x) * np.pi * np.sum( reshape(q2, (n, -1)) * np.sin(np.pi * indeksi2[1:] * x) * transpose(indeksi2[1:] * np.cos(np.pi * indeksi2[1:] * x)) )
    p3 = np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi2[1:] * np.cos(np.pi * indeksi2[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (1, n, 1)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (1, 1, n)))
    p3 += np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi2[1:] * np.cos(np.pi * indeksi2[1:] * x), (1, n, 1)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (1, 1, n)))
    p3 += np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi2[1:] * np.cos(np.pi * indeksi2[1:] * x), (1, 1, n)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi2[1:] * x), (1, n, 1)))
    return L * ro_0(x) / (L + ro_0(x) * (p1+p2+np.sum(p3)))
def dv(x, t, v):
    return np.pi * np.sum(indeksi[1:] * v * np.cos(np.pi * x*indeksi[1:]))
def w(x, t, omega):
    return np.sum(omega * np.sin(np.pi * indeksi[1:] * x))
def dw(x, t, omega):
    return np.pi * np.sum(indeksi[1:] * omega * np.cos(np.pi * x*indeksi[1:]))
def th(x, t, teta):
    return np.sum(teta * np.cos(np.pi * x * indeksi))
def dth(x, t, teta):
    return -np.pi * np.sum(indeksi * teta * sin(np.pi*x*indeksi))
def zz(x, t, z):
    return z_da_ne*np.sum(z * np.cos(np.pi*indeksi*x))
def dz(x, t, z):
    return -z_da_ne*np.pi * np.sum(indeksi * z * sin(np.pi*x*indeksi))
def lamb(k):
    if(k):
        return 2.
    return 1.

def v_podintegralna(x, t, i, v, teta, q, q2, q3):
    rn = r(x, t, q)
    drn = dr(x, t, q)
    vn = dv(x, t, v)
    ron = ro(x, t, q, q2, q3)
    thn = th(x, t, teta)
    dvn = dv(x, t, v)
    return  2 * (R/L * ron**p * thn - (lam+2*mi)/L**2 * (ron * 2 * rn * drn* vn + ron * rn**2 * dvn)) * (2 * rn * drn * sin(np.pi * i * x) + rn**2 *  np.pi * i * cos(np.pi * i * x))
def v_jednadzba(t, i, v, teta, q, q2, q3):
    zbroj = 0.
    for it in range (deg):
        zbroj = zbroj + tezine[it] * v_podintegralna(cvorovi[it], t, i, v, teta, q, q2, q3)
    return zbroj

def omega_podintegralna(x, t, i, omega, q, q2, q3):
    rn = r(x, t, q)
    drn = dr(x, t, q)
    ron = ro(x, t, q, q2, q3)
    wn = w(x, t, omega)
    dwn = dw(x, t, omega)
    return 2/j_I * ( -(c_0+2*c_d)/L**2 * ron * (2 * rn * drn * wn + rn**2 * dwn) * (2 * rn * drn * sin(np.pi * i * x) + rn**2 * np.pi * i * cos(np.pi * i * x)) - 4*mi_r * wn / ron * sin(np.pi * i * x))
def omega_jednadzba(t, i, omega, q, q2, q3):
    zbroj = 0.
    for it in range (deg):
        zbroj = zbroj + tezine[it] * omega_podintegralna(cvorovi[it], t, i, omega, q, q2, q3)
    return zbroj

def teta_podintegralna(x, t, i, v, omega, teta, z, q, q2, q3):
    rn = r(x, t, q)
    drn = dr(x, t, q)
    ron = ro(x, t, q, q2, q3)
    dthn = dth(x, t, teta)
    thn = th(x, t, teta)
    vn = dv(x, t, v)
    dvn = dv(x, t, v)
    wn = w(x, t, omega)
    dwn = dw(x, t, omega)
    zn = zz(x, t, z)
    # - sin ili + sin u prvom dijelu u redu ispod?
    return lamb(i)/c_v * ( (-kappa/L**2 * rn**4 * ron * dthn + 4/L * rn * (mi*vn**2 + c_d*wn**2)) * (-np.pi * i * sin(np.pi * i * x)) + ( -R/L * ron**p * thn * (2*rn*drn*vn + rn**2 * dvn) + (lam+2*mi)/L**2 * ron * (2*rn*drn*vn + rn**2 * dvn)**2 + 4*mi_r * wn**2 / ron + (c_0+2*c_d)/L**2 * ron * (2*rn*drn*wn + rn**2 * dwn)**2 ) * cos(np.pi * i * x) )
def teta_jednadzba(t, i, v, omega, teta, z, q, q2, q3):
    zbroj = 0.
    for it in range (deg):
        zbroj = zbroj + tezine[it] * teta_podintegralna(cvorovi[it], t, i, v, omega, teta, z, q, q2, q3)
    return zbroj

def z_podintegralna(x, t, i, teta, z, q, q2, q3):
    ron = ro(x, t, q, q2, q3)
    thn = th(x, t, teta)
    zn = zz(x, t, z)
    dzn = dz(x, t, z)
    return z_da_ne*lamb(i) * (-d/L**2 * ron**2 * dzn * (np.pi * i * sin(np.pi * i * x)) - fja_r(ron, thn, zn) * cos(np.pi * i * x))
def z_jednadzba(t, i, teta, z, q, q2, q3):
    zbroj = 0.
    for it in range (deg):
        zbroj = zbroj + tezine[it] * z_podintegralna(cvorovi[it], t, i, teta, z, q, q2, q3)
    return z_da_ne*zbroj

def q_jednadzba(t, i, v): # i = 1, ..., n
    return v[i-1]

def q2_jednadzba(t, i, q, v):
    i1 = (i-1)//n
    i2 = (i-1)%n
    return q[i1] * v[i2]

def q3_jednadzba(t, i, q, v):
    i1 = (i-1)//(n**2)
    n2 = (i-1)%(n**2)
    i2 = n2//n
    i3 = n2%n
    return q[i1] * q[i2] * v[i3]


def F(t, y):

    v = y[0 : n]
    omega = y[n : 2*n]
    teta = y[2*n : 3*n+1]
    z = y[3*n+1 : 4*n+2]
    q = y[4*n+2 : 5*n+2]
    q2 = y[5*n+2 : 5*n+2 + n**2]
    q3 = y[5*n+2 + n**2 : ]

    Y = np.empty(np.shape(y))
    for it in range(0, n):
        Y[it] = v_jednadzba(t, it+1, v, teta, q, q2, q3)
    for it in range(n, 2*n):
        Y[it] = omega_jednadzba(t, (it-n)+1, omega, q, q2, q3)
    for it in range(2*n, 3*n+1):
        Y[it] = teta_jednadzba(t, (it-2*n), v, omega, teta, z, q, q2, q3)
    for it in range(3*n+1, 4*n+2):
        Y[it] = z_jednadzba(t, (it-(3*n+1)), teta, z, q, q2, q3)
    for it in range(4*n+2, 5*n+2):
        Y[it] = q_jednadzba(t, (it-(4*n+2))+1, v)
    for it in range(5*n+2, 5*n+2 + n**2):
        Y[it] = q2_jednadzba(t, (it-(5*n+2))+1, q, v)
    for it in range(5*n+2 + n**2, np.size(Y)):
        Y[it] = q3_jednadzba(t, (it-(5*n+2 + n**2))+1, q, v)
        
    global pratim_te
    if(pratim_te%500==0):
        print('iter = ', pratim_te)
        print('t=', t, ', v=', Y[int(n/2)])
    pratim_te +=1
        
    return Y



#------------------------------------------------------------------------------
# RJESENJE py

#t_rj = np.loadtxt('t_rj_veeeliko.txt') # rjesenje.t
#y_rj = np.loadtxt('y_rj_veeeliko.txt') # rjesenje.y

start = time.time()
rjesenje = solve_ivp(F, (0,T), pocetni, method = 'BDF', dense_output = True)
print('time = ', time.time() - start)
t_rj = rjesenje.t
y_rj = rjesenje.y
#tplot = np.linspace(0,1, 500)
#yni = rjesenje.sol(tplot)
#print(np.shape(yni))

#------------------------------------------------------------------------------
# np.savetxt('t.txt', t_rj, delimiter=' ')
# np.savetxt('y.txt', y_rj, delimiter=' ')
v_rj = y_rj[0 : n, :]
omega_rj = y_rj[n : 2*n, :]
teta_rj = y_rj[2*n : 3*n+1, :]
z_rj = y_rj[3*n+1 : 4*n+2, :]
q_rj = y_rj[4*n+2 : 5*n+2, :]
q2_rj = y_rj[5*n+2 : 5*n+2+n**2, :]
q3_rj = y_rj[5*n+2+n**2 : , :]

#------------------------------------------------------------------------------
# PRIKAZ
t_ind = np.size(t_rj)-1
bla = [0., 0.1, 0.5, 1]#[0, 4, 4.4, 4.7, 4.9, 5] # vremenski trenuci koji se prikazuju
nbla = np.size(bla)

## ro -------------------------------------------------------------------------
def ro_rj_f(x, t_ind):
    n = np.size(q_rj[:, t_ind])
    zbroj = 0
    for it in range(0, n):
        zbroj += 2 * r_0(x) * dr_0(x) * np.sum(q_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)) + r_0(x)**2  * np.pi * np.sum(q_rj[it, t_ind] * (it+1) * np.cos(np.pi * (it+1) * x))
        for it2 in range(0, n):
            zbroj += 2 * dr_0(x) * np.sum( q2_rj[it*n+it2, t_ind] * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x)) + 2 * r_0(x) * np.pi * np.sum( q2_rj[it*n+it2, t_ind] * ((it+1) * np.cos(np.pi * (it+1) * x)) * np.sin(np.pi * (it2+1) * x) ) + 2 * r_0(x) * np.pi * np.sum( q2_rj[it*n+it2, t_ind] * np.sin(np.pi * (it+1) * x) * (it2+1) * np.cos(np.pi * (it2+1) * x) )
            for it3 in range(0, n):
                zbroj += np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it+1) * np.cos(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x) * np.sin(np.pi * (it3+1) * x))
                zbroj += np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it2+1) * np.cos(np.pi * (it2+1) * x) * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it3+1) * x))
                zbroj += np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it3+1) * np.cos(np.pi * (it3+1) * x) * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x))

    return L * ro_0(x) / (L + ro_0(x) * zbroj)

def interp_ro_rj_f(x, t):
    zbroj = 0.
    for it in range(0, n):
        zbroj += (2*r_0(x)*dr_0(x)*(rjesenje.sol(t))[4*n+2+it]*sin(np.pi*(it+1)*x) + r_0(x)**2 *(rjesenje.sol(t))[4*n+2+it]*np.pi*(it+1)*cos(np.pi*(it+1)*x))
    for it in range(0, n*n):
        i1 = it//n
        i2 = it%n
        zbroj += 2*dr_0(x)*(rjesenje.sol(t))[5*n+2+it]*sin(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)
        zbroj += 2*r_0(x)*(rjesenje.sol(t))[5*n+2+it]*cos(np.pi*(i1+1)*x)*np.pi*(i1+1)*sin(np.pi*(i2+1)*x)
        zbroj += 2*r_0(x)*(rjesenje.sol(t))[5*n+2+it]*sin(np.pi*(i1+1)*x)*cos(np.pi*(i2+1)*x)*np.pi*(i2+1)
    for it in range(0, n**3):
        i1 = it//(n**2)
        n2 = it%(n**2)
        i2 = n2//n
        i3 = n2%n
        zbroj += (rjesenje.sol(t))[5*n+2+n**2+it]*cos(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)*sin(np.pi*(i3+1)*x)*np.pi*(i1+1)
        zbroj += (rjesenje.sol(t))[5*n+2+n**2+it]*sin(np.pi*(i1+1)*x)*cos(np.pi*(i2+1)*x)*sin(np.pi*(i3+1)*x)*np.pi*(i2+1)
        zbroj += (rjesenje.sol(t))[5*n+2+n**2+it]*sin(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)*cos(np.pi*(i3+1)*x)*np.pi*(i3+1)
    return L * ro_0(x) / (L + ro_0(x) * zbroj)

fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
xplot = np.linspace(0,1, 200)
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_ro_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r'$\rho^n(\cdot, t)$')
plt.show()
fig.savefig('slike/ro_t_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = ro_rj_f(xplot, it)
fig = plt.figure(figsize=(19, 15))
plt.contourf(t_rj, xplot, Z, 500,  cmap='viridis')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$\rho$")
plt.colorbar();
plt.show()
fig.savefig('slike/ro_k_pr' + str(primjer) + '.png')


## v --------------------------------------------------------------------------
def v_rj_f(x, t_ind):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + v_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)
    return zbroj

def interp_v_rj_f(x, t):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + (rjesenje.sol(t))[it] * np.sin(np.pi * (it+1) * x)
    return zbroj

fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
xplot = np.linspace(0,1, 200)
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_v_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r"$v^n(\cdot, t)$")
plt.show()
fig.savefig('slike/v_t_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = v_rj_f(xplot, it)
fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
plt.contourf(t_rj, xplot, Z, 500,  cmap='PiYG')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$v$")
plt.colorbar();
plt.show()
fig.savefig('slike/v_k_pr' + str(primjer) + '.png')

## omega ----------------------------------------------------------------------
def omega_rj_f(x, t_ind):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + omega_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)
    return zbroj

def interp_omega_rj_f(x, t):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + (rjesenje.sol(t))[n+it] * np.sin(np.pi * (it+1) * x)
    return zbroj

fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
xplot = np.linspace(0,1, 200)
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_omega_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r"$\omega^n(\cdot, t)$")
plt.show()
fig.savefig('slike/omega_t_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = omega_rj_f(xplot, it)
fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
plt.contourf(t_rj, xplot, Z, 500,  cmap='plasma')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$\omega$")
plt.colorbar();
plt.show()
fig.savefig('slike/omega_k_pr' + str(primjer) + '.png')

## teta -----------------------------------------------------------------------
def teta_rj_f(x, t_ind):
    zbroj = 0.
    for it in range(0, n+1):
        zbroj = zbroj + teta_rj[it, t_ind] * np.cos(np.pi * it * x)
    return zbroj

def interp_teta_rj_f(x, t):
    zbroj = 0.
    for it in range(0, n+1):
        zbroj = zbroj + (rjesenje.sol(t))[2*n+it] * np.cos(np.pi * it * x)
    return zbroj

fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_teta_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r"$\theta^n(\cdot, t)$")
plt.show()
fig.savefig('slike/teta_t_pr' + str(primjer) + '.png')

fig = plt.figure(figsize=(19, 15))
ttplot = np.linspace(0, T, 100)
xtplot = np.linspace(0.975, 1, 9)
xtplot[-2] = 0.9795918367346939
# plt.title.set_size(20)
for i in range(np.size(xtplot)):
    plt.subplot(3, 3, i+1)
    plt.plot(ttplot, interp_teta_rj_f(xtplot[i], ttplot), 'k--', label = "x = " + str(xtplot[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('t')
plt.suptitle(r"$\theta^n(\cdot, t)$")
plt.show()
fig.savefig('slike/teta_tt_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = teta_rj_f(xplot, it)
    for ix in range(np.size(xplot)):
        if(Z[ix, it] < 0):
            print(xplot[ix], t_rj[it], Z[ix, it])
fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
plt.contourf(t_rj, xplot, Z, 500,  cmap='spring')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$\theta$")
plt.colorbar();
plt.show()
fig.savefig('slike/teta_k_pr' + str(primjer) + '.png')

## z -----------------------------------------------------------------------
def z_rj_f(x, t_ind):
    zbroj = 0.
    for it in range(0, n+1):
        zbroj = zbroj + z_rj[it, t_ind] * np.cos(np.pi * it * x)
    return zbroj

def interp_z_rj_f(x, t):
    zbroj = 0.
    for it in range(0, n+1):
        zbroj = zbroj + (rjesenje.sol(t))[3*n+1+it] * np.cos(np.pi * it * x)
    return zbroj

fig = plt.figure(figsize=(19, 15))
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_z_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r"$z^n(\cdot, t)$")
plt.show()
fig.savefig('slike/z_t_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = z_rj_f(xplot, it)
fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
plt.contourf(t_rj, xplot, Z, 500,  cmap='cool')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$z$")
plt.colorbar();
plt.show()
fig.savefig('slike/z_k_pr' + str(primjer) + '.png')

## tlak ----------------------------------------------------------------------
def interp_fja_ro(x, t):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + (rjesenje.sol(t))[4*n+2+it] * np.pi * (it+1) * np.cos(np.pi * (it+1) * x)
    return (L*ro_0(x) /(L + ro_0(x) * zbroj))

def interp_fja_teta(x, t):
    zbroj = 0.
    for it in range(0, n+1):
        zbroj = zbroj + (rjesenje.sol(t))[2*n+it] * np.cos(np.pi * it * x)
    return zbroj
def tlak(x, t):
    return R * (interp_fja_ro(x, t))**p * interp_fja_teta(x, t)

koliko = 64+1
xplot = np.linspace(0, 1, koliko)
tplot = np.linspace(0, T, koliko)

plak = np.zeros((np.size(xplot), np.size(tplot)))
for it in range(np.size(tplot)):
    for ix in range(np.size(xplot)):
        plak[ix, it] = tlak(xplot[ix], tplot[it])


fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
plt.title('P')
plt.contourf(tplot, xplot, plak, 500, cmap = 'hot')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
# plt.title(r"$P$")
cbar = plt.colorbar()
fig.savefig('slike/tlak_k_pr' + str(primjer) + '.png')
print('plak', np.min(plak), np.max(plak))

## micajuće -------------------------------------------------------------------
#fig = plt.figure(100)
#xplot = np.linspace(0, 1, 200)
#for ind in range(0, np.size(t_rj), 10):
#    plt.cla()
#    plt.ylim([-0.1, 0.1])
#    plt.xlim([0,1])
#    plt.plot(xplot, v_rj_f(xplot, ind), label = "t = " + str(t_rj[ind]))
#    plt.title(str(ind) + " / " + str(np.size(t_rj)) )
#    plt.legend()
#    plt.grid()
#    plt.pause(0.3)



#------------------------------------------------------------------------------
# v funkcija
def interp_fja_v(x, t):
    zbroj = 0.
    for it in range(0, n):
        zbroj = zbroj + (rjesenje.sol(t))[it] * np.sin(np.pi * (it+1) * x)
    return zbroj

## pomak
def integr_v(t, x):
    return interp_fja_v(x, t)
def pomak (x, t):
    cvorovip = t * 0.5 * (cvorovi + 1)
    tezinep = t * 0.5 * tezine
    return np.sum (tezinep * integr_v(cvorovip, x))

koliko = 65
xplot = np.linspace(0, 1, koliko)
tplot = np.linspace(0, T, koliko)
#Z = np.zeros((np.size(xplot), np.size(tplot)))
#for it in range(np.size(tplot)):
#    for ix in range(np.size(xplot)):
#        Z[ix, it] = pomak(xplot[ix], tplot[it])

## prikaz nekih t-ova
fig = plt.figure(figsize=(19, 15))
# plt.title.set_size(20)
ind = 0
boja = ['blue', 'red', 'green', 'yellow', 'black', 'gray', 'purple', 'chocolate', 'orange', 'pink']
for t in [T*0.00001, T*0.0001, T*0.0015, T*0.002, T*0.01, T*0.08, T*0.2, T*0.37, T*0.6, T*1]:
    yplot = np.zeros(np.size(xplot))
    for ix in range(np.size(xplot)):
        yplot[ix] = pomak(xplot[ix], t)
    plt.plot(yplot, xplot, boja[ind], label = "t = " + str(t))
    ind+=1
plt.title('Displacement from the initial position')
plt.ylabel('x')
plt.grid()
plt.legend()
plt.show()
fig.savefig('slike/pomak_pr' + str(primjer) + '.png')


## r -------------------------------------------------------------------------
def r_rj_f(x, t_ind):
    n = np.size(q_rj[:, t_ind])
    zbroj = r_0(x)
    for it in range(0, n):
        zbroj+= np.sum(q_rj[it, t_ind] * np.sin(np.pi * (it+1) * x))
    return zbroj

def interp_r_rj_f(x, t):
    zbroj = r_0(x)
    for it in range(0, n):
        zbroj += ((rjesenje.sol(t))[4*n+2+it]*sin(np.pi*(it+1)*x) )
    return zbroj

fig = plt.figure(figsize=(19, 15))
#plt.title.set_size(20)
xplot = np.linspace(0,1, 200)
for i in range(nbla):
    plt.subplot(2, 2, i+1)
    plt.plot(xplot, interp_r_rj_f(xplot, T*bla[i]), 'k--', label = "t = " + str(T*bla[i]))
    plt.legend()
    plt.grid()
    plt.xlabel('x')
plt.suptitle(r'$r^n(\cdot, t)$')
plt.show()
fig.savefig('slike/r_t_pr' + str(primjer) + '.png')

xplot = np.linspace(0, 1)
Z = np.zeros((np.size(xplot), np.size(t_rj)))
for it in range(0, np.size(t_rj)):
    Z[:, it] = r_rj_f(xplot, it)
fig = plt.figure(figsize=(19, 15))
plt.contourf(t_rj, xplot, Z, 500,  cmap='viridis')
#plt.clim(0, 0.001)
plt.xlabel('t')
plt.ylabel('x')
plt.title(r"$r$")
plt.colorbar();
plt.show()
fig.savefig('slike/r_k_pr' + str(primjer) + '.png')

# Euler
# interp_r_rj_f(x, t)
def ksi(x):
    return np.sqrt(3*L*x + a**3)

def iks(k):
    return (k**3-a**3)/(3*L)

from scipy.optimize import fsolve
def r_nult(k, t, r):
    return interp_r_rj_f(iks(k), t) - r
def iks_od_r(r, t):
    x0 = (a+b)/2
    return iks(fsolve(r_nult, x0, args =(t, r)))

def ro_kon(r, t):
    return interp_ro_rj_f(iks_od_r(r, t), t)
def v_kon(r, t):
    return interp_v_rj_f(iks_od_r(r, t), t)
def omega_kon(r, t):
    return interp_omega_rj_f(iks_od_r(r, t), t)
def teta_kon(r, t):
    return interp_teta_rj_f(iks_od_r(r, t), t)

# ro
fig = plt.figure(figsize=(19, 15))
#plt.title.set_size(20)
kovi = np.radians(np.linspace(0, 360))
rovi = np.linspace(a, b, 50)
rr, k = np.meshgrid(rovi, kovi)
Z = np.zeros(np.shape(rr) + (nbla, ))
for i in range(nbla):
    t = T*bla[i]
    for it in range(np.size(rovi)):
        Z[0, it, i] = ro_kon(rr[0, it], t)
        for it2 in range(np.size(kovi)):
            Z[it2, it, i] = Z[0, it, i]
levels = np.linspace(np.amin(Z), np.amax(Z))
for i in range(nbla):  
    t = T*bla[i]
    plt.subplot(2, 2, i+1, polar=True)
    plt.gca().title.set_text("t = " + str(t))      
    plt.gca().set_xticks([])
    plt.gca().set_rlim([0, b])
    plt.gca().set_rticks([a, b])
    plt.contourf(k, rr, Z[:, :, i], levels)
    plt.grid()
    plt.colorbar()
plt.suptitle(r'$\rho_E(\cdot, t)$')
plt.show()
fig.savefig('slike/ro_E_pr' + str(primjer) + '.png')

# v
fig = plt.figure(figsize=(19, 15))
for i in range(nbla):
    t = T*bla[i]
    plt.subplot(2, 2, i+1, polar=True)
    plt.gca().title.set_text("t = " + str(t))
    rr, k = np.meshgrid(rovi, kovi)
    Z = np.zeros(np.shape(rr))
    for it in range(np.size(rovi)):
        Z[0, it] = v_kon(rr[0, it], t)
        for it2 in range(np.size(kovi)):
            Z[it2, it] = Z[0, it]
    plt.gca().set_xticks([])
    plt.gca().set_rlim([0, b])
    plt.gca().set_rticks([a, b])
    plt.contourf(k, rr, Z)
    plt.grid()
    plt.colorbar()
plt.suptitle(r'$v_E(\cdot, t)$')
plt.show()
fig.savefig('slike/v_E_pr' + str(primjer) + '.png')


# omega
fig = plt.figure(figsize=(19, 15))
for i in range(nbla):
    t = T*bla[i]
    plt.subplot(2, 2, i+1, polar=True)
    plt.gca().title.set_text("t = " + str(t))
    rr, k = np.meshgrid(rovi, kovi)
    Z = np.zeros(np.shape(rr))
    for it in range(np.size(rovi)):
        Z[0, it] = omega_kon(rr[0, it], t)
        for it2 in range(np.size(kovi)):
            Z[it2, it] = Z[0, it]
    plt.gca().set_xticks([])
    plt.gca().set_rlim([0, b])
    plt.gca().set_rticks([a, b])
    plt.contourf(k, rr, Z)
    plt.grid()
    plt.colorbar()
plt.suptitle(r'$\omega_E(\cdot, t)$')
plt.show()
fig.savefig('slike/omega_E_pr' + str(primjer) + '.png')


# teta
fig = plt.figure(figsize=(19, 15))
#plt.title.set_size(20)
kovi = np.radians(np.linspace(0, 360))
rovi = np.linspace(a, b, 50)
rr, k = np.meshgrid(rovi, kovi)
Z = np.zeros(np.shape(rr) + (nbla, ))
for i in range(nbla):
    t = T*bla[i]
    for it in range(np.size(rovi)):
        Z[0, it, i] = teta_kon(rr[0, it], t)
        for it2 in range(np.size(kovi)):
            Z[it2, it, i] = Z[0, it, i]
levels = np.linspace(np.amin(Z), np.amax(Z))
for i in range(nbla):  
    t = T*bla[i]
    plt.subplot(2, 2, i+1, polar=True)
    plt.gca().title.set_text("t = " + str(t))      
    plt.gca().set_xticks([])
    plt.gca().set_rlim([0, b])
    plt.gca().set_rticks([a, b])
    plt.contourf(k, rr, Z[:, :, i], levels)
    plt.grid()
    plt.colorbar()
plt.suptitle(r'$\theta_E(\cdot, t)$')
plt.show()
fig.savefig('slike/teta_E_pr' + str(primjer) + '.png')

# micajuće
koliko = 5
tplot = np.linspace(0, T, koliko)
Z = np.zeros(np.shape(rr) + (koliko,))
for i in range(koliko):
    t = tplot[i]
    for it in range(np.size(rovi)):
        Z[0, it, i] = teta_kon(rr[0, it], t)
        for it2 in range(np.size(kovi)):
            Z[it2, it, i] = Z[0, it, i]
levels = np.linspace(np.amin(Z), np.amax(Z))

plt.subplots(subplot_kw={'projection': 'polar'})
t = tplot[0]
plt.gca().set_xticks([])
plt.gca().set_ylim([0, b])
plt.gca().set_yticks([a, b])
plt.contourf(k, rr, Z[:, :, i], levels,  cmap=plt.cm.bone)
plt.grid()
plt.colorbar()
plt.show()
plt.pause(0.01)
for i in range(koliko):
    t = tplot[i]
    plt.cla()
    plt.gca().title.set_text("t = " + str(t))
    plt.gca().set_xticks([])
    plt.gca().set_ylim([0, b])
    plt.gca().set_yticks([a, b])
    plt.contourf(k, rr, Z[:, :, i], levels, cmap=plt.cm.bone)
    plt.grid()
    plt.show()
    plt.pause(0.01)









