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

parametri = 1
primjer = 1
pratim_te = 0 # brojac
z_da_ne = 0 # 0 ili 1

T = 10
n = 8


if(parametri == 1):
    scale_r = 0.6
    scale_v = 0.002
    scale_o = 0.002
    scale_t = 0.1
    
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

elif(parametri == 3):
    scale_r = 0.1
    scale_v = 0.1
    scale_o = 0.1
    scale_t = 0.1
    
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




rjesenje = []

povi = [1., 1.5, 4., 10.]
for pi in range(np.size(povi)):
    p = povi[pi]
    
    indeksi = np.array(list(range(n+1)))
    
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
            return L / (ro_0(x) * (r_0(x))**2)
    
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
        p2 = 2 * dr_0(x) * np.sum( reshape(q2, (n, -1)) * np.sin(np.pi * indeksi[1:] * x) * transpose(np.sin(np.pi * indeksi[1:] * x))) + 2 * r_0(x) * np.pi * np.sum( reshape(q2, (n, -1)) * (indeksi[1:] * np.cos(np.pi * indeksi[1:] * x)) * transpose(np.sin(np.pi * indeksi[1:] * x)) ) + 2 * r_0(x) * np.pi * np.sum( reshape(q2, (n, -1)) * np.sin(np.pi * indeksi[1:] * x) * transpose(indeksi[1:] * np.cos(np.pi * indeksi[1:] * x)) )
        p3 = np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi[1:] * np.cos(np.pi * indeksi[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi[1:] * x), (1, n, 1)) * reshape(np.sin(np.pi * indeksi[1:] * x), (1, 1, n)))
        p3 += np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi[1:] * np.cos(np.pi * indeksi[1:] * x), (1, n, 1)) * reshape(np.sin(np.pi * indeksi[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi[1:] * x), (1, 1, n)))
        p3 += np.pi * reshape(q3, (n, n, -1)) * (reshape(indeksi[1:] * np.cos(np.pi * indeksi[1:] * x), (1, 1, n)) * reshape(np.sin(np.pi * indeksi[1:] * x), (n, 1, 1)) * reshape(np.sin(np.pi * indeksi[1:] * x), (1, n, 1)))
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
        return np.pi * np.sum(indeksi * teta * sin(np.pi*x*indeksi))
    def zz(x, t, z):
        return z_da_ne*np.sum(z * np.cos(np.pi*indeksi*x))
    def dz(x, t, z):
        return z_da_ne*np.pi * np.sum(indeksi * z * sin(np.pi*x*indeksi))
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
        if(p>1.):
            ronp = np.power(ron, p)
        else:
            ronp = ron
        return  2 * ((R/L * ronp * thn - (lam+2*mi)/L**2 * (ron * 2 * rn * drn* vn + ron * rn**2 * dvn)) * (2 * rn * drn * sin(np.pi * i * x) + rn**2 *  np.pi * i * cos(np.pi * i * x)))
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
        if(p>1.):
            ronp = np.power(ron, p)
        else:
            ronp = ron
        return lamb(i)/c_v * ( (-kappa/L**2 * rn**4 * ron * dthn + 4/L * rn * (mi*vn**2 + c_d*wn**2)) * (np.pi * i * sin(np.pi * i * x)) + ( -R/L * ronp * thn * (2*rn*drn*vn + rn**2 * dvn) + (lam+2*mi)/L**2 * ron * (2*rn*drn*vn + rn**2 * dvn)**2 + 4*mi_r * wn**2 / ron + (c_0+2*c_d)/L**2 * ron * (2*rn*drn*wn + rn**2 * dwn)**2 + delta * fja_r(ron, thn, zn) ) * cos(np.pi * i * x) )
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
            print(p, ' iter = ', pratim_te)
            print('t=', t, ', v=', Y[int(n/2)])
        pratim_te +=1    
        
        return Y
    
    
    
    #------------------------------------------------------------------------------
    # RJESENJE py
    
    #t_rj = np.loadtxt('t_rj_veeeliko.txt') # rjesenje.t
    #y_rj = np.loadtxt('y_rj_veeeliko.txt') # rjesenje.y
    
    start = time.time()
    rjesenje_p = solve_ivp(F, (0,T), pocetni, method = 'BDF', dense_output = True)
    rjesenje.append(rjesenje_p)
    print('time sol = ', time.time() - start, 'feval =', rjesenje_p.nfev)
    
    
    
    
    
#------------------------------------------------------------------------------
# PRIKAZ
plt.close('all')
boje = ['b', 'r', 'g', 'k', 'y']
crte = ['-', ':', '--', '-.', '-']
markeri = ["v", " ", " ", " ", " "]

tplot_size = 50
rovi_size = 50

Zr = np.zeros((tplot_size, rovi_size, 4))
Zv = np.zeros((tplot_size, rovi_size, 4))
Zo = np.zeros((tplot_size, rovi_size, 4))
Zt = np.zeros((tplot_size, rovi_size, 4))
scale_r = 0.6
scale_v = 0.002
scale_o = 0.002
scale_t = 0.1
    

for pi in range(np.size(povi)):
    start = time.time()
    p = povi[pi] 
    rjesenje_pi = rjesenje[pi]
    t_rj = rjesenje_pi.t
    y_rj = rjesenje_pi.y
    #tplot = np.linspace(0,1, 500)
    #yni = rjesenje_pi.sol(tplot)
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
    bla = [0, 0.1, 0.5, 1]#[0, 4, 4.4, 4.7, 4.9, 5] # vremenski trenuci koji se prikazuju
    nbla = np.size(bla)

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
            zbroj += ((rjesenje_pi.sol(t))[4*n+2+it]*sin(np.pi*(it+1)*x) )
        return zbroj

    ## ro -------------------------------------------------------------------------   
    def ro_rj_f(x, t_ind):
        n = np.size(q_rj[:, t_ind])
        zbroj = 0.
        for it in range(0, n):
            zbroj+= 2 * r_0(x) * dr_0(x) * np.sum(q_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)) + r_0(x)**2  * np.pi * np.sum(q_rj[it, t_ind] * (it+1) * np.cos(np.pi * (it+1) * x))
            for it2 in range(0, n):
                zbroj+= 2 * dr_0(x) * np.sum( q2_rj[it*n+it2, t_ind] * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x)) + 2 * r_0(x) * np.pi * np.sum( q2_rj[it*n+it2, t_ind] * ((it+1) * np.cos(np.pi * (it+1) * x)) * np.sin(np.pi * (it2+1) * x) ) + 2 * r_0(x) * np.pi * np.sum( q2_rj[it*n+it2, t_ind] * np.sin(np.pi * (it+1) * x) * (it2+1) * np.cos(np.pi * (it2+1) * x) )
                for it3 in range(0, n):
                    zbroj+=np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it+1) * np.cos(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x) * np.sin(np.pi * (it3+1) * x))
                    zbroj+= np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it2+1) * np.cos(np.pi * (it2+1) * x) * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it3+1) * x))
                    zbroj+= np.pi * q3_rj[it*n*n+it2*n+it3, t_ind] * ((it3+1) * np.cos(np.pi * (it3+1) * x) * np.sin(np.pi * (it+1) * x) * np.sin(np.pi * (it2+1) * x))
    
        return L * ro_0(x) / (L + ro_0(x) * zbroj)
    
    def interp_ro_rj_f(x, t):
        zbroj = 0.
        for it in range(0, n):
            zbroj += (2*r_0(x)*dr_0(x)*(rjesenje_pi.sol(t))[4*n+2+it]*sin(np.pi*(it+1)*x) + r_0(x)**2 *(rjesenje_pi.sol(t))[4*n+2+it]*np.pi*(it+1)*cos(np.pi*(it+1)*x))
        for it in range(0, n*n):
            i1 = it//n
            i2 = it%n
            zbroj += 2*dr_0(x)*(rjesenje_pi.sol(t))[5*n+2+it]*sin(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)
            zbroj += 2*r_0(x)*(rjesenje_pi.sol(t))[5*n+2+it]*cos(np.pi*(i1+1)*x)*np.pi*(i1+1)*sin(np.pi*(i2+1)*x)
            zbroj += 2*r_0(x)*(rjesenje_pi.sol(t))[5*n+2+it]*sin(np.pi*(i1+1)*x)*cos(np.pi*(i2+1)*x)*np.pi*(i2+1)
        for it in range(0, n**3):
            i1 = it//(n**2)
            n2 = it%(n**2)
            i2 = n2//n
            i3 = n2%n
            zbroj += (rjesenje_pi.sol(t))[5*n+2+n**2+it]*cos(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)*sin(np.pi*(i3+1)*x)*np.pi*(i1+1)
            zbroj += (rjesenje_pi.sol(t))[5*n+2+n**2+it]*sin(np.pi*(i1+1)*x)*cos(np.pi*(i2+1)*x)*sin(np.pi*(i3+1)*x)*np.pi*(i2+1)
            zbroj += (rjesenje_pi.sol(t))[5*n+2+n**2+it]*sin(np.pi*(i1+1)*x)*sin(np.pi*(i2+1)*x)*cos(np.pi*(i3+1)*x)*np.pi*(i3+1)
        return L * ro_0(x) / (L + ro_0(x) * zbroj)
    
    
    # Euler
    # interp_r_rj_f(x, t)
    # def ksi(x):
    #     return np.sqrt(3*L*x + a**3)
    
    # def iks(k):
    #     return (k**3-a**3)/(3*L)
    
    # from scipy.optimize import fsolve
    # def r_nult(k, t, r):
    #     return interp_r_rj_f(iks(k), t) - r
    # def iks_od_r(r, t):
    #     x0 = (a+b)/2
    #     return iks(fsolve(r_nult, x0, args =(t, r)))
    
    # def ro_kon(r, t):
    #     return interp_ro_rj_f(iks_od_r(r, t), t)
    
    ## v --------------------------------------------------------------------------
    def v_rj_f(x, t_ind):
        zbroj = 0.
        for it in range(0, n):
            zbroj = zbroj + v_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)
        return zbroj
    
    def interp_v_rj_f(x, t):
        zbroj = 0.
        for it in range(0, n):
            zbroj = zbroj + (rjesenje_pi.sol(t))[it] * np.sin(np.pi * (it+1) * x)
        return zbroj
        
    # def v_kon(r, t):
    #     return interp_v_rj_f(iks_od_r(r, t), t)
    
    ## omega ----------------------------------------------------------------------
    def omega_rj_f(x, t_ind):
        zbroj = 0.
        for it in range(0, n):
            zbroj = zbroj + omega_rj[it, t_ind] * np.sin(np.pi * (it+1) * x)
        return zbroj
    
    def interp_omega_rj_f(x, t):
        zbroj = 0.
        for it in range(0, n):
            zbroj = zbroj + (rjesenje_pi.sol(t))[n+it] * np.sin(np.pi * (it+1) * x)
        return zbroj
    
    # def omega_kon(r, t):
    #     return interp_omega_rj_f(iks_od_r(r, t), t)
    
    
    ## teta -----------------------------------------------------------------------
    def teta_rj_f(x, t_ind):
        zbroj = 0.
        for it in range(0, n+1):
            zbroj = zbroj + teta_rj[it, t_ind] * np.cos(np.pi * it * x)
        return zbroj
    
    def interp_teta_rj_f(x, t):
        zbroj = 0.
        for it in range(0, n+1):
            zbroj = zbroj + (rjesenje_pi.sol(t))[2*n+it] * np.cos(np.pi * it * x)
        return zbroj
    
    # def teta_kon(r, t):
    #     return interp_teta_rj_f(iks_od_r(r, t), t)
    
    ## z -----------------------------------------------------------------------
    def z_rj_f(x, t_ind):
        zbroj = 0.
        for it in range(0, n+1):
            zbroj = zbroj + z_rj[it, t_ind] * np.cos(np.pi * it * x)
        return zbroj
    
    def interp_z_rj_f(x, t):
        zbroj = 0.
        for it in range(0, n+1):
            zbroj = zbroj + (rjesenje_pi.sol(t))[3*n+1+it] * np.cos(np.pi * it * x)
        return zbroj

    
    
    
    # ro graf
    fig = plt.figure(1, figsize=(19, 15))
    fig.suptitle(r'$\rho(\cdot, t)$')
    rovi = np.linspace(0,1, rovi_size)
    for i in range(nbla):
        t = T*bla[i]
        plt.subplot(2, 2, i+1)
        plt.gca().title.set_text("t = " + str(t))
        Z = np.zeros(np.shape(rovi))
        for it in range(np.size(rovi)):
            Z[it] = interp_ro_rj_f(rovi[it], t)
        plt.plot(rovi, Z, label = 'p = '+str(p), color = boje[pi], linestyle = crte[pi], marker = markeri[pi])
        plt.gca().legend()
        plt.gca().grid(True)
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/ro_povi_pr' + str(primjer) + '.png')
    
    tplot = np.linspace(0, scale_r*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    for it in range(0, np.size(rovi[1])):
        for it2 in range(0, np.size(tplot[0])):
            Zr[it2, it, pi] = interp_ro_rj_f(rovi[it2, it], tplot[it2, it]) 
        
    # v graf
    fig = plt.figure(2, figsize=(19, 15))
    fig.suptitle(r'$v(\cdot, t)$')
    rovi = np.linspace(0,1, rovi_size)
    for i in range(nbla):
        t = T*bla[i]
        plt.subplot(2, 2, i+1)
        plt.gca().title.set_text("t = " + str(t))
        Z = np.zeros(np.shape(rovi))
        for it in range(np.size(rovi)):
            Z[it] = interp_v_rj_f(rovi[it], t)
        plt.plot(rovi, Z, label = 'p = '+str(p), color = boje[pi], linestyle = crte[pi], marker = markeri[pi])
        plt.gca().legend()
        plt.gca().grid(True)
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/v_povi_pr' + str(primjer) + '.png')
 
    tplot = np.linspace(0, scale_v*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    for it in range(0, np.size(rovi[1])):
        for it2 in range(0, np.size(tplot[0])):
            Zv[it2, it, pi] = interp_v_rj_f(rovi[it2, it], tplot[it2, it]) 

    # omega graf
    fig = plt.figure(3, figsize=(19, 15))
    fig.suptitle(r'$\omega(\cdot, t)$')
    rovi = np.linspace(0,1, rovi_size)
    for i in range(nbla):
        t = T*bla[i]
        plt.subplot(2, 2, i+1)
        plt.gca().title.set_text("t = " + str(t))
        Z = np.zeros(np.shape(rovi))
        for it in range(np.size(rovi)):
            Z[it] = interp_omega_rj_f(rovi[it], t)
        plt.plot(rovi, Z, label = 'p = '+str(p), color = boje[pi], linestyle = crte[pi], marker = markeri[pi])
        plt.gca().legend()
        plt.gca().grid(True)
    if(pi == np.size(povi)-1):
        plt.show()        
        fig.savefig('slike/omega_povi_pr' + str(primjer) + '.png')
        
    tplot = np.linspace(0, scale_o*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    for it in range(0, np.size(rovi[1])):
        for it2 in range(0, np.size(tplot[0])):
            Zo[it2, it, pi] = interp_omega_rj_f(rovi[it2, it], tplot[it2, it]) 
        
    # teta graf
    fig = plt.figure(4, figsize=(19, 15))
    fig.suptitle(r'$\theta(\cdot, t)$')
    rovi = np.linspace(0,1, rovi_size)
    for i in range(nbla):
        t = T*bla[i]
        plt.subplot(2, 2, i+1)
        plt.gca().title.set_text("t = " + str(t))
        Z = np.zeros(np.shape(rovi))
        for it in range(np.size(rovi)):
            Z[it] = interp_teta_rj_f(rovi[it], t)
        plt.plot(rovi, Z, label = 'p = '+str(p), color = boje[pi], linestyle = crte[pi], marker = markeri[pi])
        plt.gca().legend()
        plt.gca().grid(True)
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/teta_povi_pr' + str(primjer) + '.png')
        
    tplot = np.linspace(0, scale_t*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    for it in range(0, np.size(rovi[1])):
        for it2 in range(0, np.size(tplot[0])):
            Zt[it2, it, pi] = interp_teta_rj_f(rovi[it2, it], tplot[it2, it]) 
    print(p, 'time =', time.time() - start)
        
 
plt.rcParams.update({'font.size': 12})
for pi in range(np.size(povi)):
    p = povi[pi] 
    # ro
    rovi = np.linspace(0,1, rovi_size)
    tplot = np.linspace(0, scale_r*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    levels = np.linspace(np.amin(Zr[:, :, pi]), np.amax(Zr[:, :, pi]), 500)
    fig = plt.figure(5, figsize=(19, 15))
    fig.suptitle(r'$\rho$')
    plt.subplot(2, 2, pi+1)
    plt.gca().title.set_text("p = " + str(p))
    plt.contourf(tplot, rovi, Zr[:, :, pi], levels,  cmap='viridis')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.colorbar()
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/ro_povi_kont_pr' + str(primjer) + '.png')
        
    # v
    rovi = np.linspace(0,1, rovi_size)
    tplot = np.linspace(0, scale_v*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    levels = np.linspace(np.amin(Zv[:, :, pi]), np.amax(Zv[:, :, pi]), 500)
    fig = plt.figure(6, figsize=(19, 15))
    fig.suptitle(r'$v$')
    plt.subplot(2, 2, pi+1)
    plt.gca().title.set_text("p = " + str(p))
    plt.contourf(tplot, rovi, Zv[:, :, pi], levels,  cmap='viridis')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.colorbar()
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/v_povi_kont_pr' + str(primjer) + '.png')
        
    # omega
    rovi = np.linspace(0,1, rovi_size)
    tplot = np.linspace(0, scale_o*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    levels = np.linspace(np.amin(Zo[:, :, pi]), np.amax(Zo[:, :, pi]), 500)
    fig = plt.figure(7, figsize=(19, 15))
    fig.suptitle(r'$\omega$')
    plt.subplot(2, 2, pi+1)
    plt.gca().title.set_text("p = " + str(p))
    plt.contourf(tplot, rovi, Zo[:, :, pi], levels,  cmap='viridis')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.colorbar()
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/omega_povi_kont_pr' + str(primjer) + '.png')
        
        
    # teta
    rovi = np.linspace(0,1, rovi_size)
    tplot = np.linspace(0, scale_t*T, tplot_size)
    rovi, tplot = np.meshgrid(rovi, tplot)
    levels = np.linspace(np.amin(Zt[:, :, pi]), np.amax(Zt[:, :, pi]), 500)
    fig = plt.figure(8, figsize=(19, 15))
    fig.suptitle(r'$\theta$')
    plt.subplot(2, 2, pi+1)
    plt.gca().title.set_text("p = " + str(p))
    plt.contourf(tplot, rovi, Zt[:, :, pi], levels,  cmap='viridis')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.colorbar()
    if(pi == np.size(povi)-1):
        plt.show()
        fig.savefig('slike/sliketeta_povi_kont_pr' + str(primjer) + '.png')
        
    
