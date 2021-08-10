#!/software/anaconda3/bin

""" 

Moist Energy Balance Model with a seasonal cycle of insolation and a
thermodynamic sea ice model, as described in the paper:

Feldl, N., and T. M. Merlis (2021), Polar amplification in idealized 
climates: the role of ice, moisture, and seasons

This code is based off the Dry EBM presented in the paper:

Wagner, T.J. and I. Eisenman, 2015: How Climate Model Complexity 
Influences Sea Ice Stability. J. Climate, 28, 3998–4014, 
https://doi.org/10.1175/JCLI-D-14-00654.1

with the following modifications:
- We have added the effect of latent energy on the diffusive 
  representation of atmospheric energy transport
- We have added functionality for disabling the effect of sea-ice 
  thermodynamics, the seasonal cycle of insolation, ice-albedo
  feedback, and latent energy transport.
- We use a global, rather than single-hemisphere, domain.

"""

import sys
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()
import mixedlayer
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.rcParams['axes.titlesize'] = 10 # reset global fig properties
mpl.rcParams['legend.fontsize'] = 6 # reset global fig properties
mpl.rcParams['legend.title_fontsize'] = 6 # reset global fig properties
mpl.rcParams['xtick.labelsize'] = 8 # reset global fig properties
mpl.rcParams['ytick.labelsize'] = 8 # reset global fig properties
mpl.rcParams['ytick.right'] = True # reset global fig properties
mpl.rcParams['axes.labelsize'] = 8 # reset global fig properties

def saturation_specific_humidity(temp,press):

  """
  We assume a single liquid-to-vapor phase transition with the parameter values 
  of the Clausius Clapeyron (CC) relation given in OGorman and Schneider (2008) 
  to determine the saturation specific humidity qs(T).

  """

  es0 = 610.78 # saturation vapor pressure at t0 (Pa)
  t0 = 273.16
  Rv = 461.5
  Lv = 2.5E6
  ep = 0.622 # ratio of gas constants of dry air and water vapor
  temp = temp + 273.15 # convert to Kelvin
  es = es0 * np.exp(-(Lv/Rv) * ((1/temp)-(1/t0)))
  qs = ep * es / press

  return qs


def model(grid, T, F=0, moist = 0, albT = 0, seas = 0, thermo = 0):

  if moist==0:
    D = 0.6 # diffusivity for heat transport (W m^-2 K^-1)
  elif moist==1:
    D = 0.3 # diffusivity for heat transport (W m^-2 K^-1)
  print(f'diffusivity for heat transport is {D} W m^-2 K^-1')
  S1 = 338 # insolation seasonal dependence (W m^-2)
  A = 196 # OLR when T = 0 (W m^-2)
  B = 1.8 # OLR temperature dependence (W m^-2 K^-1)
  #cw = 9.8 # ocean mixed layer heat capacity (W yr m^-2 K^-1)
  cw = mixedlayer.heatcapacity(60) # ocean mixed layer heat capacity (W yr m^-2 K^-1)
  S0 = 420 # insolation at equator (W m^-2)
  S2 = 240 # insolation spatial dependence (W m^-2)
  a0 = 0.7 # ice-free co-albedo at equator
  a2 = 0.1 # ice=free co-albedo spatial dependence
  ai = 0.4 # co-albedo where there is sea ice
  Fb = 0 # heat flux from ocean below (W m^-2)
  k = 2 # sea ice thermal conductivity (W m^-1 K^-1)
  Lf = 9.5 # sea ice latent heat of fusion (W yr m^-3)
  cg = 0.098 #0.01*cw # ghost layer heat capacity(W yr m^-2 K^-1)
  tau = 1e-5 # ghost layer coupling timescale (yr)
  Lv = 2.5E6 # latent heat of vaporization (J kg^-1)
  cp = 1004.6 # heat capacity of air at constant pressure (J kg^-1 K^-1)
  RH = 0.8 # relative humidity
  Ps = 1E5 # surface pressure (Pa)

  # Read grid dict
  n = grid['n']; dx = grid['dx']; x = grid['x']; xb = grid['xb']
  nt = grid['nt']; dur = grid['dur']; dt = grid['dt']
  print(f'Running model for {n} grid cells and {dur} years ...')

  # Diffusion Operator (WE15, Appendix A) 
  lam = D/dx**2*(1-xb**2)
  L1=np.append(0, -lam) 
  L2=np.append(-lam, 0) 
  L3=-L1-L2
  diffop = - np.diag(L3) - np.diag(L2[:n-1],1) - np.diag(L1[1:n],-1)

  # Definitions for implicit scheme on Tg
  cg_tau = cg/tau
  dt_tau = dt/tau
  dc = dt_tau*cg_tau
  kappa = (1+dt_tau)*np.identity(n)-dt*diffop/cg

  # Seasonal forcing (WE15 eq.3)
  if seas == 0:
    S1 = 0.0

  ty = np.arange(dt/2,1+dt/2,dt)
  S = (np.tile(S0-S2*x**2,[nt,1])-
       np.tile(S1*np.cos(2*np.pi*ty),[n,1]).T*np.tile(x,[nt,1]))

  # zero out negative insolation
  S = np.where(S<0,0,S)

  # Further definitions
  M = B+cg_tau
  aw = a0-a2*x**2 # open water albedo
  kLf = k*Lf

  # Set up output arrays
  Efin = np.zeros((n,nt)) 
  Tfin = np.zeros((n,nt))
  T0fin = np.zeros((n,nt))
  ASRfin = np.zeros((n,nt))
  tfin = np.linspace(0,1,nt)

  # Initial conditions
  Tg = T
  E = cw*T

  # Integration (see WE15_NumericIntegration.pdf)
  # Loop over Years 
  for years in range(dur):
    # Loop within One Year
    for i in range(nt):

      # forcing
      if albT == 1:
        alpha = aw*(E>0) + ai*(E<0) #WE15, eq.4
      else:
        alpha = aw

      C = alpha*S[i,:] + cg_tau*Tg - A + F

      # surface temperature
      if thermo == 1:
        T0 = C/(M-kLf/E) #WE15, eq.A3
      else:
        T0 = E/cw

      # store final year 
      if years==(dur-1): 
        Efin[:,i] = E
        Tfin[:,i] = T
        T0fin[:,i] = T0
        ASRfin[:,i] = alpha*S[i,:]

      T = E/cw*(E>=0)+T0*(E<0)*(T0<0) #WE15, eq.9
      # Forward Euler on E
      E = E+dt*(C-M*T+Fb) #WE15, eq.A2

      # Implicit Euler on Tg
      if moist == 1:

        # Forward Euler on diffusion of latent heat
        q = RH * saturation_specific_humidity(Tg,Ps)
        rhs1 = np.matmul(dt*diffop/cg, Lv*q/cp)

        if thermo == 1:
        # FM21, eq. 3
          Tg = np.linalg.solve(kappa-np.diag(dc/(M-kLf/E)*(T0<0)*(E<0)),
                               Tg + rhs1 + (dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A+F)/(M-kLf/E)*(T0<0)*(E<0))))
        else:
          Tg = np.linalg.solve(kappa,
                               Tg + rhs1 + dt_tau*(E/cw) )

      elif moist == 0:
        if thermo ==1:
        #WE15, eq. A1
          Tg = np.linalg.solve(kappa-np.diag(dc/(M-kLf/E)*(T0<0)*(E<0)),
                               Tg + (dt_tau*(E/cw*(E>=0)+(ai*S[i,:]-A+F)/(M-kLf/E)*(T0<0)*(E<0))))
        else:
          Tg = np.linalg.solve(kappa,
                               Tg + dt_tau*(E/cw) )

  print(f'{np.mean(Tfin, axis=(0,1))} global mean temp')
  print(f'{np.ptp(np.mean(Tfin, axis=1))} equator-pole temp difference')
  print(f'{np.mean(S, axis=(0,1))} global mean annual mean inso')
  print(f'{np.mean(ASRfin, axis=(0,1)) - A - B*np.mean(Tfin, axis=(0,1))} energy balance')

  return tfin, Efin, Tfin, T0fin, ASRfin


def main():

  # Set up grid and time-stepping
  n = 120
  dx = 2./n #grid box width
  x = np.linspace(-1+dx/2,1-dx/2,n) #native grid
  xb = np.linspace(-1+dx,1-dx,n-1) 
  nt = 1000
  dur= 100
  dt = 1./nt
  grid = {'n': n, 'dx': dx, 'x': x, 'xb': xb, 'nt': nt, 'dur': dur, 'dt': dt} 

  # Integrate EBM
  Ti = 7.5+20*(1-2*x**2) # initial condition 
  F = 0
  tfin, Efin, Tfin, T0fin, ASRfin = model(grid, Ti, F, moist = 1, albT = 1, seas = 1, thermo = 1)

  print("--- %s seconds ---" % (time.time() - start_time))

  # WE15, Figure 2: Default Steady State Climatology 
  fig = plt.figure(figsize=(10,5))
  plt.suptitle('Annual mean global mean: {:.2f}deg C, Equator-pole difference: {:.2f}deg C, Polar seasonal amplitude: {:.2f}deg C'.format(np.mean(Tfin, axis=(0,1)), np.ptp(np.mean(Tfin, axis=1)), np.ptp(Tfin[-1,:])))

  # northern hemi only for plot
  n_2 = int(n/2)
  x_n = x[-n_2:]
  Tfin = Tfin[-n_2:,:]
  Efin = Efin[-n_2:,:]
  T0fin = T0fin[-n_2:,:]

  # seasons and ice edge
  # winter/summer occur when hemispheric T is min/max
  winter = np.argmin(np.mean(Tfin, axis=0))
  summer = np.argmax(np.mean(Tfin, axis=0))
  ice = np.where(Efin<0,np.expand_dims(x_n,1),1)
  xi = np.min(ice, axis=0)
  Lf = 9.5 # sea ice latent heat of fusion (W yr m^-3)
  icethick = -Efin/Lf*(Efin<0)

  # plot enthalpy (Fig 2a)
  plt.subplot(141)
  clevsE = np.arange(-300,301,50)
  plt.contourf(tfin,x_n,Efin,clevsE)
  plt.colorbar()
  # plot ice edge on E
  plt.contour(tfin,x_n,icethick,[0],colors='k')
  plt.xlabel('t (final year)')
  plt.ylabel('x')
  plt.ylim(0,1)
  plt.title(r'E (J m$^{-2}$)')

  # plot temperature (Fig 2b)
  plt.subplot(142)
  clevsT = np.arange(-30,31,5)
  plt.contourf(tfin,x_n,Tfin,clevsT)
  plt.colorbar()
  # plot T=0 contour (the region between ice edge and T=0 contour is the
  # region of summer ice surface melt)
  plt.contour(tfin,x_n,icethick,[0],colors='k')
  plt.contour(tfin,x_n,T0fin,[0],colors='r')
  plt.xlabel('t (final year)')
  plt.ylabel('x')
  plt.ylim(0,1)
  plt.title(r'T ($^\circ$C)')
   
  # plot the ice thickness (Fig 2c)
  plt.subplot(1,4,3)
  clevsh = np.arange(0.00001,5.5,0.5)
  plt.contourf(tfin,x_n,icethick,clevsh)
  plt.colorbar()
  # plot ice edge on h
  plt.contour(tfin,x_n,icethick,[0],colors='k')
  plt.plot([tfin[winter], tfin[winter]],[0,max(x_n)],'k')
  plt.plot([tfin[summer], tfin[summer]],[0,max(x_n)],'k--')
  plt.xlabel('t (final year)')
  plt.ylabel('x')
  plt.ylim(0,1)
  plt.title('h (m)')
   
  # plot temperature profiles (Fig 2d)
  plt.subplot(444)
  Summer, = plt.plot(x_n,Tfin[:,summer],'k--',label='summer')
  Winter, = plt.plot(x_n,Tfin[:,winter],'k',label='winter')
  plt.plot([0,1],[0,0],'k')
  plt.xlabel('x')
  plt.ylabel(r'T ($^\circ$C)')
  plt.legend(handles = [Summer,Winter],loc=0, fontsize=8)
   
  # plot ice thickness profiles (Fig 2e)
  plt.subplot(448)
  plt.plot(x_n,icethick[:,summer],'k--')
  plt.plot(x_n,icethick[:,winter],'k')
  plt.plot([0,1], [0,0],'k')
  plt.xlim([0.7,1])
  plt.xlabel('x')
  plt.ylabel('h (m)')

  # plot seasonal thickness cycle at pole (Fig 2f)
  plt.subplot(4,4,12)
  plt.plot(tfin,icethick[-1,:],'k')
  plt.xlabel('t (final year)')
  plt.ylabel(r'h$_p$ (m)')
   
  # plot ice edge seasonal cycle (Fig 2g)
  plt.subplot(4,4,16)
  xideg = np.degrees(np.arcsin(xi))
  plt.plot(tfin,xideg,'k-')
  plt.ylim([40,90])
  plt.xlabel('t (final year)')
  plt.ylabel(r'$\theta_i$ (deg)')

  plt.tight_layout() 
  plt.show()


if __name__ == '__main__':
  main()

