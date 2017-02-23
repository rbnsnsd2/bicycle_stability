#!/usr/bin/env python
# Created by Saul Robinson
"""
    The original research paper that this python code replicates can be found here:
    http://rspa.royalsocietypublishing.org/content/463/2084/1955
    Linearized dynamics equations for the balance and steer of a bicycle: a benchmark and review
    J.P Meijaard, Jim M Papadopoulos, Andy Ruina, A.L Schwab
    Published 8 August 2007.DOI: 10.1098/rspa.2007.1857
"""

from pylab import plot, axhline, plt, show
from numpy import array, sin, cos, arange, zeros, complex128, poly1d, sort

def main():

    ###################################
    # Parameters for the bicycle model#
    ###################################
    w = 1.02     # Wheelbase (m) = 1.02
    c = .08      # trail (m) = .08
    pi = 3.141592653589793238462643383279
    lam = pi/10. # Steering tilt (90-HT angle) =pi/10
    g = 9.81     # gravity
    # Rear wheel parameters
    rr = .3      # radius (m)
    ra = .03     # radius rim/tyre (m)
    mr = 2.       # mass (kg)
    #Irxx = .0603 # moment of inertia (kg m^2)
    Irxx = (.625*ra**2 + .5*rr**2)*mr*.7
    #Iryy = .12   # moment of inertia (kg m^2)
    Iryy = (.75*ra**2 + rr**2)*mr*.7
    # Body (rear body and frame) parameters
    xb = .3      # center of mass (m) =.3
    zb = -.9     # center of mass (m) =-.9
    mb = 85.      # mass (kg) = 85
    Ibxx = 19.2   # moment of inertia =9.2
    Ibxz = 2.4   # moment of inertia =2.4
    Ibyy = 11.    # moment of inertia =11
    Ibzz = 2.8   # moment of inertia
    IB = array([[Ibxx, 0, Ibxz], [0, Ibyy, 0], [Ibxz, 0, Ibzz]]) # body inertia matrix

    # Front (handlebar and fork) parameters
    xh = .9      # center of mass (m) =.9
    zh = -.7     # center of mass (m) =-.7
    mh = 4       # mass (kg)
    Ihxx = .05892   # moment of inertia =.05892
    Ihxz = -.00756  # moment of inertia
    Ihyy = .06      # moment of inertia =.06
    Ihzz = .00708   # moment of inertia =.00708
    IH = array([[Ihxx, 0, Ihxz], [0, Ihyy, 0], [Ihxz, 0, Ihzz]])

    # Front wheel parameters
    rf = .35     # radius (m)
    fa = .03     # radius rim/tyre (m)
    mf = 3       # mass (kg)
    Ifxx = .1405 # moment of inertia =.1405
    #IfxxCG = (.625*fa**2 + .5*fr**2)*mf*.7
    Ifyy = .28   # moment of inertia =.28
    #IfyyCG = (.75*fa**2 + fr**2)*mf*.7
    
    ########################################
    # End Parameters for the bicycle model #
    ########################################


    #####################################
    # Definition of combined parameters #
    #####################################
    mt = mr + mb + mh + mf                  # A1
    xt = (xb*mb + xh*mh + w*mf)/mt          # A2
    zt = (-rr*mr + zb*mb + zh*mh -rf*mf)/mt # A3

    Itxx = Irxx + Ibxx + Ihxx + Ifxx + mr*rr**2 + mb*zb**2 + mh*zh**2 + mf*rf**2 # A4
    Itxz = Ibxz + Ihxz - mb*xb*zb - mh*xh*zh + mf*w*rf                           # A5

    Irzz = Irxx # A6
    Ifzz = Ifxx # 

    Itzz = Irzz + Ibzz + Ihzz + Ifzz + mb*xb**2 + mh*xh**2 + mf*w**2 # A7
    ma = mh + mf                                                     # A8
    xa = (xh*mh + w*mf)/ma                                           # A9
    za = (zh*mh - rf*mf)/ma                                          # A9

    Iaxx = Ihxx + Ifxx + mh*(zh - za)**2 + mf*(rf + za)**2       # A10
    Iaxz = Ihxz - mh*(xh - xa)*(zh - za) + mf*(w - xa)*(rf + za) # A11
    Iazz = Ihzz + Ifzz + mh*(xh - xa)**2 + mf*(w - xa)**2        # A12

    LAM = array([[sin(lam)], [0], [cos(lam)]])

    ua = (xa - w - c)*cos(lam) - za*sin(lam) # A13

    Iall = ma*ua**2 + Iaxx*sin(lam)**2 + 2*Iaxz*sin(lam)*cos(lam) + Iazz*cos(lam)**2 # A14
    Ialx = -ma*ua*za + Iaxx*sin(lam) + Iaxz*cos(lam)                                 # A15
    Ialz = ma*ua*xa + Iaxz*sin(lam) + Iazz*cos(lam)                                  # A16

    mu = (c/w)*cos(lam)    # A17

    sr = Iryy/rr           # A18
    sf = Ifyy/rf           #
    st = sr + sf           #

    sa = ma*ua + mu*mt*xt  # A19

    Mpp = Itxx                          # A20
    Mpd = Ialx + mu*Itxz                #
    Mdp = Mpd                           #
    Mdd = Iall + 2*mu*Ialz + mu**2*Itzz #

    MM = array([[Mpp, Mpd], [Mdp, Mdd]]) # A21

    K0pp = mt*zt           # A22
    K0pd = - sa            #
    K0dp = K0pd            #
    K0dd = - sa*sin(lam)   #

    K0 = array([[K0pp, K0pd], [K0dp, K0dd]]) # A23
    KK0 = g*K0
    
    K2pp = 0                               # A24
    K2pd = ((st - mt*zt)/w)*cos(lam)       #
    K2dp = 0                               #
    K2dd = ((sa + sf*sin(lam))/w)*cos(lam) #
    
    K2 = array([[K2pp, K2pd], [K2dp, K2dd]]) # A25

    C1pp = 0                                                  # A 26
    C1pd = mu*st + sf*cos(lam) + (Itxz/w)*cos(lam) - mu*mt*zt #
    C1dp = -(mu*st + sf*cos(lam))                             #
    C1dd = (Ialz/w)*cos(lam) + mu*(sa + (Itzz/w)*cos(lam))    #
    
    C1 = array([[C1pp, C1pd], [C1dp, C1dd]]) # A27
    #print('M', MM, 'K0', K0, 'K2', K2, 'C1', C1)
    
    #########################################
    # End Definition of combined parameters #
    #########################################

    ###################
    # Main functional #
    ###################
    V = arange(0.,10.,.0125)
    #print(V)

    # Redefine the C1 matrix to a 3D matrix
    # thus losing the variable v
    CC1 = zeros((len(V), 2, 2))
    KK2 = zeros((len(V), 2, 2))
    ZZ = zeros((len(V),4)).astype(complex128)

    for H in range(len(V)):
        CC1[H,:,:] = V[H]*C1
        KK2[H,:,:] = V[H]**2*K2
        Lam4 = -MM[0,1]*MM[1,0] + MM[0,0]*MM[1,1]
        Lam3 = CC1[H,1,1]*MM[0,0] - CC1[H,1,0]*MM[0,1] - CC1[H,0,1]*MM[1,0] + CC1[H,0,0]*MM[1,1]
        Lam2 = - CC1[H,0,1]*CC1[H,1,0] + CC1[H,0,0]*CC1[H,1,1] + KK0[1,1]*MM[0,0] + KK2[H,1,1]*MM[0,0] - KK0[1,0]*MM[0,1] - KK2[H,1,0]*MM[0,1] - KK0[0,1]*MM[1,0] - KK2[H,0,1]*MM[1,0] + KK0[0,0]*MM[1,1] + KK2[H,0,0]*MM[1,1]
        Lam1 = CC1[H,1,1]*KK0[0,0] - CC1[H,1,0]*KK0[0,1] - CC1[H,0,1]*KK0[1,0] + CC1[H,0,0]*KK0[1,1] + CC1[H,1,1]*KK2[H,0,0] - CC1[H,1,0]*KK2[H,0,1] - CC1[H,0,1]*KK2[H,1,0] + CC1[H,0,0]*KK2[H,1,1]
        Lam0 = -KK0[0,1]*KK0[1,0] + KK0[0,0]*KK0[1,1] + KK0[1,1]*KK2[H,0,0] - KK0[1,0]*KK2[H,0,1] - KK0[0,1]*KK2[H,1,0] - KK2[H,0,1]*KK2[H,1,0] + KK0[0,0]*KK2[H,1,1] + KK2[H,0,0]*KK2[H,1,1]
        zz = poly1d([Lam4, Lam3, Lam2, Lam1, Lam0])
        ZZ[H,:] = zz.r # find the roots and put them in an array
        
    ss = ZZ.shape

    #######################
    # End Main functional #
    #######################

    pp = sort(ZZ, axis=1) # axis 1 is across rows
    
    plot(V, pp[:,0], V, pp[:,1], V, pp[:,2], V, pp[:,3])
    #plot( V, pp[:,2], V, pp[:,3])
    axhline(y=0)
    plt.axis([0, 10, -10, 10])
    show()
        
if __name__ == "__main__":
    main()
