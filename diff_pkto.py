import numpy as np
from classes_pkto import Vector

G = 398600.44
Rst0 = np.array([[0], [0], [0]])
Aw0 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

def dnadtn(r, v, Rst = Rst0, Aw = Aw0):
    Rm = ((Rst + r)**2)**0.5
    RVm = (Rst + r)*(Aw*(Rst + r) + v)
    Vm = ((Aw*(Rst + r) + v)**2)**0.5

    res = [r, v]
    a = -G*(Rst + r)*((Rst + r)**2)**(-1.5)
    res.append(a)
    a1 = G*(3*RVm*(Rst + r) - Rm**2*(Aw*(Rst + r) + v))/Rm**5
    res.append(a1)
    a2 =  -G*(2*G*Rm*(Rst + r) + 15*RVm**2*(Rst + r) - 6*RVm*Rm**2*(Aw*(Rst + r) + v) - 3*Rm**2*Vm**2*(Rst + r))/Rm**7
    res.append(a2)
    a3 = G*(30*G*RVm*Rm*(Rst + r) - 8*G*Rm**3*(Aw*(Rst + r) + v) + 105*RVm**3*(Rst + r) - 45*RVm**2*Rm**2*(Aw*(Rst + r) + v) - 45*RVm*Rm**2*Vm**2*(Rst + r) + 9*Rm**4*Vm**2*(Aw*(Rst + r) + v))/Rm**9
    res.append(a3)
    return res

def drdt(r, v, Rst = Rst0, Aw = Aw0):
    Rd = dnadtn(r, v, Rst, Aw)
    rd2 = -Aw**2*(Rst + r) - 2*Aw*v + Rd[2]
    rd3 = -Aw*(Aw**2*(Rst + r) + 3*Aw*v + 3*rd2) + Rd[3]
    rd4 = -Aw*(Aw**3*(Rst + r) + 4*Aw**2*v + 6*Aw*rd2 + 4*rd3) + Rd[4]
    rd5 = -Aw*(Aw**4*(Rst + r) + 5*Aw**3*v + 10*Aw**2*rd2 + 10*Aw*rd3 + 5*rd4) + Rd[5]
    res = [r, v, rd2, rd3, rd4, rd5]
    return res

def dadx(r, r_x, v, v_x, Rst = Rst0, Aw = Aw0):
    Rm = ((Rst + r)**2)**0.5
    RVm = (Rst + r)*(Aw*(Rst + r) + v)
    Vm = ((Aw*(Rst + r) + v)**2)**0.5
    RRx = (Rst + r)*r_x
    RVx = (Rst + r)*(Aw*r_x + v_x)
    VRx = (Aw*(Rst + r) + v)*r_x

    res = [r, v]
    a = -G*(r_x*Rm**2 - 3*(Rst + r)*(r_x*(Rst + r)))*Rm**-5
    res.append(a)
    a1 = G*(-3*RVm*(5*(Rst + r)*RRx - Rm**2*r_x) - Rm**4*(Aw*r_x + v_x) + 3*Rm**2*((Rst + r)*(RVx + VRx) + RRx*(Aw*(Rst + r) + v)))/Rm**7
    res.append(a1)
    a2 =  0
    res.append(a2)
    a3 = 0
    res.append(a3)
    return res

def drdx(r, r_x, v, v_x, Rst = Rst0, Aw = Aw0):
    Rd = dadx(r, r_x, v, v_x, Rst, Aw)
    rd2 = -Aw**2*r_x - 2*Aw*v_x + Rd[2]
    rd3 = -Aw*(Aw**2*r_x + 3*Aw*v_x + 3*rd2) + Rd[3]
    rd4 = Vector([0, 0, 0]) # -Aw*(Aw**3*r_x + 4*Aw**2*v_x + 6*Aw*rd2 + 4*rd3) + Rd[4]
    rd5 = Vector([0, 0, 0]) # -Aw*(Aw**4*r_x + 5*Aw**3*v_x + 10*Aw**2*rd2 + 10*Aw*rd3 + 5*rd4) + Rd[5]
    res = [r_x, v_x, rd2, rd3, rd4, rd5]
    return res

def ak(r, Rst = Rst0):
    res = -6.0*G*((Rst + r)**2)**(-2.5) * (Rst + r)
    return res

def bk(r, Rst = Rst0, Aw = Aw0):
    bk1 = Vector([1, 1, 1])
    bk2 = -2*Aw * Vector([1, 1, 1])
    bk3 = 3*Aw**2 * Vector([1, 1, 1]) + 2.0*G*((Rst + r)**2)**(-1.5) * Vector([1, 1, 1])
    bk4 = -Aw*(12.0*G*r + 4*Aw**2*((Rst + r)**2)**1.5 * (Rst + r) + 8.0*G * (Rst + r))*((Rst + r)**2)**(-2.5) * (Rst + r) * Vector([1, 1, 1])
    res = [bk1, bk2, bk3, bk4]
    return res

def ck(r, Rst = Rst0, Aw = Aw0):
    ck0 = r
    ck1 = Vector([0, 0, 0])
    ck2 = -(Aw**2*((Rst + r)**2)**1.5*(Rst + r) + G*(Rst + r))*((Rst + r)**2)**(-1.5)
    ck3 = 2*Aw**3*(Rst + r) + 3*Aw*G*Rst*((Rst + r)**2)**(-1.5) + 5.0*Aw*G*r*((Rst + r)**2)**(-1.5)
    ck4 = -1/((Rst + r)**2)**3 * (3*Aw**4*((Rst + r)**2)**3* (Rst + r) + 2*Aw**2*G*((Rst+r)**2)**0.5*(3*Rst**2 + 10*Rst*r + 10*r**2)* (Rst + r) + 2*G**2* (Rst + r)) # * (Rst + r)
    res = [ck0, ck1, ck2, ck3, ck4]
    return res

if __name__=='__main__':

    from dnxdtn import *

    bt = baltar(5)
    print('a:', 1/2 * bt[4].diff(v).diff(v))
    for i in range(4):
        print('b%i:' % i, bt[i].diff(v))
    print('b4:', (bt[4].diff(v) - bt[4].diff(v).diff(v)*v).simplify())
    for i in range(4):
        print('c%i:' % i, (bt[i] - bt[i].diff(v)*v).simplify())
    sl1 = 1/2 * (bt[4].diff(v).diff(v)*v*v).simplify()
    sl2 = (bt[4].diff(v)*v - bt[4].diff(v).diff(v)*v*v).simplify()
    print('c4:', (bt[4] - sl1 - sl2).simplify())
