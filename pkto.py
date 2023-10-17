import pandas as pd
from math import pi, sin, cos, asin, acos
from gauss import Gauss
from classes_pkto import *
from diff_pkto import *

import os
file_spFiV = 'pkto_spFiV.txt'
file_spFiR = 'pkto_spFiR.txt'
file_procMinimizeFunc = 'pkto_procMinimizeFunc.txt'
file_Bl13 = 'pkto_Bl13.txt'
file_V_borders = 'pkto_V_borders.txt'
list_files = [file_spFiV, file_spFiR, file_procMinimizeFunc, file_Bl13, file_V_borders]
for i in range(len(list_files)):
    if list_files[i] in os.listdir():
        os.remove(list_files[i])

# Глобальные перменные
# Переменные, не зависящие от станции
Rlev, Rprav = 300, 7000
minFi_dR = [1000, 500, 100]
minFi_dV = 1
beta0 = pi/4
dbeta = 1
MAXPOINTS = 100

# Константы вращения Земли и положения станции
w = Vector()
wst = Vector()
R0 = Vector()
R0st = Vector()
Aw = Matrix()

# Число опорных точек и среднее время
Num, Tc = None, None

# Общие переменные
u, v, V_u, V_v, uu = None, None, None, None, None
dbeta = None

uv = case_par()

def write_points(file):
    '''
    Функция записи опорных точек и точки сглаживания в глобальные переменные.\n
    file - переменная на названием файла
    '''
    global uv
    global u, v, V_u, V_v, uu
    global list_files, file_spFiV, file_spFiR, file_procMinimizeFunc, file_Bl13, file_V_borders

    file_spFiR = file + '/' + file_spFiR
    file_spFiV = file + '/' + file_spFiV
    file_procMinimizeFunc = file + '/' + file_procMinimizeFunc
    file_Bl13 = file + '/' + file_Bl13
    file_V_borders = file + '/' + file_V_borders
    for i in range(len(list_files)):
        # list_files[i] = file + '/' + list_files[i]
        if list_files[i] in os.listdir(file):
            os.remove(file + '/' + list_files[i])

    points = file + '/slava_points.txt'
    data = pd.read_csv(points, sep=';')#file, usecols=['t', 'x[0]', 'x[1]', 'bx[0]', 'bx[1]'])

    for i in range(Num):
        # idx = idxs[i]
        # ot = kpar['_mot'][i] [idx]
        uv.t[i] = data['t'][i] # ot._t
        uv.u[i] = data['u'][i]# cos(data['x[0]'][i]) # cos(ot._e)
        uv.v[i] = data['v'][i]# sin(data['x[1]'][i]) # sin(ot._q)

        uv.pu[i] = data['pu'][i]# data['bx[0]'][i] / sin(data['x[0]'][i])**2 # ot._be / sin(ot._e) / sin(ot._e)
        uv.pv[i] = data['pv'][i]# data['bx[1]'][i] / sin(data['x[1]'][i])**2 # ot._bq / sin(ot._q) / sin(ot._q)

    secpoints = file + '/slava_secpoints.txt'
    data = pd.read_csv(secpoints, sep=';')

    u = data['u'][0]# cos(data['x[0]'][Num // 2]) # cos(kpar['Eps'])
    v = data['v'][0]# sin(data['x[1]'][Num // 2]) # sin(kpar['Tet'])
    V_u = data['du'][0]# (data['x[0]'][Num // 2] - data['x[0]'][Num // 2 - 1]) / (data['t'][Num // 2] - data['t'][Num // 2 - 1])
    # V_u *= sin(data['x[0]'][Num // 2]) # - kpar['V_Eps'] * sin(kpar['Eps'])
    V_v = data['dv'][0]# (data['x[1]'][Num // 2] - data['x[1]'][Num // 2 - 1]) / (data['t'][Num // 2] - data['t'][Num // 2 - 1])
    # V_v *= cos(data['x[1]'][Num // 2]) # kpar['V_Tet'] * cos(kpar['Tet'])

    uu = data['uu'][0]# cos(data['x[1]'][Num // 2])**2 - u**2 # cos(kpar['Tet']) * cos(kpar['Tet']) - u * u
    
    if uu < 0:
        raise PKTO_Error('uu')
    # else:
    #     uu = uu**0.5

# Блок 1
def start_point(file):
    '''
    Функция подготовки программы. Записывает значения переменных из файла и выбирает направление движения по расстоянию.\n
    file - переменная на названием файла
    '''

    write_points(file)

    global dbeta

    if -1 * (V_u * u/uu + V_v * v/uu) * sin(beta0) + V_v * cos(beta0) >= 0:
        dbeta = -1
    else:
        dbeta = 1

    return 0

def Bl13_1(R_min, fi_min, Rpred):

    golden = 1 - (-1 + 5**0.5) / 2

    if fi_min[0] > fi_min[1]:
        r = golden * (R_min[0] - Rpred)
    else:
        r = golden * (R_min[1] - Rpred)

    with open(file_Bl13, 'a') as file:
        file.write('R1: %s\nR2: %s\nR3: 0\nRpred: %s\nY0: %s\nRr: %s\n' % (R_min[0], R_min[1], Rpred, r, Rpred + r))

    return Rpred + r

def Bl13_2(R_min, fi_min, Rpred, fi_pred):

    golden = 1 - (-1 + 5**0.5) / 2

    Y0, Y1, Z1, Z2 = None, None, None, None

    if Rpred >= 0.5 * (R_min[0] + R_min[1]):
        Y0 = golden * (R_min[0] - Rpred)
    else:
        Y0 = golden * (R_min[1] - Rpred)
    
    fi = [fi_min[i] - fi_pred for i in range(3)]
    dR = [R_min[i] - Rpred for i in range(3)]

    al = [
        dR[0] * fi[1] * fi[2] * (dR[1] - dR[2]),
        dR[1] * fi[0] * fi[2] * (dR[2] - dR[0]),
        dR[2] * fi[0] * fi[1] * (dR[0] - dR[1]),
        dR[1] * dR[2] * fi[0] * (dR[1] - dR[2]),
        dR[0] * dR[2] * fi[1] * (dR[2] - dR[0]),
        dR[0] * dR[1] * fi[2] * (dR[0] - dR[1])
    ]

    delt = []

    if al[0] + al[1] + al[3]:
        delt.append(-1/(al[0] + al[1] + al[2]))
    else:
        return Rpred
    
    if dR[0] * dR[1] * dR[2] * (dR[0] - dR[1]) * (dR[0] - dR[2]) * (dR[1] - dR[2]):
        delt.append( 1/(dR[0]*dR[1]*dR[2]*(dR[0] - dR[1])*(dR[0] - dR[2])*(dR[1] - dR[2])))
    else:
        return Rpred
    
    a = [(dR[0]*al[0] + dR[1]*al[1] + dR[2]*al[2])*delt[0], (al[3] + al[4] + al[5])*delt[1]]
    b = [-(dR[0]*al[3] + dR[1]*al[4] + dR[2]*al[5])*delt[0], -(al[3]*(dR[1] + dR[2]) + al[4]*(dR[0] + dR[2]) + al[5]*(dR[0] + dR[1]))*delt[1]]
    c = [(al[3] + al[4] + al[5])*delt[0], (dR[1]*dR[2]*al[3] + dR[0]*dR[2]*al[4] + dR[0]*dR[1]*al[5])*delt[1]]

    ac = b[1]**2 - 3*a[1]*c[1]
    if ac < 0:
        Z2 = 0
    else:
        Z2 = (- b[1] + ac**0.5) / (3 * a[1])

    abc = b[0]**2 - a[0]*b[0]*c[0]
    if abc < 0:
        Z1 = 0
    else:
        Z1 = (- b[0] + abc**0.5) / c[0]
    
    if (Rpred + Z1) > R_min[0] and (Rpred + Z1) < R_min[1] and Z1:
        Y1 = Z1
    else:
        Y1 = Z2
    
    if abs(Y1) > abs(Y0) or (Rpred + Y1) < R_min[0] or (Rpred + Y1) > R_min[1]:
        Rr = Rpred + Y0
    else:
        Rr = Rpred + Y1
    
    with open(file_Bl13, "a") as file:
        for i in range(len(R_min)):
            file.write('R%s: %s\n' % (i+1, R_min[i]))
        file.write('Rpred: %s\n' % Rpred)
        for i in range(len(fi)):
            file.write('fi%s: %s\n' % (i+1, fi[i]))
        for i in range(len(dR)):
            file.write('dR%s: %s\n' % (i+1, dR[i]))
        for i in range(len(al)):
            file.write('al%s: %s\n' % (i+1, al[i]))
        for i in range(len(delt)):
            file.write('delt%s: %s\n' % (i+1, delt[i]))
        for i in range(len(a)):
            file.write('a%s: %s\n' % (i+1, a[i]))
            file.write('b%s: %s\n' % (i+1, b[i]))
            file.write('c%s: %s\n' % (i+1, c[i]))
        file.write('Y0: %s\n' % Y0)
        file.write('Y1: %s\n' % Y1)
        file.write('Z1: %s\n' % Z1)
        file.write('Z2: %s\n' % Z2)

    return Rr

def Bl19(pr):
    if pr > 0:
        R = 3900
        V = - 3.7
    
    if uv.t[Num - 1] - uv.t[0] < 100:
    # if 1:
        print(uv.t[Num - 1])
        sigmaR0 = 700
        sigmaR1 = 2
    else:
        sigmaR0 = 500
        sigmaR1 = 1
    
    return 0

# Рассчет fi(V)
def spFiV(R, V, n, V_n, Pr):
    '''
    Функция вычисления поправок и значения критерия.\n
    R - модуль радиус-вектора\n
    V - модуль вектора скорости\n
    n- еденичный направляющий вектор\n
    V_n - аналогично скорости\n
    Pr - признаки вычисления
    '''
    
    # Радиус-вектор и скорость объекта
    Rn = R*n
    V_Rn = V*n + R*V_n

    # Их производные по соответствующим компонентам
    dRndu  = R*Vector([- u/uu, 0, 1])
    dRndv  = R*Vector([- v/uu, 1, 0])
    dRndVu = R*Vector([0, 0, 0])
    dRndVv = R*Vector([0, 0, 0])

    dV_Rndu = V*Vector([- u/uu, 0, 1]) + R*Vector([- (u*v*V_v - v**2*V_u + V_u) / uu**3, 0, 0])
    dV_Rndv = V*Vector([- v/uu, 1, 0]) + R*Vector([- (u*v*V_u - u**2*V_v + V_v) / uu**3, 0, 0])
    dV_RndV_u = V*Vector([0, 0, 0]) + R*Vector([-u / uu, 0, 1])
    dV_RndV_v = V*Vector([0, 0, 0]) + R*Vector([-v / uu, 1, 0])

    # Вычисление 5-ых производных по времени
    d5rdt5 = drdt(Rn, V_Rn, R0st, Aw)
    d6rdt5du = drdx(Rn, dRndu, V_Rn, dV_Rndu, R0st, Aw)
    d6rdt5dv = drdx(Rn, dRndv, V_Rn, dV_Rndv, R0st, Aw)
    d6rdt5dVu = drdx(Rn, dRndVu, V_Rn, dV_RndV_u, R0st, Aw)
    d6rdt5dVv = drdx(Rn, dRndVv, V_Rn, dV_RndV_v, R0st, Aw)

    const = 0
    cofs = [[0 for i in range(4)] for i in range(4)]
    x = [0 for i in range(4)]
    
    with open(file_spFiV, 'a') as file:
        file.write('R: %s; V: %s\n' % (R, V))

    for k in range(Num):
        t = uv.t[k] - Tc
        
        # Вектор координат, экстраполированный на время k-го замера
        Rk = d5rdt5[0] + d5rdt5[1]*t + 1/2*d5rdt5[2]*t**2 + 1/6*d5rdt5[3]*t**3 + 1/24*d5rdt5[4]*t**4 + 1/120*d5rdt5[5]*t**5
        with open(file_spFiV, 'a') as file:
            file.write('Rk: %s; t: %s\n' % ((Rk**2)**0.5, t))
            file.write('X0.x;X0.y;X0.z: %s;%s;%s\n' % (d5rdt5[0].x, d5rdt5[0].y, d5rdt5[0].z))
            file.write('X1.x;X1.y;X1.z: %s;%s;%s\n' % (d5rdt5[1].x, d5rdt5[1].y, d5rdt5[1].z))
            file.write('X2.x;X2.y;X2.z: %s;%s;%s\n' % (0.5*d5rdt5[2].x, 0.5*d5rdt5[2].y, 0.5*d5rdt5[2].z))
            file.write('X3.x;X3.y;X3.z: %s;%s;%s\n' % (1/6*d5rdt5[3].x, 1/6*d5rdt5[3].y, 1/6*d5rdt5[3].z))
            file.write('X4.x;X4.y;X4.z: %s;%s;%s\n' % (1/24*d5rdt5[4].x, 1/24*d5rdt5[4].y, 1/24*d5rdt5[4].z))
            file.write('X5.x;X5.y;X5.z: %s;%s;%s\n' % (1/120*d5rdt5[5].x, 1/120*d5rdt5[5].y, 1/120*d5rdt5[5].z))
        # print((Rk**2)**0.5)
        # print('X0', d5rdt5[0].x, d5rdt5[0].y, d5rdt5[0].z)
        # print('X1', d5rdt5[1].x, d5rdt5[1].y, d5rdt5[1].z)
        # print('X2', 0.5*d5rdt5[2].x, 0.5*d5rdt5[2].y, 0.5*d5rdt5[2].z)
        # print('X3', 1/6*d5rdt5[3].x, 1/6*d5rdt5[3].y, 1/6*d5rdt5[3].z)
        # print('X4', 1/24*d5rdt5[4].x, 1/24*d5rdt5[4].y, 1/24*d5rdt5[4].z)
        # print('X5', 1/120*d5rdt5[5].x, 1/120*d5rdt5[5].y, 1/120*d5rdt5[5].z)

        # Вектор измерений
        Rt = Vector([0, uv.v[k], uv.u[k]])

        # Вектор весов
        p = Vector([0, uv.pv[k], uv.pu[k]])

        # Вектор разности эекстрополяции и он же с весами
        dRk = 1/(Rk**2)**0.5 * Rk - Rt
        pdRk = Vector([p.x*dRk.x, p.y*dRk.y, p.z*dRk.z])
        # print(Rk.y)

        # Вектора производных по u, v и их скоростям от производных по времени на время k-го замера
        dxdu  = d6rdt5du[0]  + d6rdt5du[1]*t  + 1/2*d6rdt5du[2]*t**2  + 1/6*d6rdt5du[3]*t**3  #+ 1/24*d6rdt5du[4]*t**4  + 1/120*d6rdt5du[5]*t**5
        dxdv  = d6rdt5dv[0]  + d6rdt5dv[1]*t  + 1/2*d6rdt5dv[2]*t**2  + 1/6*d6rdt5dv[3]*t**3  #+ 1/24*d6rdt5dv[4]*t**4  + 1/120*d6rdt5dv[5]*t**5
        dxdVu = d6rdt5dVu[0] + d6rdt5dVu[1]*t + 1/2*d6rdt5dVu[2]*t**2 + 1/6*d6rdt5dVu[3]*t**3 #+ 1/24*d6rdt5dVu[4]*t**4 + 1/120*d6rdt5dVu[5]*t**5
        dxdVv = d6rdt5dVv[0] + d6rdt5dVv[1]*t + 1/2*d6rdt5dVv[2]*t**2 + 1/6*d6rdt5dVv[3]*t**3 #+ 1/24*d6rdt5dVv[4]*t**4 + 1/120*d6rdt5dVv[5]*t**5

        with open(file_spFiV, 'a') as file:
            file.write('dxdu: %s\n' % (dxdu**2)**0.5)
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (d6rdt5du[0].x, d6rdt5du[0].y, d6rdt5du[0].z))
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (d6rdt5du[1].x, d6rdt5du[1].y, d6rdt5du[1].z))
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (0.5*d6rdt5du[2].x, 0.5*d6rdt5du[2].y, 0.5*d6rdt5du[2].z))
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (1/6*d6rdt5du[3].x, 1/6*d6rdt5du[3].y, 1/6*d6rdt5du[3].z))
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (1/24*d6rdt5du[4].x, 1/24*d6rdt5du[4].y, 1/24*d6rdt5du[4].z))
            file.write('dxdu.x;dxdu.y;dxdu.z: %s;%s;%s\n' % (1/120*d6rdt5du[5].x, 1/120*d6rdt5du[5].y, 1/120*d6rdt5du[5].z))
            file.write('dxdv: %s\n' % (dxdv**2)**0.5)
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (d6rdt5dv[0].x, d6rdt5dv[0].y, d6rdt5dv[0].z))
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (d6rdt5dv[1].x, d6rdt5dv[1].y, d6rdt5dv[1].z))
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (0.5*d6rdt5dv[2].x, 0.5*d6rdt5dv[2].y, 0.5*d6rdt5dv[2].z))
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (1/6*d6rdt5dv[3].x, 1/6*d6rdt5dv[3].y, 1/6*d6rdt5dv[3].z))
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (1/24*d6rdt5dv[4].x, 1/24*d6rdt5dv[4].y, 1/24*d6rdt5dv[4].z))
            file.write('dxdv.x;dxdv.y;dxdv.z: %s;%s;%s\n' % (1/120*d6rdt5dv[5].x, 1/120*d6rdt5dv[5].y, 1/120*d6rdt5dv[5].z))
            file.write('dxdVu: %s\n' % (dxdVu**2)**0.5)
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (d6rdt5dVu[0].x, d6rdt5dVu[0].y, d6rdt5dVu[0].z))
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (d6rdt5dVu[1].x, d6rdt5dVu[1].y, d6rdt5dVu[1].z))
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (0.5*d6rdt5dVu[2].x, 0.5*d6rdt5dVu[2].y, 0.5*d6rdt5dVu[2].z))
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (1/6*d6rdt5dVu[3].x, 1/6*d6rdt5dVu[3].y, 1/6*d6rdt5dVu[3].z))
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (1/24*d6rdt5dVu[4].x, 1/24*d6rdt5dVu[4].y, 1/24*d6rdt5dVu[4].z))
            file.write('dxdVu.x;dxdVu.y;dxdVu.z: %s;%s;%s\n' % (1/120*d6rdt5dVu[5].x, 1/120*d6rdt5dVu[5].y, 1/120*d6rdt5dVu[5].z))
            file.write('dxdVv: %s\n' % (dxdVv**2)**0.5)
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (d6rdt5dVv[0].x, d6rdt5dVv[0].y, d6rdt5dVv[0].z))
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (d6rdt5dVv[1].x, d6rdt5dVv[1].y, d6rdt5dVv[1].z))
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (0.5*d6rdt5dVv[2].x, 0.5*d6rdt5dVv[2].y, 0.5*d6rdt5dVv[2].z))
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (1/6*d6rdt5dVv[3].x, 1/6*d6rdt5dVv[3].y, 1/6*d6rdt5dVv[3].z))
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (1/24*d6rdt5dVv[4].x, 1/24*d6rdt5dVv[4].y, 1/24*d6rdt5dVv[4].z))
            file.write('dxdVv.x;dxdVv.y;dxdVv.z: %s;%s;%s\n' % (1/120*d6rdt5dVv[5].x, 1/120*d6rdt5dVv[5].y, 1/120*d6rdt5dVv[5].z))

        # Они же нормированные на длину
        # dxdu = 1/(dxdu**2)**0.5 * dxdu
        # dxdv = 1/(dxdv**2)**0.5 * dxdv
        # dxdVu = 1/(dxdVu**2)**0.5 * dxdVu
        # dxdVv = 1/(dxdVv**2)**0.5 * dxdVv

        # Вектор разности экстополяции для производных
        dxdu = 1/(Rk**2)**0.5 * (dxdu - 1/Rk**2 * Rk*(Rk*dxdu))
        dxdv = 1/(Rk**2)**0.5 * (dxdv - 1/Rk**2 * Rk*(Rk*dxdv))
        dxdVu = 1/(Rk**2)**0.5 * (dxdVu - 1/Rk**2 * Rk*(Rk*dxdVu))
        dxdVv = 1/(Rk**2)**0.5 * (dxdVv - 1/Rk**2 * Rk*(Rk*dxdVv))

        # Вектора производных с весами
        pdxdu = Vector([p.x*dxdu.x, p.y*dxdu.y, p.z*dxdu.z])
        pdxdv = Vector([p.x*dxdv.x, p.y*dxdv.y, p.z*dxdv.z])
        pdxdVu = Vector([p.x*dxdVu.x, p.y*dxdVu.y, p.z*dxdVu.z])
        pdxdVv = Vector([p.x*dxdVv.x, p.y*dxdVv.y, p.z*dxdVv.z])

        # Рассчет коэффициентов системы уравнения
        const   += pdRk*dRk

        x[0]    += pdRk*dxdu
        x[1]    += pdRk*dxdv
        x[2]    += pdRk*dxdVu
        x[3]    += pdRk*dxdVv

        cof0 = [pdxdu*dxdu,  pdxdu*dxdv,  pdxdu*dxdVu,  pdxdu*dxdVv]
        cof1 = [pdxdv*dxdu,  pdxdv*dxdv,  pdxdv*dxdVu,  pdxdv*dxdVv]
        cof2 = [pdxdVu*dxdu, pdxdVu*dxdv, pdxdVu*dxdVu, pdxdVu*dxdVv]
        cof3 = [pdxdVv*dxdu, pdxdVv*dxdv, pdxdVv*dxdVu, pdxdVv*dxdVv]

        cofs[0] = [x+y for x,y in zip(cofs[0], cof0)]
        cofs[1] = [x+y for x,y in zip(cofs[1], cof1)]
        cofs[2] = [x+y for x,y in zip(cofs[2], cof2)]
        cofs[3] = [x+y for x,y in zip(cofs[3], cof3)]
    
    # Решение СЛАУ методом Гаусса
    mx = [-x[0], -x[1], -x[2], -x[3]]
    
    ans = Gauss(cofs, mx)
    
    # Рассчет значения fi(V)
    fi_V = const + x[0]*ans[0] + x[1]*ans[1] + x[2]*ans[2] + x[3]*ans[3]
    # print(fi_V)

    with open(file_spFiV, 'a') as file:
        file.write('const: %s\nx[0]: %s\nx[1]: %s\nx[2]: %s\nx[3]: %s\n' % (const, x[0], x[1], x[2], x[3]))
        file.write('ans[0]: %s\nans[1]: %s\nans[2]: %s\nans[3]: %s\n' % (ans[0], ans[1], ans[2], ans[3]))
        file.write('%s\n' % fi_V)

    # Запись результата
    res      = Target()
    res.R    = R
    res.V    = V
    res.fi_R = fi_V
    res.fi_V = fi_V
    res.du   = ans[0]
    res.dv   = ans[1]
    res.dV_u = ans[2]
    res.dV_v = ans[3]

    return res

# Рассчет границ скорости
def V_borders(t = 0, R = 0, n = Vector([0, 0, 0]), V_n = Vector([0, 0, 0]), borders = (0, 0)):
    '''
    Функция вычисления границ скорости. Создает радиус-вектор, находит его корни, делает обратное преобразование,
    находит пересечения со стандартными значениями.\n
    t - момент времени для вычисления\n
    R - модуль расстояния\n
    n - еденичный направляющий вектор\n
    V_n - аналогично скорость
    '''

    # Радиус-вектор к объекту и коэффициенты его неравенства x(t) > 0
    Rn = R * n

    a = ak(Rn, R0st)
    b = bk(Rn, R0st, Aw)
    c = ck(Rn, R0st, Aw)

    a = (                                                 1/24* a  *t**4 ).x
    b = (        b[0]*t + 1/2*b[1]*t**2 + 1/6*b[2]*t**3 + 1/24*b[3]*t**4 ).x
    c = ( c[0] + c[1]*t + 1/2*c[2]*t**2 + 1/6*c[3]*t**3 + 1/24*c[4]*t**4 ).x
    
    d = b**2 - 4*a*c
    
    if a == 0 or (d < 0 and a < 0):
        print('a == 0 or (d < 0 and a < 0):', a, d)
        return (3, borders)

    if d >= 0:
        # Корни радиус-вектора
        V1 = (- b - d**0.5) / (2 * a)
        V2 = (- b + d**0.5) / (2 * a)

        # Обратное преобразование для получения корней модуля
        V1 = (V1 - R * V_n.x) / n.x
        V2 = (V2 - R * V_n.x) / n.x
        
        with open(file_V_borders, 'a') as file:
            file.write('R: %s; t: %s\n' % (R, t))
            file.write('V1: %s\nV2: %s\n' % (V1, V2))

        if a >= 0:
            if V1 < borders[0] and V2 > borders[1]:
                print('V1 < borders[0] and V2 > borders[1]')
                return (3, borders)
            elif V1 > borders[0] and V2 < borders[1]:
                return (0, (V2, borders[1]))
            elif V1 > borders[0] and V2 > borders[1]:
                return (0, (borders[0], V1))
            elif V1 < borders[0] and V2 < borders[1]:
                return (0, (V2, borders[1]))
        else:
            if V1 < borders[0] or V2 > borders[1]:
                print('V1 < borders[0] or V2 > borders[1]')
                return (3, borders)
            elif V2 < borders[0] and V1 < borders[1]:
                return (0, (borders[0], V1))
            elif V2 > borders[0] and V1 > borders[1]:
                return (0, (V2, borders[1]))
            elif V2 > borders[0] and V1 < borders[1]:
                return (0, (V2, V1))
            else:
                return (0, borders)
    else:
        print('d<0')
        return (3, borders)

# Рассчет fi(R)
def spFiR(min_points, R, ch, Pr):
    '''
    Функция рассчета поправок к параметрам движения. Возвращает объект Target с подсчитанными значениями в случаи успеха и входные
    параметры с нулевой скоростью в случаи неудачи.\n
    min_points - сохраненные точки минимума\n
    R - текущее расстояние для рассчета\n
    ch - счетчик\n
    Pr - признаки
    '''

    chn = 0

    # Создаем еденичные вектора направления и скорости
    n = Vector([uu, v, u])
    V_n = Vector([-(V_v*v + V_u*u)/uu, V_v, V_u])
    
    if 100 - (R*V_n)**2 >= 0.01: # 100 - V**2 >= 0.01
        Vprav = (100 - (R*V_n)**2) ** 0.5
    else:
        Vprav = 0.1

    Vlev = -Vprav
    
    if Pr[2] == 0:

        t = [uv.t[0] - Tc, uv.t[Num - 1] - Tc]
        
        for i in range(2):
            res_bd = V_borders(t[i], R, n, V_n, (Vlev, Vprav))
            if res_bd[0]:
                Pr[0] = res_bd[0]
                break
            Vlev = res_bd[1][0]
            Vprav = res_bd[1][1]
    
    with open(file_spFiR, 'a') as file:
        file.write('V borders: %s;%s\n' % (Vlev, Vprav))

    if Pr[0]:
        return spFiV(R, 0, n, V_n, Pr[2])
    
    if ch[1] == 0:
        if abs(Vlev) < abs(Vprav):
            Vst = 0.5 * dbeta * Vlev
        else:
            Vst = 0.5 * dbeta * Vprav
    elif ch[1] == 1:
        Vst = min_points[1]
    else:
        Vst = 0
        for i in range(0, ch[0] * 3, 3):
            dd1 = 1
            for j in range(0, ch[0] * 3, 3):
                if i != j and min_points[i] != min_points[j]:
                    dd1 *= (R - min_points[j]) / (min_points[i] - min_points[j])
            Vst += min_points[i+1] * dd1
    
    if (ch[1] >= 1 and (Vst < Vlev or Vst > Vprav)):
        Pr[1] = 0
        if abs(Vlev) < abs(Vprav):
            Vst = 0.5 * dbeta * Vlev
        else:
            Vst = 0.5 * dbeta * Vprav

    dV = minFi_dV
    
    fi_a = spFiV(R, Vst, n, V_n, Pr[2])
    chn += 1
    with open(file_spFiR, 'a') as file:
        file.write('Vst: %s\n' % Vst)
        file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
    # print('V%s' % chn, fi_a.V, fi_a.fi_V)
    
    if Pr[2] != 1:
        Vk = Vprav - 0.01
        Pr[1] = 0

        fi_b = fi_a
        V = fi_b.V + dV
        fi_a = spFiV(R, V, n, V_n, Pr[2])
        chn += 1
        with open(file_spFiR, 'a') as file:
            file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
        # print('V%s' % chn, fi_a.V, fi_a.fi_V)

        if fi_a.fi_V > fi_b.fi_V:
            dV = - dV
            if dV < 0:
                Vk = Vprav - 0.01
            else:
                Vk = Vlev + 0.01
            V1 = fi_a.V
            fi_V1 = fi_a.fi_V

        #     fi_b = fi_a
        # else:
        #     V1 = fi_b.V
        #     fi_V1 = fi_b.fi_V

        V = fi_b.V + dV
        fi_a = spFiV(R, V, n, V_n, Pr[2])
        chn += 1
        with open(file_spFiR, 'a') as file:
            file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
        # print('V%s' % chn, fi_a.V, fi_a.fi_V)
        
        if fi_a.fi_V < fi_b.fi_V:
            V1 = fi_b.V
            fi_V1 = fi_b.fi_V
            fi_b = fi_a
            while V != Vk:
                V = fi_b.V + dV
                # fi_b = fi_a
                fi_a = spFiV(R, V, n, V_n, Pr[2])
                chn += 1
                with open(file_spFiR, 'a') as file:
                    file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
                # print('V%s' % chn, fi_a.V, fi_a.fi_V)

                if fi_a.fi_V > fi_b.fi_V:
                    break

                V1 = fi_b.V
                fi_V1 = fi_b.fi_V

                fi_b = fi_a
            else:
                if dV < 1:
                    Pr[1] = -1
                else:
                    Pr[1] = 1
                fi_a.R = R
                Pr[2] = 0
                print_min_points(fi_a, min_points, ch, Pr)
                return fi_a

        Vpred = fi_b.V
        
        if fi_a.V > fi_b.V:
            V2 = fi_a.V
            fi_V2 = fi_a.fi_V
        else:
            V2 = V1
            fi_V2 = fi_V1
            V1 = fi_a.V
            fi_V1 = fi_a.fi_V
        
        # print(V1, V2, Vpred, fi_b.V, fi_V1, fi_V2)
        V = Bl13_1([V1, V2], [fi_V1, fi_V2], fi_b.V)

        fi_a = spFiV(R, V, n, V_n, Pr[2])
        chn += 1
        with open(file_spFiR, 'a') as file:
            file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
        # print('V%s' % chn, fi_a.V, fi_a.fi_V)
        if fi_a.fi_V > fi_b.fi_V:
            if fi_a.V > fi_b.V:
                V3 = V2
                fi_V3 = fi_V2
                V2 = fi_a.V
                fi_V2 = fi_a.fi_V
            else:
                V3 = V1
                fi_V3 = fi_V1
                V1 = fi_a.V
                fi_V1 = fi_a.fi_V
        else:
            if fi_a.V > fi_b.V:
                V3 = V1
                fi_V3 = fi_V1
                V1 = fi_b.V
                fi_V1 = fi_b.fi_V
                fi_b = fi_a
            else:
                V3 = V2
                fi_V3 = fi_V2
                V2 = fi_b.V
                fi_V2 = fi_b.fi_V
                fi_b = fi_a
        
        # print('V1, V2, V3, fi_a.V, Vpred:', V1, V2, V3, fi_a.V, Vpred)
        # print('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
        while abs(V2 - V1) > 0.0001 and abs(Vpred - fi_a.V) > 0.0001 and chn < 100:
            Vpred = fi_a.V
            
            V = Bl13_2([V1, V2, V3], [fi_V1, fi_V2, fi_V3], fi_b.V, fi_b.fi_V)

            fi_a = spFiV(R, V, n, V_n, Pr[2])
            chn += 1
            with open(file_spFiR, 'a') as file:
                file.write('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s\n' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})
            # print('V%(i)s; fi_V%(i)s: %(V)s; %(fi)s' % {'i': chn, 'V': fi_a.V, 'fi': fi_a.fi_V})

            V21 = V2
            V31 = V3
            fi_V31 = fi_V3
            
            if fi_a.fi_V > fi_b.fi_V:
                if fi_a.V > fi_b.V:
                    V3 = V2
                    fi_V3 = fi_V2
                    V2 = fi_a.V
                    fi_V2 = fi_a.fi_V
                else:
                    V3 = V1
                    fi_V3 = fi_V1
                    V1 = fi_a.V
                    fi_V1 = fi_a.fi_V
            else:
                if fi_a.V > fi_b.V:
                    V3 = V1
                    fi_V3 = fi_V1
                    V1 = fi_b.V
                    fi_V1 = fi_b.fi_V
                    fi_b = fi_a
                else:
                    V3 = V2
                    fi_V3 = fi_V2
                    V2 = fi_b.V
                    fi_V2 = fi_b.fi_V
                    fi_b = fi_a

            if (V21 > V31 and V2 != V21) or (V21 < V31 and V2 == V21) and abs(V1 - V3) > abs(V31 - V2):
                V3 = V31
                fi_V3 = fi_V31

            # print(V1, V2, V3, fi_b.V)
            if chn > 100:
                fi_a.R = R
                fi_a.fi_R = fi_a.fi_V
                Pr[0] = 4
                return fi_a
        


        if abs(V2 - V1) <= 0.0001 or abs(Vpred - fi_a.V) < 0.0001:
            if fi_a.fi_V > fi_b.fi_V:
                fi_a = fi_b
                fi_a.R = R
            else:
                fi_a.R = R
                fi_a.fi_R = fi_a.fi_V
        
        if abs(Vst - fi_a.V) < 0.0001 and Pr[3] == 1:
            Pr[2] = 1
        else:
            Pr[2] = 0
        
        ch[1] += 1
        print_min_points(fi_a, min_points, ch, Pr)
        return fi_a
    else:
        ch[1] += 1
        print_min_points(fi_a, min_points, ch, Pr)
        return fi_a

def print_min_points(fi_a, min_points, ch, Pr):
    if ch[0] == 5 and Pr[3] != 0:
        min_points[15] = fi_a.R
        min_points[16] = fi_a.V
        min_points[17] = fi_a.fi_V
    
    if ch[0] == 5 and Pr[3] == 0:
        ss = 0.0001
        kk1 = 18
        for i in range(0,17,3):
            if ss < abs(min_points[i] - fi_a.R):
                ss = abs(min_points[i] - fi_a.R)
                kk1 = i
        if kk1 != 18:
            min_points[kk1] = fi_a.R
            min_points[kk1+1] = fi_a.V
            min_points[kk1+2] = fi_a.fi_V
    
    if ch[0] <= 4:
        min_points[ch[0]*3] = fi_a.R
        min_points[ch[0]*3+1] = fi_a.V
        min_points[ch[0]*3+2] = fi_a.fi_V
        ch[0] += 1

def procMinimizeFunc(file):
    '''
    Главный алгоритм
    '''

    print('Start main program')
    # Возвращаемое занчение; 0 в случаи неудачи, 1 в случаи выполнения
    ret_value = 0
    # Четыре признака
    Pr = [0, 0, 0, 0]
    # Счетчики числа вычислений радиальной скорости и дальности
    ch = [0, 0]

    # 4 точки минимума и значение функции в них
    R_min = [None, None, None, None]
    V_min = [None, None, None, None]
    fi_min = [None, None, None, None]
    min_points = [0 for i in range(18)]

    # Формирование констант станции и вращения Земли для рассчетов
    global R0st, Aw, wst
    data = pd.read_csv(file + '/slava_geo.txt', sep=';')
    R0st = Vector([data['R0xt'][0], data['R0yt'][0], data['R0zt'][0]])
    wst = Vector([data['wxt'][0], data['wyt'][0], data['wzt'][0]])
    Aw = Matrix([[0, -wst.z, wst.y], [wst.z, 0, -wst.x], [-wst.y, wst.x, 0]])

    # Количество опорных точек и взвешенное время
    global Num, Tc
    data = pd.read_csv(file + '/slava_secpoints.txt', sep=';')
    Num = int(data['Num'][0])
    Tc = data['Tc'][0]

    # Объекты для хранения переменных на текущем и предыдущем шагу
    fi_now, fi_last = Target(), Target()

    # Создание и запись опорных точек и точки сглаживания
    start_point(file)

    # dR - величина шага по дальности
    # Rk - граница в сторону движения
    dR = minFi_dR[0]
    Rk = Rprav

    fi_last = spFiR(min_points, 0.5 * (Rlev + Rprav), ch, Pr)
    with open(file_procMinimizeFunc, 'a') as file:
        file.write('First calculate in the middle of R value is done. Parameters:\nR: %s\nV: %s\nfi_V: %s\n' % (fi_last.R, fi_last.V, fi_last.fi_R))
    print('First calculate in the middle of R value is done. Parameters:\nR:%s\nV:%s\nfi_V:%s' % (fi_last.R, fi_last.V, fi_last.fi_R))
    
    if Pr[0]:
        Bl19(Pr[0])
        return ret_value
    
    if Pr[1] < 0:
        dR = - dR
        Rk = Rlev
    
    Rr = fi_last.R + dR
    
    fi_now = spFiR(min_points, Rr, ch, Pr)

    if Pr[0]:
        Bl19(Pr[0])
        return ret_value

    with open(file_procMinimizeFunc, 'a') as file:
        file.write('Second calculete is done. Parameters:\nR: %s\nV: %s\nfi_V: %s\n' % (fi_now.R, fi_now.V, fi_now.fi_R))
    print('Second calculete is done', fi_now.R, fi_now.fi_R)

    if fi_now.fi_R - fi_last.fi_R > 0:
        dR = - dR
        R_min[0] = fi_now.R
        V_min[0] = fi_now.V
        fi_min[0] = fi_now.fi_R
        
        if dR > 0:
            Rk = Rprav
        else:
            Rk = Rlev
    else:
        R_min[0] = fi_last.R
        V_min[0] = fi_last.V
        fi_min[0] = fi_last.fi_R
        fi_last = fi_now

    Rr = fi_last.R + dR
    
    fi_now = spFiR(min_points, Rr, ch, Pr)

    with open(file_procMinimizeFunc, 'a') as file:
        file.write('Third calculete is done. Parameters:\nR: %s\nV: %s\nfi_V: %s\n' % (fi_now.R, fi_now.V, fi_now.fi_R))
    print('Third calculate is done', fi_now.R, fi_now.fi_R)

    while fi_now.fi_R - fi_last.fi_R < 0:

        R_min[0] = fi_last.R
        V_min[0] = fi_last.V
        fi_min[0] = fi_last.fi_R
        fi_last = fi_now
        
        if Pr[0] or fi_last.R == Rk:
            print('Eto fiasco, bratan')
            Bl19(Pr[0])
            return ret_value
        
        if abs(fi_now.R - Rk) > minFi_dR[0] + minFi_dR[1]:
            Rr = fi_now.R + dR
        else:
            fi_Rk = spFiR(min_points, Rk, ch, Pr)
            fi_now = spFiR(min_points, Rk - minFi_dR[2], ch, Pr)

            if fi_Rk.fi_R < fi_now.fi_R:
                print('Min outside the sigment')
                Pr[0] = 1
                Bl19(Pr[0])
                return ret_value
            else:
                fi_last = fi_now
                fi_now = fi_Rk
                break
        
        fi_now = spFiR(min_points, Rr, ch, Pr)
        
        with open(file_procMinimizeFunc, 'a') as file:
            file.write('Third+ calculete is done. Parameters:\nR: %s\nV: %s\nfi_V: %s\n' % (fi_now.R, fi_now.V, fi_now.fi_R))
    
    Pr[3] = 1
    Rpred = fi_last.R

    if fi_now.R > fi_last.R:
        R_min[1] = fi_now.R
        V_min[1] = fi_now.V
        fi_min[1] = fi_now.fi_R
    else:
        R_min[1] = R_min[0]
        V_min[1] = V_min[0]
        fi_min[1] = fi_min[0]
        R_min[0] = fi_now.R
        V_min[0] = fi_now.V
        fi_min[0] = fi_now.fi_R

    with open(file_procMinimizeFunc, 'a') as file:
        file.write('First part is done. Parameters:\nR1: %s\nV1: %s\nfi1: %s\nR2: %s\nV2: %s\nfi2: %s\nRs: %s\nVs: %s\nfis: %s\n' % (R_min[0], V_min[0], fi_min[0], R_min[1], V_min[1], fi_min[1], fi_last.R, fi_last.V, fi_last.fi_R))
    print('First part is done. Parameters:\nR1:%s\nV1:%s\nfi1:%s\nR2:%s\nV2:%s\nfi2:%s' % (R_min[0], V_min[0], fi_min[0], R_min[1], V_min[1], fi_min[1]))

    Rr = Bl13_1(R_min, fi_min, fi_last.R)

    fi_now = spFiR(min_points, Rr, ch, Pr)
    
    if Pr[0]:
        Bl19(Pr[0])
        return ret_value
    
    R21, R31, fiR31 = 0, 0, 0

    if fi_now.fi_R > fi_last.fi_R:
        if Rr > fi_last.R:
            R_min[2] = R_min[1]
            V_min[2] = V_min[1]
            fi_min[2] = fi_min[1]
            R_min[1] = Rr
            V_min[1] = fi_now.V
            fi_min[1] = fi_now.fi_R
        else:
            R_min[2] = R_min[0]
            V_min[2] = V_min[0]
            fi_min[2] = fi_min[0]
            R_min[0] = Rr
            V_min[0] = fi_now.V
            fi_min[0] = fi_now.fi_R
    else:
        if Rr > fi_last.R:
            R_min[2] = R_min[0]
            V_min[2] = V_min[0]
            fi_min[2] = fi_min[0]
            R_min[0] = fi_last.R
            V_min[0] = fi_last.V
            fi_min[0] = fi_last.fi_R
            fi_last = fi_now
        else:
            R_min[2] = R_min[1]
            V_min[2] = V_min[1]
            fi_min[2] = fi_min[1]
            R_min[1] = fi_last.R
            V_min[1] = fi_last.V
            fi_min[1] = fi_last.fi_R
            fi_last = fi_now

    with open(file_procMinimizeFunc, 'a') as file:
        file.write('Second part is done. Parameters:\nR1: %s\nV1: %s\nfi1: %s\nR2: %s\nV2: %s\nfi2: %s\nR3: %s\nV3: %s\nfi3: %s\nRs: %s\nVs: %s\nfis: %s\n' % (R_min[0], V_min[0], fi_min[0], R_min[1], V_min[1], fi_min[1], R_min[2], V_min[2], fi_min[2], fi_last.R, fi_last.V, fi_last.fi_R))
    print('Second part is done. Parameters:\nR1:%s\nV1:%s\nfi1:%s\nR2:%s\nV2:%s\nfi2:%s\nR3:%s\nV3:%s\nfi3:%s' % (R_min[0], V_min[0], fi_min[0], R_min[1], V_min[1], fi_min[1], R_min[2], V_min[2], fi_min[2]))

    while abs(R_min[1] - R_min[0]) > 0.1 and abs(Rpred - Rr) > 0.1 and ch[1] < 100:

        Rpred = fi_now.R
        
        # print('Rr before:', Rr)
        Rr = Bl13_2(R_min, fi_min, fi_last.R, fi_last.fi_R)
        # print('Rr after:', Rr)
        
        fi_now = spFiR(min_points, Rr, ch, Pr)

        if Pr[0]:
            Bl19(Pr[0])
            return ret_value
        
        aa = min_points[2]
        ii = 18

        if min_points[15] != 0:
            for ik in range(0, 18, 3):
                if min_points[ik] < R_min[1] or min_points[ik] > R_min[2]:
                    if aa <= min_points[ik+2]:
                        aa = min_points[ik+2]
                        ii = ik
            
            if ii != 18:
                min_points[ii] = fi_now.R
                min_points[ii+1] = fi_now.V
                min_points[ii+2] = fi_now.fi_R
        
        R21 = R_min[1]
        R31 = R_min[2]
        fiR31 = fi_min[2]

        if fi_now.fi_R > fi_last.fi_R:
            if Rr > fi_last.R:
                R_min[2] = R_min[1]
                V_min[2] = V_min[1]
                fi_min[2] = fi_min[1]
                R_min[1] = Rr
                V_min[1] = fi_now.V
                fi_min[1] = fi_now.fi_R
            else:
                R_min[2] = R_min[0]
                V_min[2] = V_min[0]
                fi_min[2] = fi_min[0]
                R_min[0] = Rr
                V_min[0] = fi_now.V
                fi_min[0] = fi_now.fi_R
        else:
            if Rr > fi_last.R:
                R_min[2] = R_min[0]
                V_min[2] = V_min[0]
                fi_min[2] = fi_min[0]
                R_min[0] = fi_last.R
                V_min[0] = fi_last.V
                fi_min[0] = fi_last.fi_R
                fi_last = fi_now
            else:
                R_min[2] = R_min[1]
                V_min[2] = V_min[1]
                fi_min[2] = fi_min[1]
                R_min[1] = fi_last.R
                V_min[1] = fi_last.V
                fi_min[1] = fi_last.fi_R
                fi_last = fi_now

        if R21 > R31 and R_min[1] != R21 and abs(R_min[0] - R_min[2]) > abs(R31 - R_min[1]):
            R_min[2] = R31
            fi_min[2] = fiR31
        if R21 < R31 and R_min[1] == R21 and abs(R_min[0] - R_min[2]) > abs(R31 - R_min[1]):
            R_min[2] = R31
            fi_min[2] = fiR31
        
        with open(file_procMinimizeFunc, 'a') as file:
            file.write('Second+ part is done. Parameters:\nR1: %s\nV1: %s\nfi1: %s\nR2: %s\nV2: %s\nfi2: %s\nR3: %s\nV3: %s\nfi3: %s\nRs: %s\nVs: %s\nfis: %s\n' % (R_min[0], V_min[0], fi_min[0], R_min[1], V_min[1], fi_min[1], R_min[2], V_min[2], fi_min[2], fi_last.R, fi_last.V, fi_last.fi_R))
        
    else:
        if ch[1] > 100:
            Pr[0] = 2
            Bl19(Pr[0])
            return ret_value
    
    if fi_now.fi_R > fi_last.fi_R:
        fi_now = fi_last

    if abs(u + fi_now.du) <= 1 and abs(v + fi_now.dv) <= 1:
        R = fi_now.R
        V = fi_now.V
        fi = fi_now.fi_R
        Eps = acos(u + fi_now.du)
        Tet = asin(v + fi_now.dv)
        V_Eps = - (V_u + fi_now.dV_u) / sin(Eps)
        V_Tet = (V_v + fi_now.dV_v) / cos(Tet)
        with open(file_procMinimizeFunc, 'a') as file:
            file.write('Final par:\n')
            file.write('R: %s\n' % R)
            file.write('V: %s\n' % V)
            file.write('fi: %s\n' % fi)
            file.write('u+du: %s\n' % (u + fi_now.du))
            file.write('v+dv: %s\n' % (v + fi_now.dv))
            file.write('V_u + dV_u: %s\n' % (V_u + fi_now.dV_u))
            file.write('V_v + dV_v: %s\n' % (V_v + fi_now.dV_v))
        print('Final par:')
        print('R:', R)
        print('V:', V)
        print('fi:', fi)
        print('u+du:', u + fi_now.du)
        print('v+dv:', v + fi_now.dv)
        print('V_u + dV_u:', V_u + fi_now.dV_u)
        print('V_v + dV_v:', V_v + fi_now.dV_v)
    else:
        return ret_value
    
    if R != 7000:
        ret_value = 1
    
    Bl19(Pr[0])

    return ret_value