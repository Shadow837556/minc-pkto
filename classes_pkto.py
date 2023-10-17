import numpy as np

MAXPOINTS = 100

class Vector:
    def __init__(self, array = [None, None, None]) -> None:
        self.array = np.array(array)
        self.x = self.array[0]
        self.y = self.array[1]
        self.z = self.array[2]
    
    def __str__(self):
        return 'Vector({}, {}, {})'.format(self.x, self.y, self.z)

    def __add__(self, other):
        return Vector(self.array + other.array)
    
    def __sub__(self, other):
        return Vector(self.array - other.array)
    
    def __mul__(self, other):
        if type(other) is Vector:
            return other.array @ self.array
        elif type(other) is Matrix:
            return Vector(other.array @ self.array)
        else:
            return Vector(other * self.array)
        
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other):
        return 1/other * self
    
    def __rtruediv__(self, other):
        return 1/other * self

    def __pow__(self, other):
        if other % 2 == 0:
            return (self.array @ self.array) ** (other // 2)
        else:
            return Vector((self.array @ self.array) ** (other // 2) * self.array)
    
    def __neg__(self):
        return (-1) * self

class Matrix:
    def __init__(self, array = [[None, None, None], [None, None, None], [None, None, None]]):
        self.array = np.array(array)

    def __str__(self):
        return str(self.array)

    def __mul__(self, other):
        if type(other) is Matrix:
            return Matrix(self.array @ other.array)
        elif type(other) is Vector:
            return Vector(self.array @ other.array)
        else:
            return Matrix(other * self.array)
    
    def __rmul__(self, other):
        return self * other
    
    def __pow__(self, other):
        res = self.array
        for i in range(other - 1):
            res = res @ self.array
        return Matrix(res)

    def __neg__(self):
        return (-1) * self

class Target:
    def __init__(self, R = None, V = None, fi_R = None, fi_V = None, du = None, dv = None, dV_u = None, dV_v = None) -> None:
        self.R = R
        self.V = V
        self.fi_R = fi_R
        self.fi_V = fi_V
        self.du = du
        self.dv = dv
        self.dV_u = dV_u
        self.dV_v = dV_v

class PKTO_Error(Exception):
    def __init__(self, name):
        self.name = name
    def __str__(self) -> str:
        print_error = 'This is a ' + self.name + ' error'
        return print_error

# Переменные, в которые происходит запись измерений опорных точек
class case_par:
    t = [None for i in range(MAXPOINTS)]
    u = [None for i in range(MAXPOINTS)]
    v = [None for i in range(MAXPOINTS)]
    d = [None for i in range(MAXPOINTS)]
    dd = [None for i in range(MAXPOINTS)]
    pu = [None for i in range(MAXPOINTS)]
    pv = [None for i in range(MAXPOINTS)]
    pd = [None for i in range(MAXPOINTS)]
    pdd = [None for i in range(MAXPOINTS)]

if __name__ == '__main__':
    print('Start pkto classes test')
    vec1 = Vector([2, 2, 2])
    vec2 = Vector([2, 2, 2])
    print('Make two Vector(2, 2, 2):', vec1, vec2, sep='\n')
    mat1 = Matrix([[2, 2, 2], [3, 3, 3], [4, 4, 4]])
    mat2 = Matrix([[2, 2, 2], [3, 3, 3], [4, 4, 4]])
    print('Make two Matrix([2, 2, 2], [3, 3, 3], [4, 4, 4])', mat1, mat2, sep='\n')
    print('Vector operations:')
    print('Add:', vec1 + vec2)
    print('Sub:', vec1 - vec2)
    print('Mul on number 2:', 2 * vec1)
    print('Mul on Vector:', vec1 * vec2)
    print('Mul on Matrix:', mat1 * vec1)
    print('div 2:', vec1/2)
    print('Pow 2:', vec1**2)
    print('Pow 3:', vec1**3)
    print('Neg:', -vec1)
    print('Matrix opertions:')
    print('Mul on number:', 3*mat1, sep='\n')
    print('Mul on Vector:', mat1*vec1)
    print('Mul on Matrix:', mat1*mat2, sep='\n')
    print('Pow 2:', mat1**2, sep='\n')
    print('Neg:', -mat1, sep='\n')

    # a = ak(R)
    # b = bk(R)
    # c = ck(R)

    # if Pr[2] == 0:
        # for i in range(2):
        #     at = 1/24*a.x*t[i]**4
        #     bt = b[0].x*t + 1/2*b[1].x*t**2 + 1/6*b[2].x*t**3 + 1/24*b[3].x*t**4
        #     ct = c[0].x + c[1].x*t + 1/2*c[2]*t*2 + 1/6*c[3].x*t**3 + 1/24*c[4].x*t**4

        #     d = bt**2 - 4*at*ct

        #     if at == 0:
        #         Pr[0] = 3
        #         break
            
        #     if d < 0 and at < 0:
        #         Pr[0] = 3
        #         break

        #     if d >= 0:
        #         V1 = (- bt - d ** 0.5) / (2 * at)
        #         V2 = (- bt + d ** 0.5) / (2 * at)
            
        #         if at >= 0:
        #             if V2 > Vprav and V1 < Vlev:
        #                 Pr[0] = 3
        #                 break
        #             elif V2 < Vprav and V1 > Vlev:
        #                 Vlev = V2
        #             elif V2 > Vprav and V1 > Vlev:
        #                 Vprav = V1
        #             elif V1 < Vlev and V2 < Vprav:
        #                 Vlev = V2
        #         else:
        #             if V2 > Vprav or V1 < Vlev:
        #                 Pr[0] = 3
        #                 break
        #             elif V2 < Vlev and V1 < Vprav:
        #                 Vprav = V1
        #             elif V2 > Vlev and V1 > Vprav:
        #                 Vlev = V2
        #             elif V2 > Vlev and V1 < Vprav:
        #                 Vlev = V2
        #                 Vprav = V1