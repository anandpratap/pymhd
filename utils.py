from pylab import sqrt
#abs(x) returns mod x

#returns distance between two points
def distance(a, b):
    tmp = pow((a.x-b.x),2) + pow((a.y-b.y),2)
    return sqrt(tmp) 


#class of point in two dimension
class point:
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y


#vector containing 8 components
class vector:
    def __init__(self):
        self.x = zeros(8)

#line containing two points a and b and length        
class line:
    def __init__(self,a=point(),b=point()):
        self.a = a
        self.b = b
        self.length = distance(a, b)

class state:
    def __init__(self):
        self.rho = 1.0
        self.u = 1.0
        self.v = 1.0
        self.w = 1.0
        self.Bx = 1.0
        self.By = 1.0
        self.Bz = 1.0
        self.p = 1.0

class cell:
    def __init__(self):
        self.u = vector()
        self.area = 0.0
        self.corner = [point(), point(), point(), point()]
        self.face = [line(), line(), line(), line()]
        self.center = point()

class matrix:
    def __init__(self):
        self.Y = zeros([8,8])

def calc_tri_area(a, b, c):
    A = distance(b, c)
    B = distance(c, a)
    C = distance(a, b)
    s = (A + B + C)/2.0
    return sqrt(s*(s-A)*(s-B)*(s-C))

def calc_quad_area(a, b, c, d):
    area = calc_tri_area(a, b, c) + calc_tri_area(a, c, d)
    return area

#min(a,b) gives minimum of a,b
def min3(a, b, c):
    return min(min(a,b),min(a,c))

def norm(a):
    n = len(a.X)
    nom = 0.0
    for i in range(n):
        nom += a.X[i]*a.X[i]
    return sqrt(nom)

def sgn(x):
    if x < 0:
        return -1.0
    else:
        return 1.0
    
#dot(a,b) gives dot product of two vectors

def matrixmul(A,b):
    "AdasdASD"
    c = vector()
    for i in range(8):
        c.X[i] = 0.0
        for j in range(8):
            c.X[i] += A.Y[i][j]*b.X[j]
            
    return c

