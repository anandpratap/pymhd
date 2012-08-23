from pylab import sqrt, zeros


def distance(a, b):
    """
    Return distance between two points.
  
    input : point a, point b
    output: float
    """
    tmp = pow((a.x-b.x),2) + pow((a.y-b.y),2)
    return sqrt(tmp) 



class point:
    """
    Define class of point in two dimension

    Attributes: x,y [Default (0,0)]
    """
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y



class vector:
    """
    Define vector containing 8 components

    Attributes: array x of length 8  [default (0)]
    """
    
    def __init__(self):
        self.X = zeros(8)


class line:
    """
    Define line joining two points a and b

    Attributes: point a, point b, float length
    """
    def __init__(self,a=point(),b=point()):
        self.a = a
        self.b = b
        self.length = distance(a, b)

class state:
    """
    Define a class state containing primitive variables.

    Attributes: rho, u, v, w, Bx, By, Bz, p [default 1.0]
    """
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
    """
    Define a cell

    Attributes: vector u, list of 4 'corner' points, list of 4 'face' lines, 
    point center, float area
    """
    def __init__(self):
        self.u = vector()
        self.area = 0.0
        self.corner = [point(), point(), point(), point()]
        self.face = [line(), line(), line(), line()]
        self.center = point()

class matrix:
    """
    Define a 8x8 array

    Attributes: Y
    """
    def __init__(self):
        self.Y = zeros([8,8])

def calc_tri_area(a, b, c):
    """
    Calculate area of a Triangle.
    
    input : point a,b,c
    return: float
    """
    A = distance(b, c)
    B = distance(c, a)
    C = distance(a, b)
    s = (A + B + C)/2.0
    return sqrt(s*(s-A)*(s-B)*(s-C))

def calc_quad_area(a, b, c, d):
    """
    Calculate area of a Quadrilateral

    input : point a,b,c,d
    return: float
    """
    area = calc_tri_area(a, b, c) + calc_tri_area(a, c, d)
    return area

#min(a,b) gives minimum of a,b
def min3(a, b, c):
    """
    Return minimum of three numbers
    """
    return min(min(a,b),min(a,c))

#max(a,b) return maximum of a,b
def max3(a, b, c):
    """
    Return maximum of three numbers
    """
    return max(max(a,b),max(a,c))


def norm(a):
    """
    Return norm of a vector
    """
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
    

def vectormul(a, b):
    """
    Returns dot product of two vector a and b
    """
    tmp_sum = 0.0
    for i in range(8):
        tmp_sum += a.X[i]*b.X[i]
    return tmp_sum

def matrixmul_mv(A,b):
    """
    Multiply a matrix A with column b

    return: column vector
    """
    c = vector()
    for i in range(8):
        c.X[i] = 0.0
        for j in range(8):
            c.X[i] += A.Y[i][j]*b.X[j]
    return c

def matrixmul_vm(a, B):
    """
    Multiply a row vector a with matrix B

    return: row vector
    """
    c = vector()
    for i in range(8):
        c.X[i] = 0.0
        for j in range(8):
            c.X[i] += a.X[j]*B.Y[j][i]
    return c

def matrixmul_mm(A, B):
    """
    Multiply two matrix A, B

    return: matrix
    """
    c = matrix()
    for i in range(8):
        for j in range(8):
            for k in range(8):
                c.Y[i][j] +=  A.Y[i][k]*B.Y[k][j]
    return c


def belongs(x, a, b):
    """
    Function to check if x belongs to [a,b]
    
    return: bool
    """
    if (x-a)*(x-b) <= 0:
        return True
    else:
        return False

def minmod_2(a, b):
    """
    Minmod function
    """
    if a*b <= 0:
        return 0.0
    else:
        if abs(a) < abs(b):
            return a
        else:
            return b

def minmod_4(a, b, c, d):
    """
    Minmod function
    """
    if a > 0 and b > 0 and c > 0 and d > 0:
        mn = a
        if b < mn:
            mn = b
        if c < mn:
            mn = c
        if d < mn:
            mn = d
        return mn;
    elif a < 0 and b < 0 and c < 0 and d < 0:
        mx = a
        if b > mx:
            mx = b
        if c > mx:
            mx = c
        if d > mx: 
            mx = d
        return mx
    else:
        return 0.0

def minmod_6(a, b, c, d, e, f):
    """
    Minmod function
    """
    if a > 0 and b > 0 and c > 0 and d > 0 and e > 0 and f > 0:
        mn = a
        if b < mn:
            mn = b
        if c < mn:
            mn = c
        if d < mn:
            mn = d
        if e < mn:
            mn = e
        if f < mn:
            mn = f
        return mn
    elif a < 0 and b < 0 and c < 0 and d < 0 and e < 0 and f < 0:
        mx = a
        if b > mx:
            mx = b
        if c > mx:
            mx = c
        if d > mx:
            mx = d
        if e > mx: 
            mx = e
        if f > mx:
            mx = f
        return mx
    else:
        return 0.0

def median(x, y, z):
    """
    Median function as given in balasara and shu
    """
    return x + minmod_2(y-x, z-x)

def reflect(UU, face, gamma):
    """
    Function to find the reflection of state vector UU across
    a given face (Solid wall boundary condition states that normal
    component of the velocity to the wall should be zero -> thus the 
    ghost cell on the other side should be a mirror image of the cell inside)
    pressure, density, magnetic field, tangential velocity remain same
    The normal velocity is reflected.

    input: UU vector, face line, gamma
    output: VV vector 
    """
    
    #Extracting the flow variables from the state vector
    rho1 = UU.X[0]
    u1 = UU.X[1] / rho1
    v1 = UU.X[2] / rho1
    w1 = UU.X[3] / rho1
    Bx1 = UU.X[4]
    By1 = UU.X[5]
    Bz1 = UU.X[6]
    E1 = UU.X[7]
    p1 = (GAMMA - 1.0) * (E1 - rho1*(u1*u1 + v1*v1 + w1*w1)/2 - (Bx1*Bx1 + By1*By1 + Bz1*Bz1)/2)

    #Finding the unit normal to the given face
    normal = point()
    if abs(face.a.x - face.b.x) < 1e-15:
        normal.x = 1
        normal.y = 0
    elif abs(face.a.y == face.b.y) < 1e-15:
        normal.x = 0
        normal.y = 1
    else:
        p = point()
        p.x = face.a.x + mod(face.a.x - face.b.x)/2
        p.y = face.a.y - (p.x - face.a.x) * (face.b.x - face.a.x) / (face.b.y - face.a.y)
        normal.x = p.x - face.a.x
        normal.y = p.y - face.a.y
        magn = sqrt(normal.x*normal.x + normal.y*normal.y)
        normal.x /= magn
        normal.y /= magn
	
    #Unit Normal found.

    #Finding the normal and tangential components of velocity
    un1 = u1 * normal.x + v1 * normal.y
    ut1 = v1 * normal.x - u1 * normal.y


    un2 = - un1	#Reflecting the normal component
    ut2 = ut1	#Keeping the tangential component same

    #Casting the reflected velocity in cartesian coordinates
    u2 = un2*normal.x - ut2*normal.y
    v2 = ut2*normal.x + un2*normal.y
	
    #All other terms are same...
    rho2 = rho1
    w2 = w1
    Bx2 = Bx1
    By2 = By1
    Bz2 = Bz1
    p2 = p1
    E2 = p2/(GAMMA-1) + rho2 * (u2*u2 + v2*v2 + w2*w2) / 2 + (Bx2*Bx2 + By2*By2 + Bz2*Bz2) / 2

    #Calculating the reflected state
    VV = vector()
    VV.X[0] = rho2
    VV.X[1] = rho2 * u2
    VV.X[2] = rho2 * v2
    VV.X[3] = rho2 * w2
    VV.X[4] = Bx2
    VV.X[5] = By2
    VV.X[6] = Bz2
    VV.X[7] = E2
    return VV
    

