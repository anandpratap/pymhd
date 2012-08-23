from utils import *
#SMALL should be defined
def average(U, V, gamma):
    """
    Calculate average of the two states

    return <vector>
    """
    rho1 = U.X[0]
    u1 = U.X[1] / rho1
    v1 = U.X[2] / rho1
    w1 = U.X[3] / rho1
    Bx1 = U.X[4]
    By1 = U.X[5]
    Bz1 = U.X[6]
    E1 = U.X[7]
    p1 = (gamma-1.0) * (E1 - rho1*(u1*u1 + v1*v1 + w1*w1)/2.0 - (Bx1*Bx1 + By1*By1 + Bz1*Bz1)/2.0)
    p_full1 = p1 + (Bx1*Bx1 + By1*By1 + Bz1*Bz1) / 2.0
	
    #Extracting flow variables of the second state
    rho2 = V.X[0]
    u2 = V.X[1] / rho2
    v2 = V.X[2] / rho2
    w2 = V.X[3] / rho2
    Bx2 = V.X[4]
    By2 = V.X[5]
    Bz2 = V.X[6]
    E2 = V.X[7]
    p2 = (gamma-1.0) * (E2 - rho2*(u2*u2 + v2*v2 + w2*w2)/2.0 - (Bx2*Bx2 + By2*By2 + Bz2*Bz2)/2.0)
    p_full2 = p2 + (Bx2*Bx2 + By2*By2 + Bz2*Bz2) / 2.0
    
    #Finding the averaged state by arithmetic averaging
    rho = (rho1 + rho2) / 2
    u = (u1 + u2) / 2
    v = (v1 + v2) / 2
    w = (w1 + w2) / 2
    Bx = (Bx1 + Bx2) / 2
    By = (By1 + By2) / 2
    Bz = (Bz1 + Bz2) / 2
    p_full = (p_full1 + p_full2) / 2
    p = p_full - (Bx*Bx + By*By + Bz*Bz) / 2
    E = p / (gamma-1.0) + rho * (u*u + v*v + w*w) / 2 + (Bx*Bx + By*By + Bz*Bz) / 2
    
    U_avg = vector()
    U_avg.X[0] = rho
    U_avg.X[1] = rho * u
    U_avg.X[2] = rho * v
    U_avg.X[3] = rho * w
    U_avg.X[4] = Bx
    U_avg.X[5] = By
    U_avg.X[6] = Bz
    U_avg.X[7] = E
    
    return U_avg

#The eigenvalues and eigenvectors have been derived for a general cell 
#face, not necessarily normal to the x- or y- axis.They have been der-
#-ived from the eigenvectors given in Powell,Roe's paper which are def-
#-ined for a cell face normal to x-axis. The eigenvectors given there 
#are in primitive variable form and have been transformed to the cons-
#-erved variable form using the transformations given in that paper

#In the following functions, "un" is the velocity component normal to
#the cell face and "ut" is the velocity component tangential to it.

def eigenvalue(k, U, gamma, normal):
    """
    Function which returns the k-th eigenvalue for Jacobian at state U .
    Ordering of eigenvalues is un-cf, un-ca, un-cs, un, un, un+cs, un+ca, un+cf
    """
    
    rho = U.X[0]
    u = U.X[1] / rho
    v = U.X[2] / rho
    w = U.X[3] / rho
    Bx = U.X[4]
    By = U.X[5]
    Bz = U.X[6]
    E = U.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    
    #Finding normal and tangential components of velocity and magnetic fields
    
    Bn = Bx * normal.x + By * normal.y
    Bt = By * normal.x - Bx * normal.y
    un = u * normal.x + v * normal.y
    ut = v * normal.x - u * normal.y
    
    A = (gamma*p + (Bx*Bx + By*By + Bz*Bz)) / rho
    a = sqrt (gamma * p / rho) #speed of sound
    
    bn = Bn / sqrt(rho)
    bt = Bt / sqrt(rho)
    bz = Bz / sqrt(rho)
	
    ca = bn 
    
    #Normalization procedure as outlined in Roe & Balsara's paper:-
    
    if ((((bn*bn) / (a*a)) < SMALL) and (((bt*bt + bz*bz) / (a*a)) > SMALL)):
        cf = sqrt(a*a + (bt*bt + bz*bz))
        cs = sqrt((a*a * bn*bn) / (a*a + (bt*bt + bz*bz)))
    elif ((((bn*bn) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = 0
    elif ((mod(bn) < (a - SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = bn
    elif ((mod(bn) > (a + SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = bn
        cs = a
    elif (((mod(bn*bn - a*a) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = a
    elif ((((a*a) / (bn*bn)) <= SMALL) and (((a*a) / (bt*bt + bz*bz)) <= SMALL)):
        cf = sqrt(bn*bn + bt*bt + bz*bz)
        cs = 0
    else:
        cs = sqrt((A - sqrt(A*A - 4*a*a*bn*bn)) / 2)
        cf = sqrt((A + sqrt(A*A - 4*a*a*bn*bn)) / 2)
	
    if k == 1:
        e_value = (un - cf) #left-going fast wave
    elif k == 2:
        e_value = (un - ca) #left-going Alfven wave
    elif k == 3:
        e_value = (un - cs) #left-going slow wave
    elif k == 4:
        e_value = (un) #divergence wave
    elif k == 5:
        e_value = (un) #entropy wave
    elif k == 6:
        e_value = (un + cs) #right-going slow wave
    elif k == 7:
        e_value = (un + ca) #right-going Alfven wave
    elif k == 8:
        e_value = (un + cf) #right-going fast wave
	

    return e_value

def lefteigenvector(k, U, gamma, normal):
    """
    Calculates left eigenvector
    """
    	
    #extracting flow variables from U
    
    rho = U.X[0]
    u = U.X[1] / rho
    v = U.X[2] / rho
    w = U.X[3] / rho
    Bx = U.X[4]
    By = U.X[5]
    Bz = U.X[6]
    E = U.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    
    #Finding normal and tangential components of velocity and magnetic fields
    
    Bn = Bx * normal.x + By * normal.y
    Bt = By * normal.x - Bx * normal.y
    un = u * normal.x + v * normal.y
    ut = v * normal.x - u * normal.y
    
    A = (gamma*p + (Bx*Bx + By*By + Bz*Bz)) / rho
    a = sqrt (gamma * p / rho) #speed of sound
    
    bn = Bn / sqrt(rho)
    bt = Bt / sqrt(rho)
    bz = Bz / sqrt(rho)
    
	#Normalization procedure as outlined in Roe & Balsara's paper:-
	
    if (((bt*bt+bz*bz)/(a*a)) < SMALL):
        beta_t = 1.0 / sqrt(2.0)
        beta_z = 1.0 / sqrt(2.0)
    else:
        beta_t = bt / sqrt(bt*bt + bz*bz)
        beta_z = bz / sqrt(bt*bt + bz*bz)
	
	
    if ((((bn*bn) / (a*a)) < SMALL) and (((bt*bt + bz*bz) / (a*a)) > SMALL)):
        cf = sqrt(a*a + (bt*bt + bz*bz))
        cs = sqrt((a*a * bn*bn) / (a*a + (bt*bt + bz*bz)))
        alpha_f = sqrt ((a*a) / (a*a + (bt*bt + bz*bz)))
        alpha_s = sqrt ((bt*bt + bz*bz) / (a*a + (bt*bt+bz*bz)))
    elif ((((bn*bn) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = 0
        alpha_f = 1
        alpha_s = 0
    elif ((mod(bn) < (a - SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = bn
        alpha_f = 1
        alpha_s = 0
    elif ((mod(bn) > (a + SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = bn
        cs = a
        alpha_f = 0
        alpha_s = 1
    elif (((mod(bn*bn - a*a) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = a
        
        if (mod((mod(bn) - a)/a) < SMALL):
            phi = 3.1416 / 2.0
        else:
            phi = atan((sqrt(bt*bt + bz*bz)) / (mod(bn) - a))
            
        alpha_f = sin (phi/2.0)
        alpha_s = cos (phi/2.0)
    elif ((((a*a) / (bn*bn)) <= SMALL) and (((a*a) / (bt*bt + bz*bz)) <= SMALL)):
        cf = sqrt(bn*bn + bt*bt + bz*bz)
        cs = 0
        alpha_f = 0
        alpha_s = 1
    else:
        cs = sqrt((A - sqrt(A*A - 4*a*a*bn*bn)) / 2)
        cf = sqrt((A + sqrt(A*A - 4*a*a*bn*bn)) / 2)
        alpha_f = sqrt((a*a - cs*cs) / (cf*cf - cs*cs))
        alpha_s = sqrt((cf*cf - a*a) / (cf*cf - cs*cs))
	

    L = vector() #left-eigenvectors in conserved variable form

    if k == 1:
		#u-cf
        L.X[0] = (alpha_f*cf*un/rho - alpha_s*cs*beta_t*sgn(bn)*ut/rho - alpha_s*cs*beta_z*sgn(bn)*w/rho + alpha_f*(gamma-1.0)*(un*un+ut*ut+w*w)/(2*rho)) / (2*a*a)
        L.X[1] = (-alpha_f*cf/rho + alpha_f*(1.0-gamma)*un/rho) / (2*a*a)
        L.X[2] = (alpha_s*cs*beta_t*sgn(bn)/rho + alpha_f*(1.0-gamma)*ut/rho) / (2*a*a)
        L.X[3] = (alpha_s*cs*beta_z*sgn(bn)/rho + alpha_f*(1.0-gamma)*w/rho) / (2*a*a)
        L.X[4] = (1.0-gamma) * Bn  * alpha_f / (2*rho*a*a)
        L.X[5] = (alpha_s*a*beta_t/sqrt(rho) + alpha_f*(1.0-gamma)*Bt/rho) / (2*a*a)
        L.X[6] = (alpha_s*a*beta_z/sqrt(rho) + alpha_f*(1.0-gamma)*Bz/rho) / (2*a*a)
        L.X[7] = (alpha_f*(gamma-1.0)/rho) / (2*a*a)
    elif k == 2:
		#u-ca
        L.X[0] = (beta_z*ut/rho - beta_t*w/rho) / sqrt(2)
        L.X[1] = 0
        L.X[2] = (-beta_z/rho) / sqrt(2)
        L.X[3] = (beta_t/rho) / sqrt(2)
        L.X[4] = 0
        L.X[5] = -beta_z/ sqrt(2*rho)
        L.X[6] = beta_t / sqrt(2*rho)
        L.X[7] = 0
    elif k == 3:
		#u-cs
        L.X[0] = (alpha_s*cs*un/rho + alpha_f*cf*beta_t*sgn(bn)*ut/rho + alpha_f*cf*beta_z*sgn(bn)*w/rho + alpha_s*(gamma-1.0)*(un*un+ut*ut+w*w)/(2*rho)) / (2*a*a)
        L.X[1] = (-alpha_s*cs/rho + alpha_s*(1.0-gamma)*un/rho) / (2*a*a)
        L.X[2] = (-alpha_f*cf*beta_t*sgn(bn)/rho + alpha_s*(1.0-gamma)*ut/rho) / (2*a*a)
        L.X[3] = (-alpha_f*cf*beta_z*sgn(bn)/rho + alpha_s*(1.0-gamma)*w/rho) / (2*a*a)
        L.X[4] = (1.0-gamma) * Bn  * alpha_s / (2*rho*a*a)
        L.X[5] = (-alpha_f*a*beta_t/sqrt(rho) + alpha_s*(1.0-gamma)*Bt/rho) / (2*a*a)
        L.X[6] = (-alpha_f*a*beta_z/sqrt(rho) + alpha_s*(1.0-gamma)*Bz/rho) / (2*a*a)
        L.X[7] = (alpha_s*(gamma-1.0)/rho) / (2*a*a)
    elif k == 4:
		#u
        L.X[0] = 0
        L.X[1] = 0
        L.X[2] = 0
        L.X[3] = 0
        L.X[4] = 1
        L.X[5] = 0
        L.X[6] = 0
        L.X[7] = 0
    elif k == 5:
		#u
        L.X[0] = 1 - (gamma-1.0)*(un*un+ut*ut+w*w)/(2*a*a)
        L.X[1] = (gamma-1.0)*un / (a*a)
        L.X[2] = (gamma-1.0)*ut / (a*a)
        L.X[3] = (gamma-1.0)*w / (a*a)
        L.X[4] = (gamma-1.0)*Bn / (a*a)
        L.X[5] = (gamma-1.0)*Bt / (a*a)
        L.X[6] = (gamma-1.0)*Bz / (a*a)
        L.X[7] = (1.0-gamma) / (a*a)
    elif k == 6:
		#u+cs
        L.X[0] = (-alpha_s*cs*un/rho - alpha_f*cf*beta_t*sgn(bn)*ut/rho - alpha_f*cf*beta_z*sgn(bn)*w/rho + alpha_s*(gamma-1.0)*(un*un+ut*ut+w*w)/(2*rho)) / (2*a*a)
        L.X[1] = (alpha_s*cs/rho + alpha_s*(1.0-gamma)*un/rho) / (2*a*a)
        L.X[2] = (alpha_f*cf*beta_t*sgn(bn)/rho + alpha_s*(1.0-gamma)*ut/rho) / (2*a*a)
        L.X[3] = (alpha_f*cf*beta_z*sgn(bn)/rho + alpha_s*(1.0-gamma)*w/rho) / (2*a*a)
        L.X[4] = (1.0-gamma) * Bn  * alpha_s / (2*rho*a*a)
        L.X[5] = (-alpha_f*a*beta_t/sqrt(rho) + alpha_s*(1.0-gamma)*Bt/rho) / (2*a*a)
        L.X[6] = (-alpha_f*a*beta_z/sqrt(rho) + alpha_s*(1.0-gamma)*Bz/rho) / (2*a*a)
        L.X[7] = (alpha_s*(gamma-1.0)/rho) / (2*a*a)
    elif k == 7:
		#u+ca
        L.X[0] = (beta_z*ut/rho - beta_t*w/rho) / sqrt(2)
        L.X[1] = 0
        L.X[2] = -(beta_z/rho) / sqrt(2)
        L.X[3] = (beta_t/rho) / sqrt(2)
        L.X[4] = 0
        L.X[5] = beta_z / sqrt(2*rho)
        L.X[6] = - beta_t / sqrt(2*rho)
        L.X[7] = 0
    elif k == 8:
		#u+cf
        L.X[0] = (-alpha_f*cf*un/rho + alpha_s*cs*beta_t*sgn(bn)*ut/rho + alpha_s*cs*beta_z*sgn(bn)*w/rho + alpha_f*(gamma-1.0)*(un*un+ut*ut+w*w)/(2*rho)) / (2*a*a)
        L.X[1] = (alpha_f*cf/rho + alpha_f*(1.0-gamma)*un/rho) / (2*a*a)
        L.X[2] = (-alpha_s*cs*beta_t*sgn(bn)/rho + alpha_f*(1.0-gamma)*ut/rho) / (2*a*a)
        L.X[3] = (-alpha_s*cs*beta_z*sgn(bn)/rho + alpha_f*(1.0-gamma)*w/rho) / (2*a*a)
        L.X[4] = (1.0-gamma) * Bn  * alpha_f / (2*rho*a*a)
        L.X[5] = (alpha_s*a*beta_t/sqrt(rho) + alpha_f*(1.0-gamma)*Bt/rho) / (2*a*a)
        L.X[6] = (alpha_s*a*beta_z/sqrt(rho) + alpha_f*(1.0-gamma)*Bz/rho) / (2*a*a)
        L.X[7] = (alpha_f*(gamma-1.0)/rho) / (2*a*a)
	

    LL = vector()
    LL.X = L.X
    LL.X[1] = L.X[1]*normal.x - L.X[2]*normal.y
    LL.X[2] = L.X[1]*normal.y + L.X[2]*normal.x
    LL.X[4] = L.X[4]*normal.x - L.X[5]*normal.y
    LL.X[5] = L.X[4]*normal.y + L.X[5]*normal.x
    return LL
	

def righteigenvector(k, U, gamma, normal):
    """
    Calculates left eigenvector
    """
    
    
    rho = U.X[0]
    u = U.X[1] / rho
    v = U.X[2] / rho
    w = U.X[3] / rho
    Bx = U.X[4]
    By = U.X[5]
    Bz = U.X[6]
    E = U.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    
    #Finding normal and tangential components of velocity and magnetic fields
    
    Bn = Bx * normal.x + By * normal.y
    Bt = By * normal.x - Bx * normal.y
    un = u * normal.x + v * normal.y
    ut = v * normal.x - u * normal.y
    
    A = (gamma*p + (Bx*Bx + By*By + Bz*Bz)) / rho
    a = sqrt (gamma * p / rho) #speed of sound

    bn = Bn / sqrt(rho)
    bt = Bt / sqrt(rho)
    bz = Bz / sqrt(rho)
    
    #Normalization procedure as outlined in Roe & Balsara's paper:-

    if (((bt*bt+bz*bz)/(a*a)) < SMALL):
        beta_t = 1.0 / sqrt(2.0)
        beta_z = 1.0 / sqrt(2.0)
    else:
        beta_t = bt / sqrt(bt*bt + bz*bz)
        beta_z = bz / sqrt(bt*bt + bz*bz)
	
	
    if ((((bn*bn) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) > SMALL)):
        cf = sqrt(a*a + (bt*bt + bz*bz))
        cs = sqrt((a*a * bn*bn) / (a*a + (bt*bt + bz*bz)))
        alpha_f = sqrt ((a*a) / (a*a + (bt*bt + bz*bz)))
        alpha_s = sqrt ((bt*bt + bz*bz) / (a*a + (bt*bt+bz*bz)))
    elif ((((bn*bn) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = 0
        alpha_f = 1
        alpha_s = 0
    elif ((mod(bn) < (a - SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = bn
        alpha_f = 1
        alpha_s = 0
    elif ((mod(bn) > (a + SMALL)) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = bn
        cs = a
        alpha_f = 0
        alpha_s = 1
    elif (((mod(bn*bn - a*a) / (a*a)) <= SMALL) and (((bt*bt + bz*bz) / (a*a)) <= SMALL)):
        cf = a
        cs = a
        
        if (mod((mod(bn) - a)/a) < SMALL):
            phi = 3.1416 / 2.0
        else:
            phi = atan((sqrt(bt*bt + bz*bz)) / (mod(bn) - a))
            
        alpha_f = sin (phi/2.0)
        alpha_s = cos (phi/2.0)
    elif ((((a*a) / (bn*bn)) <= SMALL) and (((a*a) / (bt*bt + bz*bz)) <= SMALL)):
        cf = sqrt(bn*bn + bt*bt + bz*bz)
        cs = 0
        alpha_f = 0
        alpha_s = 1
    else:
        cs = sqrt((A - sqrt(A*A - 4*a*a*bn*bn)) / 2)
        cf = sqrt((A + sqrt(A*A - 4*a*a*bn*bn)) / 2)
        alpha_f = sqrt((a*a - cs*cs) / (cf*cf - cs*cs))
        alpha_s = sqrt((cf*cf - a*a) / (cf*cf - cs*cs))
	
	
    R = vector() #right-eigenvector in conserved variable form

    if k == 1:
		#u-cf
        R.X[0] = alpha_f*rho
        R.X[1] = un*alpha_f*rho - rho*alpha_f*cf
        R.X[2] = ut*alpha_f*rho + rho*alpha_s*cs*beta_t*sgn(bn)
        R.X[3] = w*alpha_f*rho + rho*alpha_s*cs*beta_z*sgn(bn)
        R.X[4] = 0
        R.X[5] = alpha_s*sqrt(rho)*a*beta_t
        R.X[6] = alpha_s*sqrt(rho)*a*beta_z
        R.X[7] = alpha_f*rho*(un*un+ut*ut+w*w)/2 - rho*un*alpha_f*cf + rho*ut*alpha_s*cs*beta_t*sgn(bn) + rho*w*alpha_s*cs*beta_z*sgn(bn) + Bt*alpha_s*sqrt(rho)*a*beta_t + Bz*alpha_s*sqrt(rho)*a*beta_z + alpha_f*rho*a*a/(gamma-1.0)
    elif k == 2:
		#u-ca
        R.X[0] = 0
        R.X[1] = 0
        R.X[2] = -rho * beta_z / sqrt(2.0)
        R.X[3] = rho * beta_t / sqrt(2.0)
        R.X[4] = 0
        R.X[5] = - sqrt(rho/2.0) * beta_z
        R.X[6] = sqrt(rho/2.0) * beta_t
        R.X[7] = (-rho*ut*beta_z + rho*w*beta_t - Bt*sqrt(rho)*beta_z + Bz*sqrt(rho)*beta_t) / sqrt(2.0)
    elif k == 3:
		#u-cs
        R.X[0] = alpha_s*rho
        R.X[1] = un*alpha_s*rho - alpha_s*cs*rho
        R.X[2] = ut*alpha_s*rho - rho*alpha_f*cf*beta_t*sgn(bn)
        R.X[3] = w*alpha_s*rho - rho*alpha_f*cf*beta_z*sgn(bn)
        R.X[4] = 0
        R.X[5] = -alpha_f*sqrt(rho)*a*beta_t
        R.X[6] = -alpha_f*sqrt(rho)*a*beta_z
        R.X[7] = alpha_s*rho*(un*un+ut*ut+w*w)/2 - rho*un*alpha_s*cs - rho*ut*alpha_f*cf*beta_t*sgn(bn) - rho*w*alpha_f*cf*beta_z*sgn(bn) - Bt*alpha_f*sqrt(rho)*a*beta_t - Bz*alpha_f*sqrt(rho)*a*beta_z + alpha_s*rho*a*a/(gamma-1.0)
    elif k == 4:
		#u
        R.X[0] = 0
        R.X[1] = 0
        R.X[2] = 0
        R.X[3] = 0
        R.X[4] = 1
        R.X[5] = 0
        R.X[6] = 0
        R.X[7] = Bn
    elif k == 5:
		#u
        R.X[0] = 1
        R.X[1] = un
        R.X[2] = ut
        R.X[3] = w
        R.X[4] = 0
        R.X[5] = 0
        R.X[6] = 0
        R.X[7] = (un*un+ut*ut+w*w)/2
    elif k == 6:
		#u+cs
        R.X[0] = alpha_s*rho
        R.X[1] = un*alpha_s*rho + alpha_s*cs*rho
        R.X[2] = ut*alpha_s*rho + rho*alpha_f*cf*beta_t*sgn(bn)
        R.X[3] = w*alpha_s*rho + rho*alpha_f*cf*beta_z*sgn(bn)
        R.X[4] = 0
        R.X[5] = -alpha_f*sqrt(rho)*a*beta_t
        R.X[6] = -alpha_f*sqrt(rho)*a*beta_z
        R.X[7] = alpha_s*rho*(un*un+ut*ut+w*w)/2 + rho*un*alpha_s*cs + rho*ut*alpha_f*cf*beta_t*sgn(bn) + rho*w*alpha_f*cf*beta_z*sgn(bn) - Bt*alpha_f*sqrt(rho)*a*beta_t - Bz*alpha_f*sqrt(rho)*a*beta_z + alpha_s*rho*a*a/(gamma-1.0)
    elif k == 7:
		#u+ca
        R.X[0] = 0
        R.X[1] = 0
        R.X[2] = -rho*beta_z/sqrt(2.0)
        R.X[3] = rho*beta_t/sqrt(2.0)
        R.X[4] = 0
        R.X[5] = sqrt(rho/2.0)*beta_z
        R.X[6] = -sqrt(rho/2.0)*beta_t
        R.X[7] = (-rho*ut*beta_z + rho*w*beta_t + Bt*sqrt(rho)*beta_z - Bz*sqrt(rho)*beta_t) / sqrt(2.0)
    elif k == 8:
		#u+cf
        R.X[0] = alpha_f*rho
        R.X[1] = un*alpha_f*rho + rho*alpha_f*cf
        R.X[2] = ut*alpha_f*rho - rho*alpha_s*cs*beta_t*sgn(bn)
        R.X[3] = w*alpha_f*rho - rho*alpha_s*cs*beta_z*sgn(bn)
        R.X[4] = 0
        R.X[5] = alpha_s*sqrt(rho)*a*beta_t
        R.X[6] = alpha_s*sqrt(rho)*a*beta_z
        R.X[7] = alpha_f*rho*(un*un+ut*ut+w*w)/2 + rho*un*alpha_f*cf - rho*ut*alpha_s*cs*beta_t*sgn(bn) - rho*w*alpha_s*cs*beta_z*sgn(bn) + Bt*alpha_s*sqrt(rho)*a*beta_t + Bz*alpha_s*sqrt(rho)*a*beta_z + alpha_f*rho*a*a/(gamma-1.0)
	
        
    RR = vector()
    RR.X = R.X
    RR.X[1] = R.X[1]*normal.x - R.X[2]*normal.y
    RR.X[2] = R.X[1]*normal.y + R.X[2]*normal.x
    RR.X[4] = R.X[4]*normal.x - R.X[5]*normal.y
    RR.X[5] = R.X[4]*normal.y + R.X[5]*normal.x
    
    return RR
    
def f(U, gamma):
    """
    Return flux along x direction
    """
    rho = U.X[0]
    u = U.X[1] / rho
    v = U.X[2] / rho
    w = U.X[3] / rho
    Bx = U.X[4]
    By = U.X[5]
    Bz = U.X[6]
    E = U.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    
    F = vector()
    F.X[0] = rho * u
    F.X[1] = rho*u*u + p + (Bx*Bx + By*By + Bz*Bz)/2 - Bx*Bx
    F.X[2] = rho*u*v - By*Bx
    F.X[3] = rho*u*w - Bz*Bx
    F.X[4] = 0
    F.X[5] = u*By - v*Bx
    F.X[6] = u*Bz - w*Bx
    F.X[7] = (E + p + (Bx*Bx + By*By + Bz*Bz)/2) * u - Bx * (u*Bx + v*By + w*Bz)
    return F

def g(U, gamma):
    """
    Return flux along y direction
    """
    rho = U.X[0]
    u = U.X[1] / rho
    v = U.X[2] / rho
    w = U.X[3] / rho
    Bx = U.X[4]
    By = U.X[5]
    Bz = U.X[6]
    E = U.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    
    G = vector()
    G.X[0] = rho * v
    G.X[1] = rho*u*v - By*Bx
    G.X[2] = rho*v*v + p + (Bx*Bx + By*By + Bz*Bz)/2 - By*By
    G.X[3] = rho*v*w - Bz*By
    G.X[4] = v*Bx - u*By
    G.X[5] = 0
    G.X[6] = v*Bz - w*By
    G.X[7] = (E + p + (Bx*Bx + By*By + Bz*Bz)/2) * v - By * (u*Bx + v*By + w*Bz)
    return G

def flux(U, normal, gamma):
    """
    Function to return the flux along a given normal for a given state U 
    f(U).nx + g(U).ny
    """
    flx = vector()
    F = vector()
    G = vector()
    
    F = f(U, gamma)
    G = g(U, gamma)
    for k in range(8):
        flx.X[k] = normal.x * F.X[k] + normal.y * G.X[k];
	
    return flx

def maxwavespeed(U, gamma):
    """
    Function to return the maximum wavespeed possible for a given state U
    """
    rho = U.u.X[0]
    u = U.u.X[1] / rho
    v = U.u.X[2] / rho
    w = U.u.X[3] / rho
    Bx = U.u.X[4]
    By = U.u.X[5]
    Bz = U.u.X[6]
    E = U.u.X[7]
    p = (gamma - 1.0) * (E - rho*(u*u + v*v + w*w)/2 - (Bx*Bx + By*By + Bz*Bz)/2)
    max_speed = sqrt(u*u + v*v + w*w) + sqrt((gamma*p + (Bx*Bx + By*By + Bz*Bz)) / rho)
    return max_speed


