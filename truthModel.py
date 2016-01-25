from sympy import *

# --------- define symbol ------------
# forces and mement 
F_x,F_y,T_z= symbols('F_x F_y T_z')
# mass
m_tot, m_s, m_u= symbols('m_tot m_s m_u')
# gravity
g= symbols('g')

# inertia
I_z,I_x= symbols('I_z I_x')

# pose
x,y,theta= symbols('x y theta')
x_global,y_global= symbols('x_global y_global')

# 1st derivative
x_d,y_d,theta_d= symbols('x_d y_d theta_d')
x_global_d,y_global_d= symbols('x_global_d y_global_d')

# 2nd derivative
x_dd,y_dd,theta_dd= symbols('x_dd y_dd theta_dd')

# roll angle
phi,phi_d,phi_dd= symbols('phi phi_d phi_dd')

# roll damping and stiffness
k_roll, c_roll= symbols('k_roll c_roll') 

# cg height
h_s= symbols('h_s')

# wheel base dimensions
a,b= symbols('a b')

# wheel forces
F_x1,F_x2,F_x3,F_x4= symbols('F_x1 F_x2 F_x3 F_x4')
F_y1,F_y2,F_y3,F_y4= symbols('F_y1 F_y2 F_y3 F_y4')

# -------- Force Equilibruim Equation ----------
eqFx= m_tot*(x_dd - y_d*theta_d) - F_x
eqFy= m_tot*(y_dd + x_d*theta_d) + m_s*h_s*phi_dd - F_y
eqTz= I_z*theta_dd - T_z
eqTx= I_x*phi_dd + c_roll*phi_d + k_roll*phi - m_s*g*h_s*phi \
+ m_s*(y_dd+x_d*theta_d)*h_s

# build matrix
eqs= expand(Matrix([eqFx,eqFy,eqTz,eqTx]))


av= Matrix([x_dd,y_dd,theta_dd,phi_dd])  # accelration vetor
fv= Matrix([F_x,F_y,T_z,0])              # force vector

""" Mass Matrix 
the equation is :
massMat*av + stateMat= fv
"""
massMat= eqs.jacobian(av)
stateMat= eqs- massMat*av + fv

""" right hand side
==> av= massMat**-1 *(-stateMat + fv)
==> rhs= massMat**-1 *(-stateMat + fv)
    ==> rhs1= massMat**-1* (-stateMat)
    ==> rhs2= massMat**-1* fv
"""
rhs1= simplify((massMat**-1)*(-stateMat))
rhs2= simplify((massMat**-1)*fv)

""" Rearrange to full-state truth model
state: x, x_d, y, y_d, theta, theta_d, phi, phi_d, x_global, y_global
"""

true_sv= Matrix([x, x_d, y, y_d, theta, theta_d, phi, phi_d, x_global, y_global])

true_sv_d= Matrix([x_d, rhs1[0]+rhs2[0], y_d, rhs1[1]+rhs2[1], \
theta_d, rhs1[2]+rhs2[2], phi_d, rhs1[3]+rhs2[3], \
x_d*cos(theta)-y_d*sin(theta), -x_d*sin(theta)-y_d*cos(theta)])

true_rhs= simplify(true_sv_d)




