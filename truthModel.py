from sympy import *

# --------- define symbol ------------
# forces and mement 
F_x,F_y,T_z= symbols('F_x F_y T_z')
# mass
m_tot, m_s, m_u= symbols('m_tot m_s m_u')

# inertia
I_z,I_x= symbols('I_z I_x')

# pose
x,y,theta= symbols('x y theta')
# 1st derivative
x_d,y_d,theta_d= symbols('x_d y_d theta_d')
# 2nd derivative
x_dd,y_dd,theta_dd= symbols('x_dd y_dd theta_dd)

# roll angle
psi,psi_d,psi_dd= symbols('psi psi_d psi_dd')

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
eqFy= m_tot*(y_dd + x_d*theta_d) + m_s*h_s*psi_dd - F_y
eqTz= I_z*theta_dd - T_z
eqTx= I_x*psi_dd + c_roll*psi_d + k_roll*psi - m_s*g*h_s*psi \
+ m_s*(y_dd+x_d*theta_d)*h_s

