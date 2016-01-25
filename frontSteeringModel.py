# -*- coding: utf-8 -*-
"""
Front Steering Model --WITHOUT-- Roll Dynamics

Assumptions:
1. h_s = 0
2. no braking ==> delta_b= 0
3. longitudinal slip =0
4. moving at constant speed x_d= const
5. front steering ==> rear steering angle =0
6. small angle approximation

State variables:
y_d, theta, theta_d, y_global

Created on Mon Jan 25 10:06:15 2016

@author: Josh Lo
"""

# run the full state truth model
import truthModel
from sympy import *

""" slip angle definition """
alpha_f,alpha_r= symbols('alpha_f alpha_r')    # slip angle
delta_f= symbols('delta_f')                    # front steering angle 

alpha_f= delta_f - (y_d + a*theta_d)/x_d
alpha_r= (b*theta_d - y_d)/x_d

""" tire lateral force
assume left and right are the same
and small angle
"""

C_f, C_r= symbols('C_f C_r')   # cornering stiffness
L_f, L_r= symbols('L_f L_r')   # tire lateral force

L_f= 2*C_f*alpha_f*1   # small angle approximation (2*C_f*alpha_f*cos(delta_f))
L_r= 2*C_r*alpha_r

""" Forcing terms expressions """

exprF_y= L_f + L_r
exprT_z= a*L_f - b*L_r

exprRhs= Matrix([true_rhs[3],true_rhs[4],true_rhs[5],true_rhs[9]])

# No Roll dynamics
exprRhs= exprRhs.subs([(psi,0),(psi_d,0),(psi_dd,0)])
exprRhs= exprRhs.subs(h_s,0)

# small angle approximation
exprRhs= exprRhs.subs([(sin(theta), theta),(cos(theta),1)])

# sub the forcing terms
exprRhs= exprRhs.subs(F_y, exprF_y)
exprRhs= exprRhs.subs(T_z, exprT_z)

exprRhs= simplify(exprRhs)

""" add steering actuation dynamics 
ANVEL steering input is used :
tau*delta_f_d + delta_f = u_s
"""
tau, u_s = symbols('tau u_s')
exprSteerAct= Matrix([(u_s-delta_f)/tau])

# append equation
exprRhs= Matrix([exprRhs, exprSteerAct])



""" obtain A, B matrix """
sv= Matrix([y_d,theta,theta_d,y_global,delta_f])  # state vector
uv= Matrix([u_s])                     # input vector

Amat= simplify(exprRhs.jacobian(sv))
Bmat= simplify(exprRhs.jacobian(uv))

# make sure no residual terms
simplify(exprRhs-Amat*sv-Bmat*uv)








