# -*- coding: utf-8 -*-
"""
Truth Model Normal Force

Created on Tue Jan 26 14:42:01 2016

@author: Josh Lo
"""

from sympy import *

""" Define Symbols """

# dimensions
a,b= symbols('a b')

# kinematics
x_dd,y_dd= symbols('x_dd y_dd')
phi= symbols('phi')

# gravity
g= symbols('g')

# sprung and unsprung mass
m_s, m_uf, m_ur= symbols('m_s m_uf m_ur')
m= symbols('m')  # total mass

# sprung mass height
h_s,h_f,h_r= symbols('h_s h_f h_r') 

# track width
t_f,t_r= symbols('t_f t_r')

# roll stiffness front and rear
k_f, k_r ,k_roll= symbols('k_f k_r k_roll')



""" 1. Load Trasfer Due to Roll """
F_rf= k_f*h_s*m_s/k_roll/t_f*(y_dd*cos(phi) + g*sin(phi))
F_rr= k_r*h_s*m_s/k_roll/t_r*(y_dd*cos(phi) + g*sin(phi))

""" 2. Load Transfer due to Roll center height"""
F_lf= m_s*b*h_f*y_dd/t_f/(a+b)
F_lr= m_s*a*h_r*y_dd/t_r/(a+b)

""" 3. Load Transfer due to Unsprung mass"""
F_uf= m_uf*y_dd*h_f/t_f
F_ur= m_ur*y_dd*h_r/t_r

""" 4. Load Transfer due to Longitudinal weight"""
# cg total height
h_cg= (h_f + h_r)/2 + h_s

F_lgf= (m_uf*h_f + m_s*h_cg + m_ur*h_r)*x_dd/(a+b)
F_lgr= -(m_uf*h_f + m_s*h_cg + m_ur*h_r)*x_dd/(a+b)

""" 5. Load Transfer due to Static Weight """
F_nsf= b/(a+b)*m*g
F_nsr= a/(a+b)*m*g


""" Normal Forces on Tires 
1. front left
2. front right
3. rear left
4. rear right
"""
F_n1= F_rf + F_lf + F_uf - 0.5*F_lgf + 0.5*F_nsf
F_n2= -F_rf - F_lf - F_uf - 0.5*F_lgf + 0.5*F_nsf
F_n3= F_rr + F_lr + F_ur + 0.5*F_lgf +0.5*F_nsr
F_n4= -F_rr -F_lr - F_ur + 0.5*F_lgf +0.5*F_nsr


