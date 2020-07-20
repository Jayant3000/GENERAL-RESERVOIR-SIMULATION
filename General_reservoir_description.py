# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:39:24 2020

@author: Dell
"""


import numpy as np
import matplotlib.pyplot as plt 
import math
import pandas as pd

print(   "\t ########################################################################### \t"  )
print(   "\t #############                General solution                           ######## \t"  )
print(   "\t #############    Reservoir Simulation implicit solution to 1D            ############### \t"  )
print(   "\t ######## closed boundary or zero flow rate Neumann boundary condition both side   ############ \t"  )
print(   "\t ##################           injection in first block(flow rate id given)  ############### \t"  )
print(   "\t ##################       production in third block  (prerssure given)   ############### \t"  )
print(   "\t ########################################################################### \t"  )


###############################################################################
"""                   importing data from csv file                          """

heterogeneous_reservoir_data = pd.read_csv('./General_reservoir_description.csv')
print(heterogeneous_reservoir_data.count())
print("\n")
#print(heterogeneous_reservoir_data)
k = np.array(heterogeneous_reservoir_data['permeability(md)'])
dx = np.array(heterogeneous_reservoir_data['gridblock_lenght(ft)'])
q = np.array(heterogeneous_reservoir_data['flow_rate()'])

print(k)

print(dx)

print(q)

###############################################################################



L = 10000
print("\n\tThe lenght of the reservoir is ", str(L))

number_nodes = len(k)-2
print("\n\tBlock node in the reservoir is ", str(number_nodes))

P0 = 1000
print("\n\tThe intial pressure of the reservoir is ", str(P0))

P_left = 0
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_left))

P_right = 0
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_right))

porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))


viscosity = 1
print("\n\tthe viscosity value is ", str(viscosity))

area = 200000
print("\n\tCross sectional area of the reservoir ", str(area))

compressibility = 1*10**(-6)
print("\n\tcompressibility of the reservoir is ", str(compressibility))

Bw = 1
print("\n\t water formation volume factor is " +str(Bw)+ " rb/stb" )

print("\n\tthe permeability(md) value of the reservoir is \n", str(k))

print("\n\tthe gridblock_lenght(ft) value of the reservoir is \n", str(dx))

print("\n\tthe flow rate value of the reservoir is \n", str(q))


###############################################################################
###############################################################################
#################### final time for simulation is  ############################

t_final = 3
print("\n\t the reservoir simulation time in days is  " +  str(t_final) + "days")

#################### time step  ###############################################

dt_increment = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt_increment)+ "day")

###############################################################################
###############################################################################
###############################################################################
############            pressure and  boundary condition          #############

pressure_previous = np.ones([number_nodes,1])*P0
print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is\n", str(pressure_previous))


###########################################################################
###############################################################################

print("########################")

def permeability(i,k,dx):
    perm = (( dx[i-1] + dx[i] )/( dx[i-1]/k[i-1] + dx[i]/k[i] ))
    return perm

transmissibility_2 = np.ones([number_nodes+1,1])
print(transmissibility_2)


print("########################")

for i in range(1,number_nodes+2):    
    print("i vaue ",i)
    print(k[i])
    print(k[i-1])
    print(dx[i])
    print(dx[i-1])
    print(permeability(i,k,dx))

    transmissibility_2[i-1] =  ( permeability(i,k,dx)*area )/( Bw*viscosity*((dx[i-1] + dx[i])/2) )

print("\nlenght = ",len(transmissibility_2))
if P_left == 0:
    print("\nhahahahaha")
    transmissibility_2[0] = 0*transmissibility_2[0] 
else:
    transmissibility_2[0] = transmissibility_2[0] 
    print("hola")    


if P_right == 0:
    print("\nhahahahaha")
    transmissibility_2[number_nodes] = 0*transmissibility_2[number_nodes] 
else:
    transmissibility_2[number_nodes] = transmissibility_2[number_nodes] 
    print("hola")    

    
print("\n transibillity matrix = \n",transmissibility_2)





         
###############################################################################
###############################################################################
###############################################################################
print("\n############## transmisibility_matrix ################\n")

transmisibility_matrix = np.zeros([number_nodes , number_nodes]) 
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))

for i in range(1,number_nodes,1):
    print("i value" , i)
    transmisibility_matrix[i][i] = transmissibility_2[i] + transmissibility_2[i+1]  
    transmisibility_matrix[i][i-1] = - transmissibility_2[i]
    transmisibility_matrix[i-1][i] = - transmissibility_2[i]

transmisibility_matrix[0][0] =  transmissibility_2[0] +  transmissibility_2[1] 
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))



###############################################################################
print("\n############## B_matrix ################\n")
    
B_matrix = np.zeros([number_nodes , number_nodes])
#print("\n B_matrix is\n", str(B_matrix))
B = np.ones([number_nodes ,1])
for i in range (0,number_nodes):
#    print(i)
    B[i] = area*dx[i+1]*porosity*compressibility
    B_matrix[i][i] = B[i]
print("\n B_matrix is\n", str(B_matrix))

    

###############################################################################
print("\n############## Q_matrix ################\n")
Q= np.zeros([number_nodes , 1])
Q[0] = transmissibility_2[0]*P_left*6.33*10**(-3)
Q[number_nodes-1] = transmissibility_2[number_nodes-1]*2*P_right*6.33*10**(-3)
#print("Q is kjghhterterg",Q)
      
Q_matrix = np.zeros([number_nodes , 1])
#print("\n Q_matrix is\n", str(Q_matrix))
for i in range (0,number_nodes):
    print(i)
    Q_matrix[i] = Q[i] + q[i+1]
    
print("\n Q_matrix is\n", str(Q_matrix))



###############################################################################
print("\n############## N_plus_1_matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N_plus_1_matrix = ( 6.33*10**(-3)*transmisibility_matrix + inverse_dt*B_matrix )

print("\n N_plus_1_matrix is\n", str(N_plus_1_matrix))

###############################################################################
print("\n############## N__matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N__matrix = (inverse_dt*B_matrix )

print("\n N__matrix is\n", str(N__matrix))

print("\n",pressure_previous)
for k in range(0, t_final, 1):
    
    #print("\n time step value",k)
#    boundary_condition_array[0][0] = 2*neta*P_left
    pressure_previous = np.dot(np.linalg.inv(N_plus_1_matrix), (np.dot(N__matrix , pressure_previous) + Q_matrix))
    print("\n",pressure_previous)

print("\nfinal pressure \n", pressure_previous )

print(N__matrix.shape)
print(pressure_previous.shape)
  
     
      
      
      
      
      
      
      
      
      
      
      
      
      
      