#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:20:26 2022

@author: sharif-al-mahmud
"""

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
from copy import deepcopy
rng = np.random.default_rng()
import random


def Update_Res(R,GG,PP):
    for s in range(len(GG)): #Updating Residual within the while loop
        R=R-np.dot(GG[s],PP[s])
    return R

def Length(g_i_j,Y_b,U,e):
    l_i = e+ np.dot(g_i_j,Y_b)/U
    return l_i


# ## H1 Function Modified from Tsao et al. 2020  (-suggested by Dr. Weyers) 

def H1(Q,Y,Pm,U,l_max,e):
    ''' 
    P= [P_min,P_max]
    '''
    if len(Q)!=len(Y):
        raise ValueError('number of sizes and number of fabric consumption does not match')
        
    r_b = deepcopy(Q)
    Y_b = Y
    
    G_a_b=[]
    P_a=[]
    P_min,P_max= Pm
    Beta=len(Y_b)
    K = range(Beta)
    while max(r_b)>0:

        #print(f'\nresidual demand:{r_b}')
        t_r_b=[a for a in r_b if a>0]
        p_i = max(P_min,min(min(t_r_b),P_max))

        g_i_j= list (np.zeros(Beta,dtype=int))   # empty list containg Beta values, Beta=Sizes
        
        l_i=0                       #inital length = 0

        K= list(range(Beta))
        random.shuffle(K)  # Random size

        for j in K:
            if r_b[j] >= P_min:

                if l_i <= l_max:
                    temp_g_i_j= g_i_j.copy()
                    temp_g_i_j[j] = math.floor(r_b[j]/p_i)
                    temp_l_i= Length(temp_g_i_j,Y_b,U,e)
                    # print (temp_l_i)

                    if temp_l_i <= l_max:
                        g_i_j=temp_g_i_j
                        l_i= temp_l_i
                    else:
                        g_i_j[j] = math.floor((l_max-l_i)/(Y_b[j]/U))
                        l_i = Length(g_i_j,Y_b,U,e)

            elif r_b[j] > 0:
                if l_i <= l_max:
                    temp_g_i_j= g_i_j.copy()
                    temp_g_i_j[j]= math.ceil(r_b[j]/p_i)
                    #print (temp_g_i_j)
                    temp_l_i= Length(temp_g_i_j,Y_b,U,e)
                    #print (temp_l_i)
                    
                    if temp_l_i <= l_max:
                        g_i_j=temp_g_i_j
                        l_i= temp_l_i
                    else:
                        g_i_j[j]= math.floor((l_max-l_i)/(Y_b[j]/U))
                        l_i = Length(g_i_j,Y_b,U,e)

            else:
                pass

        if l_i<l_max and l_i>e:
            G_a_b.append (g_i_j)
            P_a.append(p_i)
            r_b= r_b-np.dot(g_i_j,p_i)

    return dict(G=G_a_b,P=P_a)


# ## H3 Function Modified from M&B 2016  (-suggested by Dr. Weyers) 

def H3 (Q,Y,Pm,U,l_max,e):

    if len(Q)!=len(Y):
        raise ValueError('number of sizes and number of fabric consumption does not match')
    
    P_min,P_max= Pm # Unpacking
    r = deepcopy(Q)
    Y_b=Y
    Beta=len(Y_b)

    G=[]
    P=[]
    while max(r)>0:
        
        rr=[i for i in r if i >0] 
        p_i = min(P_max,max (np.gcd.reduce(rr),P_min))  #GCD
        
        g= list (np.zeros(Beta,dtype=int))   # empty list containg Beta values, Beta=Sizes, betas
        
        l=0                       #inital length = 0

        K= list(range(Beta))
        random.shuffle(K)  # Random size

        for j in K:
            if r[j] >0:
                g[j]=max(0,min(math.ceil(r[j]/p_i),math.floor((l_max-l)/(Y_b[j]/U))))
                l= Length(g,Y_b,U,e)
        
        G.append (g)
        P.append(p_i)
        r= r-np.dot(g,p_i)
    
    return dict(G=G,P=P)


# ## H5 Function - Using Random sizes ( suggested by Dr. Weyers)

def H5(Q,Y,Pm,U,l_max,e):

    if len(Q)!=len(Y):
        raise ValueError('number of sizes and number of fabric consumption does not match')
    
    P_min,P_max = Pm  #Unpack
    
    i=1
    r = deepcopy(Q)
    Y_b=Y
    Beta= len(Y_b)
    m_max = 2*Beta
    G=[]
    P=[]
    
    while max(r)>0:
        for j in range(Beta):
            if r[j]>0 and r[j]<P_min:
                r[j]=P_min

        p_i = rng.integers(max(P_min,1),min(P_max,max(P_min,max(r))),endpoint=True)
        #print('p_i=',p_i)

        g_i_j= list(np.zeros(Beta,dtype=int))   # empty list containg Beta values, 
        l_i=0                       #inital length = 0

        K= list(range(Beta))
        random.shuffle(K)  # Random size

        for j in K:
            if r[j]<=0 or l_i>= l_max:
                g_i_j[j]=0
                continue
            
            temp_g_i_j= g_i_j.copy()
            temp_g_i_j[j]=rng.integers(0,math.floor(r[j]/p_i)+1)
            temp_l_i= Length(temp_g_i_j,Y_b,U,e)
            
            if temp_l_i <= l_max:
                g_i_j= temp_g_i_j
                l_i= temp_l_i
            else:
                g_i_j[j]= rng.integers(0, math.floor((l_max-l_i)/(Y_b[j]/U))+1)
                l_i= Length(g_i_j,Y_b,U,e)

        r= r-np.dot(g_i_j,p_i)
        G.append (g_i_j)
        P.append(p_i)

        i+=1
        if i >= m_max and max(r)>0:
            #print ('H3')
            h3=H3(Q=r,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)  # Use H3 algorithm to pack all sizes
            g=h3['G']
            p= h3['P']    
            
            G=G+g
            P=P+p
            i=len(P)+1
            break

    return dict(G=G,P=P)

def H3_MB2016(Q,Y,Pm,U,l_max,e):
    
    P_min,P_max = Pm
    r_b = deepcopy(Q)
    Beta=len(r_b)
    G=[]
    P=[]
    while max(r_b)>0:

        p_i =min(P_max,max(np.gcd.reduce([x for x in r_b if x>0]),P_min))

        g_i_j= list (np.zeros(Beta,dtype=int))   # empty list containg Beta values, Beta=Sizes
        l_i=0                       #inital length = 0

        temp_g_i_j= g_i_j.copy()
        for j in range (Beta):       
            temp_g_i_j [j] =  math.ceil(r_b[j]/p_i)
        temp_l_i = Length(temp_g_i_j,Y,U,e)
       

        while temp_l_i >= l_max:    #loop until length < max lenght
            idx= temp_g_i_j.index(max(temp_g_i_j))
            temp_g_i_j[idx] = temp_g_i_j[idx]-1   #reducing number of occurences from the maximum occurences in the section 
    
            temp_l_i = Length(temp_g_i_j,Y,U,e)

        g_i_j= temp_g_i_j
        l_i=temp_l_i

        r_b= r_b-np.dot(g_i_j,p_i)
        G.append (g_i_j)
        P.append(p_i)

    return dict(G=G,P=P)

def length(j,g_i_j,Y_b,U,e):
    l_i = e+np.dot(g_i_j,Y_b[:j+1])/U
    return l_i

def H1_T2020(Q,Y,Pm,U,l_max,e):
    
    P_min,P_max = Pm   #unpack
    Y_b=Y
    r_b = deepcopy(Q)
    Beta=len(r_b)
    
    G=[]
    P=[]
    while max(r_b)>0:
        p_i = max(P_min,min(min([x for x in r_b if x>0]),P_max))

        g_i_j= []
        l_i=0
        for j in range(Beta):
            if r_b[j] >= P_min:
                if l_i <= l_max:
                    temp_g_i_j= g_i_j.copy()
                    temp_g_i_j.append(math.floor(r_b[j]/p_i))

                    temp_l_i= length(j,temp_g_i_j,Y_b,U,e)

                    if temp_l_i <= l_max:
                        g_i_j=temp_g_i_j
                        l_i= temp_l_i
                    else:
                        g_i_j.append (math.floor((l_max-l_i)/(Y_b[j]/U)))
                        l_i = length(j,g_i_j,Y_b,U,e)
                else:
                    g_i_j.append (0)

            elif r_b[j] > 0:
                if l_i <= l_max:
                    temp_g_i_j= g_i_j.copy()
                    temp_g_i_j.append(math.ceil(r_b[j]/p_i))
                    temp_l_i= length(j,temp_g_i_j,Y_b,U,e)

                    if temp_l_i <= l_max:
                        g_i_j=temp_g_i_j
                        l_i= temp_l_i
                    else:
                        g_i_j.append(math.floor((l_max-l_i)/(Y_b[j]/U)))
                        l_i = length(j,g_i_j,Y_b,U,e)
                else:
                    g_i_j.append (0)

            else:
                g_i_j.append (0)

        r_b= r_b-np.dot(g_i_j,p_i)
        G.append (g_i_j)
        P.append(p_i)

    return dict(G=G,P=P)

def main():
    print("main function")
    ##Data
    Q_b = [18,172,214,254,227,187,187,79]   # b=beta
    Y_b = [1.14,1.18,1.215,1.225,1.247,1.26,1.275,1.285] # Fabric yield rate per garment of size 𝛽
    U = 0.85 # Marker Utilization Rate
    l_max= 20 #max marker length
    e= .083   # Fabric end allowance  
    P_min,P_max= 10,400 #min and max number of ply in one section
    
    Sol_1=H1(Q=Q_b,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)
    Sol_2=H3(Q=Q_b,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)
    Sol_3=H5(Q=Q_b,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)
    Sol_4=H3_MB2016(Q=Q_b,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)
    Sol_5=H1_T2020(Q=Q_b,Y=Y_b,Pm=[P_min,P_max],U=U,l_max=l_max,e=e)
    print(Sol_1)
    print(Sol_2)
    print(Sol_3)
    print(Sol_4)
    print(Sol_5)

if __name__ == "__main__":
    main()
