#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 19:55:20 2019

"""

from docplex.mp.model import *
import numpy as np
import pandas as pd
import time
import itertools as it
from random import randint

start = time.time()
################### for different sizes just change this part rest is the same ###################

n=5 # of nodes

ig=[1,2,5] #generator nodes
not_ig=[] #non_generator nodes
for i in range(n):
    if i+1 not in ig:
        not_ig.append(i+1)

slack_bus=[3]

input_file = 'Small_size_case2.xlsx'
output_file="total_enumaration_med_case_1.csv"

#########################################################################################

df = pd.read_excel (input_file,header = None, sheet_name='Cost') #header=0 
de = pd.read_excel (input_file,header = None, sheet_name='Demand') 
dw = pd.read_excel (input_file,header = None, sheet_name='Bidset')
dq = pd.read_excel (input_file,header = None, sheet_name='Pmax')
dg = pd.read_excel (input_file,header = None, sheet_name='Fmax')
dr = pd.read_excel (input_file,header = None, sheet_name='Y')

kkk=0
fmax_set={}
pmax_set={}
Search_Iterations=0
Number=2
def DCOPF(Bid): #Bid'i tuple olarak gir 
        
    m=Model('DC-OPF')
        
    LMP=np.zeros(n)
    Bid=list(Bid) 
    
    p=m.continuous_var_list(n,lb=0, name="p")
    theta=m.continuous_var_list(n,lb=-(1e+20), name="theta") ####lower bound -infinity 
    m.add_constraints(-p[i-1]+pmax[i-1]>=0 for i in ig)
    MarkClears = m.add_constraints(p[i]==m.sum(y[i][j]*(theta[i]-theta[j]) for j in range(n))+d[i] for i in range(n))
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
                m.add_constraint(y[i][j]*(theta[i]-theta[j]) <= fmax[i][j])
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:  
                m.add_constraint(y[i][j]*(theta[i]-theta[j])>= -fmax[i][j] )
    m.add_constraints(theta[i-1]==0 for i in slack_bus )
    m.add_constraints(p[i-1]==0 for i in not_ig)
    obj=m.sum(Bid[i-1]*p[i-1] for i in ig)
    m.minimize(obj)
    m.solve() 
    
    status=m.solve_details.status_code  
    # print(status)
    if status == 3:
        return status
    # if m.get_solve_status==
    else:
        LMP=m.dual_values(MarkClears)
    
        Profit_dcopf=np.zeros(n)
        power=np.zeros(n)
        for i in ig:
            power[i-1]=m.get_var_by_name(str(p[i-1])).solution_value*100
        
        for i in ig:
            Profit_dcopf[i-1]= (LMP[i-1]-c[i-1])*power[i-1]
            
        Obj[tuple(Bid)]=obj.solution_value
        DCOPF_matrix[tuple(Bid)]=[k for k in Profit_dcopf ]
        Power_matrix[tuple(Bid)]=[k for k in power]
        
        return (DCOPF_matrix, Power_matrix, Obj)

def Nash(Bid): #Bid'i tuple olarak gir 
#    global DCOPF_matrix
    max_r_array=np.zeros(n)
    bool_nash=0
#    Neighbor_DCOPF_matrix={}
               
    Bid_original=tuple(list(Bid).copy())     
    
#    if Bid not in Nash_set  :
    bid_index=np.zeros(n)
    for i in ig:
        for j in range(len(bidset[i])):
            if Bid[i-1]==bidset[i][j] : 
                bid_index[i-1]=j
                break
    for i in ig:
        bid_index_new=bid_index.copy()
        for j in range(len(bidset[i])):
            if j!=bid_index[i-1] : 
                bid_index_new[i-1]=int(j)
                Bid=list(Bid)
                for l in ig:
                    Bid[l-1]=bidset[l][int(bid_index_new[l-1])]
                Bid=tuple(Bid)
                if DCOPF_matrix[Bid_original][i-1] < DCOPF_matrix[Bid][i-1]:
                    return (Nash_set,max_r_array ,bool_nash)
        
    Nash_set[tuple(Bid_original)]=DCOPF_matrix[tuple(Bid_original)]
    bool_nash=1
                
    return ( Nash_set,max_r_array,bool_nash)

def collusive(Bid):
    if tuple(Bid) not in Nash_set:
        for i in ig:
            if r_star[i-1]>=DCOPF_matrix[tuple(Bid)][i-1]:
                return(collusive_set)
#                continue
#            else:                
        collusive_set.append(Bid)
    return(collusive_set)

while kkk<=Number:
    Search_Iterations+=1
    dgg=np.zeros((len(dq),len(dq)))
    for i in range(len(dq)):
        for j in range(len(dq)):
            if dg[i][j]>0:
                dgg[i][j]=randint(10, 500)
                dgg[j][i]=dgg[i][j]
    
    dqq=np.zeros(len(dq))
    dq=np.array(dq)
    for i in range(len(dq)):
        if dq[i]>0:
            dqq[i]=randint(10, 500)
            # dqq[i]=dq[i]
    # dqq=[356, 309, 0, 0, 384, 491, 0]
    dqq=np.array(dqq)
    y=np.array(dr)
    c = list(df[0]) #cost
    pmax = dqq/100
    d = [a/100 for a in list(de[0])]
    fmax = np.array(dgg)/100
    # fmax=dgg
    Flow=np.zeros((n,n))
    bid=np.array(dw)  #bidset
    Theta=np.zeros(n)
    
    bidset={}
    for i in ig:
        bidset[i]=[x for x in list(bid[i-1]) if x!=0]
    
    
    num_col_bidset=bid.shape[1] # of columns
    
    Obj={}
    DCOPF_matrix={}
    Power_matrix={}
    
    
    
    #(DCOPF_matrix, Power_matrix, Obj)=DCOPF((21,22,0,0,33,34,0))

    ### intertools package   
    #bidset=sorted(bidset)           
    combinations = it.product(*(bidset[xx] for xx in bidset  ))
    strategy_combinations=list(combinations)
    strategy_combinations2=[]
    
    for  Bid in strategy_combinations: ##adding non-generators 0
        Bid=list(Bid)
        if len(Bid)<n:
            for i in range(1,n+1):
                if i in not_ig:
                   Bid.insert(i-1,0)  
        strategy_combinations2.append(Bid)
    ssss=0
    for Bid in strategy_combinations2:
        if DCOPF(Bid)!= 3:
            ssss=1
            (DCOPF_matrix, Power_matrix, Obj)=DCOPF(Bid)
            
    if ssss==1:
        
        Nash_set={}
         
        
        
        for Bid in strategy_combinations2:
            ( Nash_set,max_r_array,bool_nash)=Nash(Bid)
        
        
        #Bid=(21,22,0,0,33,34,0)
            
        r_star=[0]*n
        for Bid in Nash_set:
            for i in ig:
                r_star[i-1]=max(r_star[i-1],Nash_set[Bid][i-1])
        
        #collect collusive if len(Nash)>0
        collusive_set=[]
        if len(Nash_set)>0:
             
            for Bid in strategy_combinations2:
                (collusive_set)=collusive(Bid)    
            
        #Bid=(51,52,0,0,33,29,0)
        if len(collusive_set)!=0:
            fmax_set[kkk]=fmax
            pmax_set[kkk]=pmax
            kkk+=1
        

end=time.time()
time=end-start
    
with open(output_file,"a") as f:
    for sets in pmax_set:
            f.write(" "+str(pmax_set[kkk])+" ")            
    for i in range(n):
      for j in range(n):
        f.write(" "+str(fmax_set[kkk])+" ")

                