    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 19:55:20 2019

@author: elifyilmaz
"""

from gurobipy import *
import numpy as np
import pandas as pd
import time
import itertools as it

# start = time.time()
################### for different sizes just change this part rest is the same ###################

n=7 # of nodes

ig=[1,2,5,6] #generator nodes
not_ig=[] #non_generator nodes
for i in range(n):
    if i+1 not in ig:
        not_ig.append(i+1)

slack_bus=[3]

input_file = 'Medium_6.xlsx'
output_file="000total_enu_medplus1.csv"
output_file2="000total_enu_nash_medplus1.csv"
output_file3="000total_enu_coll_medplus1.csv"

#########################################################################################

def total_enu(input_file,output_file,n ,ig,not_ig,slack_bus ):

    import time
    start=time.time()
    df = pd.read_excel (input_file,header = None, sheet_name='Cost') #header=0 
    de = pd.read_excel (input_file,header = None, sheet_name='Demand') 
    dw = pd.read_excel (input_file,header = None, sheet_name='Bidset')
    dq = pd.read_excel (input_file,header = None, sheet_name='Pmax')
    dg = pd.read_excel (input_file,header = None, sheet_name='Fmax')
    dr = pd.read_excel (input_file,header = None, sheet_name='Y')
    
    
    y=np.array(dr)
    c = list(df[0]) #cost
    pmax = [a/100 for a in list(dq[0])]
    d = [a/100 for a in list(de[0])]
    fmax = np.array(dg)/100
    Flow=np.zeros((n,n))
    bid=np.array(dw)  #bidset
    Theta=np.zeros(n)
    
    Profit_dcopf=np.zeros(n)
    power=np.zeros(n)
    theta_list=np.zeros(n)
    #    flow=np.zeros((n,n))
    
    
    bidset={}
    for i in ig:
        bidset[i]=[x for x in list(bid[i-1]) if x!=0]
    
    
    num_col_bidset=bid.shape[1] # of columns
    
    Obj={}
    Obj1={}
    Obj2={}
    Obj3={}
    Obj4={}
    DCOPF_matrix={}
    Power_matrix={}
    LMP_matrix={}
    Theta_matrix={}
    Flow_matrix={}
    
    
    def DCOPF(Bid): #Bid'i tuple olarak gir 
        global flow_list
        m=Model('DC-OPF')
        LMP=np.zeros(n)
        Bid=list(Bid) 
        Profit_dcopf=np.zeros(n)
        power=np.zeros(n)
        
        p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p")
        theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta") 
        m.addConstrs(-p[i-1]+pmax[i-1]>=0 for i in ig)
        MarkClears = m.addConstrs(p[i]==quicksum(y[i][j]*(theta[i]-theta[j]) for j in range(n))+d[i] for i in range(n))
        for i in range(n):
            for j in range(n): 
                if y[i][j] >0:
                    m.addConstr(y[i][j]*(theta[i]-theta[j]) <= fmax[i][j])
        
        for i in range(n):
            for j in range(n):
                if y[i][j] >0:  
                    m.addConstr(y[i][j]*(theta[i]-theta[j])>= -fmax[i][j] )
        m.addConstrs(theta[i-1]==0 for i in slack_bus )
        m.addConstrs(p[i-1]==0 for i in not_ig)
    #    m.addConstr(obj=quicksum(Bid[i-1]*p[i-1] for i in ig))
        m.setObjective(quicksum(Bid[i-1]*p[i-1] for i in ig), GRB.MINIMIZE)
        m.optimize()
        
        
        for i in range(n):
            LMP[i]=MarkClears[i].pi
        
        for i in ig:
            power[i-1]= p[i-1].x *100 ##############
        
        for i in ig:
            Profit_dcopf[i-1]=(LMP[i-1]-c[i-1])*power[i-1]
        
        for i in range(n):
            theta_list[i]=theta[i].x
            
        flow_list=[]
        for i in range(n):
            for j in range(n):
                if y[i][j]>0:
                    flow_list.append(y[i][j]*(theta_list[i]-theta_list[j]))
    #                flow[i][j]=y[i][j]*(theta_list[i]-theta_list[j])
                    
            
        Obj[tuple(Bid)]=m.objVal
        min_r=10e+100
        for i in ig:
            min_r=min(min_r,p[i-1].x*100*(LMP[i-1]-c[i-1]))
        min_b=10e+100
        for i in ig:
            min_b=min(min_b,Bid[i-1])
    #    print(min_b)
        max_r=0
        for i in ig:
            max_r=max(max_r,p[i-1].x*100*(LMP[i-1]-c[i-1]))
        Obj1[tuple(Bid)]=min_r
        Obj2[tuple(Bid)]=sum(p[i-1].x*(LMP[i-1]-c[i-1]) for i in ig)
        Obj3[tuple(Bid)]=max_r-min_r
        min_profit=10e+100
        for i in ig:
            min_profit=min(min_profit,p[i-1].x*(Bid[i-1]-c[i-1]) )
        Obj4[tuple(Bid)]=0.3*min_b+0.7*min_profit
    #    print(0.3*min_b+0.7*min_profit)
        DCOPF_matrix[tuple(Bid)]=[k for k in Profit_dcopf ]
        Power_matrix[tuple(Bid)]=[k for k in power]
        LMP_matrix[tuple(Bid)]=[k for k in LMP]
        Theta_matrix[tuple(Bid)]=[k for k in theta_list]
        Flow_matrix[tuple(Bid)]=[k for k in flow_list]
        
    #    with open(output_file, "a") as f:
    #        for i in ig:
    #            f.write(" "+str(Bid[i-1])+" ")
    #        f.write(" "+str(m.objVal)+" ")   
    #        f.write(" "+str(min_r)+" ") 
    #        f.write(" "+str(sum(p[i-1].x*(LMP[i-1]-c[i-1]) for i in ig))+" ") 
    #        f.write(" "+str(max_r-min_r)+" ") 
    #        f.write(" "+str(0.3*min_b+0.7*min_profit)+" ") 
    #        for i in ig:    
    #            f.write(" "+str(Profit_dcopf[i-1])+" ")
    #        for i in ig:
    #            f.write(" "+str(power[i-1])+" ")
    #        for i in range(n):
    #            f.write(" "+str(LMP[i])+" ")
    #        for i in range(n):
    #            f.write(" "+str(theta_list[i])+" ")
    #        for i in range(len(flow_list)):
    #            f.write(" "+str(flow_list[i])+" ")            
    #        f.write("\n")
        
        
        
        return (DCOPF_matrix, Power_matrix, Obj)
    
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
        
    for Bid in strategy_combinations2:
        (DCOPF_matrix, Power_matrix, Obj)=DCOPF(Bid)
    
    Nash_set={}
     
    
    def Nash(Bid): #Bid'i tuple olarak gir 
    #    global DCOPF_matrix
    #    max_r_array=np.zeros(n)
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
                    if DCOPF_matrix[Bid_original][i-1] +0.01 < DCOPF_matrix[Bid][i-1]:
                        return ( Nash_set ,bool_nash)
            
        Nash_set[tuple(Bid_original)]=DCOPF_matrix[tuple(Bid_original)]
        bool_nash=1
                    
        return ( Nash_set,bool_nash)
    
    for Bid in strategy_combinations2:
        ( Nash_set,bool_nash)=Nash(Bid)
    
    
        
        
    r_star=[0]*n
    for Bid in Nash_set:
        for i in ig:
            r_star[i-1]=max(r_star[i-1],Nash_set[Bid][i-1])
    
    
    collusive_set=[]
    
    def collusive(Bid):
        if tuple(Bid) not in Nash_set:
            for i in ig:
                if r_star[i-1]+0.01>DCOPF_matrix[tuple(Bid)][i-1]:
                    return(collusive_set)
    #                continue
    #            else:                
            collusive_set.append(Bid)
        return(collusive_set)
        
    if len(Nash_set)>0:  
    	for Bid in strategy_combinations2:
    		(collusive_set)=collusive(Bid)
        
        
    #Bid=(51,52,0,0,33,29,0)   
    end=time.time()
    time=end-start
    
    with open(output_file,"a") as f:
        f.write(" "+str(len(collusive_set))+" ")
        f.write(" "+str(len(Nash_set))+" ")
        f.write(" "+str(round(time,2))+" ")
        f.write("\n")
    
           
    # for Bid in Nash_set:
    #     with open(output_file2, "a") as f:
    #         for i in ig:
    #             f.write(" "+str(Bid[i-1])+" ")
    #         f.write(" "+str(Obj[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj1[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj2[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj3[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj4[tuple(Bid)])+" ")   
    #         for i in ig:    
    #             f.write(" "+str(DCOPF_matrix[tuple(Bid)][i-1])+" ")
    #         for i in ig:
    #             f.write(" "+str(Power_matrix[tuple(Bid)][i-1])+" ")
    #         for i in range(n):
    #             f.write(" "+str(LMP_matrix[tuple(Bid)][i])+" ")
    #         for i in range(n):
    #             f.write(" "+str(Theta_matrix[tuple(Bid)][i])+" ")
    #         for i in range(len(Flow_matrix[tuple(Bid)])):
    #             f.write(" "+str(Flow_matrix[tuple(Bid)][i])+" ")            
    #         f.write("\n")
        
        
        
    # for Bid in collusive_set:
    #     with open(output_file3, "a") as f:
    #         for i in ig:
    #             f.write(" "+str(Bid[i-1])+" ")
    #         f.write(" "+str(Obj[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj1[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj2[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj3[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj4[tuple(Bid)])+" ")   
    #         for i in ig:    
    #             f.write(" "+str(DCOPF_matrix[tuple(Bid)][i-1])+" ")
    #         for i in ig:
    #             f.write(" "+str(Power_matrix[tuple(Bid)][i-1])+" ")
    #         for i in range(n):
    #             f.write(" "+str(LMP_matrix[tuple(Bid)][i])+" ")
    #         for i in range(n):
    #             f.write(" "+str(Theta_matrix[tuple(Bid)][i])+" ")
    #         for i in range(len(Flow_matrix[tuple(Bid)])):
    #             f.write(" "+str(Flow_matrix[tuple(Bid)][i])+" ")            
    #         f.write("\n")
       
    # for Bid in strategy_combinations2:
    #     with open(output_file, "a") as f:
    #         for i in ig:
    #             f.write(" "+str(Bid[i-1])+" ")
    #         f.write(" "+str(Obj[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj1[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj2[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj3[tuple(Bid)])+" ")   
    #         f.write(" "+str(Obj4[tuple(Bid)])+" ")   
    #         for i in ig:    
    #             f.write(" "+str(DCOPF_matrix[tuple(Bid)][i-1])+" ")
    #         for i in ig:
    #             f.write(" "+str(Power_matrix[tuple(Bid)][i-1])+" ")
    #         for i in range(n):
    #             f.write(" "+str(LMP_matrix[tuple(Bid)][i])+" ")
    #         for i in range(n):
    #             f.write(" "+str(Theta_matrix[tuple(Bid)][i])+" ")
    #         for i in range(len(Flow_matrix[tuple(Bid)])):
    #             f.write(" "+str(Flow_matrix[tuple(Bid)][i])+" ")            
    #         f.write("\n")
        

    return 

total_enu(input_file, output_file, n, ig, not_ig, slack_bus)


