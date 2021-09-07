#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 23:44:56 2019

@author: elifyilmaz
"""

from gurobipy import *
import numpy as np
import pandas as pd
import itertools as it

#start = time.time()
################### for different sizes just change this part rest is the same ###################

#n=7 # of nodes
#
#ig=[1,2,5,6] #generator nodes
#not_ig=[] #non_generator nodes
#for i in range(n):
#    if i+1 not in ig:
#        not_ig.append(i+1)
#
#slack_bus=[3]
#
#input_file = 'Medium_2.xlsx'
#output_file="pp_swpn_nons_ref4_small.csv"
#dtt=pd.read_excel('Medium_2.xlsx', header=None, sheet_name='CollusiveStrategies')

#########################################################################################

def pp_swpn_nons_ref4(input_file,output_file,dtt,n ,ig,not_ig,slack_bus,s):
    import time
    start = time.time()
    global suspicious
    suspicious=[]
    global dcopf_counter
    dcopf_counter=0

    collusive_set=np.array(dtt)
    
    df = pd.read_excel (input_file,header = None, sheet_name='Cost') #header=0 
    de = pd.read_excel (input_file,header = None, sheet_name='Demand') 
    dw = pd.read_excel (input_file,header = None, sheet_name='Bidset')
    dq = pd.read_excel (input_file,header = None, sheet_name='Pmax')
    dg = pd.read_excel (input_file,header = None, sheet_name='Fmax')
    dr = pd.read_excel (input_file,header = None, sheet_name='Y')
    
    y=np.array(dr)
    fmax1=np.array(dg)
    bid=np.array(dw)  #bidset
    c = list(df[0]) #cost
    d1=list(de[0])  #demand
    pmax1=list(dq[0])  #pmax
    
    Theta=np.zeros(n)
    Flow=np.zeros((n,n))
    
    bidset={}
    for i in ig:
        bidset[i]=[x for x in list(bid[i-1]) if x!=0]
    
    a=max(len(bidset[i]) for i in ig ) #for column number
    
    pmax = [a/100 for a in pmax1]
    d = [a/100 for a in d1]
    fmax = fmax1/100
    
    ee=1 #to count number of alternative solutions
    
    num_col_bidset=bid.shape[1] # of columns
    
    m=Model("Ref2")
    
    lambdaa = m.addVar(vtype=GRB.CONTINUOUS, name="lambdaa")
    theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta")
    LMP=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="LMP")
    Bidd=m.addVars(n,vtype=GRB.INTEGER, name="Bid")
    p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p")
    phi=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="phi")
    zz=m.addVars(n,a,lb=0,vtype=GRB.CONTINUOUS, name="zz")
    psineg=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psineg")
    psipos=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psipos")
    b=m.addVars(n,a,vtype=GRB.BINARY, name="b") #attention to # of columns (ig,k)
    
    #redundants
    #w=m.addVars(n,vtype=GRB.BINARY, name="w")
    #x=m.addVars(n,vtype=GRB.BINARY, name="x")
    #v=m.addVars(n,vtype=GRB.BINARY, name="v")
    #z=m.addVars(n,n,vtype=GRB.BINARY, name="z")
    #u=m.addVars(n,n,vtype=GRB.BINARY, name="u")
    #t=m.addVars(n,n,vtype=GRB.BINARY, name="t")
    
    v1=m.addVars(n,s,vtype=GRB.CONTINUOUS, name="v1")
    v2=m.addVars(n,s,vtype=GRB.CONTINUOUS, name="v2")
    v3=m.addVars(n,n,s,vtype=GRB.CONTINUOUS, name="v3")
    v4=m.addVars(n,n,s,vtype=GRB.CONTINUOUS, name="v4")
    
    
    m.setObjective(lambdaa, GRB.MAXIMIZE)
    for i in range(n):
        m.addSOS(GRB.SOS_TYPE1, [v1[i,j] for j in range(s)]) 
    for i in range(n):
        m.addSOS(GRB.SOS_TYPE1, [v2[i,j] for j in range(s)])
    for i in range(n):
        for j in range(n):
            m.addSOS(GRB.SOS_TYPE1, [v3[i,j,ddd] for ddd in range(s)]) 
    for i in range(n):
        for j in range(n):
            m.addSOS(GRB.SOS_TYPE1, [v4[i,j,ddd] for ddd in range(s)])
    
    #Elif will replace "range(a)" with her range
    m.addConstrs(lambdaa<= -p[i-1]*c[i-1]+ quicksum( bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i]))) for i in ig)
    
    m.addConstrs(quicksum(b[i-1,k] for k in range(len(bidset[i])))== 1 for i in ig)
    
    m.addConstrs(zz[i-1,k]<= pmax[i-1]*b[i-1,k] for i in ig for k in range(len(bidset[i])))
    
    m.addConstrs(zz[i-1,k]<= p[i-1] for i in ig for k in range(len(bidset[i])) )
    
    m.addConstrs(zz[i-1,k]>= p[i-1]-pmax[i-1]*(1-b[i-1,k]) for i in ig for k in range(len(bidset[i])) )
    
    m.addConstrs(Bidd[i-1]== quicksum(bid[i-1][k]*b[i-1,k]for k in range(len(bidset[i]))) for i in ig)
    
    m.addConstrs(p[i]== quicksum(y[i][j]*(theta[i]-theta[j]) for j in range(n))+d[i] for i in range(n))
    
    m.addConstrs(quicksum(y[i][j]*(LMP[j]-LMP[i]) for j in range(n))
    
    +quicksum(y[i][j]*(psineg[i,j]-psipos[i,j])for j in range(n)) 
    
    + quicksum(y[j][i]*(psipos[j,i]-psineg[j,i]) for j in range(n))== 0 for i in range(n))
    
    m.addConstrs(v1[i-1,1] == p[i-1] for i in ig)
    m.addConstrs(v1[i-1,0] == Bidd[i-1]-LMP[i-1]+phi[i-1]  for i in ig)
    
    m.addConstrs(v2[i-1,0] == pmax[i-1] - p[i-1] for i in ig)
    m.addConstrs(v2[i-1,1] == phi[i-1] for i in ig)
    
    m.addConstrs(Bidd[i-1]-LMP[i-1]+phi[i-1]>=0 for i in ig) #new #22.5.20
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
                m.addConstr(y[i][j]*(theta[i]-theta[j]) <= fmax[i][j]) #new #22.5.20
     
    for i in range(n):
        for j in range(n):
            if y[i][j] >0:  
                m.addConstr(y[i][j]*(theta[i]-theta[j])>= -fmax[i][j]  ) #new #22.5.20
    
    m.addConstrs(-p[i-1]+pmax[i-1]>=0 for i in ig) #new #22.5.20
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
               m.addConstr(v3[i,j,0] == fmax[i][j] - y[i][j]*(theta[i]-theta[j])) #new #22.5.20
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
               m.addConstr(v3[i,j,0] == fmax[i][j] - y[i][j]*(theta[i]-theta[j]))
    
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
              m.addConstr(v3[i,j,1] == psipos[i,j])
              
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
              m.addConstr(v4[i,j,0] == fmax[i][j] + y[i][j]*(theta[i]-theta[j]))
    
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
              m.addConstr(v4[i,j,1] == psineg[i,j])             
     
         
    m.addConstrs(theta[i-1]==0 for i in slack_bus )
    m.addConstrs(p[i-1]==0 for i in not_ig)
    
    m.Params.Heuristics=0
    
#    m.optimize()
    
    def DCOPF(Bid): #Bid'i tuple olarak gir 
        
        global dcopf_counter
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
        m.setObjective(quicksum(Bid[i-1]*p[i-1] for i in ig), GRB.MINIMIZE)
        m.optimize()
        
        for i in range(n):
            LMP[i]=MarkClears[i].pi
    #    print(LMP)   
        
        for i in ig:
            power[i-1]= p[i-1].x *100 ##############
        
        for i in ig:
            Profit_dcopf[i-1]=(LMP[i-1]-c[i-1])*power[i-1] #####round 4
            
        Obj=m.objVal
        
        dcopf_counter+=1
        
        return (Profit_dcopf, power, Obj)
    
    
    DCOPF_matrix={}  #(Bids:Profits)
    Nash_set={}
    suspicious=[]
    primer_suspicious=[] #suspicious strategies before nash is found
    
    
    
    def primer_suspicious_check(Bid):
        added_to_primer_suspicious=0
        for ii in ig:
            if DCOPF_matrix[Bid][ii-1]==0:
                return (primer_suspicious,added_to_primer_suspicious)
        primer_suspicious.append(list(Bid)) 
        added_to_primer_suspicious=1           
        return (primer_suspicious,added_to_primer_suspicious)
    
    
        
    def neighbor_search(Bid):
        global suspicious
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
                    if Bid not in DCOPF_matrix:
                        (Profit_dcopf, power, obj)=DCOPF(Bid)
                        DCOPF_matrix[Bid]=[k for k in Profit_dcopf ]                    
                        bidzero=False                   
                        for ii in ig:
                            if DCOPF_matrix[Bid][ii-1]==0:
                                bidzero=True
                                break
                        if bidzero == False:
                            add_cut(Bid)
                            suspicious_check(Bid)
    				
        return(suspicious)
    
    def primer_neighbor_search(Bid):
        global suspicious
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
                    if Bid not in DCOPF_matrix:
                        (Profit_dcopf, power, obj)=DCOPF(Bid)
                        DCOPF_matrix[Bid]=[k for k in Profit_dcopf ]
                        (primer_suspicious,added_to_primer_suspicious)=primer_suspicious_check(Bid)
                        if added_to_primer_suspicious==1:
                            add_cut(Bid)
        return(suspicious)
    
    
    
    def Nash(Bid): #Bid'i tuple olarak gir 
        global primer_suspicious
        bool_nash=0
                   
        if Bid not in DCOPF_matrix:
            (Profit_dcopf, power, obj)=DCOPF(Bid)
            DCOPF_matrix[Bid]=[i for i in Profit_dcopf ]
            (primer_suspicious,added_to_primer_suspicious)=primer_suspicious_check(Bid)
            if added_to_primer_suspicious==1:
                add_cut(Bid)
        Bid_original=tuple(list(Bid).copy())     
        
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
                    if Bid not in DCOPF_matrix:
                        (Profit_dcopf, power, obj)=DCOPF(Bid)
                        DCOPF_matrix[Bid]=[k for k in Profit_dcopf ]
                        (primer_suspicious,added_to_primer_suspicious)=primer_suspicious_check(Bid)
                        if added_to_primer_suspicious==1:
                            add_cut(Bid)          
                            primer_neighbor_search(Bid)
                    if DCOPF_matrix[Bid_original][i-1]+0.01<DCOPF_matrix[Bid][i-1]:
                        return ( Nash_set ,bool_nash, primer_suspicious)
                        
                
        Nash_set[Bid_original]=DCOPF_matrix[Bid_original]
        bool_nash=1
                    
        return ( Nash_set,bool_nash,primer_suspicious )
     
    ### intertools package   
    #bidset=sorted(bidset)           
    combinations = it.product(*(bidset[xx] for xx in bidset  ))
    strategy_combinations=list(combinations)
    #print(list(combinations))       
            
            
    def finding_nash():
        global primer_suspicious
        for  Bid in strategy_combinations:
            Bid=list(Bid)
            if len(Bid)<n:
                for i in range(1,n+1):
                    if i not in ig:
                       Bid.insert(i-1,0)                    
            Bid=tuple(Bid)   
            if Bid not in DCOPF_matrix:
                (Profit_dcopf, power, obj)=DCOPF(Bid)
                DCOPF_matrix[Bid]=[k for k in Profit_dcopf ]
            else: 
                if list(Bid) in primer_suspicious: 
                    continue
            (primer_suspicious,added_to_primer_suspicious)=primer_suspicious_check(Bid)
            if added_to_primer_suspicious==1:
                add_cut(Bid)
                continue
            else:
                (Nash_set,bool_nash,primer_suspicious )=Nash(Bid)
                if bool_nash==1:
                    add_cut(Bid) ###adding nash as a cut                
                    return(Nash_set)
            
        return (Nash_set)            
            
    def add_cut(Bid):
        bid_index=np.zeros(n)
        for i in ig:
            for j in range(len(bidset[i])):
                if Bid[i-1]==bidset[i][j] : 
                    bid_index[i-1]=j
                    break
        m.addConstr(quicksum(b[i-1,bid_index[i-1]]for i in ig)<=len(ig)-1) 
        return    
    def suspicious_check(Bid): #according to Nash
        global suspicious
        added_to_suspicious=0
        for ii in ig:
            if DCOPF_matrix[tuple(Bid)][ii-1]==0:
                return (suspicious,added_to_suspicious)
        
        for k in Nash_set:
            for j in ig:
                if Nash_set[k][j-1]+0.01>DCOPF_matrix[tuple(Bid)][j-1]: 
                    return (suspicious,added_to_suspicious)
            added_to_suspicious=1
            suspicious.append(list(Bid))
            
        return (suspicious,added_to_suspicious)
                
    finding_nash()
    
    
    m.optimize()
    
    Profit_dcopf=np.zeros(n)
        
    Bid=[0]*n
    for i in range(n):
        Bid[i]=round(Bidd[i].x) 
    for i in ig:
        Profit_dcopf[i-1]=(LMP[i-1].x-c[i-1])*p[i-1].x*100 
        
        
#    (Profit_dcopf, power, Obj)=DCOPF(tuple(Bid))
    DCOPF_matrix[tuple(Bid)]=[k for k in Profit_dcopf ]
          
    if lambdaa.x>0.00001:    
        (suspicious,added_to_suspicious)=suspicious_check(tuple(Bid))
    neighbor_search(Bid)
            
    while lambdaa.x>0.00001:
        add_cut(Bid)
        m.optimize()
        #to escape from infeasiblity
        status=m.status
        if status == 3 or status == 4:
                break
        if lambdaa.x>0.00001:
            ee+=1
            Bid=[0]*n
            Profit_dcopf=np.zeros(n)

            for i in range(n):
                Bid[i]=round(Bidd[i].x)      
#            (Profit_dcopf, power, Obj)=DCOPF(tuple(Bid))
            for i in ig:
                Profit_dcopf[i-1]=(LMP[i-1].x-c[i-1])*p[i-1].x*100 
            
            DCOPF_matrix[tuple(Bid)]=[k for k in Profit_dcopf ]
            (suspicious,added_to_suspicious)=suspicious_check(Bid)            
                
            neighbor_search(Bid)
        else:
            break
    
            
    for i in range(len(primer_suspicious)):
        (suspicious,added_to_suspicious)=suspicious_check(primer_suspicious[i])    
        
    
    end=time.time()        
    time=end-start
            
    num_collusive=0 #to count number of found collusives
    
    for Bid in suspicious:
        if Bid in collusive_set.tolist():
            num_collusive+=1
    
    with open(output_file, "a") as f:
        f.write(" "+str(ee)+" ") #number of solved bilevel
        f.write(" "+str(len(suspicious))+" ")
        f.write(" "+str(len(collusive_set) / len(suspicious))+" ")
        f.write(" "+str(num_collusive / len(collusive_set))+" ")
        f.write(" "+str(round(time,2))+" ")
        f.write(" "+str(dcopf_counter)+" ") 
        f.write(" "+str(len(DCOPF_matrix)-ee)+" ") 
        f.write("\n")

    return


#pp_swpn_nons_ref4(input_file,output_file,dtt,n ,ig,not_ig,slack_bus,s)






