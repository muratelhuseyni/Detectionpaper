#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 23:44:56 2019

@author: elifyilmaz
"""

#binary exp. is used in lambda<=... constraint


from gurobipy import *
import numpy as np
import pandas as pd
import itertools as it

################### for different sizes just change this part rest is the same ###################

#n=5 # of nodes
#
#ig=[1,2,5] #generator nodes
#not_ig=[] #non_generator nodes
#for i in range(n):
#    if i+1 not in ig:
#        not_ig.append(i+1)
#
#slack_bus=[3]
#
#input_file = 'Small_5.xlsx'
#output_file="rs_nashcut2_woutobjs_small.csv"
#dtt=pd.read_excel('Small_7.xlsx', header=None, sheet_name='CollusiveStrategies')

#########################################################################################

def nashcut2_ref2_woutobjs(input_file,output_file,dtt,n ,ig,not_ig,slack_bus):
    import time

    start = time.time()

    df = pd.read_excel (input_file,header = None, sheet_name='Cost') #header=0 
    de = pd.read_excel (input_file,header = None, sheet_name='Demand') 
    dw = pd.read_excel (input_file,header = None, sheet_name='Bidset')
    dq = pd.read_excel (input_file,header = None, sheet_name='Pmax')
    dg = pd.read_excel (input_file,header = None, sheet_name='Fmax')
    dr = pd.read_excel (input_file,header = None, sheet_name='Y')
    
    collusive_set=np.array(dtt)
    y=np.array(dr)
    fmax1=np.array(dg)
    bid=np.array(dw)  #bidset
    c = list(df[0]) #cost
    d1=list(de[0])  #demand
    pmax1=list(dq[0])  #pmax
    
    Theta=np.zeros(n)
    Flow=np.zeros((n,n))
    
    DCOPF_matrix={}#new
    primer_suspicious=[] #new
    #Nash_set={} #new
    
    bidset={}
    for i in ig:
        bidset[i]=[x for x in list(bid[i-1]) if x!=0]
    
    a=len(bid[0])  #for column number
    
    pmax = [a/100 for a in pmax1]
    d = [a/100 for a in d1]
    fmax = fmax1/100
    
    ee=1 #to count number of alternative solutions
    
    num_col_bidset=bid.shape[1] # of columns
    
    m=Model("Ref1")
    
    def add_cut(Bid):
        bid_index=np.zeros(n)
        for i in ig:
            for j in range(len(bidset[i])):
                if Bid[i-1]==bidset[i][j] : 
                    bid_index[i-1]=j
                    break
        m.addConstr(quicksum(b[i-1,bid_index[i-1]]for i in ig)<=len(ig)-1) 
        return 
    
    #def BuildModel(m):
    
    theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta")
    LMP=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="LMP")
    Bidd=m.addVars(n, vtype=GRB.INTEGER, name="Bid") 
    p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p")    
    phi=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="phi")
    b=m.addVars(n,a,vtype=GRB.BINARY, name="b") #attention to # of columns (ig,k)
    psineg=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psineg")
    psipos=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psipos")
    lambdaa = m.addVar(vtype=GRB.CONTINUOUS, name="lambdaa")
    zz=m.addVars(n,a,vtype=GRB.CONTINUOUS, name="zz") ########## attention to a,n  
    
    m.addConstrs(zz[i-1,k]<= pmax[i-1]*b[i-1,k] for i in ig for k in range(len(bidset[i])))
    m.addConstrs(zz[i-1,k]<= p[i-1] for i in ig for k in range(len(bidset[i])) )
    m.addConstrs(zz[i-1,k]>= p[i-1]-pmax[i-1]*(1-b[i-1,k]) for i in ig for k in range(len(bidset[i])) )
    
    m.addConstrs(quicksum(b[i-1,k] for k in range(len(bidset[i])))==1 for i in ig)
    m.addConstrs(Bidd[i-1]==quicksum(bidset[i][k]*b[i-1,k]for k in range(len(bidset[i]))) for i in ig) 
    
    m.addConstrs(Bidd[i-1]-LMP[i-1]+phi[i-1]>=0 for i in ig)
    
    m.addConstrs(quicksum(y[i][j]*(LMP[j]-LMP[i]) for j in range(n))
    +quicksum(y[i][j]*(psineg[i,j]-psipos[i,j])for j in range(n)) 
    + quicksum(y[j][i]*(psipos[j,i]-psineg[j,i]) for j in range(n))==0 for i in range(n))
    
    
    for i in range(n):
        for j in range(n): 
            if y[i][j] >0:
                m.addConstr(y[i][j]*(theta[i]-theta[j]) <= fmax[i][j])
    
    for i in range(n):
        for j in range(n):
            if y[i][j] >0:  
                m.addConstr(y[i][j]*(theta[i]-theta[j])>= -fmax[i][j]  )
    
    
    m.addConstrs(-p[i-1]+pmax[i-1]>=0 for i in ig)
    
    m.addConstrs(p[i]-d[i]==quicksum(y[i][j]*(theta[i]-theta[j]) for j in range(n)) for i in range(n))
    
    #old
    #m.addConstr(quicksum(Bidd[i-1]*p[i-1] for i in ig)==quicksum(d[i]*LMP[i] for i in range(n))- quicksum(pmax[i]*phi[i] for i in range(n)) -quicksum(fmax[i][j]*(psipos[i,j]+psineg[i,j])for i in range(n) for j in range(n)))    
    #new
    m.addConstr(quicksum(quicksum(bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i])) ) for i in ig)==quicksum(d[i]*LMP[i] for i in range(n))- quicksum(pmax[i]*phi[i] for i in range(n)) -quicksum(fmax[i][j]*(psipos[i,j]+psineg[i,j])for i in range(n) for j in range(n))      )
    
    m.addConstrs(theta[i-1]==0 for i in slack_bus )
    m.addConstrs(p[i-1]==0 for i in not_ig)
    
    #    m.params.NonConvex = 2 
    
    def RHSmodels(m,fac,bidindex,obj,Bidd):
        
        theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta"+str(fac)+str(bidindex))
        LMP=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="LMP"+str(fac)+str(bidindex))
        #Bidd=m.addVars(n, vtype=GRB.INTEGER, name="Bid") 
        p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p"+str(fac)+str(bidindex))    
        phi=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="phi"+str(fac)+str(bidindex))
        b=m.addVars(n,a,vtype=GRB.BINARY, name="b"+str(fac)+str(bidindex)) #attention to # of columns (ig,k)
        psineg=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psineg"+str(fac)+str(bidindex))
        psipos=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psipos"+str(fac)+str(bidindex))
       # lambdaa = m.addVar(vtype=GRB.CONTINUOUS, name="lambdaa")
        zz=m.addVars(n,a,vtype=GRB.CONTINUOUS, name="zz"+str(fac)+str(bidindex)) ########## attention to a,n  
        
        #m.setObjective(lambdaa, GRB.MAXIMIZE)
        #m.addConstrs(lambdaa<= -p[i-1]*c[i-1]+ quicksum( bidset[i][k]*zz[i-1,k] for k in range(len(bidset[i]))) for i in ig) #old
        #m.addConstrs(lambdaa<= -p[i-1]*c[i-1]+ p[i-1]* LMP[i-1] for i in ig) #new
        
        m.addConstrs(zz[i-1,k]<= pmax[i-1]*b[i-1,k] for i in ig for k in range(len(bidset[i])))
        m.addConstrs(zz[i-1,k]<= p[i-1] for i in ig for k in range(len(bidset[i])) )
        m.addConstrs(zz[i-1,k]>= p[i-1]-pmax[i-1]*(1-b[i-1,k]) for i in ig for k in range(len(bidset[i])) )
        
        m.addConstrs(quicksum(b[i-1,k] for k in range(len(bidset[i])))==1 for i in ig)
       # m.addConstrs(Bidd[i-1]==quicksum(bidset[i][k]*b[i-1,k]for k in range(len(bidset[i]))) for i in ig) 
        
        #m.addConstrs(Bidd[i-1]-LMP[i-1]+phi[i-1]>=0 for i in ig)
        
        m.addConstrs(quicksum(y[i][j]*(LMP[j]-LMP[i]) for j in range(n))
        +quicksum(y[i][j]*(psineg[i,j]-psipos[i,j])for j in range(n)) 
        + quicksum(y[j][i]*(psipos[j,i]-psineg[j,i]) for j in range(n))==0 for i in range(n))
        
        
        for i in range(n):
            for j in range(n): 
                if y[i][j] >0:
                    m.addConstr(y[i][j]*(theta[i]-theta[j]) <= fmax[i][j])
        
        for i in range(n):
            for j in range(n):
                if y[i][j] >0:  
                    m.addConstr(y[i][j]*(theta[i]-theta[j])>= -fmax[i][j]  )
        
        
        m.addConstrs(-p[i-1]+pmax[i-1]>=0 for i in ig)
        
        m.addConstrs(p[i]-d[i]==quicksum(y[i][j]*(theta[i]-theta[j]) for j in range(n)) for i in range(n))
        
        #old
        #m.addConstr(quicksum(Bidd[i-1]*p[i-1] for i in ig)==quicksum(d[i]*LMP[i] for i in range(n))- quicksum(pmax[i]*phi[i] for i in range(n)) -quicksum(fmax[i][j]*(psipos[i,j]+psineg[i,j])for i in range(n) for j in range(n)))    
        #new
        m.addConstr(quicksum(quicksum(bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i])) ) for i in ig)==quicksum(d[i]*LMP[i] for i in range(n))- quicksum(pmax[i]*phi[i] for i in range(n)) -quicksum(fmax[i][j]*(psipos[i,j]+psineg[i,j])for i in range(n) for j in range(n))      )
        
        m.addConstrs(theta[i-1]==0 for i in slack_bus )
        m.addConstrs(p[i-1]==0 for i in not_ig)
        
        
        m.addConstr(b[fac-1,bidindex] == 1) #fixleme
        
        for i in ig:
            if i == fac:
                m.addConstr(bidset[fac][bidindex]-LMP[i-1]+phi[i-1]>=0)
            else: 
                m.addConstr(Bidd[i-1]-LMP[i-1]+phi[i-1]>=0)
                m.addConstr(Bidd[i-1]==quicksum(bid[i-1][k]*b[i-1,k] for k in range(len(bidset[i]))))
            
        #nash constraint
        #objrhs = -p[fac-1]*c[fac-1]+ m.sum(2**r * vir[fac-1,r] for r in mccormick[i])
        m.addConstr(obj >= -p[fac-1]*c[fac-1]+ quicksum( bid[fac-1][k]*zz[fac-1,k] for k in range(len(bidset[fac])))+pmax[fac-1]*phi[fac-1])
        
    #    return obj
    
    #BuildModel(m)
    
    for i in ig: 
       obj = -p[i-1]*c[i-1]+ quicksum( bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i])))+pmax[i-1]*phi[i-1]
       #m.add_constraint(obj >= objrhs) #new
       for k in range(len(bidset[i])):
           RHSmodels(m,i,k,obj,Bidd)


    m.Params.Heuristics=0        
    m.optimize()
    
    
    nash = []
       
    Bid=[0]*n
    for i in range(n): 
        Bid[i]=round(Bidd[i].x)  
    
    nash.append(Bid)
    Profit_dcopf=np.zeros(n)
    for i in ig:
            Profit_dcopf[i-1]=(LMP[i-1].x-c[i-1])*p[i-1].x
    
    m.dispose()     
    
    #new model for collusives
    
    m=Model("Ref2") 
    
#    lambdaa = m.addVar(vtype=GRB.CONTINUOUS, name="lambdaa")
    theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta")
    LMP=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="LMP")
    Bidd=m.addVars(n,vtype=GRB.INTEGER, name="Bid")
    p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p")
    phi=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="phi")
    zz=m.addVars(n,a,lb=0,vtype=GRB.CONTINUOUS, name="zz")
    psineg=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psineg")
    psipos=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psipos")
    b=m.addVars(n,a,vtype=GRB.BINARY, name="b") #attention to # of columns (ig,k)    
    
    v1=m.addVars(n,s,vtype=GRB.CONTINUOUS, name="v1")
    v2=m.addVars(n,s,vtype=GRB.CONTINUOUS, name="v2")
    v3=m.addVars(n,n,s,vtype=GRB.CONTINUOUS, name="v3")
    v4=m.addVars(n,n,s,vtype=GRB.CONTINUOUS, name="v4")
    
    
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
    
    m.addConstrs(LMP[i-1]>= Bidd[i-1] for i in ig) #new   
    
    #collusive cut
    for i in ig:
            m.addConstr( -p[i-1]*c[i-1] + quicksum( bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i])))+pmax[i-1]*phi[i-1] >= Profit_dcopf[i-1]+0.0001)

    
    suspicious=[]

    m.Params.Heuristics=0     
    m.optimize()

    
    Bid=[0]*n
    for i in range(n):
        Bid[i]=round(Bidd[i].x)
    
    first_sol=Bid
    suspicious.append(Bid)
            
    #while lambdaa.x>0.00001:
    while True:
        add_cut(Bid)
        m.optimize()
        status = m.status                
        if status == 3 or status == 4:
                break
        #if lambdaa.x>0.00001:
        ee+=1
        Bid=[0]*n
        for i in range(n):
            Bid[i]=round(Bidd[i].x)
        suspicious.append(Bid)
    #    else:
    #        break
                        
    end=time.time()        
    time=end-start
            
    
    num_collusive=0 #to count number of found collusives

    for Bid in suspicious:
        if Bid in collusive_set.tolist():
            num_collusive+=1
    

    with open(output_file, "a") as f:
        f.write(" "+str(ee)+" ")
        f.write(" "+str(len(suspicious))+" ")
        f.write(" "+str(len(collusive_set) / len(suspicious))+" ")
        f.write(" "+str(num_collusive / len(collusive_set))+" ")
        f.write(" "+str(round(time,2))+" ")
        f.write("\n")
    return