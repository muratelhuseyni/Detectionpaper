#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 23:44:56 2019

@author: elifyilmaz
"""

# Linear, Gurobi

from gurobipy import *
import numpy as np
import pandas as pd
import itertools as it

################### for different sizes just change this part rest is the same ###################

n=5 # of nodes

ig=[1,2,5] #generator nodes
not_ig=[] #non_generator nodes
for i in range(n):
    if i+1 not in ig:
        not_ig.append(i+1)

slack_bus=[3]

input_file = 'Small_5.xlsx'
output_file="detect4_gurobi_small.csv"
dtt=pd.read_excel('Small_7.xlsx', header=None, sheet_name='CollusiveStrategies')

#########################################################################################

def linear(input_file,output_file,dtt,n ,ig,not_ig,slack_bus):

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
    
    bidset={}
    for i in ig:
        bidset[i]=[x for x in list(bid[i-1]) if x!=0]
    
    a=max(len(bidset[i]) for i in ig ) #for column number
    
    pmax = [a/100 for a in pmax1]
    d = [a/100 for a in d1]
    fmax = fmax1/100
    
    ee=1 #to count number of alternative solutions
    
    num_col_bidset=bid.shape[1] # of columns
    
    m=Model("Ref1")
    
    theta=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="theta")
    LMP=m.addVars(n,lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name="LMP")
    Bidd=m.addVars(n, vtype=GRB.INTEGER, name="Bid") 
    p=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="p")    
    phi=m.addVars(n,lb=0,vtype=GRB.CONTINUOUS, name="phi")
    b=m.addVars(n,a,vtype=GRB.BINARY, name="b") #attention to # of columns (ig,k)
    psineg=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psineg")
    psipos=m.addVars(n,n,lb=0,vtype=GRB.CONTINUOUS, name="psipos")
    lambdaa = m.addVar(vtype=GRB.CONTINUOUS, name="lambdaa")
    zz=m.addVars(n,a, vtype=GRB.CONTINUOUS, name="zz")
    
    
    
    m.setObjective(lambdaa, GRB.MAXIMIZE)
    
    m.addConstrs(quicksum(b[i-1,k] for k in range(len(bidset[i])))==1 for i in ig)
    m.addConstrs(Bidd[i-1]==quicksum(bidset[i][k]*b[i-1,k]for k in range(len(bidset[i]))) for i in ig) 
    
    m.addConstrs(lambdaa<= -p[i-1]*c[i-1]+ quicksum( bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i]))) for i in ig)
    m.addConstrs(Bidd[i-1]-LMP[i-1]+phi[i-1]>=0 for i in ig)
    
    m.addConstrs(zz[i-1,k]<= pmax[i-1]*b[i-1,k] for i in ig for k in range(len(bidset[i])))
    m.addConstrs(zz[i-1,k]<= p[i-1] for i in ig for k in range(len(bidset[i])) )
    m.addConstrs(zz[i-1,k]>= p[i-1]-pmax[i-1]*(1-b[i-1,k]) for i in ig for k in range(len(bidset[i])) )
    
    
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
    m.addConstr(quicksum(quicksum(bid[i-1][k]*zz[i-1,k] for k in range(len(bidset[i]))) for i in ig)==quicksum(d[i]*LMP[i] for i in range(n))- quicksum(pmax[i]*phi[i] for i in range(n)) -quicksum(fmax[i][j]*(psipos[i,j]+psineg[i,j])for i in range(n) for j in range(n)))    
    
    m.addConstrs(theta[i-1]==0 for i in slack_bus )
    m.addConstrs(p[i-1]==0 for i in not_ig)
    
    m.Params.Heuristics=0
    
    m.optimize()
    
    
    
    def add_cut(Bid):
        bid_index=np.zeros(n)
        for i in ig:
            for j in range(len(bidset[i])):
                if Bid[i-1]==bidset[i][j] : 
                    bid_index[i-1]=j
                    break
        m.addConstr(quicksum(b[i-1,bid_index[i-1]]for i in ig)<=len(ig)-1) 
        return  
      
    suspicious=[]
    
    Bid=[0]*n
    for i in range(n):
        Bid[i]=round(Bidd[i].x)
    
    first_sol=Bid
    suspicious.append(Bid)
    
            
    while lambdaa.x>0.00001:
        add_cut(Bid)
        m.optimize()
        if lambdaa.x>0.00001:
            ee+=1
            Bid=[0]*n
            for i in range(n):
                Bid[i]=round(Bidd[i].x)
            suspicious.append(Bid)
        else:
            break
    
    
    end=time.time()        
    time=end-start
            
    
       
    num_collusive=0 #to count number of found collusives
    
    for Bid in suspicious:
        if Bid in collusive_set.tolist():
            num_collusive+=1
    
    with open(output_file, "a") as f:
        f.write(" "+str(input_file)+" ")
        f.write(" "+str(ee)+" ")
        f.write(" "+str(len(suspicious))+" ")
        f.write(" "+str(len(collusive_set) / len(suspicious))+" ")
        f.write(" "+str(num_collusive / len(collusive_set))+" ")
        f.write(" "+str(round(time,2))+" ")
        f.write("\n")
    return

linear(input_file,output_file,dtt,n ,ig,not_ig,slack_bus)








