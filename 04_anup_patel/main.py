#!/usr/bin/env python
# coding: utf-8

#### Author :: Anup Patel

####################################################################################################

import pandas as pd
import numpy as np
from scipy import stats

####################################################################################################


##### Read File

with open('../data/Raw Data_GeneSpring.txt') as fp:
    # 1. iterate over file line-by-line
    # 2. strip line of newline symbols
    # 3. split line by spaces into list (of number strings)
    # 4. convert number substrings to int values
    # 5. convert map object to list
    data = [list(map(str, line.strip().split("\t"))) for line in fp]

data=np.array(data)
print("Code Running.... Please wait\n")
#Gene Symbol
#for Gene Symbol
gene_symbol=[]
for i in range(1,len(data)):
        if(len(data[i])>49):
            gene_symbol.append(data[i][49])
        else:
            gene_symbol.append('NA')

for i in range(1,len(data)):
    for j in range(len(data[1])):
        if(j!=0):
            data[i][j]=float(data[i][j])

tmp=[]
for i in range(1,len(data)):
    tmp.append(data[i][3])

# Male Non Smokers - (106-117)
# Male Smokers - (118-129)
# FeMale Non Smokers - (130-141)
# FeMale Smokers - (142-153)

#non_smokers_male
non_smokers_male=[[0.0]*12 for p in range(1,len(data))]
for i in range(1,len(data)):
    for j in range(1,13):
        non_smokers_male[i-1][j-1]=data[i][j]

#smokers_male
smokers_male=[[0.0]*12 for p in range(1,len(data))]
x=-1
for i in range(1,len(data)):
    x=x+1
    y=0
    for j in range(13,25):
        smokers_male[x][y]=data[i][j]
        y=y+1
        
#non_smokers_female
non_smokers_female=[[0.0]*12 for p in range(1,len(data))]
x=-1
for i in range(1,len(data)):
    x=x+1
    y=0
    for j in range(25,37):
        non_smokers_female[x][y]=data[i][j]
        y=y+1
        

#smokers_female
smokers_female=[[0.0]*12 for p in range(1,len(data))]
x=-1
for i in range(1,len(data)):
    x=x+1
    y=0
    for j in range(37,49):
        smokers_female[x][y]=data[i][j]
        y=y+1

male=np.hstack((non_smokers_male,smokers_male))
female=np.hstack((non_smokers_female,smokers_female))

all_data=np.hstack((male,female))

####################################################################################################

#### Null Hypothesis

# Male Non Smokers - (106-117)
# Male Smokers - (118-129)
# FeMale Non Smokers - (130-141)
# FeMale Smokers - (142-153)
#matrix=[male_mean,female_mean,non_smoker_mean,smoker_mean]

A_null=[[0]*4 for p in range(48)]
for i in range(12):
    A_null[i][0]=1
    A_null[i][2]=1
    
for i in range(12,24):
    A_null[i][0]=1
    A_null[i][3]=1

for i in range(24,36):
    A_null[i][1]=1
    A_null[i][2]=1
    
for i in range(36,48):
    A_null[i][1]=1
    A_null[i][3]=1
    


#### Alternative Hypothesis

#matrix=[male_ns_mean,male_s_mean,female_ns_mean,female_s_mean]

A=[[0]*4 for p in range(48)]
for i in range(12):
    A[i][0]=1
    
for i in range(12,24):
    A[i][1]=1

for i in range(24,36):
    A[i][2]=1
    
for i in range(36,48):
    A[i][3]=1


A=np.matrix(A)
A_null=np.matrix(A_null)
#Rank Calculation
r1=np.linalg.matrix_rank(A_null) 
r2=np.linalg.matrix_rank(A)

####################################################################################################

#Calculation of F statistics
n=48
I=np.identity(48)
F=[]
for i in range(len(all_data)):
   
    h=all_data[i]
    tmp1=np.matmul(h.T, I - (np.matmul(np.matmul(A_null,np.linalg.pinv(np.matmul(A_null.T,A_null))),A_null.T)))
    tmp2=np.matmul(h.T, I - (np.matmul(np.matmul(A,np.linalg.pinv(np.matmul(A.T,A))),A.T)))
    if(np.matmul(tmp2,h)==0):
        f_stats= ((np.matmul(tmp1,h)/0.0000000000000000001)-1)*(n-r2)/(r2-r1)
        F.append(f_stats)
    else:
        f_stats= ((np.matmul(tmp1,h)/np.matmul(tmp2,h))-1)*(n-r2)/(r2-r1)
        F.append(f_stats)
    
F=np.array(F)
F=F.tolist()

for i in range(len(F)):
    F[i]=F[i][0]

##### p-value calculation

dfd=48-r2
dfn=r2-r1
p_score=1-stats.f.cdf(F,dfn,dfd)
print("p-Values Generated")

###### p-value histogram

import matplotlib.pyplot as plt
plt.hist(p_score,bins=20)
plt.savefig("histogram.png")
print("Histogram for p-values saved")

####################################################################################################
count=0
dict_genes={}
genes_list=[]
genes_symbol_list=[]
for i in range(len(p_score)):
    if(p_score[i]<0.05):
        count=count+1
        dict_genes[i]=[data[i][0],p_score[i][0]]
        genes_list.append(data[i][0])
        if(gene_symbol[i]!='NA'):
            genes_symbol_list.append(gene_symbol[i])
#print(count)
    
####################################################################################################


with open("genes_probe_list.txt", "w") as txt_file:
    for line in genes_list:
        txt_file.write("".join(line) + "\n")



with open("genes_symbol_list.txt", "w") as txt_file:
    for line in genes_symbol_list:
        txt_file.write("".join(line) + "\n")
print("Genes Symbol list file generated")

####################################################################################################

import operator
sorted_dict = sorted(dict_genes.values(), key=operator.itemgetter(1), reverse=False)

#### FDR 

q=0.05
n0=811*0.2
ps=[]
for i in range(len(p_score)):
    ps.append(p_score[i][0])

ps=np.sort(ps)
c=0
for i in range(0,len(ps)):
    if((n0*ps[i]/(i+1))<=q):
        c=c+1

####################################################################################################


#### Find Intersection  (File 1)

with open('../data/XenobioticMetabolism1.txt') as fp:
    next(fp)
    next(fp)
    data_f = [list(map(str, line.strip().split("\t"))) for line in fp]

check=[]
for i in range(len(data_f)):
    check.append(data_f[i][0])

count=0
intersection1=[]
for p in genes_symbol_list:
    for q in check:
        if(p==q):
            intersection1.append(p)
            count=count+1

print("Genes Symbol Intersetion with Xenobiotic: ", intersection1)
print("Count ::", len(intersection1) )

####################################################################################################

# ### Find Intersection ( File 2) 

with open('../data/FreeRadicalResponse.txt') as fp:
    next(fp)
    next(fp)
    data_f = [list(map(str, line.strip().split("\t"))) for line in fp]

check=[]
for i in range(len(data_f)):
    check.append(data_f[i][0])
    

count=0
intersection2=[]
for p in genes_symbol_list:
    for q in check:
        if(p==q):
            intersection2.append(p)
            count=count+1
#### No Intersection

####################################################################################################

#### Find Intersection (File 3) 

with open('../data/DNARepair1.txt') as fp:
    next(fp)
    next(fp)
    data_f = [list(map(str, line.strip().split("\t"))) for line in fp]


check=[]
for i in range(len(data_f)):
    check.append(data_f[i][0])

count=0
intersection3=[]
for p in genes_symbol_list:
    for q in check:
        if(p==q):
            intersection3.append(p)
            count=count+1

print("Genes Symbol Intersetion with DNA Repair:", intersection3)
print("Count ::", len(intersection3) )

####################################################################################################

#### Find Intersection (File 4) 

with open('../data/NKCellCytotoxicity.txt') as fp:
    next(fp)
    next(fp)
    data_f = [list(map(str, line.strip().split("\t"))) for line in fp]

check=[]
for i in range(len(data_f)):
    check.append(data_f[i][0])


count=0
intersection4=[]
for p in genes_symbol_list:
    for q in check:
        if(p==q):
            intersection4.append(p)
            count=count+1

print("Genes Symbol Intersetion with NKCellCytotoxicity : ", intersection4)
print("Count ::", len(intersection4) )

####################################################################################################

#### Grouping lists 

female_s_up=[]
female_s_down=[]
male_s_up=[]
male_s_down=[]


#### Grouping (1)

rows=[]
tmp=[]
for i in range(len(gene_symbol)):
    for j in range(len(intersection1)):
        if(gene_symbol[i]==intersection1[j]):
            #print(intersection1[j])
            tmp.append(intersection1[j])
            rows.append(i)
            break 

print("\nGrouping for Xenobiotic")
print("Probe-Name","Gene-Symbol","Male-NS","Male-S","Female-NS","Female-S")
for i in rows:
    male_ns=all_data[i,0:12]
    male_s=all_data[i,12:24]
    female_ns=all_data[i,24:36]
    female_s=all_data[i,36:48]
    male_ns_mean=np.mean(male_ns)
    male_s_mean=np.mean(male_s)
    female_ns_mean=np.mean(female_ns)
    female_s_mean=np.mean(female_s)
    if(male_ns_mean<male_s_mean):
        male_s_up.append(gene_symbol[i])
    if(male_ns_mean>male_s_mean):
        male_s_down.append(gene_symbol[i])
    if(female_ns_mean<female_s_mean):
        female_s_up.append(gene_symbol[i])
    if(female_ns_mean>female_s_mean):
        female_s_down.append(gene_symbol[i])
    
    print(data[i+1][0],gene_symbol[i],male_ns_mean,male_s_mean,female_ns_mean,female_s_mean)


#### Grouping (2)

rows=[]
tmp=[]
for i in range(len(gene_symbol)):
    for j in range(len(intersection2)):
        if(gene_symbol[i]==intersection2[j]):
            #print(intersection1[j])
            tmp.append(intersection2[j])
            rows.append(i)
            break 


#### Grouping (3) 

rows=[]
tmp=[]
for i in range(len(gene_symbol)):
    for j in range(len(intersection3)):
        if(gene_symbol[i]==intersection3[j]):
            #print(intersection1[j])
            tmp.append(intersection3[j])
            rows.append(i)
            break 

print("\nGrouping for DNA Repair")
print("Probe-Name","Gene-Symbol","Male-NS","Male-S","Female-NS","Female-S")
for i in rows:
    male_ns=all_data[i,0:12]
    male_s=all_data[i,12:24]
    female_ns=all_data[i,24:36]
    female_s=all_data[i,36:48]
    male_ns_mean=np.mean(male_ns)
    male_s_mean=np.mean(male_s)
    female_ns_mean=np.mean(female_ns)
    female_s_mean=np.mean(female_s)
    if(male_ns_mean<male_s_mean):
        male_s_up.append(gene_symbol[i])
    if(male_ns_mean>male_s_mean):
        male_s_down.append(gene_symbol[i])
    if(female_ns_mean<female_s_mean):
        female_s_up.append(gene_symbol[i])
    if(female_ns_mean>female_s_mean):
        female_s_down.append(gene_symbol[i])
    print(data[i+1][0],gene_symbol[i],male_ns_mean,male_s_mean,female_ns_mean,female_s_mean)


#### Grouping (4) 

rows=[]
tmp=[]
for i in range(len(gene_symbol)):
    for j in range(len(intersection4)):
        if(gene_symbol[i]==intersection4[j]):
            #print(intersection1[j])
            tmp.append(intersection4[j])
            rows.append(i)
            break 

print("\nGrouping for Natural Killer Cell")
print("Probe-Name","Gene-Symbol","Male-NS","Male-S","Female-NS","Female-S")
for i in rows:
    male_ns=all_data[i,0:12]
    male_s=all_data[i,12:24]
    female_ns=all_data[i,24:36]
    female_s=all_data[i,36:48]
    male_ns_mean=np.mean(male_ns)
    male_s_mean=np.mean(male_s)
    female_ns_mean=np.mean(female_ns)
    female_s_mean=np.mean(female_s)
    if(male_ns_mean<male_s_mean):
        male_s_up.append(gene_symbol[i])
    if(male_ns_mean>male_s_mean):
        male_s_down.append(gene_symbol[i])
    if(female_ns_mean<female_s_mean):
        female_s_up.append(gene_symbol[i])
    if(female_ns_mean>female_s_mean):
        female_s_down.append(gene_symbol[i])
    print(data[i+1][0],gene_symbol[i],male_ns_mean,male_s_mean,female_ns_mean,female_s_mean)

####################################################################################################


print("\n****************** Final Grouping ***********************")
print("Female Smokers up genes" ,set(female_s_up))
print("Female Smokers Down genes", set(female_s_down))
print("Male Smokers up genes" ,set(male_s_up))
print("Male Smokers down genes" ,set(male_s_down))

####################################################################################################




