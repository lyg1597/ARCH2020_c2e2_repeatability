from math import *
import numpy as np
import numpy.linalg as la
import sys
import time
from sympy import Derivative
def jcalc(listvalue,curstate):
    ret = []
    if curstate ==1:
        vx= listvalue[0]
        vy= listvalue[1]
        px= listvalue[2]
        py= listvalue[3]
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=1
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=1
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        return ret
    if curstate ==2:
        vx= listvalue[0]
        vy= listvalue[1]
        px= listvalue[2]
        py= listvalue[3]
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=1
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=1
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        Entry=0
        ret.append(Entry)
        return ret

if int(state) == 1:
    notbloating=[4,5,6]
    thin_indicate=[]
    bloating=[0,1,2,3]
    not_thin=[0,1,2,3]
if int(state) == 2:
    notbloating=[4,5,6]
    thin_indicate=[]
    bloating=[0,1,2,3]
    not_thin=[0,1,2,3]
thin_Jac_ind = []
for not_thin_ele in (not_thin):
    new_non_ele = not_thin_ele
    for not_bloat_ele in notbloating:
        if not_bloat_ele < not_thin_ele:
            new_non_ele-=1
    thin_Jac_ind.append(new_non_ele)

for i in range (len(thin_indicate)):
   thin_indicate[i]+=1
for i in range (len(bloating)):
   bloating[i]+=1
for i in range (len(not_thin)):
   not_thin[i]+=1

for i in range (len(notbloating)):
   #delta.pop(notbloating[i])
   notbloating[i]+=1

f = open('../work-dir/SimuOutput', 'r')
x = f.readlines()
f.close()
start = time.time()
numofvar = len(x[0].rstrip().split())
Simulation_data = [[] for i in range(len(x))]
curline = 0
for line in x:
    if line.rstrip():
        dataLine = line.rstrip().split()
    for i in range (numofvar):
        Simulation_data[curline].append(float(dataLine[i]))
    curline+=1
Simulation_data = np.array(Simulation_data)
CT_step = min(CT_step,len(Simulation_data))
blowting = np.zeros((int(len(x)),len(not_thin)))
deltacopy = list(delta)
blowting[0,:] = [delta[tmp] for tmp in [tmp_list-1 for tmp_list in not_thin]]
numofvar-=(1+len(notbloating))
Reach_tube = np.zeros((len(Simulation_data),len(Simulation_data[0])))
if numofvar>0:
   for i in range (0,len(Simulation_data),CT_step):
       mean_array = np.mean(Simulation_data[i:i+CT_step,bloating],axis=0)
       if len(mean_array)!=numofvar:
           sys.exit('Error in calculating the mean matrix for coordinated transformation')
       Jacobian_CT = jcalc(mean_array,int(state))
       Jacobian_CT = np.reshape(Jacobian_CT,(numofvar,numofvar))
       Jacobian_CT = Jacobian_CT[:,thin_Jac_ind]
       Jacobian_CT = Jacobian_CT[thin_Jac_ind,:]
       Eigenvalues,CT_matrix = la.eig(Jacobian_CT)
       Local_Lipschitz = la.norm(Jacobian_CT)
       for j in range (i,min(i+CT_step,len(x)),2):
           if j+1 > len(Simulation_data):
               sys.exit('Error accessing the last rectangle')
           Delta_time = Simulation_data[j+1,0] - Simulation_data[j,0]
           Simulation_box = Simulation_data[j:j+2,bloating]
           Current_Jacobian = jcalc(np.mean(Simulation_box,axis=0),int(state))
           Current_Jacobian = np.reshape(Current_Jacobian,(numofvar,numofvar))
           Current_Jacobian = Current_Jacobian[:,thin_Jac_ind]
           Current_Jacobian = Current_Jacobian[thin_Jac_ind,:]
           Current_Jordan = np.dot(la.inv(CT_matrix),np.dot(Current_Jacobian,CT_matrix))
           Current_lambda = max(la.eigvalsh(np.transpose(Current_Jordan)+Current_Jordan))/2
           Disturb_matrix = np.reshape(jcalc((np.amax(Simulation_box,axis=0)+max(blowting[j,:])/2),int(state)),(numofvar,numofvar))-np.reshape(jcalc((np.amin(Simulation_box,axis=0)-max(blowting[j,:])/2),int(state)),(numofvar,numofvar))
           Disturbance = la.norm(Disturb_matrix,ord=2)*exp(Local_Lipschitz*Delta_time)
           Current_lambda = Current_lambda + Disturbance
           if j == 0:
               blowting[1,:] = blowting[0,:]*exp(Current_lambda * Delta_time)
           else:
               blowting[j,:] = blowting[j-1,:]
               blowting[j+1,:] = blowting[j,:] * exp(Current_lambda * Delta_time)
       for cnt in range (min(CT_step,len(Simulation_data)-i)):
           blowting[i+cnt,:] = blowting[i+cnt,:] * la.cond(CT_matrix)
   for i in range (0,len(Simulation_data),2):
       Reach_tube[i,not_thin] = np.amin(Simulation_data[i:i+2,not_thin],axis=0)-blowting[i,:]
       Reach_tube[i+1,not_thin] = np.amax(Simulation_data[i:i+2,not_thin],axis=0)+blowting[i,:]
       Reach_tube[i,notbloating] = Simulation_data[i,notbloating] - [deltacopy[tmp] for tmp in [tmp_list-1 for tmp_list in notbloating]]
       Reach_tube[i+1,notbloating] = Simulation_data[i+1,notbloating] + [deltacopy[tmp] for tmp in [tmp_list-1 for tmp_list in notbloating]]
       Reach_tube[i,thin_indicate] = Simulation_data[i,thin_indicate]
       Reach_tube[i+1,thin_indicate] = Simulation_data[i+1,thin_indicate]
for i in range (0,len(Simulation_data),2):
   Reach_tube[i,0] = Simulation_data[i,0]
   Reach_tube[i+1,0] = Simulation_data[i+1,0]

end = time.time()
print('ComputeLDF Time: ' + str(end-start) + 's')
print('Final Bloating Size:')
print(blowting[-1,:])
f = open('../work-dir/reachtube.dat', 'w')
bloatstring = ''
bloatstring += state
bloatstring += '\n'
for i in range (Reach_tube.shape[0]):
    for j in range (Reach_tube.shape[1]):
        bloatstring+=str(Reach_tube[i][j])
        bloatstring+=' '
    bloatstring+='\n'
f.write(bloatstring)
f.close()
