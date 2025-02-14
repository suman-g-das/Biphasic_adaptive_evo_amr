# Wright-Fisher model on TIL landscape

import math
import numpy as np
import sys
import itertools
from numpy.random import seed
from numpy.random import randn
from numpy.random import exponential
np.set_printoptions(threshold=sys.maxsize)

seed(568612)

def decimal_to_binary(decimal_num,L):
    binary_arr = [int(x) for x in list('{0:0b}'.format(decimal_num))]
    binary_arr=[0]*(L-len(binary_arr))+binary_arr
    return binary_arr
    
def binary_to_decimal(binary_arr,L):
	dec_num=0
	for i in range(L):
		dec_num=dec_num+binary_arr[L-1-i]*2**i
	return dec_num

	
def geno_landscape_til():
	# additive effects
	global u,v,ugeno,vgeno,mgeno
	mgeno=np.zeros(2**L)
	for i in range(2**L):
		mgeno[i]=float(np.sum(decimal_to_binary(i,L)))
	v=np.zeros(L)
	u=np.zeros(L)
	vgeno=np.zeros(2**L)
	ugeno=np.zeros(2**L)
	fitness=np.zeros(2**L)

	v=np.random.normal(0.0,1.0,L)
	v=np.absolute(v)
	for i in range(L):
		flg=0
		while flg==0:
			z=np.abs(1.*(np.random.normal(0.0,1.0)))
			if (z<2.*v[i]):
				u[i]=-z
				flg=1
	for i in range(2**L):
		geno=decimal_to_binary(i,L)
		ugeno[i]=np.sum(u*geno)
		vgeno[i]=np.sum(v*geno)
	for i in range(2**L):
		geno=decimal_to_binary(i,L)
		fitness[i]=ugeno[i]-np.log(1+np.exp(2.0*(conc-vgeno[i])))
	ftns=np.exp(fitness)
	return(ftns)

def ftns_uv(ul,vl,log_conc):
	fitness=ul-np.log(1+np.exp(2.0*(log_conc-vl)))
	return(np.exp(fitness))

def hamming_distance(g1,g2,L):
	g1_binary=decimal_to_binary(g1,L)
	g2_binary=decimal_to_binary(g2,L)
	distance=0.0
	for i in range(L):
		distance=distance+float(abs(g1_binary[i]-g2_binary[i]))
	return(distance)


def one_neighbor(L,genotype,k):
	geno_binary=[int(x) for x in list('{0:0b}'.format(genotype))]
	geno_binary=[0]*(L-len(geno_binary))+geno_binary
	geno_binary[k]=1-geno_binary[k]
	neighbor_geno_dec=0
	for i in range(L):
		neighbor_geno_dec=neighbor_geno_dec+2**(L-i-1)*geno_binary[i]
	return(neighbor_geno_dec)
	
	
def neighbors_all(L,genotype):
	nbrs_set=[]
	for k in range(L):
		nbr_one=one_neighbor(L,genotype,k)
		nbrs_set.append(nbr_one)
	return(nbrs_set)


def selected_pops(pops,pop_list):
	probs=np.zeros(len(pops))
	for i in pop_list:
		ftns[int(i)]=ftns_uv(ugeno[i],vgeno[i],log_conc)
	for i in range(len(pops)):
		probs[i]=pops[i]*ftns[int(pop_list[i])]
	probs=probs/(np.sum(probs))
	pops_offspring=np.random.multinomial(np.sum(pops), probs)
	return(pops_offspring.tolist())

def mutated_pops(pops,pop_list,mu):
	type_num=len(pops)
	p_pops=np.full(len(pops),mu)
	mut_num=np.zeros(len(pops))
	for i in range(len(pops)):
		mut_num[i]=np.random.binomial(pops[i],mu)
		if mut_num[i]>pops[i]:
			mut_num[i]=pops[i]
	pops_new=pops	
	for i in range(len(pops)):
		pops_new[i]=pops_new[i]-mut_num[i]
		nbrs=neighbors_all(L,pop_list[i])
		nbrs_muts=np.random.multinomial(mut_num[i], [1./float(L)]*L)
		for j in range(len(nbrs)):
			if nbrs[j] not in pop_list:
				if nbrs_muts[j]>0:
					pop_list.append(nbrs[j])
					pops_new.append(nbrs_muts[j])
			else:
				index=pop_list.index(nbrs[j])
				pops_new[index]=pops_new[index]+nbrs_muts[j]
	
	return([pops_new,pop_list])



def trajectory(pops_ini,mu,t_lim):
	global log_conc
	pop_list=[]
	pops=[]
	type_counter=0
	pop_mean=np.zeros((t_lim,4))
	for i in range(2**L):
		if pops_ini[i]>0:
			pops.append(pops_ini[i])
			type_counter=type_counter+1
			pop_list.append(i)
			
	h_seq=[]
	upop=np.zeros(t_lim)
	vpop=np.zeros(t_lim)
	fpop=np.zeros(t_lim)
	
	for i in pop_list:
		pop_mean[0]=pop_mean[0]+np.array([mgeno[i]*pops[pop_list.index(i)],ugeno[i]*pops[pop_list.index(i)],vgeno[i]*pops[pop_list.index(i)],ftns[i]*pops[pop_list.index(i)]])
	pop_mean[0]=pop_mean[0]/float(np.sum(pops))
	for k in range(1,t_lim):
		log_conc=np.random.normal(conc,fluc*conc)
		pops=selected_pops(pops,pop_list)
		pops2=1*pops
		pop_list2=1*pop_list
		del_counter=0
		for i in range(len(pops2)):
			if pops2[i]==0:
				pop_list.remove(pop_list2[i])
				del pops[i-del_counter]
				del_counter=del_counter+1
		xy=mutated_pops(pops,pop_list,mu)
		pops=xy[0]
		pop_list=xy[1]
		
		pops2=1*pops
		pop_list2=1*pop_list
		del_counter=0
		for i in range(len(pops2)):
			if pops2[i]==0:
				pop_list.remove(pop_list2[i])
				del pops[i-del_counter]
				del_counter=del_counter+1
		for i in range(len(pops)):
			if pops[i]==0:
				pop_list.remove(pop_list[i])
				del pops[i]

		for i in pop_list:
			pop_mean[k]=pop_mean[k]+np.array([mgeno[i]*pops[pop_list.index(i)],ugeno[i]*pops[pop_list.index(i)],vgeno[i]*pops[pop_list.index(i)],ftns[i]*pops[pop_list.index(i)]])
		pop_mean[k]=pop_mean[k]/float(np.sum(pops))


	return(pop_mean)

def one_landscape(pops_ini,mu,trajectory_num):
	traj_pop_mean=np.zeros((t_lim,4))
	for i in range(trajectory_num):
		traj_pop_mean=traj_pop_mean+trajectory(pops_ini,mu,t_lim)
	traj_pop_mean=traj_pop_mean/float(trajectory_num)
	return(traj_pop_mean)
		
		
def landscape_ensemble(eps,sigma,mu,landscape_num,trajectory_num):
	global ftns
	ls_traj_pop_mean=np.zeros((t_lim,4))

	for j in range(landscape_num):
		pops_ini=np.zeros(2**L)
		pops_ini[0]=N
		freqs_ini=pops_ini/(np.sum(pops_ini))
		ftns=geno_landscape_til()
		xy=one_landscape(pops_ini,mu,trajectory_num)
		ls_traj_pop_mean=ls_traj_pop_mean+xy
	ls_traj_pop_mean=ls_traj_pop_mean/(float(landscape_num))
	print(' ',file=f)
	for i in range(t_lim):
		print(i,*ls_traj_pop_mean[i],file=f)

	return(ls_traj_pop_mean)


L=20
sigma=.1
eps=.014
N=1000000
vavg=0.7974940
conc=.5*vavg*float(L)
fluc=.1
mu=.1/float(N)


landscape_num=10000
trajectory_num=1
t_lim=100000
f=open("file.d","w+")
output=landscape_ensemble(eps,sigma,mu,landscape_num,trajectory_num)



				
	

	
			
	

