from neuron import h,gui
from matplotlib import pyplot as p
import numpy as np
from tempCorr import *
from nernst import *
from math import pi
import random


tstop = .05*1000
h.tstop=tstop
h.dt=0.1
h.celsius = (T-273)



def define_geometry(icc):

	icc.diam = 4.4			#7.136*2
	icc.L = 65.92					
	icc.nseg = 1
	icc.Ra = 50		#50
	icc.cm = .025*(10**5)/h.area(0.5,sec=icc)
	

def insert_mechanisms(icc):

	icc.insert('pas')
	icc.g_pas = .01
	icc(0.5).g_pas=0

	icc.insert('Na')
	icc.nai = nai_
	icc.nao = nao_
	icc.ena = ena_
	icc.G_Na_Na= 0		#20
	icc.tau_f_Na_Na = tau_f_Na
	icc.tau_d_Na_Na = tau_d_Na

	icc.insert('nscc')
	icc.G_NSCC_nscc = 0	#12.15
	icc.tau_NSCC_nscc = tau_NSCC

	icc.insert('ERG')
	icc.G_ERG_ERG = 0	#2.5
	icc.tau_ERG_ERG = tau_ERG 
	icc.ki= ki_	#23
	icc.ko= ko_	#23
	icc.ek= ek_	#23

	icc.insert('bk')
	icc.G_bk_bk= 0 #23+T_correction_BK

	icc.insert('Kb')
	icc.G_Kb_Kb=0	#.15

	icc.insert('kv11')
	icc.G_Kv11_kv11=0   #6.3
	icc.tau_d_kv11_kv11 = tau_d_kv11
	icc.tau_f_kv11_kv11 = tau_f_kv11

	icc.insert('vddr')
	icc.G_Ca_VDDR_vddr=0	#3
	icc.cao = 2.5
	icc.tau_d_VDDR_vddr = tau_d_VDDR
	icc.tau_f_VDDR_vddr = tau_f_VDDR

	icc.insert('ltype')
	icc.G_Ca_Ltype_ltype=0 #2
	icc.tau_f_Ltype_ltype = tau_f_Ltype
	icc.tau_d_Ltype_ltype = tau_d_Ltype
	icc.tau_f_Ca_Ltype_ltype = tau_f_Ca_Ltype

	icc.insert('pmca')
	icc.J_max_PMCA_pmca= 0 #0.088464
	icc.Vol_pmca = pi*((icc.diam/2)**2)*icc.L
	icc.a_pmca = h.area(0.5,sec=icc)		#area

	time_corr = 10**(-3)
	icc.insert('concyto')
	J_max_leak=0.01*time_corr
	icc.J_max_leak_concyto= 0	#0.01(1/s) 
	icc.Vol_concyto = pi*((icc.diam/2)**2)*icc.L
	icc.a_concyto = h.area(0.5,sec=icc)		#area



	icc.insert('conpu')
	icc.Vol_conpu = pi*((icc.diam/2)**2)*icc.L
	icc.J_max_leak_conpu= 0 #0.01(mM/s) 
	icc.Jmax_serca_conpu = 0		# 1.8333	
	icc.J_ERleak_conpu = 0		# 1.666667	

	icc.Jmax_IP3_conpu = 	0		#50000  (1/s)
	icc.Jmax_NaCa_conpu =  0 	#0.05(mM/s)
	icc.Jmax_uni_conpu = 0 		#5000 mM/s

	icc.insert('ano1')
	icc.g_Ano1_ano1 = 0	#20
	icc.cli = cli_
	icc.clo = clo_
	icc.ecl = ecl_

	icc.insert('cacl')
	icc.G_Cacl_cacl = 0	#10.1
	icc.tau_act_Cacl_cacl = tau_act_Cacl


	active_sites = [0.5]
	a=1/(10*h.area(0.5))			# conversion factor for conductances - page 11
	

	for i in active_sites:

		
		icc(i).G_Na_Na= a*20			#20		

		icc(i).G_NSCC_nscc = a*12.15	#12.15

		icc(i).G_ERG_ERG = a*2.5	#2.5

		icc(i).G_bk_bk= a*(23+T_correction_BK) #23+T_correction_BK

		icc(i).G_Kb_Kb=a*0.15	#.15

		icc(i).G_Kv11_kv11=a*6.3   #6.3

		icc(i).G_Ca_VDDR_vddr=a*3	#3

		icc(i).G_Ca_Ltype_ltype=a*2 #2


		icc(i).J_max_PMCA_pmca= 0.088464 #0.088464

		icc(i).J_max_leak_concyto= J_max_leak	#0.01(mM/s) 

		icc(i).J_max_leak_conpu= J_max_leak #0.01(mM/s) 

		icc(i).Jmax_serca_conpu = 1.8333*time_corr		# 1.8333	

		icc(i).J_ERleak_conpu = 1.666667*time_corr		# 1.666667	

		icc(i).Jmax_IP3_conpu = 	50000*time_corr		#50000  (1/s)

		icc(i).Jmax_NaCa_conpu = 0.05*time_corr 	#0.05(mM/s)

		icc(i).Jmax_uni_conpu = 5000*time_corr 		#5000 mM/s



		icc(i).g_Ano1_ano1 = a*20	#20

		icc(i).G_Cacl_cacl = a*10.1	#10.1

 

def initialize_network(ncells): 
	cells =  [ [[[] for i in range(ncells)] for i in range(ncells)] for i in range(ncells) ]
	gap_junctions_f_x=[]
	gap_junctions_b_x=[]
	gap_junctions_f_y=[]
	gap_junctions_b_y=[]
	gap_junctions_f_z=[]
	gap_junctions_b_z=[]
	
	for i in range(ncells):
		for j in range(ncells):
			for k in range(ncells):
			 	cells[i][j][k] = h.Section()
			 	define_geometry(cells[i][j][k])
			 	insert_mechanisms(cells[i][j][k])


	for k in range(ncells):	
		for i in range(ncells):
			gj_fx=[]
			gj_bx=[]
			gj_fy=[]
			gj_by=[]
			gj_fz=[]
			gj_bz=[]
			for j in range(ncells-1):
			
				gj_fx.append(h.gap())
				gj_bx.append(h.gap())
				gj_fx[j].loc(cells[j][i][k](1))
				gj_bx[j].loc(cells[j+1][i][k](0))
				h.setpointer(cells[j+1][i][k](0)._ref_v,'vgap',gj_fx[j])  
				h.setpointer(cells[j][i][k](1)._ref_v,'vgap',gj_bx[j])

				gj_fy.append(h.gap())
				gj_by.append(h.gap())
				gj_fy[j].loc(cells[i][j][k](1))
				gj_by[j].loc(cells[i][j+1][k](0))
				h.setpointer(cells[i][j+1][k](0)._ref_v,'vgap',gj_fy[j])  
				h.setpointer(cells[i][j][k](1)._ref_v,'vgap',gj_by[j])
				
				gj_fz.append(h.gap())
				gj_bz.append(h.gap())
				gj_fz[j].loc(cells[i][k][j](1))
				gj_bz[j].loc(cells[i][k][j+1](0))
				h.setpointer(cells[i][k][j+1](0)._ref_v,'vgap',gj_fz[j])  
				h.setpointer(cells[i][k][j](1)._ref_v,'vgap',gj_bz[j])


			gap_junctions_f_x.append(gj_fx)
			gap_junctions_f_y.append(gj_fy)
			gap_junctions_b_x.append(gj_bx)
			gap_junctions_b_y.append(gj_by)
			gap_junctions_f_z.append(gj_fz)
			gap_junctions_b_z.append(gj_bz)
	gap_junctions = [gap_junctions_f_x,gap_junctions_b_x,gap_junctions_f_y,gap_junctions_b_y,gap_junctions_f_z,gap_junctions_b_z]
	return gap_junctions, cells
	# return gap_junctions_f_x,gap_junctions_f_y,gap_junctions_b_x,gap_junctions_b_y, cells

##------- Voltage clamps : 

# vclamp1 = h.SEClamp(cells[0](0.5))
# vclamp1.dur1 = 4000
# vclamp1.amp1=-70.0
# vclamp1.rs=.00001

# vclamp2 = h.SEClamp(cells[1](0.5))
# vclamp2.dur1 = 8000
# vclamp2.amp1=-70.0
# vclamp2.rs=.00001

# vclamp3 = h.SEClamp(cells[2](0.5))
# vclamp3.dur1 = 15000
# vclamp3.amp1=-70.0
# vclamp3.rs=.00001


#Recording vectors
def  set_recording_vectors(cells):
	
	t= h.Vector()
	t.record(h._ref_t)

	
	v = [ [[[] for i in range(len(cells))] for i in range(len(cells))] for i in range(len(cells)) ]
	for i in range(len(cells)):
		for j in range(len(cells)):
			for k in range(len(cells)):

				v[i][j][k] = h.Vector()
				v[i][j][k].record(cells[i][j][k](0.5)._ref_v)
	return v,t


h.v_init = -70

ncells = [2,3,4]

for N_cells in ncells:
	gap_junctions,cells = initialize_network(N_cells)
	# gap_junctions_f_x,gap_junctions_f_y,gap_junctions_b_x,gap_junctions_b_y,cells = initialize_network(N_cells)
	
	v,t = set_recording_vectors(cells)
	# cells[1][1][1].IP3_conpu= 0.00058 
	# print(cells[1][0].IP3_conpu) 
	# print((cells)) 
	# print(len(gap_junctions[0]))
	# cells[0][1].IP3_conpu =0.00066
	# cells[1][0].IP3_conpu =0.00069
	# cells[1][1].IP3_conpu =0.00063

	# ip3 = [0.00006,0.000068,0.00006]	
	#ip3 = [(0.00058 + random.randrange(10)*0.00001) for x in range(N_cells)]
	#for i in range(N_cells):
	#	cells[i].IP3_conpu = ip3[i]

	ip3 = [ [[[] for i in range(N_cells)] for i in range(N_cells)] for i in range(N_cells) ]
	for i in range(N_cells):
		for j in range(N_cells):
			for k in range(N_cells):
				ip3[i][j][k] = (0.00058 + random.randrange(10)*0.00001)
				cells[i][j][k].IP3_conpu = ip3[i][j][k]

	# cells[0][0][0].IP3_conpu = 0
	# print(ip3)			

	h.run()
	
	for i in range(N_cells):
		for j in range(N_cells):
			# v[i][j] = v[i][j].to_python()
			pass
	# t = t.to_python()
	# t = np.array(t.to_python())
	# t = t/1000
	# v = v.to_python()
	np.save('v_ncells_3d_'+str(N_cells)+'.npy', v)
	np.save('t.npy', t)
	np.save('ip3_v_ncells_3d_'+str(N_cells)+'.npy', ip3)


	p.figure()

	v_n = np.load('v_ncells_3d_'+str(N_cells)+'.npy')
	t = np.load('t.npy')
	for i in range(N_cells):
		for j in range(N_cells):
			for k in range(N_cells):
				p.plot(t,v_n[i][j][k] , label = 'cell '+ str(i)+ str(j)+str(k),color = (i/(N_cells),j/(N_cells)*1.1,k/(N_cells)*1.1))
				# p.plot(t,v[1], color = 'blue', label = 'cell 1')
				# p.plot(t,v[2], color = 'black', label = 'cell 2')

				p.xlabel('time(s)')
				p.ylabel('Membrane Voltage(mV)')

	p.legend(loc = 'upper center', ncol = 3*N_cells, bbox_to_anchor=(0.5,1.1))
p.show()








