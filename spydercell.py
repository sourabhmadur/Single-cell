from neuron import h, gui
import matplotlib.pyplot as plt
import numpy as np
from tempCorr import *
from nernst import *
from matplotlib.collections import LineCollection
from plots import *
from math import pi
# from run import *
import numpy as np

tstop = 500
h.tstop = tstop
h.dt=0.1
h.celsius = (T-273)


#geometry
#area = 4*100*8+2*20*8



def define_geometry(icc):

	icc.diam = 4.4			#7.136*2

	icc.L = 65.92					
	icc.nseg = 1
	
	icc.Ra = 50		#50

	area = 4*100*8+2*20*8
	icc.cm = .025*(10**5)/h.area(0.5,sec=icc)
	print('capacitance=',icc.cm)

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

	time_corr = 10**(-3)
	icc.insert('concyto')
	J_max_leak=0.01*time_corr
	icc.J_max_leak_concyto= 0	#0.01(mM/s) 
	icc.Vol_concyto = pi*((icc.diam/2)**2)*icc.L 

	icc.insert('conpu')
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
	a=1/(10*h.area(0.5))
	

	for i in active_sites:

		
		icc(i).G_Na_Na= a*20			#.0022		#this is the tuned value  #20 - total conductance

		icc(i).G_NSCC_nscc = a*0	#12.15

		icc(i).G_ERG_ERG = a*2.5	#2.5

		icc(i).G_bk_bk= a*(23+T_correction_BK) #23+T_correction_BK

		icc(i).G_Kb_Kb=a*0.15	#.15

		icc(i).G_Kv11_kv11=a*6.3   #6.3

		icc(i).G_Ca_VDDR_vddr=a*3	#3

		icc(i).G_Ca_Ltype_ltype=a*2 #2


		icc(i).J_max_PMCA_pmca= 0.088464*time_corr #0.088464

		icc(i).J_max_leak_concyto= J_max_leak	#0.01(mM/s) 

		icc(i).J_max_leak_conpu= J_max_leak #0.01(mM/s) 

		icc(i).Jmax_serca_conpu = 1.8333*time_corr		# 1.8333	

		icc(i).J_ERleak_conpu = 1.666667*time_corr		# 1.666667	

		icc(i).Jmax_IP3_conpu = 	50000*time_corr		#50000  (1/s)

		icc(i).Jmax_NaCa_conpu = 0.05*time_corr 	#0.05(mM/s)

		icc(i).Jmax_uni_conpu = 5000*time_corr 		#5000 mM/s



		icc(i).g_Ano1_ano1 = a*20	#20

		icc(i).G_Cacl_cacl = a*10.1	#10.1


def run_and_record(icc,v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG):



	v.record(icc(0.5)._ref_v)
	t.record(h._ref_t)
	ina.record(icc(.5)._ref_ina)  
	ik.record(icc(.5)._ref_ik)
	# icl.record(icc(.5)._ref_icl)
	# ica.record(icc(.5)._ref_ica)
	# ins.record(icc(.5)._ref_i_nscc)
	# eca.record(icc(.5)._ref_eca)
	cai.record(icc(.5)._ref_cai)
	# cao.record(icc(.5)._ref_cao)
	# ek.record(icc(.5)._ref_ek)
	capui.record(icc(.5)._ref_capui)
	caeri.record(icc(.5)._ref_caeri)
	cami.record(icc(.5)._ref_cami)
	# cli.record(icc(.5)._ref_cli)
	# clo.record(icc(.5)._ref_clo)
	# oano.record(icc(.5)._ref_O_Ano1_ano1)
	# ecl.record(icc(.5)._ref_ecl)
	# ena.record(icc(.5)._ref_ena)
	
	ica_vddr.record(icc(.5)._ref_ica_vddr)
	icl_ano1.record(icc(.5)._ref_icl_ano1)
	icl_cacl.record(icc(.5)._ref_icl_cacl)
	# ina_Na.record(icc(.5)._ref_ina_Na)
	ik_kv11.record(icc(.5)._ref_ik_kv11)
	ik_bk.record(icc(.5)._ref_ik_bk)
	ica_ltype.record(icc(.5)._ref_ica_ltype)
	# i_nscc.record(icc(.5)._ref_i_nscc)
	ik_Kb.record(icc(.5)._ref_ik_Kb)
	# ica_pmca.record(icc(.5)._ref_ica_pmca)
	ik_ERG.record(icc(.5)._ref_ik_ERG)

	h.run()




v = h.Vector()             # Membrane potential vector
t = h.Vector()             # Time stamp vector

cai =h.Vector()
cao =h.Vector()
capui =h.Vector()
caeri = h.Vector()
cami = h.Vector()
cli = h.Vector()
clo = h.Vector()


eca =h.Vector()
ek = h.Vector()
ecl = h.Vector()
ena = h.Vector()

oano = h.Vector()

ina = h.Vector()			#uA/c
ik = h.Vector()
ica = h.Vector()
ins = h.Vector()
icl = h.Vector()

ica_vddr = h.Vector()
icl_ano1 = h.Vector()
icl_cacl= h.Vector()
ina_Na = h.Vector() 
ik_kv11= h.Vector() 
ik_bk= h.Vector()
ica_ltype= h.Vector()
i_nscc= h.Vector() 
ik_Kb= h.Vector() 
ica_pmca= h.Vector() 
ik_ERG= h.Vector()

variables = [v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG]
# v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG = runThisAlready(icc,v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG)



icc = h.Section(name='icc')
define_geometry(icc)
insert_mechanisms(icc)

# icc1 = h.Section(name='icc1')
# icc1.nseg=1001
# icc1.insert('pas')
# icc1.L = 5000
# icc1.diam=1
# icc1.connect(icc(1))
# icc1.cm=25
v_e = h.Vector()
v_e1 = h.Vector()  
v_e2 = h.Vector()
v_e3 = h.Vector()
v_e4 = h.Vector()
v_e.record(icc(0.6)._ref_v)
v_e1.record(icc(0.7)._ref_v)
v_e2.record(icc(.8)._ref_v)
v_e3.record(icc(.9)._ref_v)
v_e4.record(icc(1)._ref_v)



vclamp = h.SEClamp(icc(0.5))
vclamp.dur1 = tstop
vclamp.amp1=-50.0
vclamp.rs=.00001
h.v_init = -80



run_and_record(icc,*variables)

t = np.array(t.to_python())
ina = h.area(0.5)*10*np.array(ina.to_python())
ik = h.area(0.5)*10*np.array(ik.to_python())
ik_ERG = h.area(0.5)*10*np.array(ik_ERG.to_python())
ik_Kb = h.area(0.5)*10*np.array(ik_Kb.to_python())
ik_kv11 = h.area(0.5)*10*np.array(ik_kv11.to_python())
ica_vddr = h.area(0.5)*10*np.array(ica_vddr.to_python())
ik_bk = h.area(0.5)*10*np.array(ik_bk.to_python())
ica_ltype = h.area(0.5)*10*np.array(ica_ltype.to_python())
icl_cacl = h.area(0.5)*10*np.array(icl_cacl.to_python())
icl_ano1 = h.area(0.5)*10*np.array(icl_ano1.to_python())


cai_vddr = np.array(cai.to_python())
capui = np.array(capui.to_python())
# ik = ik.to_python()

# ina_new = [i * (h.area(0.5)*10) for i in ina]
# ik_new = [j * (h.area(0.5)*10) for j in ik]

# plt.figure()
# plt.plot(t,v , label = 'v(0.5)', color= 'red')
# # plt.plot(t,v_e , label = 'v(0.6)', color= 'blue')
# # plt.plot(t,v_e1 , label = 'v(0.7)', color= 'green')
# # plt.plot(t,v_e2 , label = 'v(0.8)', color= 'yellow')
# # plt.plot(t,v_e3 , label = 'v(0.9)', color= 'black')
# # plt.plot(t,v_e4 , label = 'v(1)', color= 'brown')

# plt.legend(loc = 'upper right')

# plt.figure(2)
# plt.plot(t,ina*10)


# plt.show()
#print(icc.tau_f_Na_Na, tau_f_Na_Na)

# plt.plot(t,ina)
# plt.plot(t,ik_ERG)
# plt.plot(t,ik_Kb)
# plt.plot(t,ik_kv11)
# plt.plot(t,ica_vddr)

# plt.plot(t,cai)
# plt.plot(t,capui)


# plt.plot(t,ik_bk)
# plt.plot(t,ica_ltype)
# plt.plot(t,icl_cacl)
plt.plot(t,icl_ano1)




plt.show()


print(h.area(0.5))

