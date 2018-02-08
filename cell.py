from neuron import h, gui
import matplotlib.pyplot as plt
import numpy as np
from tempCorr import *
from nernst import *
from matplotlib.collections import LineCollection
from plots import *
from run import *

tstop = 10
h.tstop = tstop
h.dt=0.0001
h.celsius = (T-273)




def define_geometry(icc):
	icc.diam = 50
	icc.L=500
	icc.nseg=10
	icc.cm = 25

def insert_mechanisms(icc):

	icc.insert('pas')

	icc.insert('Na')
	icc.nai = nai
	icc.nao = nao
	icc.G_Na_Na= 20		#20
	icc.tau_f_Na_Na = tau_f_Na
	icc.tau_d_Na_Na = tau_d_Na

	icc.insert('nscc')
	icc.G_NSCC_nscc = 12.15	#12.15
	icc.tau_NSCC_nscc = tau_NSCC

	icc.insert('ERG')
	icc.G_ERG_ERG = 2.5	#2.5
	icc.tau_ERG_ERG = tau_ERG 
	icc.ki= ki	#23
	icc.ko= ko	#23

	icc.insert('bk')
	icc.G_bk_bk= 23+T_correction_BK #23+T_correction_BK

	icc.insert('Kb')
	icc.G_Kb_Kb=0.15	#.15

	icc.insert('kv11')
	icc.G_Kv11_kv11=6.3   #6.3
	icc.tau_d_kv11_kv11 = tau_d_kv11
	icc.tau_f_kv11_kv11 = tau_f_kv11

	icc.insert('vddr')
	icc.G_Ca_VDDR_vddr=3	#3
	icc.cao = 2.5
	icc.tau_d_VDDR_vddr = tau_d_VDDR
	icc.tau_f_VDDR_vddr = tau_f_VDDR

	icc.insert('ltype')
	icc.G_Ca_Ltype_ltype=2 #2
	icc.tau_f_Ltype_ltype = tau_f_Ltype
	icc.tau_d_Ltype_ltype = tau_d_Ltype
	icc.tau_f_Ca_Ltype_ltype = tau_f_Ca_Ltype

	icc.insert('pmca')
	icc.J_max_PMCA_pmca= 0.088464 #0.088464


	icc.insert('concyto')
	J_max_leak=0.01
	icc.J_max_leak_concyto= J_max_leak	#0.01(mM/s) 

	icc.insert('conpu')
	icc.J_max_leak_conpu= J_max_leak #0.01(mM/s) 
	icc.Jmax_serca_conpu = 1.8333		# 1.8333	
	icc.J_ERleak_conpu = 1.666667		# 1.666667	

	icc.Jmax_IP3_conpu = 	50000		#50000  (1/s)
	icc.Jmax_NaCa_conpu =  0.05 	#0.05(mM/s)
	icc.Jmax_uni_conpu = 5000 		#5000 mM/s

	icc.insert('ano1')
	icc.g_Ano1_ano1 = 20	#20
	icc.cli = cli
	icc.clo = clo
	icc.ecl = ecl

	icc.insert('cacl')
	icc.G_Cacl_cacl = 10.1	#10.1
	icc.tau_act_Cacl_cacl = tau_act_Cacl



icc = h.Section(name='icc')
icc1 = h.Section(name='icc1')


define_geometry(icc)
define_geometry(icc1)
icc1.nseg = 1
insert_mechanisms(icc)

icc1.insert('pas')
# icc.connect(icc1(0), 1)
  
#icc.Jmax_serca_conpu =Jmax_serca

#vclamp = h.SEClamp(icc(0.5))
 #vclamp.dur1 = tstop
# vclamp.amp1=-50
# vclamp.rs=.00001

h.v_init = -70


v = h.Vector()             # Membrane potential vector
v_e = h.Vector()


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

#print(icc.cao)
#print(icc.cai)
#h.finitialize(-50)
# h.run()
variables = [v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG]
# v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG = runThisAlready(icc,v,ina,t,ik,ica,icl,eca,ins,cai,cao,ek,capui,caeri,cami,cli,clo,oano,ecl,ena,ica_vddr,icl_ano1, icl_cacl, ina_Na, ik_kv11, ik_bk, ica_ltype, i_nscc, ik_Kb, ica_pmca, ik_ERG)

# variables = run_and_record(icc,*variables)
v_e.record(icc1(0.9)._ref_v)
run_and_record(icc,*variables)

#print(eca[1])




lwidths = []
segments = []
lc = []
points=[]
y=[]


print(type(t))
x=t
y.append(-15*np.ones(np.shape(t)))
y.append(-20*np.ones(np.shape(t)))
y.append(-25*np.ones(np.shape(t)))
y.append(-30*np.ones(np.shape(t)))
y.append(-35*np.ones(np.shape(t)))

y.append(-40*np.ones(np.shape(t)))
y.append(-45*np.ones(np.shape(t)))
y.append(-50*np.ones(np.shape(t)))
y.append(-55*np.ones(np.shape(t)))



# lwidths.append(1+ 10*np.abs(ina)/(max(np.abs(ina))) )
# lwidths.append(1+ 10*np.abs(icl)/(max(np.abs(icl))) )
# lwidths.append(1+ 10*np.abs(ica)/(max(np.abs(ica))) )
# lwidths.append( 1+ 10*np.abs(ik)/(max(np.abs(ik))) )

currents=['i_cacl','i_kb','i_bk','i_ERG','na','i_nscc','i_vddr','i_ltype','i_ano1']

lwidths.append(1+ 5*np.abs(icl_cacl)/(max(np.abs(icl_cacl))) )

lwidths.append(1+ 5*np.abs(ik_Kb)/(max(np.abs(ik_Kb))) )

lwidths.append(1+ 5*np.abs(ik_bk)/(max(np.abs(ik_bk))) )
lwidths.append(1+ 5*np.abs(ik_ERG)/(max(np.abs(ik_ERG))) )
lwidths.append(1+ 5*np.abs(ina_Na)/(max(np.abs(ina_Na))) )

lwidths.append(1+ 5*np.abs(i_nscc)/(max(np.abs(i_nscc))) )
lwidths.append(1+ 5*np.abs(ica_vddr)/(max(np.abs(ica_vddr))) )
lwidths.append(1+ 5*np.abs(ica_ltype)/(max(np.abs(ica_ltype))) )
lwidths.append(1+ 5*np.abs(icl_ano1)/(max(np.abs(icl_ano1))) )

# lwidths.append(1+ 10*np.abs(icl)/(max(np.abs(icl))) )
# lwidths.append(1+ 10*np.abs(icl)/(max(np.abs(icl))) )
# lwidths.append(1+ 10*np.abs(ica)/(max(np.abs(ica))) )
# lwidths.append( 1+ 10*np.abs(ik)/(max(np.abs(ik))) )





fig, a = plt.subplots()

# plotting currents
# color = ['xkcd:chartreuse','xkcd:purple','xkcd:teal','xkcd:violet','xkcd:yellow','xkcd:black','xkcd:fuchsia','xkcd:magenta','xkcd:lime']


# for i in range(len(currents)):

# 		points.append( np.array([x, y[i]]).T.reshape(-1, 1, 2))
# 		segments.append(np.concatenate([points[i][:-1], points[i][1:]], axis=1))
# 		lc.append(LineCollection(segments[i], linewidths=lwidths[i],color=color[i],label = currents[i] ))

# 		a.add_collection(lc[i])
# 		# plt.legend(str(i))

# # plt.plot(v)
a.plot(t,v , label = 'v', color= 'red')
a.plot(t,v_e , label = 'v', color= 'blue')
a.set_xlim(0,tstop)
a.set_ylim(-100,0)
a.legend()
# plt.plot(t,v)
fig.show()


plt.show()


