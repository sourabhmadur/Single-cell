import matplotlib.pyplot as plt


def plotCalcium(t,cai,capui,caeri,cami):

	plt.subplot(4,1,1)
	plt.plot(t,cai)
	plt.ylabel('[Ca] in cytosol(mM)')
	
	plt.subplot(4,1,2)
	plt.plot(t,cami)
	plt.ylabel('[Ca] in Mitochondria(mM)')


	plt.subplot(4,1,3)
	plt.plot(t,capui)
	plt.ylabel('[Ca] in PU(mM)')

	plt.subplot(4,1,4)
	plt.plot(t,caeri)
	plt.ylabel('[Ca] in ER(mM)')

	plt.xlabel('time(s)')

def plotcaconlim(cai,capui,caeri,cami):

	plt.subplot(2,3,1)
	plt.plot(cai,capui)
	plt.xlabel('cai(mM)')
	plt.ylabel('capui(mM)')
	
	plt.subplot(2,3,2)
	plt.plot(cai,caeri)
	plt.xlabel('cai(mM)')
	plt.ylabel('caeri(mM)')


	plt.subplot(2,3,3)
	plt.plot(cai,cami)
	plt.xlabel('cai(mM)')
	plt.ylabel('cami(mM)')



	plt.subplot(2,3,4)
	plt.plot(capui,caeri)
	plt.xlabel('capui(mM)')
	plt.ylabel('caeri(mM)')


	plt.subplot(2,3,5)
	plt.plot(capui,cami)
	plt.xlabel('capui(mM)')
	plt.ylabel('cami(mM)')

	plt.subplot(2,3,6)
	plt.plot(caeri,cami)
	plt.xlabel('caeri(mM)')
	plt.ylabel('cami(mM)')

def mempot(t,v):

	plt.plot(t,v)
	plt.xlabel('time(s)')
	plt.ylabel('Membrane Potential(mV)')



	




