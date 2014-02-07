	Program two dimensional Bidomain Model with LRD model

c	Rodriguez et al.'s implementation on Luo-Rudy Dynamics (LRD) Model 

c	Author: Sunil M Kandel
c	Date: 12/20/2013
c	Purpose: to study the role on the ischemic increase of Ke during acute 
c	myocardial ischemia

c 	Detailed list of equations and model description are providedd in
c 	Circ Res 1991;68:1501-1526
c	Circ Res 1994;74:1071-1096
c	Circ Res 1994;74:1097-1113
c	Circ Res 1995;77:140-152
c	Am J Physiol Heart Circ Physiol 2002; 283:H490-H500
c	Circ Res 1996;79:208-221

c********************************************************************
c       variables to be set by parameter file (TODO)
c********************************************************************
c       pacing interval (ms)
        integer pintvl
c********************************************************************

        integer, parameter :: VMFILE=21
        integer, parameter :: APDFILE=20

        character :: sep=","
	integer k,idt,m,n,nn
	parameter(nn = 101)
	double precision Vcell,Ageo,Acap,Vmyo,Vmito,Vsr, Vnsr,Vjsr,Vcleft
	double precision vcleftinitial, vcleftfinal,vdotmax(nn,nn)
	double precision Vm(nn,nn),vmnew(nn,nn),dvmdt(nn,nn)
	double precision dvmdtnew(nn,nn),dvmdtmax(nn,nn), dvmdtold(nn,nn)
	double precision dt, t, tstim, stimtime,tjsrol(nn,nn),clock(nn,nn)
	integer iclock(nn,nn),iclock1(nn,nn)
	double precision stim(nn,nn),it(nn,nn),istim(nn,nn),S1
	double precision Nai(nn,nn),Nae(nn,nn),Ki(nn,nn),Ke(nn,nn),
     &                   Cai(nn,nn),cae(nn,nn),cmdn(nn,nn),trpn(nn,nn),
     &                   csqn(nn,nn)
	double precision dNai(nn,nn),dNae(nn,nn),dki(nn,nn),dke(nn,nn),
     &				dcai(nn,nn),dcae(nn,nn)
	double precision Naiont(nn,nn),kiont(nn,nn),caiont(nn,nn),
     &                    kiont1(nn,nn)
	common/const/    r1,length, pi, R, frdy, temp,RTF,zNa,Zk,zca

	double precision INa(nn,nn),GNa, ENa(nn,nn),INaH
	double precision am(nn,nn),bm(nn,nn),ah(nn,nn),bh(nn,nn)
	double precision aj(nn,nn),bj(nn,nn),tm0(nn,nn)
	double precision th0(nn,nn),tj0(nn,nn),m0(nn,nn)
	double precision h0(nn,nn),jj0(nn,nn)
	double precision m1(nn,nn),h(nn,nn),jj(nn,nn)

	double precision INas(nn,nn), Inasfinal

	double precision Ilca(nn,nn),Ilcak(nn,nn),Ilcana(nn,nn)
	double precision fca(nn,nn),Ilcatot(nn,nn)
	double precision Ibarca(nn,nn),Ibarcak(nn,nn),Ibarcana(nn,nn)
	double precision ad(nn,nn),bd(nn,nn),af(nn,nn),bf(nn,nn)
	double precision tf0(nn,nn),f0(nn,nn),td0(nn,nn)
	double precision d0(nn,nn),d(nn,nn),f(nn,nn),km2ca,km1ca
	common/const/ km1ca,pca,pna,pk,ycai,ycae,ynai,ynae,yki,yke

	double precision Icat(nn,nn),gcat,eca(nn,nn)
	double precision b(nn,nn),b0(nn,nn),tb0(nn,nn)
	double precision g(nn,nn),g0(nn,nn),tg0(nn,nn)

	double precision Ikr(nn,nn), gkr(nn,nn),ekr(nn,nn),xr(nn,nn) 
	double precision txr0(nn,nn),rr(nn,nn),xr0(nn,nn)

	double precision Iks(nn,nn),Gbarks(nn,nn),Pca1(nn,nn),Eks(nn,nn)
	double precision Xs0(nn,nn),Xs(nn,nn),xs1(nn,nn),xs2(nn,nn)
	double precision xs10(nn,nn),xs20(nn,nn),gbarks1(nn,nn)
	double precision tXs0(nn,nn),txs10(nn,nn),txs20(nn,nn)
	
	double precision Ek(nn,nn),Gk(nn,nn),Ik(nn,nn),pnak
	double precision ax(nn,nn),bx(nn,nn),tx0(nn,nn),x0(nn,nn)
	double precision tx10(nn,nn),xi(nn,nn),x1(nn,nn)

	double precision Gk1(nn,nn),Ek1(nn,nn),ak1(nn,nn),bk1(nn,nn)
	double precision K10(nn,nn),tk10(nn,nn),Ik1(nn,nn)

	double precision Gkp,Ekp(nn,nn),Kp(nn,nn),Ikp(nn,nn)

	double precision Ikatp(nn,nn),gkatp(nn,nn),fatp(nn,nn)
	double precision ekatp(nn,nn),gamakatp(nn,nn),fmatp(nn,nn)
	double precision fnatp(nn,nn)
	double precision khmg(nn,nn),khna(nn,nn),kmatp(nn,nn),hatp(nn,nn)
	double precision mg(nn,nn),atp(nn,nn),adp(nn,nn)
	double precision sigmakatp, rhokatp,ft,fatpfactor,fatpfinal

	double precision INaca(nn,nn)
	double precision knaca,kmna,kmca,ksat,eta

	double precision sigma(nn,nn), sigma1(nn,nn),fNak(nn,nn)
	double precision INak(nn,nn), Ibarnak, kmnai,kmke
	double precision INakmax, finhib,finhibfinal

	double precision Ibarnsna(nn,nn),Ibarnsk(nn,nn),Ensca(nn,nn)
	double precision Insna(nn,nn),Insk(nn,nn),Insca(nn,nn)
	double precision kmnsca, pnsca

	double precision Ipca(nn,nn)
	double precision kmpca,ibarpca

	double precision Gcab,Ecan(nn,nn),Icab(nn,nn)

	double precision Gnab,ENan(nn,nn),INab(nn,nn)
	double precision Iv(nn,nn)

	double precision grelcicr(nn,nn), grelcicrbar(nn,nn)
	double precision tauon, tauoff
	double precision caizero(nn,nn), caitwo(nn,nn), dcaitwo(nn,nn) 
	double precision dcaith,Irelcicr(nn,nn),dcaidt(nn,nn)
	double precision bjsr(nn,nn),cjsr(nn,nn),jsr(nn,nn), djsr(nn,nn)
	double precision nsr(nn,nn), dnsr(nn,nn)
	double precision greljsrol(nn,nn), grelbarjsrol(nn,nn)
	double precision Ireljsrol(nn,nn)
	double precision Ileak(nn,nn), Iup(nn,nn), Itr(nn,nn), tautr
	double precision kmrel,csqnbar,kmcsqn,Iupbar,nsrbar,kmup,kleak
	double precision kmtrpn,kmcmdn,trpnbar,cmdnbar
	double precision catotal(nn,nn),bmyo(nn,nn),cmyo(nn,nn)
	double precision dmyo(nn,nn), gpig(nn,nn),gpig1(nn,nn)

	double precision tdvdtmax,tmax(nn,nn),vpeak(nn,nn)
	double precision iflag4(nn,nn),taudiff,tdotmax(nn,nn)

	nnn = 2

c       write(6,*) ' start program'
c	Cell Geometry
c 	Dimensions:Length(L)=100 micrometer,Radius(r)=11 micrometer
	Length=0.0001
	r1=.000011
	Pi=3.141592

c	Vcell...Cell Volume
c	Vcell = 10000*pi*r1*r1*length
	Vcell = 38e-6
c 	Ageo...Geometric membrane area: (Ageo=2*pi*r**2+2*pi*r*L) (cm**2)
	Ageo=2*pi*r1**2+2*pi*r1*Length
c 	Acap... Capacitive membrane area (cm**2)
c	Acap = Ageo*2
	Acap = 1.534e-4

c 	Vmyo... Myoplasm volume:(Vmyo=Vcell*68%)(ul)
	Vmyo=Vcell*0.68
c	Vmito...Mitochondria volume:(Vmito=Vcell*26%) (uL)
	Vmito=Vcell*0.26
c	Vsr... SR volume: (Vsr=Vcell*6%) (uL)
	Vsr=vcell*0.06
c	Vnsr...NSR volume:(Vnsr=Vcell*5.52%) (uL)
	Vnsr=Vcell*0.0552
c	Vjsr...JSR volume:(Vjsr=Vcell*0.48%) (uL)
	Vjsr=Vcell*0.0048
c	Vcleft... Cleft volume:(Vcleft=(Vcell/88%)*12%) (uL)
c	Vcleft=Vcell*0.12/0.88
	vcleftinitial = 5.182e-6
	vcleftfinal   = 4.094e-6

c**********************************
c	Voltages
c**********************************
c	Vm... Membrane voltage (mV)
c	Vmnew... New voltage(mV)
c	dVmdt... Change in voltage/change in time (mV/ms)
c	dvmdtnew... New dvmdt (mV/ms)
c	dvmdtold... Old dvmdt (mV/ms)
c	dvmdtmax... Max dvmdt
c	vdotmax...  Max dvmdt

c***********************************
c	Time steps
c***********************************
c       idt... steps per time unit (1/ms)
c	dt... time steps (ms)
c	t... total time (ms)
c	taudiff... time constant for diffusion of ions from the interstitial
c                  clefts to the bulk extracellular medium

c***********************************
c	Total current and stimulus
c***********************************
c	stim... Stimulus current (uA/cm^2)
c	tstim... time stimulus is applied (ms)
c	stimtime... time period during which stimulus is applied
c	it... total current (uA/cm^2)
c	is1s... starting time of first stimulus 
c	is1e... ending time of first stimulus
c       pintvl... pacing interval (ms)

c***********************************************************
c 	Standard ionic Concentration
c***********************************************************
c 	Ki... intracellualr K concentration (mM/L)
c	Ke... extracellualr K concentration (mM/L)
c 	Nai... intracellular Na concentration (mM/L)
c	Nae... extracellular Na concentration (mM/L)
c	Cai... intracellualr Ca concentration (mM/L)
c	Cae... extracellular Ca concentration (mM/L)
c	cmdn... Calmodulin buffered Ca concentration (mM/L)
c	trpn...Troponin Buffered Ca concentration (mM/L)
c 	jsr... JSR Ca concentration (mM/L)
c	nsr... NSR Ca concentration (mM/L)
c	csqn... Calsequestrin Buffered Ca Concentration (mM/L)
c
c****************************************************************
c	Some constants
c****************************************************************
c	R...Universal gas constant(j/kmol*K)
	R = 8314.
c	frdy...Faraday's constant(C/mol)
	frdy = 96485.
c	temp...Temperature(K)
	temp = 310.
c	RT/F = RTF
	RTF  = 26.7

c	Ion Valance
c	zca...Ca valence electron
c	zna...Na valence electron
c	zk... K valence electron
	zCa = 2.
	zNa = 1.
	zK  = 1.

c************************************************************************
c	Fast Sodium Current
c
c 	ina... Fast sodium current (uA/uF)
c	gna... Max. conductance of the Na channel (millisiemens/uF)
c	ena... Reversal potential of Na (mV)
c	am...  Na alpha-m opening rate constant (ms^-1)
c	bm...  Na beta-m closing rate constant (ms^-1)
c	ah...  Na alpha-h opening rate constant (ms^-1)
c	bh...  Na beta-h closing rate constant (ms^-1)
c	aj...  Na alpha-j opening rate constant (ms^-1)
c	bj...  Na beta-j closing rate constant (ms^-1)
c	m1...  Na activation
c	h...   Na inactivation
c	jj...  Na inactivation

c***********************************************************************
c	Ischemic-activated Na inward current

c	INas... New ischemia-activated Na inward current
c	Inasfinal...Final value of INas
c**********************************************************************
c	Current through L-type Ca Channel
c
c	ilca...    Ca current through L-type Ca channel (uA/uF)
c	ilcana...  Na current through L-tye Ca channel  (uA/uF)
c	ilcak...   K current through L-type Ca channel (uA/uF)
c	ilcatot... Total current through L-tye Ca channel (uA/uF)
c	ibarca...  Max.Ca current through Ca channel (uA/uF)
c	ibarcak... Max. K current through Ca channel (uA/uF)
c	ibarcana...Max. Na current through Ca channel (uA/uF)
c	fca...     Ca dependent inactivation gate of the L-type Ca channel
c	d...       Voltage dependent activation gate
c	d0...      Steady-state value of activation gate d
c	td0...     Time constant of gate d (ms^-1)
c	f...       Voltage dependent inactivation gate
c	f0...      Steady-state value of inactivation gate f
c	tf0...     Time constant of gate f (ms^-1)
c	Km2ca...   Half saturation concentration of Ca channel (mM/L)
c 	pca...     permiability of membrane to Ca (cm/s)
c	ycai...    Activity coefficient of intracellular Ca ion
c	ycae...    Activity coefficient of extracellular Ca ion
c	pna...     permiability of membrane to Na (cm/s)
c	ynai...    Activity coefficient of intracellular Na ion
c	ynae...    Activity coefficient of extracellular Na ion
c	pk...      permiability of membrane to K (cm/s)
c	yki...     Activity coefficient of intracellular K ion
c	yke...     Activity coefficient of extracellular K ion

c	constants for current through L-type Ca channel
	Km2Ca = 0.0006	
	PCa   = 0.00054
	PNa   = 0.000000675
	PK    = 0.000000193

c	y represents gamma here	
	ycai = 1.
	yCae = 0.341
	yNai = 0.75
	yKi  = 0.75
	yNae = 0.75
	yke  = 0.75

c******************************************************************
c	Current through T-type Ca Channel
c
c	Icat... Ca current through T-type Ca channel (uA/uF)
c	gcat... Max. conductance of the T-type Ca channel (uS/uF)
c	b0...   Voltage dependent activation gate
c	b...    steady state value of activation gate b0.
c	g0...   Voltage dependent inactivation gate
c	g...    Steady state value of inactivation gate g0
c	tb0...  Time constant of gate b (ms^-1)
c	tg0...  Time constant of gate g (ms^-1)

c******************************************************************
c	Rapidly Activating Potassium Current (Ikr)
c
c	Ikr... Rapidly activating K current (uA/uF)
c	gkr... Channel conductance of Rapidly Activating K Current (uS/uF)
c	ekr... Reversal potential of rapidly activating K current (mV)
c	Xr...  Rapidly activating K time-dependent activation gate
c	xr0... Steady state value of inactivation gate xr
c	txr0...Time constant of gate Xr (ms^-1)
c	rr...  K time-independent inactivation

c*******************************************************************
c	Slow Component of the Delayed Rectifier K current (Iks)
c
c	Iks... Slow component of the delayed rectifier K current
c	Gbarks... Max. Conductance of slowly activating K current
c	pca1...
c	Xs...  Slow activating k time-dependent inactivation gate
c	Xs0... Steady-state value of inactivation gate Xs
c	txs0...time constant of gate Xs
c	Eks... Reversal potential of slowly activating K current (mV)

c*******************************************************************
c	Time-dependent Potassium current
c
c	Ek...Reversal potential of K ion (mV)
c	Gk...Max. conductance of K ion channel (mS/uF)
c	Ik...time dependent K current (uA/uF)
c	ax...K alpha-x rate constant (ms^-1)
c	bx...K beta-x rate consant (ms^-1)
c	x... activation gate of IK
c	xi...inactivation gate of Ik

	Pnak = 0.01833

c******************************************************************
c	Time-independent Potassium current
c
c	Ek1... Reversal potential of K1 ion (mV)
c	Gk1... Max. conductance of K1 ion channel (mS/uF)
c	Ik1... time independent K current (uA/uF)
c	ak1... K alpha-K1 rate constant (ms^-1)
c	bk1... K beta-k1 rate constant (ms^-1)
c	k10... steady state value of K1 inactivation
c	tk10...time constant for gate K1 (ms^-1)
		
c******************************************************************
c	Plateau Potassium current: Ikp
c	
c	Gkp... channel conductance of plateau K current (mS/uF)
c	Ekp... Reversal potential of plateau K current (mV)
c	Kp...  K plateau factor
c	Ikp... Plateau K current (uA/uF)

c******************************************************************
c	ATP sensitive K channel Current (Ikatp) 
c
c	****Used by Shaw-Rudy (1997)******
c
c	Ikatp... ATP-sensitive K current (uA/uF)
c	Ekatp... K reversal potential (mV)
c	gkbaratp... Conductance of the ATP-sensitive K Channel (mS/uF)
c	gkatp... Maximum conductance of the ATP-sensitive K channel (mS/uF)
c	patp...  Percentage availibility of open channels
c	natp...  K dependence of ATP-sensitive K current
c	Nicholsarea... Nichol's area (cm^2)
c	ATPi...  Intracellular ATP concentration (mM)
c	hatp...  Hill coefficient
c	katp...  Half-maximum saturation point of ATP-sensitive K current (mM)
c******************************************************************************
c	ATP-Sensitive Potassium Current (IKatp)
c
c	****(used by Rodriguez et al (2002)*****
c
c	gkatp... conductance of a single fully open ATP channel
c	gamakatp... Unitary conductance in the abscence of Nai and Mgi
c	fmatp... inward rectification caused by Mgi
c	khmg... half-maximum saturation constant
c	fnatp... inward rectification caused by Nai
c	khna... half-maximum saturation constant
c	ft... temperature effect
c	atp... intracellular ATP concentration
c	adp... intracellular ADP concentration
c	fatp... fraction of open Katp channel
c	fatpfinal... final value of fatp after 14 min linear increase
c 	rhokatp... Maximum channel open probability
c	sigmakatp... channel density
c	ekatp... Nernst potential of ATP-sensitive K channels
c	IKatp... Current through ATP sensitive K channels
c*************************************************************************
c	Sodium-Calcium Exchanger: Inaca
c
c	Inaca... sodium-calcium exchange current (uA/uF) 
c	knaca... scaling factor of Inaca (uA/uF)
c	kmna...  Half saturation concentration of Na-Ca channel
c	ksat...  Saturation factor of Inaca at very negative potentials
c	eta...   position of the energy barrier controlling voltage 
c	         dependence of Inaca

	KNaCa = 2000.
	KmNa  = 87.5
	KmCa  = 1.38
	Ksat  = 0.1
	eta   = 0.35
c******************************************************************
c	Sodium Potassium Pump current: Inak

c	finhib... Degree of inihibition of NaK pump
c	finhibfinal... Final value of finhib	
c	sigma...  Extracellular Sodium ion dependence factor of fnak
c	fnak...   Voltage dependance parameter of inak
c	inak...   Nak pump current (uA/uF)
c	Ibarnak...Max. current through Na-K pump (uA/uF)
c	kmnai...  Half saturation concentration of NaK pump (mM)
c	kmke...   Half saturation concentration of Nak pump (mM)

c	Ibarnak = 1.50 (used in LRD 1994)
c	Ibarnak = 2.75 (Used in Rodriguez et al. 2002)
	Ibarnak = 2.75
	KmNai   = 10.
	KmKe    = 1.5
c******************************************************************
c	Sarcolemmal Ca pump: Ipca
	
c	ipca...    Sarcolemmal Ca pump current (uA/uF) 
c	ibarpca... Max Ca current through sarcolemmal Ca pump (uA/uF)
c	kmpca...   Half saturation concentration of Sarcolemmal ca pump
	
	IbarpCa = 1.15
	KmpCa   = 0.0005
c*********************************************************************
c	Ca Background current: Icab
	
c	icab... Ca background current (uA/uF)
c	gcab... Max. conductance of Ca background (mS/uF) 
c	ecan... Nernst potential for Ca (mV)
c*********************************************************************
c	Na Background current: Inab
	
c	inab... Na background current (uA/uF)
c	gnab... Max. conductance of Na background (mS/uF)
c	enan... Nernst potential for Na (mV)
c*********************************************************************
c	Myoplasmic Na ion concentration changes
c	Naiont...total Na Ion flow (uA/uF)
c	dnai...  change in Intracellular Na concentration (mM/L)
c	dnae...  change in extracellular Na concentration (mM/L)

c	Myoplasmic K ion concentration changes
c	Kiont... total K ion flow (uA/uF)
c	dki...   change in intracellular K concentration (mM/L)
c	dke...   change in extracellular K concentration (mM/L)

c	Myoplasmic Ca ion concentration changes
c	Caiont... total Ca ion flow (uA/uF)
c	dcai...   Change in myoplasmic Ca ion concentration (mM/L)
c	dcae...   change in extracellular Ca concentration (mM/L)
c	Catotal...total myoplasmic Ca concentration (mM/L)
c
c**************************************************************************
c	NSR Ca ion Concentration Changes
c
c	dnsr...  Change in Ca ion in the NSR (mM/L)
c	Iup...   Ca uptake from myoplasm to NSR (mM/LmS)
c	Ileak... Ca leakage from NSR to myoplasm (mM/Lms)
c	kleak... Rate constant of Ca leakage from NSR (1/ms)
c	Kmup...  Half saturation concentration of Iup	
c	Iupbar...Max. current throught Iup channel (mM/Lms)
c	nsrbar...Max. Ca ion concentration in NSR (mM/L)

c	Ca uptake and leakage of NSR
	Iupbar = 0.005
	nsrbar = 15.
	kmup   = 0.00092
c***************************************************************************
c	JSR Ca ion concentration changes
c
c	djsr...     Change in Ca ion in the JSR (mM/L)
c	tauon...    time constant of activation of Ca release from JSR (ms)
c	tauoff...   time constant of deactivation of Ca release from JSR (ms)
c	clock...    t = 0 at time of CICR (ms)
c	caizero...  Total myoplasmic Ca at time dvm/dt max. (mM/L)
c	caitwo...   Total myoplasmic ca 2 ms. after dv/dt max. (mM/L)
c	dcaith...   Threshold for external triggering of Ca release from JSR (mM/L)
c 	grelcicrbar...Max. rate constant of Ca release from JSR due to CICR (ms^-1)
c	grelcicr... Rate constant of Ca release from JSR due to CICR (ms^-1)
c	irelcicr... Ca release from jsr to myo. due to CICR (mM/L ms)
c	ireljsrol...Ca release from jsr to myo. due to JSR overload (mM/L ms)
c	bjsr...     b variable for analytical computation of ca in jsr (mM/L)
c 	cjsr...     c variable for analytical computation of ca in jsr (mM/L)

c	buffering parameters
c	cmdnbar...Max. Ca buffered in CMDN (mM/L)
c	trpnbar...Max. Ca buffered in TRPN (mM/L)
c	csqnbar...Max. Ca buffered in CSQN (mM/L)
c	kmcmdn... Equilibrium constant of buffering for CMDN (mM/L)
c	kmtrpn... Equilibrium constant of buffering for TRPN (mM/L)
c	kmcsqn... Equilibrium constant of buffering for CSQN (mM/L)

c	catotal...Total myoplasmic ca concentration (mM/L)
c 	bmyo...   b variable for analytical computation of ca in myoplasm (mM/L)
c	cmyo...   c variable for analytical computation of ca in myoplasm (mM/L)
c	dmyo...   d variable for analytical computation of ca in myoplasm (mM/L)

	trpnbar = 0.07
	cmdnbar = 0.05
	kmtrpn = 0.0005
	kmcmdn = 0.00238

c 	CICR of JSR
	dcaith =  0.00018
	kmrel  =  0.0008
	tauon  = 2.
	tauoff = 2.

c	Ca buffer in JSR and csqn
	csqnbar = 10.
	kmcsqn  = 0.8
c********************************************************************************
c 	translocation of Ca ions from NSR to JSR
c
c	Itr...   tranlocation current of Ca ions from NSR to JSR (mM/Lms)
c	tautr... time const. of Ca transfer from NSR to JSR (ms)
c*****************************************************************************		

c 	Time loop conditions
        idt = 100
	dt  = 1./idt
	nt  = 300*3600
	tstim = 10

	is1s = 0
	is1e = 500
        pintvl = 500

        S1     =  -14
        S2     =  500*0

c	Output parameters

	ioutput = 100
        icount = ioutput-1        	
        iout = 22
	
c 	Initial Ion concentrations (in mM)
	i = 1
	j = 1

	      nai(i,j) = 10.5
	      nae(i,j) = 140.0
	      ki(i,j)  = 144.5
	      ke(i,j)  = 5.4
	      cai(i,j) = 0.00012
	      cae(i,j) = 1.8
	      mg(i,j)  = 0.5
	      atp(i,j) = 3.0
	      adp(i,j) = 0.000036
	      INaH     = -0.14

c 	Initial Gate conditions

	      Vm(i,j) = -87.7
	      m1(i,j) = 0.0023
	      h(i,j)  = 0.9917
	      jj(i,j) = 0.9945
	      d(i,j)  = 0.
	      f(i,j)  = 0.9997
	      b(i,j)  = 0.0011
	      g(i,j)  = 0.9929
	      x1(i,j) = 0.0001
	      xi(i,j) = 0.9889
	      xr(i,j) = 0.00015
	      xs(i,j) = 0.0048
	      

c 	Initial Conditions
	      jsr(i,j)  = 1.74
	      nsr(i,j)  = 1.74
	      trpn(i,j) = 0.014
	      cmdn(i,j) = 0.0024
	      csqn(i,j) = 6.86
	      iclock(i,j) = 0
	      clock(i,j)  = 0.
	      tjsrol(i,j) = 0	  
	      dvmdtold(i,j)=0
	      iclock1(i,j) = 0
	      	    
c       output files
        open(APDFILE,file='apd.csv',status='unknown')
        open(VMFILE,file='vm.csv',status='unknown')
c        write(VMFILE, '(A15,A1,A15,A1,A15)') 'time_ms',sep,'voltage_mV',
c     &         sep,'dVdt_mVms'

        iflag2 = 0
        iflag3 = 0

c***********************************************************************
c	Time loop
	dO  k = 1,nT

	j = 1
	i = 1

c********************************************************************
c       stimulus code
c********************************************************************
	istim(i,j) = 0

c       at the start of stimulus
        if(MOD((k-is1s),pintvl*idt).eq.0) then
          vrest = vmnew(i,j)
          iflag2 = 0
          iflag3 = 0
        end if

c       during the stimulus
        if(k.gt.is1s.and.MOD(k-is1s,pintvl*idt).lt.(is1e-is1s)) then
          clock(i,j) = 0
          iclock1(i,j) = 0
          istim(i,j) = S1
        endif
c********************************************************************

	if (k.lt.60*50000) then
	    taudiff = 1000
	    fatpfactor = 0
	    finhib = 0
	    vcleft = vcleftinitial

	else if (k.gt.60*50000.and.k.lt.(60*50000+14*6e6)) then
	    taudiff = 100e6
	    fatpfactor =((k-60*50000)/(14*6e6-60*50000))*fatpfinal
	    finhib = ((k-60*50000)/(14*6e6-60*50000))*finhibfinal
	    vcleft = vcleftinitial-((k-60*50000)/(14*6e6-60*50000))*
     &                    (vcleftinitial-vcleftfinal)
	    
	else 
	    taudiff = 100e10
	    fatpfactor = fatpfinal
	    finhib = finhibfinal
	    vcleft = vcleftfinal
	   
	end if

	if (k.lt.300*50000) then
	    Inasfinal = 0
	else if (k.gt.300*50000.and.k.lt.(300*50000+12*6e6)) then
 	    Inasfinal = ((k-300*50000)/(14.5*6e6-300*50000))*(-1.2)
	else 
	    Inasfinal = -1.2
 	end if

c****************************************************************************
c	Calculation of different Ionic currents in the sarcolemma
c****************************************************************************
c	Fast sodium current: INa
	GNa = 16.
	ENa(i,j)= RTF*log(nae(i,j)/nai(i,j))
	am(i,j) = 0.32*(Vm(i,j)+47.13)/(1.-exp(-0.1*(Vm(i,j)+47.13)))
        bm(i,j) = 0.08*exp(-Vm(i,j)/11.)

	if(Vm(i,j).ge.-40.) then
         ah(i,j) = 0
	 bh(i,j) = 1./(0.13*(1.+exp((Vm(i,j)+10.66)/(-11.1))))
	 aj(i,j) = 0
	 bj(i,j) = 0.3*exp(-2.535e-7*Vm(i,j))/(1.+exp(-0.1*
     &          (Vm(i,j)+32.)))
	else
	 ah(i,j) = 0.135*exp((80.+Vm(i,j))/(-6.8))
	 bh(i,j) = 3.56*exp(0.079*Vm(i,j))+3.1e5*exp(0.35*Vm(i,j))
	 aj(i,j) = (-1.2714e5*exp(0.2444*Vm(i,j))-3.474e-5*exp(-0.04391*
     &            Vm(i,j)))*(Vm(i,j)+37.78)/(1.+exp(0.311*
     &            (Vm(i,j)+79.23)))
	 bj(i,j) = 0.1212*exp(-0.01052*Vm(i,j))/(1.+exp(-0.1378*
     &           (Vm(i,j)+40.14)))
         end if
	
	tm0(i,j) = 1./(am(i,j)+bm(i,j))
	th0(i,j) = 1./(ah(i,j)+bh(i,j))
	tj0(i,j) = 1./(aj(i,j)+bj(i,j))

	if (tm0(i,j).le.dt) tm0(i,j) = dt
	if (th0(i,j).le.dt) th0(i,j) = dt
	if (tj0(i,j).le.dt) tj0(i,j) = dt	
        
	m0(i,j) = am(i,j)*tm0(i,j)
	h0(i,j) = ah(i,j)*th0(i,j)
	jj0(i,j)= aj(i,j)*tj0(i,j)

	m1(i,j) = m1(i,j)  + (dt/tm0(i,j))*(m0(i,j)- m1(i,j))
	h(i,j)  = h(i,j)   + (dt/th0(i,j))*(h0(i,j)- h(i,j))
	jj(i,j) = jj(i,j)  + (dt/tj0(i,j))*(jj0(i,j)-jj(i,j))


	INa(i,j) = GNa*m1(i,j)**3*h(i,j)*jj(i,j)*(Vm(i,j)-ENa(i,j))

c 	Calculation of INas

	INaS(i,j) = Inasfinal
c	INaS(i,j) = 0

c	currents through the L-type Calcium Channel

	d0(i,j) = 1./(1.+exp(-(Vm(i,j)+10.)/6.24))

	td0(i,j) = d0(i,j)*(1.-exp(-(Vm(i,j)+10.)/6.24))/
     &	     (0.035*(Vm(i,j)+10.))
	  

	f0(i,j) = (1./(1.+exp((Vm(i,j)+32.)/8.))) + (0.6/
     &              (1.+exp((50.-Vm(i,j))/20.)))
	tf0(i,j) = 1./(0.0197*exp(-(0.0337*(Vm(i,j)+10.))**2)+0.02)


	d(i,j) = d(i,j) + (dt/td0(i,j))*(d0(i,j) - d(i,j))
	f(i,j) = f(i,j) + (dt/tf0(i,j))*(f0(i,j) - f(i,j))
                 

	IbarCa(i,j) =  PCa*ZCa**2*Vm(i,j)*frdy/RTF*(ycai*cai(i,j)*
     &            exp(zca*Vm(i,j)/RTF)-yCae*cae(i,j))/
     &            (exp(zca*Vm(i,j)/RTF)-1.)

	IbarCaK(i,j) = PK* ZK**2 *Vm(i,j)*frdy/RTF*(yKi *ki(i,j) *
     &            exp(zK*Vm(i,j)/RTF)-yKe*ke(i,j))/
     &            (exp(zK*Vm(i,j)/RTF)-1.)

	Ibarcana(i,j) = PNa*ZNa**2*Vm(i,j)*frdy/RTF*(yNai*Nai(i,j)*
     &            exp(zNa*Vm(i,j)/RTF)-yNae*nae(i,j))/
     &            (exp(zNa*Vm(i,j)/RTF)-1.)

c	km2ca = 0.0006
	fca(i,j) = 1./(1.+(cai(i,j)/km2ca))
 
	Ilca(i,j)   = d(i,j)*f(i,j)*fca(i,j)*Ibarca(i,j)
	Ilcana(i,j) = d(i,j)*f(i,j)*fca(i,j)*Ibarcana(i,j)
	Ilcak(i,j)  = d(i,j)*f(i,j)*fca(i,j)*Ibarcak(i,j)

	IlCatot(i,j) = Ilca(i,j) + Ilcana(i,j) + Ilcak(i,j)

c	Current through T-type Ca Channel

	b0(i,j) = 1/(1 + exp(-(vm(i,j)+14.)/10.8))
	tb0(i,j) = 3.7 + 6.1/(1+exp((vm(i,j)+25.)/4.5))

	g0(i,j) = 1/(1 + exp((vm(i,j)+60)/5.6))
	If (vm(i,j).le.0) then
	   tg0(i,j) = -0.875*vm(i,j) + 12.
	else
	   tg0(i,j) = 12.
	end if

	b(i,j) = b(i,j) + (dt/tb0(i,j))*(b0(i,j)-b(i,j))
	g(i,j) = g(i,j) + (dt/tg0(i,j))*(g0(i,j)-g(i,j))

	gcat = 0.05
	eca(i,j)  = (RTF/2)*log(cae(i,j)/cai(i,j))

	icat(i,j) = gcat*b(i,j)**2*g(i,j)*(Vm(i,j) - eca(i,j))
	
c	Rapidly Activating Potassium Current (Ikr)

	gkr(i,j) = 0.02614*sqrt(ke(i,j)/5.4)
	ekr(i,j) = RTF*log(ke(i,j)/ki(i,j))

	Xr0(i,j) =1/(1 + exp(-(Vm(i,j) + 21.5)/7.5))
	txr0(i,j)=1/(0.00138*(Vm(i,j)+14.2)/(1-exp(-0.123*(vm(i,j)+14.2)))
     &         +  0.00061*(vm(i,j)+38.9)/(exp(0.145*(vm(i,j)+38.9))-1))

	Xr(i,j) = Xr(i,j) + (dt/txr0(i,j))*(Xr0(i,j)-Xr(i,j))

	rr(i,j) = 1/(1 + exp((Vm(i,j)+9)/22.4))


	Ikr(i,j) = gkr(i,j)*Xr(i,j)*rr(i,j)*(vm(i,j) - ekr(i,j))

c	Slow Component of the Delayed Rectifier K current (Iks)

	Pca1(i,j)   = -log10(cai(i,j)) + 3.
	Gbarks(i,j) = 0.057 + 0.19/(1 + exp((-7.2 + pca1(i,j))/0.6))

	Eks(i,j)=RTF*log((Ke(i,j)+0.01833*Nae(i,j))/
     &                             (Ki(i,j)+0.01833*Nai(i,j)))

	Xs0(i,j)  = 1./(1 + exp(-(vm(i,j) - 1.5)/16.7))

	txs10(i,j)=1./(0.0000719*(vm(i,j)+30)/
     &                       (1-exp(-0.148*(vm(i,j)+30))) +
     &           0.000131*(vm(i,j)+30)/(exp(0.0687*(vm(i,j)+30))-1))
	Xs(i,j) = Xs(i,j) + (dt/txs10(i,j))*(Xs0(i,j)-Xs(i,j))

	Iks(i,j) = Gbarks(i,j)*Xs(i,j)*xs(i,j)*(Vm(i,j) - Eks(i,j))

c	Time independent Potassium current: Ik1 (time independent)

	Gk1(i,j) = 0.75*sqrt(ke(i,j)/5.4)
	EK1(i,j) = RTF*log(ke(i,j)/Ki(i,j))

	aK1(i,j) = 1.02/(1.+exp(0.2385*(Vm(i,j)-EK1(i,j)-59.215)))

	bK1(i,j) = (0.49124*exp(0.08032*(Vm(i,j)-EK1(i,j)+5.476))+
     &          exp(0.06175*(Vm(i,j)-EK1(i,j)-594.31)))/
     &          (1.+exp(-0.5143*(Vm(i,j)-EK1(i,j)+4.753)))

	tK10(i,j) = 1./(aK1(i,j)+bK1(i,j))
	K10(i,j)  = aK1(i,j)*tK10(i,j)

	IK1(i,j) = Gk1(i,j)*k10(i,j)*(Vm(i,j)-EK1(i,j))

c	Plateau Potassium current: IKp (time independent)

c	GKp      = 0.0183
	Gkp	 = 0.00552
	EKp(i,j) = EK1(i,j)
	Kp(i,j)  = 1./(1.+exp((7.488-Vm(i,j))/5.98))
	IKp(i,j) = GKp*Kp(i,j)*(Vm(i,j)-EKp(i,j))


c 	Formulation of Ikatp proposed by Ferrero et al. 

	sigmakatp = 3.8
	rhokatp   = 0.91
	fatpfinal = 0.008
	
	gkatp(i,j) = gamakatp(i,j)*fmatp(i,j)*fnatp(i,j)*ft
	
	gamakatp(i,j) = 35.375*(ke(i,j)/5.4)**0.24

	fmatp(i,j) = 1/(1 + (Mg(i,j)/khmg(i,j)))
	khmg(i,j) = (0.65/sqrt(ke(i,j)+5))*exp((-2*0.32/RTF)*vm(i,j))

	fnatp(i,j) = 1/(1 + (Nai(i,j)/khna(i,j))**2)
	khna(i,j) = 25.9*exp((-0.35/RTF)*vm(i,j))

	ft = 1.3**(Temp-309)/10

	
	kmatp(i,j) = (35.8 + 17.9*adp(i,j)**0.256)
	hatp(i,j) = 1.3 + 0.74*exp(-0.09*adp(i,j))
	fatp(i,j) =fatpfactor/(1 + (atp(i,j)/kmatp(i,j))**hatp(i,j))

	ekatp(i,j) = RTF*log(ke(i,j)/ki(i,j))

	Ikatp(i,j) = sigmakatp*gkatp(i,j)*rhokatp*fatp(i,j)*
     &			(vm(i,j)-ekatp(i,j))

c	Ikatp(i,j) = 0

c	Na-Ca exchange current:INaCa

	INaCa(i,j)= KNaCa*(1./(KmNa**3+nae(i,j)**3))*(1./(KmCa+cae(i,j)))*
     &             (1./(1.+Ksat*exp((eta-1.)*Vm(i,j)*(1./RTF))))*
     &             (exp(eta*Vm(i,j)*(1./RTF))*nai(i,j)**3*cae(i,j) - 
     &             exp((eta-1.)*Vm(i,j)*(1./RTF))*nae(i,j)**3*cai(i,j))
c	INaCa(i,j) = 0

c	Na-K Pump current: INaK (time independent)

	finhibfinal = 0.35
	Inakmax = IbarNak*(1 - finhib)

	sigma(i,j) = 0.14286*(exp(nae(i,j)/67.3)-1.)
	Sigma1(i,j)=(1./(1.+(KmNai/Nai(i,j))**1.5))*ke(i,j)/(ke(i,j)+KmKe)

	fNak(i,j) = 1./(1.+0.1245*exp(-0.1*Vm(i,j)/RTF)+
     &           0.0365*sigma(i,j)*exp(-Vm(i,j)/RTF))
	INaK(i,j) = Inakmax*fNaK(i,j)*Sigma1(i,j)

c	Sarcolemmal Calcium pump current: IpCa (time independent)
	
	IpCa(i,j) = IbarpCa*(cai(i,j)/(KmpCa+cai(i,j)))

c	Calcium Background Current:ICab (time independent)

	Gcab = 0.003016
	EcaN(i,j) = (RTF/2.)*log(cae(i,j)/cai(i,j))
	ICab(i,j) = GCab*(Vm(i,j)- ECaN(i,j))

c	Sodium Background Current: INab (time independent)

	GNab = 0.00141
	ENan(i,j) = ENa(i,j)
	INab(i,j) = GNab*(Vm(i,j)- ENan(i,j))

c	Total time independent current: Iv
	Iv(i,j) = IK1(i,j) + IKp(i,j) + IpCa(i,j) + INab(i,j)
     &            + ICab(i,j) + INaK(i,j)

c	Calculation of total independent currents

	Naiont(i,j) = INa(i,j) + INab(i,j) + Ilcana(i,j) + 3.*INak(i,j) 
     &                + Insna(i,j) + 3.*INaca(i,j) + Inas(i,j)
	Kiont(i,j)  = Ikr(i,j) + Iks(i,j) + Ik1(i,j) + IKp(i,j) +  
     &                IlcaK(i,j)+ InsK(i,j) - 2.*INak(i,j) + Ikatp(i,j)
	caiont(i,j) = ilca(i,j) + Ipca(i,j) + Icab(i,j) - 2.*INaca(i,j)
     &                + Icat(i,j)

	it(i,j) = Naiont(i,j) + Kiont(i,j) + Caiont(i,j)


c********************************************************************
c	intracellualr & extracellular ion concentrations 
c********************************************************************

	dNai(i,j) = -dt*((Naiont(i,j)+INaH)*Acap)/(Vmyo*zNa*frdy)
	nai(i,j)  =  dNai(i,j) + nai(i,j)

	dKi(i,j) = -dt*(Kiont(i,j)*Acap)/(Vmyo*zk*frdy)
	ki(i,j)  =  dki(i,j) + ki(i,j)

	dNae(i,j) = dt*((Naiont(i,j)+INaH)*Acap)/(Vcleft*frdy)
     &              -dt*(Nae(i,j)-140)/taudiff
	nae(i,j)  =  dnae(i,j)+nae(i,j)

	dKe(i,j) = dt*(Kiont(i,j)*Acap)/(Vcleft*frdy)
     &              -dt*(Ke(i,j)-5.4)/taudiff
	Ke(i,j)  =  dKe(i,j)+Ke(i,j)

	dCae(i,j) = dt*(caiont(i,j)*Acap)/(2*Vcleft*frdy)
     &              - dt*(Cae(i,j)-1.8)/taudiff
	cae(i,j)  =  dCae(i,j)+Cae(i,j)

c*********************************************************************
c	Main formula for Vm
c*********************************************************************
	Vmnew(i,j) = vm(i,j) - (it(i,j) + istim(i,j))*dt
	dvmdtnew(i,j) = (vmnew(i,j) - vm(i,j))/dt	
c*********************************************************************

c*********************************************************************
c	APD 90
c*********************************************************************
	if(iflag2.eq.0) then
	  if(k.gt.is1s) then
	    if(vmnew(i,j).lt.vm(i,j)) then
              tmax(i,j) = k*dt
              vpeak(i,j)= vm(i,j)
	      iflag2=1
	    end if
	  end if
	end if
 	
	
	if((iflag2.eq.1).and.(iflag3.eq.0)) then
	  if(k.gt.is1e+400) then
            if(vmnew(i,j).lt.(vrest + (vpeak(i,j)-vrest)*0.1)) then
              Apd=k*dt-tmax(i,j)
	      iflag3=1
	      write(APDFILE,'(F15.2,A1,ES15.3E3,A1,ES15.3E3)')
     &              tmax(i,j),sep,apd,sep,vdotmax(i,j)
	    end if
	  end if
	end if
c*********************************************************************
	
	
c*************************************************************************
c
c      iclock(i,j) = 0 means before the 2 ms integrating time
c      iclock(i,j) = 1 means during the 2 ms integrating time
c      iclock(i,j) = 2 means after the 2 ms integrating time
c
c      note: do not even look for Vdotmax if before stimulus (artifact)
c
	If(k.lt.is1s) then
         continue
        else
        If(dvmdtnew(i,j).lt.dvmdtold(i,j).and.dvmdtnew(i,j)
     &                .gt.100.and.iclock(i,j).eq.0) then
            caizero(i,j) = cai(i,j)+ cmdn(i,j) + trpn(i,j) + dcai(i,j)
     &                     + csqn(i,j) + jsr(i,j) + nsr(i,j) 

	    vdotmax(i,j) = dvmdtnew(i,j)
	  
	    tdvdtmax = k*dt
            iclock(i,j) = 1
            clock(i,j) = 0.
c           write(6,*) ' dvdtmax now ',k*dt,i,j

          end if
         end if
	
	 if (clock(i,j).gt.2.0.and.iclock(i,j).eq.1) then
	   caitwo(i,j) = cai(i,j)+ cmdn(i,j) + trpn(i,j) + dcai(i,j)
     &                   + csqn(i,j) + jsr(i,j) + nsr(i,j)

           dcaitwo(i,j) = caitwo(i,j) - caizero(i,j)
	
	    iclock(i,j) = 0
	
           If (dcaitwo(i,j).gt.dcaith) then
              grelcicrbar(i,j)= 60.
c             write(6,*) ' cicr!!!!!!', k*dt,i,j

           else 
             grelcicrbar(i,j)=0.
           end if
	   
	      dvmdtold(i,j) = 0    
	      clock(i,j) = 0
	
           end if

              clock(i,j)=clock(i,j)+dt

	 grelcicr(i,j)=grelcicrbar(i,j)*((dcaitwo(i,j) - 0.00018)/
     &                        (kmrel + dcaitwo(i,j) - 0.00018))*
     &        (1.-exp(-clock(i,j)/tauon))*exp(-clock(i,j)/tauoff)

	 Irelcicr(i,j) = grelcicr(i,j)*(jsr(i,j)-cai(i,j))

c	 Irelcicr(i,j) = 0
	 
c	Ca uptake and leakage of NSR
	
	kleak = Iupbar / nsrbar
	Ileak(i,j) = kleak * nsr(i,j)
	Iup(i,j) = Iupbar * cai(i,j) / (cai(i,j)+kmup)

c     	Iup(i,j) = 0

c 	Translocation of Ca ions from NSR to JSR
	tautr = 180. 
	Itr(i,j) = (nsr(i,j)-jsr(i,j))/tautr

c	Calculation of calcium concentration in NSR
	dnsr(i,j) = dt*(Iup(i,j) - Ileak(i,j) - Itr(i,j)*vjsr/vnsr)
	nsr(i,j) = nsr(i,j) + dnsr(i,j)

c	Calculation of calcium concentration in JSR
c 	Analytical expression for Ca2+ buffering in the junctional 
c	sarcoplasmic reticulum (JSR) and in the cytosol under steady 
c	state conditions taken from appendix 1 of Zeng et al 1995 paper

	csqn(i,j) = csqnbar * (jsr(i,j)/(jsr(i,j) + kmcsqn))
	djsr(i,j) = dt*(itr(i,j) - Irelcicr(i,j))


	bjsr(i,j) = csqnbar - csqn(i,j) -djsr(i,j) - jsr(i,j) + kmcsqn
	cjsr(i,j) = kmcsqn * (csqn(i,j) + djsr(i,j) + jsr(i,j))
	jsr(i,j) = (sqrt(bjsr(i,j)**2 + 4. * cjsr(i,j)) - bjsr(i,j))/2.


	dcai(i,j) = -dt*(((caiont(i,j)*Acap)/(Vmyo*zca*frdy))
     &	   - (Irelcicr(i,j)*vjsr/vmyo) 
     &     + ((Iup(i,j) - Ileak(i,j))*vnsr/vmyo))


	dcaidt(i,j) = dcai(i,j)/dt

	trpn(i,j) = trpnbar * (cai(i,j)/(cai(i,j) + kmtrpn))
	cmdn(i,j) = cmdnbar * (cai(i,j)/(cai(i,j) + kmcmdn))

	catotal(i,j) = trpn(i,j) + cmdn(i,j) + dcai(i,j) + cai(i,j)
 	bmyo(i,j) = cmdnbar + trpnbar - catotal(i,j) + kmtrpn + kmcmdn
	cmyo(i,j) = (kmcmdn*kmtrpn) - (catotal(i,j)*(kmtrpn + kmcmdn))
     &              + (trpnbar*kmcmdn) + (cmdnbar*kmtrpn)
	dmyo(i,j) = -kmtrpn * kmcmdn * catotal(i,j)

	gpig(i,j) = sqrt(bmyo(i,j)*bmyo(i,j)-3.*cmyo(i,j))


	gpig1(i,j) = cos(acos((9.*bmyo(i,j)*cmyo(i,j)-2.*bmyo(i,j)**3
     &     -27.*dmyo(i,j))/(2.*(bmyo(i,j)**2 - 3.*cmyo(i,j))**1.5))/3.)

	cai(i,j) = (2.*gpig(i,j)*gpig1(i,j)/3.) - (bmyo(i,j)/3.)


	vm(i,j) =  vmnew(i,j)
	dvmdtold(i,j) = dvmdtnew(i,j)

	 icount=icount+1
           if(icount.eq.ioutput) then
              iout=iout+1
              icount=0
              time=k*dt
	      i = 1
	      j = 1

c            write(26,*) time,Vmnew(i,j),INa(i,j),Ilca(i,j),Ik(i,j),
c    &              Iv(i,j),Inaca(i,j),Ik1(i,j), Ikp(i,j),
c    &              Inak(i,j),Ipca(i,j),INab(i,j),Icab(i,j),
c    &              cai(i,j),jsr(i,j),nsr(i,j),trpn(i,j),cmdn(i,j),
c    &              csqn(i,j),Irelcicr(i,j),Iup(i,j),Ileak(i,j),
c    &              itr(i,j),fca(i,j),dvmdtnew(i,j),It(i,j)

c             write(VMFILE,'(F15.2,A1,ES15.3E3,A1,ES15.3E3)') time,sep,
c     &              vmnew(i,j),sep,dvmdtnew(i,j)
	   end if

c	End of time loop
	end do

        close(VMFILE)
        close(APDFILE)

	Stop 
	End


	

	
	
