	Program two dimensional Bidomain Model with LRD model

c	Terkildsen et al.'s implementation on Luo-Rudy Dynamics (LRD) Model 

c	Author: Sunil M Kandel
c	Date: 12/20/2013
c	Purpose: to study the role on the ischemic increase of Ke during acute 
c	myocardial ischemia

c 	Detailed list of equations and model description are providedd in
c 	Circ Res 1991;68:1501-1526
c	Circ Res 1994;74:1071-1096
c	Circ Res 1994;74:1097-1113
c	Circ Res 1995;77:140-152
c	Am J Physiol Heart Circ Physiol 2007; 293:3036-3045 (Terkildsen model)
c	Circ Res 1996;79:208-221


	integer m,n,nn
	parameter(nn = 101)
	double precision Vcell,Ageo,Acap,Vmyo,Vmito,Vsr, Vnsr,Vjsr,Vcleft
	double precision vcleftinitial, vcleftfinal,vdotmax
	double precision Vm,vmnew,dvmdt
	double precision dvmdtnew,dvmdtmax, dvmdtold
	double precision dt, t, tstim, stimtime,tjsrol,clock
	integer iclock,iclock1
	double precision stim,it,istim,S1
	double precision Nai,Nae,Ki,Ke,
     &                   Cai,cae,cmdn,trpn,
     &                   csqn
	double precision dNai,dNae,dki,dke,
     &				dcai,dcae
	double precision Naiont,kiont,caiont,
     &                    kiont1
	common/const/    r1,length, pi, R, frdy, temp,RTF,zNa,Zk,zca

	double precision INa,GNa, ENa,INaH
	double precision am,bm,ah,bh
	double precision aj,bj,tm0
	double precision th0,tj0,m0
	double precision h0,jj0
	double precision m1,h,jj

c	double precision INas, Inasfinal

	double precision Ilca,Ilcak,Ilcana
	double precision fca,Ilcatot
	double precision Ibarca,Ibarcak,Ibarcana
	double precision ad,bd,af,bf
	double precision tf0,f0,td0
	double precision d0,d,f,km2ca,km1ca
	common/const/ km1ca,pca,pna,pk,ycai,ycae,ynai,ynae,yki,yke

	double precision Icat,gcat,eca
	double precision b,b0,tb0
	double precision g,g0,tg0

	double precision Ikr, gkr,ekr,xr 
	double precision txr0,rr,xr0

	double precision Iks,Gbarks,Pca1,Eks
	double precision Xs0,Xs,xs1,xs2
	double precision xs10,xs20,gbarks1
	double precision tXs0,txs10,txs20
	
	double precision Ek,Gk,Ik,pnak
	double precision ax,bx,tx0,x0
	double precision tx10,xi,x1

	double precision Gk1,Ek1,ak1,bk1
	double precision K10,tk10,Ik1

	double precision Gkp,Ekp,Kp,Ikp

c	double precision Ikatp,gkatp,fatp
c	double precision ekatp,gamakatp,fmatp
c	double precision fnatp
c	double precision khmg,khna,kmatp,hatp
c	double precision mg,atp,adp
c	double precision sigmakatp, rhokatp,ft,fatpfactor,fatpfinal

	double precision ATP, ADP,fatp,f_bound_adp,f_close_adp,
     &                   fadp,fkatp, gkatp, ekatp,Ikatp

	double precision INaca
	double precision knaca,kmna,kmca,ksat,eta

c	double precision sigma, sigma1,fNak
c	double precision INak1, Ibarnak, kmnai,kmke
c	double precision INakmax, finhib,finhibfinal

	double precision Ibarnsna,Ibarnsk,Ensca
	double precision Insna,Insk,Insca
	double precision kmnsca, pnsca

	double precision Ipca
	double precision kmpca,ibarpca

	double precision Gcab,Ecan,Icab

	double precision Gnab,ENan,INab
	double precision Iv

	double precision grelcicr, grelcicrbar
	double precision tauon, tauoff
	double precision caizero, caitwo, dcaitwo 
	double precision dcaith,Irelcicr,dcaidt
	double precision bjsr,cjsr,jsr, djsr
	double precision nsr, dnsr
	double precision greljsrol, grelbarjsrol
	double precision Ireljsrol
	double precision Ileak, Iup, Itr, tautr
	double precision kmrel,csqnbar,kmcsqn,Iupbar,nsrbar,kmup,kleak
	double precision kmtrpn,kmcmdn,trpnbar,cmdnbar
	double precision catotal,bmyo,cmyo
	double precision dmyo, gpig,gpig1

	double precision tdvdtmax,tmax,vpeak
	double precision iflag4,taudiff,tdotmax

	double precision  CrP,cr_tot,Cr_nak,PH,protons,AMP,
     &                   Pi_max

	double precision Nai_number,ki_number,Cai_number,Nae_number,
     &	                 ke_number,cae_number
	double precision osmol,imperm_i,imperm_e,water_flux,dvmyo

	double precision Inak, Kd_nak_nae1,  Kd_nak_nai1, 
     &                   Nae_bound1,Nai_bound1,Nae_bound2,Nai_bound2,
     &                   Ke_bound, Ki_bound, MgATP_bound,alpha_f_nak_1,
     &                   alpha_f_nak_2, alpha_f_nak_3, alpha_f_nak_4, 
     &			 alpha_r_nak_1, alpha_r_nak_2, alpha_r_nak_3,
     &			 alpha_r_nak_4, sigma_nak, v_cyc_nak,pi_nak,
     &			 pi_nak_tot

	double precision K0d_nak_nae, K0d_nak_nai, Kd_nak_nae2,	Kd_nak_nai2,
     &    	       Kd_nak_ke, Kd_nak_ki, Kd_nak_MgATP, rate_f_nak_1, 
     &                 rate_f_nak_2, rate_f_nak_3, rate_f_nak_4,rate_r_nak_1,
     &		       rate_r_nak_2, rate_r_nak_3, rate_r_nak_4, delta_nak,
     &                 density_nak, Kd_nak_kpi, kd_nak_hpi, Kd_nak_napi

c	nnn = 2

c        write(6,*) ' start program'
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
c	Vcleft=Vcell*(0.12/0.88)
	vcleft = 5.182e-6
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
	dt = 0.01
	nt = 30000*3000
	tstim = 10
c	taudiff = 1000

	is1s = 0
	is1e = 500

        S1     =  -14
        S2     =  500*0

c	Output parameters

	ioutput = 10000
        icount = ioutput-1        	
        iout = 22
	
c 	Initial Ion concentrations (in mM)
	
	      nai = 10.0
	      nae = 140.0
	      ki  = 145
	      ke  = 5.4
	      cai = 0.00013
	      cae = 1.8
c	      mg  = 0.5
              atp = 6.9485
	      amp = 1.8425 e-4
	      crp = 13.2946
	      cr_nak = 8.9054
	      ph  = 7.0953
	      adp = 0.034919
	      pi_nak_tot = 0.8
	      protons    = 8.029e-5

c 	Initial Gate conditions

	      Vm = -85.5
	      m1 = 0.0028
	      h  = 0.9851
	      jj = 0.9906
	      d  = 0.
	      f  = 0.9993
	      b  = 0.0014
	      g  = 0.9887
	      x1 = 0.0001
	      xi = 0.9889
	      xr = 0.00021
	      xs = 0.00558

	
c 	Initial Conditions
	      jsr  = 2.10
	      nsr  = 2.10
	      trpn = 0.016
	      cmdn = 0.0029
	      csqn = 7.25
	      iclock = 0
	      clock  = 0.
	      tjsrol = 0	  
	      dvmdtold=0
	      iclock1 = 0

c	Initial Value
	K0d_nak_nae  = 141.271240234375000
	K0d_nak_nai  = 0.120986891616325920E-07
	Kd_nak_nae2  = 356.009094238281250 
	Kd_nak_nai2  = 154.389923095703125
 	Kd_nak_ke    = 6.80234527587890625
	Kd_nak_ki    = 255.130677771665574
	Kd_nak_MgATP = 2.92630076408386230  
	rate_f_nak_1 = 1664.12402343750000
 	rate_f_nak_2 = 110.425056457519531 
	rate_f_nak_3 = 2314.25756835937500
	rate_f_nak_4 = 462.377563476562500
	rate_r_nak_1 = 264.146662389765481
	rate_r_nak_2 = 14.0371952056884766 
	rate_r_nak_3 = 7599999.50000000000
	rate_r_nak_4 = 1702.52624511718750
	delta_nak    = -0.179529249668121338
	density_nak  = 7450

c	Pi_nak_tot   = 0.8

	Kd_nak_kpi   = 292.0
	kd_nak_hpi   = 6.77
	KD_nak_napi  = 224.0
	      	    
c       output files
        open(20,file='time.dat',status='unknown')
        open(21,file='vm.dat',status='unknown')
        open(22,file='ca.dat',status='unknown')

c***********************************************************************
c	Time loop
	dO  k = 1,nT

	Istim = 0
	  
c***************************************************************	
	  do n = 1, 3500

            if(k.eq.(is1s+(n-1)*30000)) then
                vrest = vmnew
                end if

            if(k.gt.(is1s+(n-1)*30000).and.
     &         k.lt.(is1e+(n-1)*30000)) then

	    iflag1 = 0
	    iflag2 = 0
	    iflag3 = 0
            clock = 0
	    iclock1 = 0

    	    Istim = S1
	        
           end if
	  end do
c*************************************************************
c	Ischemic Metabolite Concentration Timecourses

	if (k.gt.200*30000) then

	t1 = (k-200*30000)*dt/60000
	t2 = t1*t1
	t3 = t2*t1

	ATP = 6.549e-4*t3 - 0.02305*t2 - 0.104837*t1 + 6.9485
                
        CrP = -0.01259 + 12.339*EXP(-0.92559*t1)+
     &                  0.96819*EXP(-0.078496*t1)

c       Cr_tot = 22.2
        Cr_nak = 22.2 - CrP

	pH = 6.12457 - 0.56697*EXP(-0.19015*t1) +
     &                  1.5377*EXP(-0.18462*t1)
        protons = 1.0e3*10.0**(-pH)

	ADP = (ATP*Cr_nak)/(CrP*protons*1.66e6)
        AMP = ADP*ADP*1.05/ATP


	Pi_nak_tot = 35.0055 - (3*ATP + 2*ADP + AMP + CrP)

	end if

c*****************************************************************************
c	VOLUME REGULATION
	
	Nai_number = Nai*vmyo
	Ki_number  = Ki*Vmyo
	cai_number = Cai*Vmyo

	Nae_number = Nae*Vcleft
	Ke_number  = Ke*Vcleft
	Cae_number = Cae*Vcleft

        if (k.gt. 200*30000) then
c	    imperm_i = 4.0e-3
            osmol = 310.0 + ((k-6e6)/(15*6e6))*15
            imperm_i = vmyo*(osmol-(Nai+Ki+Cai))
c            imperm_e = vcleft*(310.0 -(Nae+Ke+Cae))
	    imperm_e = 8.4363e-4
  	

c        Note - multiply by 0.1 for units consistency
c	 Lp... Hydraulic conductivity of the membrane = 1.2e-10 L/(Ns)
         water_flux = 1.2e-10*8.314*310*((Nai+Ki+Cai+
     &           (imperm_i/vmyo))-(Nae+Ke+Cae+
     &             (imperm_e/vcleft)))

         dvmyo = dt*1.534e-4*water_flux*0.1
	 vmyo  = vmyo + dvmyo
	 Vcleft = 43.182e-6 - vmyo/0.68
	
	end if

	Nai = Nai_number/Vmyo
	Ki  = Ki_number/Vmyo
	Cai = Cai_number/vmyo

	Nae = Nae_number/Vcleft
	ke  = Ke_number/vcleft
	Cae = Cae_number/vcleft

c****************************************************************************
c	Calculation of different Ionic currents in the sarcolemma
c****************************************************************************
c	Fast sodium current: INa
	GNa = 16.
	ENa= RTF*log(nae/nai)
	am = 0.32*(Vm+47.13)/(1.-exp(-0.1*(Vm+47.13)))
        bm = 0.08*exp(-Vm/11.)

	if(Vm.ge.-40.) then
         ah = 0
	 bh = 1./(0.13*(1.+exp((Vm+10.66)/(-11.1))))
	 aj = 0
	 bj = 0.3*exp(-2.535e-7*Vm)/(1.+exp(-0.1*
     &          (Vm+32.)))
	else
	 ah = 0.135*exp((80.+Vm)/(-6.8))
	 bh = 3.56*exp(0.079*Vm)+3.1e5*exp(0.35*Vm)
	 aj = (-1.2714e5*exp(0.2444*Vm)-3.474e-5*exp(-0.04391*
     &            Vm))*(Vm+37.78)/(1.+exp(0.311*
     &            (Vm+79.23)))
	 bj = 0.1212*exp(-0.01052*Vm)/(1.+exp(-0.1378*
     &           (Vm+40.14)))
         end if
	
	tm0 = 1./(am+bm)
	th0 = 1./(ah+bh)
	tj0 = 1./(aj+bj)

	if (tm0.le.dt) tm0 = dt
	if (th0.le.dt) th0 = dt
	if (tj0.le.dt) tj0 = dt	
        
	m0 = am*tm0
	h0 = ah*th0
	jj0= aj*tj0

	m1 = m1  + (dt/tm0)*(m0- m1)
	h  = h   + (dt/th0)*(h0- h)
	jj = jj  + (dt/tj0)*(jj0-jj)


	INa = GNa*m1**3*h*jj*(Vm-ENa)

c	currents through the L-type Calcium Channel

	d0 = 1./(1.+exp(-(Vm+10.)/6.24))

	td0 = d0*(1.-exp(-(Vm+10.)/6.24))/
     &	     (0.035*(Vm+10.))
	  

	f0 = (1./(1.+exp((Vm+32.)/8.))) + (0.6/
     &              (1.+exp((50.-Vm)/20.)))
	tf0 = 1./(0.0197*exp(-(0.0337*(Vm+10.))**2)+0.02)


	d = d + (dt/td0)*(d0 - d)
	f = f + (dt/tf0)*(f0 - f)
                 

	IbarCa =  PCa*ZCa**2*Vm*frdy/RTF*(ycai*cai*
     &            exp(zca*Vm/RTF)-yCae*cae)/
     &            (exp(zca*Vm/RTF)-1.)

	IbarCaK = PK* ZK**2 *Vm*frdy/RTF*(yKi *ki *
     &            exp(zK*Vm/RTF)-yKe*ke)/
     &            (exp(zK*Vm/RTF)-1.)

	Ibarcana = PNa*ZNa**2*Vm*frdy/RTF*(yNai*Nai*
     &            exp(zNa*Vm/RTF)-yNae*nae)/
     &            (exp(zNa*Vm/RTF)-1.)

c	km2ca = 0.0006
	fca = 1./(1.+(cai/km2ca))
 
	Ilca   = d*f*fca*Ibarca
	Ilcana = d*f*fca*Ibarcana
	Ilcak  = d*f*fca*Ibarcak

	IlCatot = Ilca + Ilcana + Ilcak

c	Current through T-type Ca Channel

	b0 = 1/(1 + exp(-(vm+14.)/10.8))
	tb0 = 3.7 + 6.1/(1+exp((vm+25.)/4.5))

	g0 = 1/(1 + exp((vm+60)/5.6))
	If (vm.le.0) then
	   tg0 = -0.875*vm + 12.
	else
	   tg0 = 12.
	end if

	b = b + (dt/tb0)*(b0-b)
	g = g + (dt/tg0)*(g0-g)

	gcat = 0.05
	eca  = (RTF/2)*log(cae/cai)

	icat = gcat*b**2*g*(Vm - eca)
	
c	Rapidly Activating Potassium Current (Ikr)

	gkr = 0.02614*sqrt(ke/5.4)
	ekr = RTF*log(ke/ki)

	Xr0 =1/(1 + exp(-(Vm + 21.5)/7.5))
	txr0=1/(0.00138*(Vm+14.2)/(1-exp(-0.123*(vm+14.2)))
     &         +  0.00061*(vm+38.9)/(exp(0.145*(vm+38.9))-1))

	Xr = Xr + (dt/txr0)*(Xr0-Xr)

	rr = 1/(1 + exp((Vm+9)/22.4))


	Ikr = gkr*Xr*rr*(vm - ekr)

c	Slow Component of the Delayed Rectifier K current (Iks)

	Pca1   = -log10(cai) + 3.
	Gbarks = 0.057 + 0.19/(1 + exp((-7.2 + pca1)/0.6))

	Eks=RTF*log((Ke+0.01833*Nae)/
     &                             (Ki+0.01833*Nai))

	Xs0  = 1./(1 + exp(-(vm - 1.5)/16.7))

	txs10=1./(0.0000719*(vm+30)/
     &                       (1-exp(-0.148*(vm+30))) +
     &           0.000131*(vm+30)/(exp(0.0687*(vm+30))-1))
	Xs = Xs + (dt/txs10)*(Xs0-Xs)

	Iks = Gbarks*Xs*xs*(Vm - Eks)

c	Time independent Potassium current: Ik1 (time independent)

	Gk1 = 0.75*sqrt(ke/5.4)
	EK1 = RTF*log(ke/Ki)

	aK1 = 1.02/(1.+exp(0.2385*(Vm-EK1-59.215)))

	bK1 = (0.49124*exp(0.08032*(Vm-EK1+5.476))+
     &          exp(0.06175*(Vm-EK1-594.31)))/
     &          (1.+exp(-0.5143*(Vm-EK1+4.753)))

	tK10 = 1./(aK1+bK1)
	K10  = aK1*tK10

	IK1 = Gk1*k10*(Vm-EK1)

c	Plateau Potassium current: IKp (time independent)

c	GKp      = 0.0183
	Gkp	 = 0.00552
	EKp = EK1
	Kp  = 1./(1.+exp((7.488-Vm)/5.98))
	IKp = GKp*Kp*(Vm-EKp)


c 	Formulation of Ikatp proposed by Terkildsen et al. 

	fatp = (1.0 - (0.05*ATP/(0.05*ATP + 12)))**4.0

	f_bound_adp = (0.95*ADP/(0.95*ADP + 0.42105))**2.0

	f_close_adp = 1.0 - f_bound_adp

	fadp = 1.0 - (f_close_adp)**4

	fkatp = fatp*(0.08*(1.0 - fadp) + 0.89*fadp)

	gkatp = 0.4*fkatp*(5.4/4.5)**0.24

	ekatp = RTF*log(ke/ki)

	Ikatp = gkatp*(Vm - EKatp)
c	Ikatp = 0

c	Na-Ca exchange current:INaCa

	INaCa= KNaCa*(1./(KmNa**3+nae**3))*(1./(KmCa+cae))*
     &             (1./(1.+Ksat*exp((eta-1.)*Vm*(1./RTF))))*
     &             (exp(eta*Vm*(1./RTF))*nai**3*cae - 
     &             exp((eta-1.)*Vm*(1./RTF))*nae**3*cai)
c	INaCa = 0
c**************************************************************************
c	Na-K Pump current: INaK (time independent)
	

c	Sodium Voltage Partition

         Kd_nak_nae1 = K0d_nak_nae*EXP((1.0 + delta_nak)*Vm*(1.0/RTF))
 

         Kd_nak_nai1 = K0d_nak_nai*EXP(delta_nak*Vm*(1.0/RTF))

c	Apparent Concentrations of enzyme-bound ions:

 	 Nae_bound1 = Nae/Kd_nak_nae1

         Nai_bound1 = Nai/Kd_nak_nai1

         Nae_bound2 = Nae/Kd_nak_nae2

         Nai_bound2 = Nai/Kd_nak_nai2

         Ke_bound = Ke/Kd_nak_ke

         Ki_bound = Ki/Kd_nak_ki

         MgATP_bound = ATP/Kd_nak_MgATP

c	Concentration of inorganic Phosphate:
	Pi_nak = Pi_nak_tot/(1.0 +(Ki/Kd_nak_kpi)+
     &           ((8.0297126e-5)/(10.0**(3.0-Kd_nak_hpi)))+
     &           (Nai/Kd_nak_napi))

c	Forward reaction rates:

         alpha_f_nak_1 = rate_f_nak_1*Nai_bound1*
     &             (Nai_bound2**2.0)/(((1.0 + Nai_bound1)*
     &             (1.0 + Nai_bound2)**2.0) +
     &             ((1.0 + Ki_bound)**2.0) - 1.0)

         alpha_f_nak_2 = rate_f_nak_2

         alpha_f_nak_3 = rate_f_nak_3*(Ke_bound**2.0)/
     &             (((1.0 + Nae_bound1)*(1.0 + Nae_bound2)**2.0)+
     &             ((1.0 + Ke_bound)**2.0) - 1.0)

         alpha_f_nak_4 = rate_f_nak_4*MgATP_bound/
     &             (1.0 + MgATP_bound)

c	Reverse rate constants
         
	 alpha_r_nak_1 = rate_r_nak_1*ADP

         alpha_r_nak_2 = rate_r_nak_2 * Nae_bound1*
     &            ((Nae_bound2)**2.0)/
     &            (((1.0 + Nae_bound1)*(1.0 + Nae_bound2)**2.0)+
     &            ((1.0 + Ke_bound)**2.0) - 1.0)

         alpha_r_nak_3 = rate_r_nak_3 * Pi_nak*protons/
     &             (1.0 + MgATP_bound)

         alpha_r_nak_4 = rate_r_nak_4*(Ki_bound**2.0)/
     &             (((1.0 + Nai_bound1)*(1.0 + Nai_bound2)**2.0)+
     &             ((1.0 + Ki_bound)**2.0)-1.0)

c	Clockwise steady-state cycle rate:

	sigma_nak =alpha_r_nak_1*alpha_r_nak_2*alpha_r_nak_3 + 
     &             alpha_r_nak_1*alpha_r_nak_2*alpha_f_nak_4 + 
     &             alpha_r_nak_1*alpha_f_nak_3*alpha_f_nak_4 +
     &             alpha_f_nak_2*alpha_f_nak_3*alpha_f_nak_4 +
     &             alpha_r_nak_2*alpha_r_nak_3*alpha_r_nak_4 +
     &             alpha_f_nak_1*alpha_r_nak_2*alpha_r_nak_3 + 
     &             alpha_f_nak_1*alpha_r_nak_2*alpha_f_nak_4 + 
     &             alpha_f_nak_1*alpha_f_nak_3*alpha_f_nak_4 +
     &             alpha_r_nak_1*alpha_r_nak_3*alpha_r_nak_4 +
     &             alpha_f_nak_2*alpha_r_nak_3*alpha_r_nak_4 +
     &             alpha_f_nak_1*alpha_f_nak_2*alpha_r_nak_3 + 
     &             alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_4 + 
     &             alpha_r_nak_1*alpha_r_nak_2*alpha_r_nak_4 + 
     &             alpha_r_nak_1*alpha_f_nak_3*alpha_r_nak_4 +
     &             alpha_f_nak_2*alpha_f_nak_3*alpha_r_nak_4 +
     &             alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_3

         v_cyc_nak = (alpha_f_nak_1*alpha_f_nak_2*alpha_f_nak_3*
     &             alpha_f_nak_4 - alpha_r_nak_1*alpha_r_nak_2*
     &             alpha_r_nak_3*alpha_r_nak_4)/sigma_nak

c	Na-K Pump Current:
	inak = v_cyc_nak*density_nak*1.6e-5

c**************************************************************************

c	Sarcolemmal Calcium pump current: IpCa (time independent)
	
	IpCa = IbarpCa*(cai/(KmpCa+cai))

c	Calcium Background Current:ICab (time independent)

	Gcab = 0.003016
	EcaN = (RTF/2.)*log(cae/cai)
	ICab = GCab*(Vm- ECaN)

c	Sodium Background Current: INab (time independent)

	GNab = 0.00141
	ENan = ENa
	INab = GNab*(Vm- ENan)

c	Total time independent current: Iv
	Iv = IK1 + IKp + IpCa + INab
     &            + ICab + INaK

c	Calculation of total independent currents

	Naiont = INa + INab + Ilcana + 3.*INak 
     &                + Insna + 3.*INaca + Inas
	Kiont  = Ikr + Iks + Ik1 + IKp + IlcaK 
     &               + InsK - 2.*INak + Ikatp + istim
	caiont = ilca + Ipca + Icab - 2.*INaca
     &                + Icat

	it = Naiont + Kiont + Caiont


c************************************************************************
c	Hund's Modification

c	Naiont = INa + INab + Ilcana 

c	Kiont  = Ikr + Iks + Ik1 + IKp + IlcaK 
c     &                + istim

c	caiont = ilca + Ipca + Icab + Icat 

c	it = Naiont+Kiont+Caiont+INak+INaca

c	dNai = -dt*((Naiont+3*INak+3*INaca)*Acap)/
c     &                           (Vmyo*zNa*frdy)
c	nai  =  dNai + nai

c	dKi = -dt*((Kiont-2*INaK)*Acap)/(Vmyo*zk*frdy)
c	ki  =  dki + ki

c	dNae = dt*((Naiont + 3*INak+3*INaca)*Acap)/
c     &              (Vcleft*frdy) 
c	nae  =  dnae+nae

c	dKe = dt*((Kiont-2*INaK)*Acap)/(Vcleft*frdy)
c     &              
c	Ke  =  dKe+Ke

c	dCae = dt*(caiont*Acap)/(2*Vcleft*frdy)
c     &              
c	cae  =  dCae+Cae

c****************************************************************************

c********************************************************************************
c	Functions that calculate intracellualr & extracellular ion concentrations 
c*********************************************************************************

	dNai = -dt*((Naiont)*Acap)/(Vmyo*zNa*frdy)
	nai  =  dNai + nai

	dKi = -dt*(Kiont*Acap)/(Vmyo*zk*frdy)
	ki  =  dki + ki


	if (k.gt.200*30000) then

	dNae = dt*((Naiont)*Acap)/(Vcleft*frdy)
	nae  =  dnae+nae

	
	dKe = dt*(Kiont*Acap)/(Vcleft*frdy)
	Ke  =  dKe+Ke

	dCae = dt*(caiont*Acap)/(2*Vcleft*frdy)	
	cae  =  dCae+Cae

	end if


c**************************************************************************************
c	Main formula for Vm
c**************************************************************************************
	Vmnew = vm - (it)*dt

	dvmdtnew = (vmnew - vm)/dt	
		
*************************************************************************
	
	If (apdflag.eq.0) then
	 If (k.gt.is1s) then
	   If(vmnew.lt.vm.and.iflag2.eq.0) then
            tmax = k*dt
            vpeak=vmnew

	  iflag2=1
	  apdflag = 1

	   end if
	  end if
	 end if
 	
	
	if (apdflag.eq.1) then
	if(k.gt.is1e+400) then
        if(vmnew.lt.(vrest + (vpeak-vrest)*0.1).and.
     &            iflag3.eq.0) then
            Apd=k*dt-tmax
	
	  iflag1=1
	  iflag3=1
	  apdflag = 0

	Write(23,*) tmax/(1000*60),apd,nai,ki,ke,Inak,
     &              vm,Ikatp,crp,pi_nak_tot,atp,adp,ph,
     &		    vmyo, vcleft
     
	write(33,*) tmax/(60*1000),Kd_nak_nae1,Kd_nak_nai1,Nae_bound1,
     &              Nai_bound1,Nae_bound2,Nai_bound2,Ke_bound,Ki_bound
            
	   end if
	  end if
	end if
	
	
c*************************************************************************
c
c      iclock = 0 means before the 2 ms integrating time
c      iclock = 1 means during the 2 ms integrating time
c      iclock = 2 means after the 2 ms integrating time
c
c      note: do not even look for Vdotmax if before stimulus (artifact)
c
	If(k.lt.is1s) then
         continue
        else
        If(dvmdtnew.lt.dvmdtold.and.dvmdtnew
     &                .gt.25.and.iclock.eq.0) then
            caizero = cai+ cmdn + trpn + dcai
     &                     + csqn + jsr + nsr 

	    vdotmax = dvmdtnew
	  
	    tdvdtmax = k*dt
            iclock = 1
            clock = 0.
            write(6,*) ' dvdtmax now ',k*dt,i,j

          end if
         end if
	
	 if (clock.gt.2.0.and.iclock.eq.1) then
	   caitwo = cai+ cmdn + trpn + dcai
     &                   + csqn + jsr + nsr

           dcaitwo = caitwo - caizero
	
	    iclock = 0
	
           If (dcaitwo.gt.dcaith) then
              grelcicrbar= 60.
              write(6,*) ' cicr!!!!!!', k*dt,i,j

           else 
             grelcicrbar=0.
           end if
	   
	      dvmdtold = 0    
	      clock = 0
	
           end if

              clock=clock+dt

	 grelcicr=grelcicrbar*((dcaitwo - 0.00018)/
     &                        (kmrel + dcaitwo - 0.00018))*
     &        (1.-exp(-clock/tauon))*exp(-clock/tauoff)

	 Irelcicr = grelcicr*(jsr-cai)

c	 Irelcicr = 0
	 
c	Ca uptake and leakage of NSR
	
	kleak = Iupbar / nsrbar
	Ileak = kleak * nsr
	Iup = Iupbar * cai / (cai+kmup)

c     	Iup = 0

c 	Translocation of Ca ions from NSR to JSR
	tautr = 180. 
	Itr = (nsr-jsr)/tautr

c	Calculation of calcium concentration in NSR
	dnsr = dt*(Iup - Ileak - Itr*vjsr/vnsr)
	nsr = nsr + dnsr

c	Calculation of calcium concentration in JSR
c 	Analytical expression for Ca2+ buffering in the junctional 
c	sarcoplasmic reticulum (JSR) and in the cytosol under steady 
c	state conditions taken from appendix 1 of Zeng et al 1995 paper

	csqn = csqnbar * (jsr/(jsr + kmcsqn))
	djsr = dt*(itr - Irelcicr)


	bjsr = csqnbar - csqn -djsr - jsr + kmcsqn
	cjsr = kmcsqn * (csqn + djsr + jsr)
	jsr = (sqrt(bjsr**2 + 4. * cjsr) - bjsr)/2.


	dcai = -dt*(((caiont*Acap)/(Vmyo*zca*frdy))
     &	   - (Irelcicr*vjsr/vmyo) 
     &     + ((Iup - Ileak)*vnsr/vmyo))


	dcaidt = dcai/dt

	trpn = trpnbar * (cai/(cai + kmtrpn))
	cmdn = cmdnbar * (cai/(cai + kmcmdn))

	catotal = trpn + cmdn + dcai + cai
 	bmyo = cmdnbar + trpnbar - catotal + kmtrpn + kmcmdn
	cmyo = (kmcmdn*kmtrpn) - (catotal*(kmtrpn + kmcmdn))
     &              + (trpnbar*kmcmdn) + (cmdnbar*kmtrpn)
	dmyo = -kmtrpn * kmcmdn * catotal

	gpig = sqrt(bmyo*bmyo-3.*cmyo)


	gpig1 = cos(acos((9.*bmyo*cmyo-2.*bmyo**3
     &     -27.*dmyo)/(2.*(bmyo**2 - 3.*cmyo)**1.5))/3.)

	cai = (2.*gpig*gpig1/3.) - (bmyo/3.)


	vm =  vmnew
	dvmdtold = dvmdtnew

	 icount=icount+1
           if(icount.eq.ioutput) then
              iout=iout+1
              icount=0
         time=k*dt
	 i = 1
	 j = 1
 
c	 write(6,*) time, vmnew
	 write(6,*)time,vmnew, nai,nae,ki,ke,vmyo,vcleft,imperm_e
c         write(6,*)vmnew, m1, jj, h,d, f
c	 write(6,*) 5,nai, ki, cai, nsr,jsr,Irelcicr
c        write(6,1000) time

c1000       format(4096(f6.1,' '))

	 
c	Write(iout,*)time,Vmnew,INa, Ilca, Ik,Iv,Inaca,Ik1, Ikp,Inak,
c     &              Ipca,INab,Icab,cai,jsr,nsr,trpn,cmdn,csqn,Irelcicr,
c     &               Iup,Ileak,itr,fca
        write(26,*) time,Vmnew,INa,Ilca,Ik,
     &              Iv,Inaca,Ik1, Ikp,
     &              Inak,Ipca,INab,Icab,
     &              cai,jsr,nsr,trpn,cmdn,
     &              csqn,Irelcicr,Iup,Ileak,
     &              itr,fca,dvmdtnew,It
        write(21,*) Vmnew
	write(20,*) time
c	write(33,*) vdotmax
	
	Write(24,*) time,vmnew,Ikatp,Inak,protons
        write(25,*) time, vm, nai, nae,ki, Ke, cai,cae,
     &          atp,adp,nsr,jsr,trpn,cmdn,csqn,m1, h,jj,
     &		d, f, b, g,x1,xi,xr,xs,imperm_e
        write(27,*) K0d_nak_nae, K0d_nak_nai, Kd_nak_nae2,Kd_nak_nai2,
     &    	Kd_nak_ke, Kd_nak_ki, Kd_nak_MgATP, rate_f_nak_1, 
     &          rate_f_nak_2, rate_f_nak_3, rate_f_nak_4,rate_r_nak_1,
     &	        rate_r_nak_2, rate_r_nak_3, rate_r_nak_4         
               
	write(30,*)time,vmnew,Ikatp,inak

	  end if

c	End of time loop
	end do

c       close files
        close(20)
        close(21)
        close(22)

	Stop 
	End


	

	
	
