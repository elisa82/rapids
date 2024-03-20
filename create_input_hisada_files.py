def create_input_grflt12s(folder, computational_param, layers, fault, sites):
	import numpy as np
	from pbs.conversions import read_slip_xta


	if fault['slip_mode'] == 'file_xta':
		slip_distribution, rake_distribution, rise_time, trup = read_slip_xta(folder, fault, layers)
		NF = fault['number_subfaults_strike'] * fault['number_subfaults_dip']

		#x consistenza con SPEED
		slip_srcmod = slip_distribution.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
		rupt_srcmod = trup.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
		rise_srcmod = rise_time.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
		rake_srcmod = rake_distribution.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])

		slip_distribution = np.reshape(np.flipud(slip_srcmod),(1, NF))[0]
		rake_distribution = np.reshape(np.flipud(rake_srcmod),(1, NF))[0]

	if fault['fault_type'] == 'point':
		slip_distribution = fault['slip']
		rake_distribution = fault['rake']

	NF = fault['number_subfaults_strike'] * fault['number_subfaults_dip']

	fid = open(folder+'/grflt12s_slip_vr.in', 'w')
	fid.write('  *Data for Dt, Duration, and Period *\n')
	fid.write('{:8.6f}{:8d}{:s}\n'.format(computational_param['dt_hisada'], computational_param['npts_hisada'],
											' : Delta Time(sec), and Number of Time(must be Power of 2 (f8.0, i8)'))
	fid.write('{:8.6f}{:8.3f}{:s}\n'.format(1/computational_param['fmax_hisada'], 0.005, ' : Minimum Period(sec), '
														  'and Imaginary Omega for Phinney (2f8.0)'))
	fid.write('  *MEDIUM DATA *\n')
	nlayers = len(layers['thk'])
	thk_hisada = layers['thk']
	thk_hisada[nlayers - 1] = 0
	fid.write('{:10d}{:s}\n'.format(nlayers, '    : NL(NUMBER OF LAYERS, max 20) t/m3, m/s, , m/s, , m'))
	for i in range(nlayers):
		layers['thk'][nlayers-1] = 0
		fid.write('{:7.3f}{:7.1f}{:7.1f}{:4.1f}{:7.1f}{:7.1f}{:4.1f}{:8.1f} '
				  '{:s}\n'.format(layers['rho'][i], layers['vp'][i]*1000, layers['qp'][i], layers['fqp'][i],
								  layers['vs'][i]*1000, layers['qs'][i], layers['fqs'][i],
								  thk_hisada[i]*1000, ' : DNS, VP, Qp, fqp, VS, Qs, fqs, Thick'))
	fid.write('  *FAULT DATA *\n')
	fid.write('{:10.1f}{:10.1f}{:s}\n'.format(fault['length']*1000, fault['width']*1000,
												' : Total Length(STR) and Width(DIP)(m) (2f10.0),'))
	ngauss = 6
	fid.write('{:5d}{:5d}{:5d}{:s}\n'.format(fault['number_subfaults_strike'], fault['number_subfaults_dip'], ngauss,
											  ' : Num.of Sub - Faults(NF=NSTR * NDIP, max 1000), Gauss Points per '
											  'Wavelength(max 6, ma li aggiusta lui aumentandoli con aumento della '
											  'frequenza di calcolo e usandone il max indicato qui (3I5)'))
	#if fault['IDx'] == 'exp':
	STFnum = 2 #per ora si può considerare solo esponenziale
	fid.write('{:10.2f}{:5d}{:s}\n'.format(0., STFnum, ': Time Delay(sec)(per hypo!!), '
													'Source Type( = 0:Ramp; =1: Smoothed; =2: exp) (F10.0, I5)'))
	fid.write('{:15.2f}{:15.2f}{:15.2f}{:s}\n'.format(fault['hypo_utm']['X'], fault['hypo_utm']['Y'],
														 fault['hypo_utm']['Z'], ' :  Location of Hypocenter(Xnord, '
														'Yest and Zdown: m) (pto di nucleazione)'))
	fid.write('{:15.2f}{:15.2f}{:15.2f}{:s}\n'.format(fault['origin']['X'], fault['origin']['Y'],
														 fault['origin']['Z'], ' :  Location of Fault '
														'Origin(Xnord, Yest and Zdown: m) (pto piu profondo a 0 di '
														'strike e dip sulla faglia.)'))
	fid.write('{:10.1f}{:10.1f}{:s}\n'.format(fault['strike'], fault['dip'], ' : Strike, Dip(deg)'))
	v_rup = fault['rupture_velocity']*1000.
	fid.write('{:10.1f}{:10.1f}{:s}\n'.format(0, v_rup, ' : SEED FOR RANDOM GENERATION(for slip and Vr) '
															'( if 0: no random, prog uses given vr); VR(2 f10 .0)'))
	TW = 1
	TWI = 0.
	fid.write('{:8d}{:8.2f}{:s}\n'.format(TW, TWI, ' : Number of Time Windows(TW, max 15)' 
													  ' and Time Interval(sec)(tra le TW)  (i8, f8.0)'))
	tau = fault['rise_time']
	if fault['IDx'] == 'exp':
		tau = tau/(2*np.pi)
	half_rise_time = tau/2
	for i in range(TW):
		fid.write('{:8.4f}{:8.4f}{:s}\n'.format(half_rise_time, half_rise_time, ' : Half - Rise time 1 & 2 for 1st ' 
																	'Time Window (sec) (tanti quanti le TW) (2f8.0)'))


	fid.write('  * dislocations (m) on sub-faults for 6 time windows (tanti quanti le NF)\n')

	for j in range(NF):
		fid.write('{:8.3f}\n'.format(slip_distribution[j]))
	fid.write('  * rake angles (degree) on sub-faults for 6 time windows (tanti quanti le NF)\n')
	for j in range(NF):
		fid.write('{:8.1f}\n'.format(rake_distribution[j]))
	fid.write('  * Static Wavenumber Integ.Data for Greenfields quadrature *\n')
	fid.write('       2.0 : The first corner w * k on real axis(ex. 2.0)\n')
	fid.write('        16 : Initial Num.of Intg.Points for Adaptive Newton - Cotes Quadrature\n')
	fid.write('      10.0 : The second corner w * k on imag.axis(ex. 10.0)\n')
	fid.write('        16 : Initial Num.of Intg.Points for Adaptive Newton - Cotes Quadrature\n')
	fid.write('  * WAVENUMBER INTEGRATION DATA FOR SIMPSONS OR FILONS FORMULAS *\n')
	fid.write('       200 : Number of Integration Points from 0 to om / Ryleigh(min)[200 / 400]\n')
	fid.write('        50 : Number of Integration Points from om / Ryl(min) to om / c(final)[50]\n')
	fid.write('      10.0 : Factor for c(final): c(final) = Ryl(min) / Factor(Use 10, usually)\n')
	fid.write('  * CHANGE OF SIGNS OF IMAGINARY PARTS OF FINAL RESULTS(FOR FFT) *\n')
	fid.write('         1 : CHANGE SIGN( = 1), NOT CHANGE SIGN( = 0)\n')
	fid.write('  * OBSERVATION POINT DATA *\n')
	#nobs = len(sites['Z'])
	nobs = 50
	fid.write('{:10d} {:s}\n'.format(nobs, ': NUMBER OF OP'))
	for k in range(nobs):
		fid.write('{:15.2f}{:15.2f}{:10.4f}\n'.format(sites['X'][k], sites['Y'][k], sites['Z'][k]))
	fid.close()
	return


def create_answers_grflt12s(folder):
	fid = open(folder +'/gr.dat', 'w')
	fid.write('0\n')  # Skip computation of static terms?
	# Enter 0 for NO  (Exact, but slow. Use for surface fault)
	# Enter 1 for YES(Fast, but approximate.Effective for buried fault)
	fid.write('0\n')  # Output Static Integrands?
	# Enter 0 for No(Regular)
	# Enter 1 for YES(Checking Purpose Only
	fid.write('0\n')  # Skip computation of dynamic terms?
	# Enter 0 for NO (Use for normal condition)
	# Enter 1 for YES (Use, when output static terms only)
	fid.write('0\n')  # Please choose the following combinations of integrands
	# 0: Regular( = Dynamic Integrand - Static Integrand)
	# 1: Dynamic( = Dynamic Integrand)
	# 2: Static( = Static Integrand, checking purpose only)
	# (Enter 1, if the fault is deep.It is faster.)
	fid.write('0\n')  # Output Dynamic Integrands?
	# Enter 0 for No(Regular)
	# Enter 1 for YES(Checking Purpose Only)
	fid.write('0\n')  # Use Dynamic Distribution of Gaussian Points on Fault?
	# Enter 0 for No  (Regular Num. of Points are fixed)
	# Enter 1 for YES (Faster, but may lose accuracy, only with dseed=0!!)
	fid.close()
	return


def create_answers_grfftspmm(folder):
	import sys
	fid = open(folder +'/grspm.dat', 'w')
	fid.write('1\n') # Choose the input datafile for FFT (Enter Number 1 to 6):
	# 1) grfault.dxyz (x,y,z-comp. of disp. from seismic sources)
	# 2) grfault.dqvw (u,v,w-comp. of disp. from seismic sources)
	# 3) grpoint.dxyz (x,y,z-comp. of disp. from point sources)
	# 4) grpoint.sxyz (x,y,z-comp. of stress from point sources)
	# 5) grpoint.dqvw (u,v,w-comp. of disp.  from point sources)
	# 6) grpoint.sqvw (u,v,w-comp. of stress from point sources)
	fid.write('0\n')  # Change signs of amplitudes? (Yes=1 or No=0) :
	# if computational_param['output_type'] == 'vel':
	# 	transfer = 0
	# if computational_param['output_type'] == 'acc':
	# 	transfer = 1
	transfer = 0 #ho messo di default hisada in velocità
	fid.write('{}\n'.format(str(transfer)))  # You can transfer the velocities to the accelerations;
	# Do not transfer (Output Velocities)       (Enter 0):
	# Transfer from velocities to accelerations (Enter 1):
	fid.write('0\n')  # Band Pass Filter? (Yes=1 or No=0) :
	fid.write('0\n')  # Increase J for smoother seismograms? (Yes=1 or No=0)
	fid.write('0\n')  # Change from (m) to (cm)? (yes=0)
	fid.close()
	return


def create_answers_phs3sQ(folder):
	fid = open(folder +'/ph.dat', 'w')
	fid.write('0               Check secular functions and eigen vectors (0=no, 1=yes)\n')
	fid.write('0.001           Tolerance\n')
	fid.write('100             Max number of iterations\n')
	fid.write('200             Number of Partitions  from CMIN to CMAX\n')
	fid.close()
	return


def create_script_hisada(folder , path_code_hisada):
	fid = open(folder +'/do_all.sh', 'w')
	fid.write('#!/bin/bash\n')
	fid.write('rm ?_dat.* *_dat\n')
	fid.write('rm fort* *disper grfault*\n')
	fid.write('"'+path_code_hisada+'/phs3sQ_slip_vr.out" < ph.dat\n')
	fid.write('"'+path_code_hisada+'/grflt12s_slip_vr_subf.out" < gr.dat\n')
	fid.write('"'+path_code_hisada+'/grfftspmm.out" < grspm.dat\n')
	fid.close()
	return