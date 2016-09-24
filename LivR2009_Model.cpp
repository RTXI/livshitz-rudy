/*
 * Copyright (C) 2011 Weill Medical College of Cornell University
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 675 Mass
 * Ave, Cambridge, MA 02139, USA.
 */

/*** INTRO
 *
 * 2009 Livshitz Rudy Model of a ventricular guinea pig myocte Biophysical
 * Journal, 2009
 * 
 * LivR2009_Model.cpp, v1.1
 *
 * Author: Francis Ortega (v1.0-v1.1)(2011)
 *
 *** NOTES
 *
 * v1.0 - initial version
 * v1.1 - added ability to change model rate
 *
 * Model equations taken from Matlab model created by Dr. Eric Sobie.  Plugin
 * constructed using only DefaultGUIModel.
 *
 * math.h function exp(x) replaced with fastEXP(x) which uses PowFast.hpp due
 * to spikes in computation time of math.h function.
 *
 ***/

#include <iostream>
#include <math.h>
#include "LivR2009_Model.h"
#include <powfast.hpp>

using namespace std;

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new LivR2009_Model();
}

static DefaultGUIModel::variable_t vars[] = {
	{ "Istim","A",DefaultGUIModel::INPUT, },
	{ "Vm","V (mv)",DefaultGUIModel::OUTPUT, },
	{ "IKs","IKs (A/F)",DefaultGUIModel::STATE, },
	{ "INa","INa (A/F)",DefaultGUIModel::STATE, },
	{ "INab","INab (A/F)",DefaultGUIModel::STATE, },
	{ "ICaL","ICaL_Na (A/F)",DefaultGUIModel::STATE, },
	{ "ICaL_K","ICaL_K (A/F)",DefaultGUIModel::STATE, },
	{ "ICab","ICab (A/F)",DefaultGUIModel::STATE, },
	{ "ICaT","ICaT (A/F)",DefaultGUIModel::STATE, },
	{ "IpCa","IpCa (A/F)",DefaultGUIModel::STATE, },
	{ "IKr","IKr (A/F)",DefaultGUIModel::STATE, },
	{ "IK1","IK1 (A/F)",DefaultGUIModel::STATE, },
	{ "IKp","IKp (A/F)",DefaultGUIModel::STATE, },
	{ "INCX","INCX (A/F)",DefaultGUIModel::STATE, },
	{ "INaK","INaK (A/F)",DefaultGUIModel::STATE, },
	{ "Itotal","Itotal (A/F)",DefaultGUIModel::STATE, },
	{ "Rate","Model Integrate Rate (Hz)",
	   DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
};

static size_t num_vars = sizeof(vars)/sizeof(DefaultGUIModel::variable_t);

/*** Model constructor using DefaultGUIModel ***/
LivR2009_Model::LivR2009_Model(void): DefaultGUIModel("LivR2009_Model",::vars,::num_vars) {

	DefaultGUIModel::createGUI( vars, num_vars );

	Cm = 100e-12; // 100pF

	// Initial conditions at 0 beats
	/*
		V = -84.7;
		Cai = 0.0822e-3;
		CaNSR = 1.25;
		CaJSR = 1.25;
		Nai = 9.71;
		Ki = 142.82;
		m = 2.46e-4;
		h = 0.99869;
		j = 0.99887;
		d = 1e-4;
		f = 0.983;
		b = 1e-4;
		g = 0.983;
		xKr = 0.229;
		xs1 = 1e-4;
		xs2 = 1e-4;
		Jrel = 1e-4;
	*/

	//Initial conditions after conc change - 800 beats
	V = -8.185872e+01;
	Cai = 1.797384e-04; 
	CaNSR = 2.463960e+00; 
	CaJSR = 1.524945e+00;
	Nai = 1.357382e+01; 
	Ki = 1.239044e+02; 
	m = 2.601169e-03; 
	h = 9.697101e-01; 
	j = 9.806867e-01; 
	d = 5.928788e-322; 
	f = 9.981991e-01; 
	b = 1.866354e-03;
	g = 9.650771e-01;
	xKr = 4.291283e-04; 
	xs1 = 2.954278e-02; 
	xs2 = 8.283927e-02;
	Jrel = 1.101473e-38;

	// Integration Variables, 100khz rate default (10khz minimum needed)
	// Model will always run at this rate, regardless of thread rate
	rate = 100000;
	dt = (1/rate)*1000;

	// Initialize states to 0 (aesthetic only)
	IKs = 0;
	INa = 0;
	INab = 0;
	ICaL = 0;
	ICaL_Na = 0;
	ICaL_K = 0;
	ICab = 0;
	ICaT = 0;
	IpCa = 0;
	IKr = 0;
	IK1 = 0;
	IKp = 0;
	INCX = 0;
	INaK = 0;
	Iion = 0;

	// Set states
	setState("IKs",IKs);
	setState("INa",INa);
	setState("INab",INab);
	setState("ICaL",ICaL);
	setState("ICaL_Na",ICaL_Na);
	setState("ICaL_K",ICaL_K);
	setState("ICab",ICab);
	setState("ICaT",ICaT);
	setState("IpCa",IpCa);
	setState("IKr",IKr);
	setState("IK1",IK1);
	setState("IKp",IKp);
	setState("INCX",INCX);
	setState("INaK",INaK);
	setState("Itotal",Iion);

	// Set parameters
	setParameter("Rate",rate);

	period = RT::System::getInstance()->getPeriod()*1e-6;
	steps = static_cast<int>(ceil(period*rate/1000.0));

	lkup= new double[20000][20];
	Vx = 0; // Voltage placeholder for lookup table
	V_min = -1000;

	// Lookup Table
	for(z=0; z<20000; z++){
	Vx = V_min+0.1*z;

	/* Voltage */
	lkup[z][0] = V_min+0.1*z;

	/* H-gate */
	lambda_na = 1 - 1/(1+exp(-(Vx+40)/0.024));
	ah = lambda_na * 0.135 * exp(-(80+Vx)/6.8);
	bh = (1-lambda_na) / (0.13*(1+exp((Vx+10.66)/(-11.1)))) + 
	     lambda_na * (3.56*exp(0.079*Vx) + 3.1 * 1e5 * exp(0.35*Vx));
	// hinf
	lkup[z][1] = ah/(ah + bh);
	// tauh
	lkup[z][2] = 1/(ah + bh);

	/* J-gate */
	aj = lambda_na * (-1.2714e5*exp(0.2444*Vx)-3.474e-5*exp(-0.04391*Vx)) * 
	     (Vx+37.78) / (1+exp(0.311*(Vx+79.23)));
	bj = (1-lambda_na) * (0.3*exp(-2.535e-7*Vx) / (1+exp(-0.1*(Vx+32)))) + 
	lambda_na * (0.1212*exp(-0.01052*Vx) / (1+exp(-0.1378*(Vx+40.14))));
	// tauj
	lkup[z][3] = 1/(aj + bj);
	// jinf
	lkup[z][4] = aj/(aj + bj);

	/* M-gate */
	if( Vx > -47.14 && Vx < -47.12)// if V = -47.13, divide by 0 error
	am = 3.41333;
	else
	am = 0.32 * (Vx+47.13) / (1-exp(-0.1*(Vx+47.13)));
	bm = 0.08 * exp(-Vx/11.0);

	// minf
	lkup[z][5] = am / (am + bm);
	// taum
	lkup[z][6] = 1/(am + bm);

	/* D-gate */
	dinf_0 = 1/(1+exp(-(Vx+10)/6.24));
	dinf_1 = 1/(1+exp(-(Vx+60)/0.024));
	// dinf
	lkup[z][7] = dinf_0 * dinf_1;
	// taud
	if( Vx > -10.01 && Vx < -9.99)// if V = -10, divide by 0 error
		lkup[z][8] = 2.30655;
	else
		lkup[z][8] =  1/(1+exp(-(Vx+10)/6.24)) * (1-exp(-(Vx+10)/6.24))/(0.035*(Vx+10));

	/* F-gate */
	// finf
	lkup[z][9] = 1/(1+exp((Vx+32)/8.0))+(0.6)/(1+exp((50-Vx)/20.0));
	// tauf
	lkup[z][10] = 1/(0.0197 * exp(-(0.0337*(Vx+10))*(0.0337*(Vx+10))) + 0.02);

	/* B-gate */
	// binf
	lkup[z][11] = 1/(1+exp(-(Vx+14.0)/10.8));
	// taub
	lkup[z][12] = (3.7 + 6.1/(1+exp((Vx+25.0)/4.5)));

	/* G-gate */
	lambda_g = 1 - 1/(1+exp(-Vx/0.0024));
	// ginf
	lkup[z][13] = 1/(1+exp((Vx+60.0)/5.6));
	// taug
	lkup[z][14] = 1 * (lambda_g*(-0.875*Vx+12.0) + 12.0 * (1-lambda_g));

	/* IKr */
	if( Vx > -30.01 && Vx < -29.99 )// if V = -30, divide by 0 error
		tau_xs1 = 411.501;
	else
		tau_xs1 = 10000/(0.719*(Vx+30) / (1-exp(-0.148*(Vx+30))) + 1.31 * (Vx+30) / (exp(0.0687*(Vx+30))-1));
	// xKrinf
	lkup[z][15] = 1/(1+exp(-(Vx+21.5)/7.5));
	// tauxKr
	if( Vx > -14.21 && Vx < -14.19 ) // if V = -14.2, divide by 0 error
		lkup[z][16] = 87.1735;
	lkup[z][16] = ( 1/(0.00138*(Vx+14.2) / (1-exp(-0.123*(Vx+14.2))) + 0.00061 * (Vx+38.9) / (exp(0.145*(Vx+38.9))-1)) );
	// tauxs1
	lkup[z][17] = tau_xs1;
	// tauxs2
	lkup[z][18] = 4 * tau_xs1;
	// xsinf
	lkup[z][19] = 1/(1+exp(-(Vx-1.5)/16.7));
	}

	// PowFast object
	powFast = new PowFast(18);

	refresh();
	resizeMe();
}

LivR2009_Model::~LivR2009_Model(void) {
	delete []lkup;
}

/*** Current and gating solver ***/
void LivR2009_Model::solve(double DT) {

	/* Reversal potentials */
	ENa = RTF * log(Nao/Nai);
	EK = RTF * log(Ko/Ki);
	EKs = RTF * log((Ko + pKNa*Nao)/(Ki + pKNa*Nai));
	ECa = 0.5 * RTF * log(Cao/Cai);

	/* Na currents */
	INa = GNa_ * m*m*m * h * j * (V - ENa);
	INab = GNab * (V - ENa);

	/* L-type Ca current */
	ICa_ = PCa * 4 * F * FRT * V * (gamma_Cai*Cai*fastEXP(2*V*FRT) - gamma_Cao*Cao)/(fastEXP(2*V*FRT) - 1);
	ICaK_ = PCa_K * F * FRT * V * (gamma_Ki*Ki*fastEXP(V*FRT) - gamma_Ko*Ko)/(fastEXP(V*FRT) - 1);
	ICaNa_ = PCa_Na * F * FRT * V * (gamma_Nai*Nai*fastEXP(V*FRT) - gamma_Nao*Nao) / (fastEXP(V*FRT) - 1) ;
	fCa = 1/(Cai/KmCa + 1);
	ICaL = ICa_ * d * f * fCa;
	ICaL_K = ICaK_* d * f * fCa;
	ICaL_Na = ICaNa_ * d * f * fCa;

	// Background calcium current
	ICab = GCab * (V - ECa);

	// Sarcolemmal calcium pump
	IpCa = IpCa_ * Cai / (Cai + KmpCa);

	/* T-type Ca current */
	ICaT = GCaT * b*b * g * (V-ECa);

	/* K currents */
	// Time independent K current
	// fastExp results in segmentation error for exp(0.6987*(V-EK+11.724)) for *unknown reason*
	xK1 = 0.004 * (1 + fastEXP(0.6987*(V-EK+11.724))) / (1+fastEXP(0.6168*(V-EK+4.872)));

	IK1 = GK1_ * sqrt(Ko/5.4) * (V-EK) / (1+xK1);

	// Fast component of the delayed rectifier K current
	RKr = 1/(fastEXP((V+9)/22.4) + 1);
	IKr = GKr_ * sqrt(Ko/5.4) * xKr * RKr * (V - EK);

	// Fast component of the delayed rectifier K current
	IKs = GKs_ * (1 + 0.6/(fastEXP(1.4*log(3.8e-5/Cai))+1)) * xs1 * xs2 * (V - EKs);

	// Plateau K current
	Kp = 1/(1+fastEXP((7.488-V)/5.98));
	IKp = GKp_ * Kp * (V - EK);

	/* Pumps and transporters */
	// Na-K pump
	sigma_NaK = (fastEXP(Nao/67.3) - 1) / 7.0; 
	fNaK = 1/(1 + 0.1245*fastEXP(-0.1*V*FRT) + 0.0365 * sigma_NaK * fastEXP(-V*FRT));
	INaK = INaK_ * fNaK * Ko / ( (Ko + KmK_NaK) * ( 1 + ((KmNa_NaK/Nai)*(KmNa_NaK/Nai))) );

	// Na-Ca exchanger
	INCX = kNCX * fastEXP((eta-1)*V*FRT) * ((Nai*Nai*Nai)*Cao*fastEXP(V*FRT)-(Nao*Nao*Nao)*Cai) / 
	(1+ksat*fastEXP((eta-1)*V*FRT) * ((Nai*Nai*Nai)*Cao*fastEXP(V*FRT)+(Nao*Nao*Nao)*Cai) );

	/* Intracellular Ca fluxes */
	// SR Ca release, uptake, and leak
	Jrelinf = alpha_rel * beta_tau * ICaL / (fastEXP(hrel*log(Krel_inf/CaJSR)) + 1);
	tau_rel = beta_tau / (Krel_tau/CaJSR + 1);
	dJreldt = - (Jrelinf + Jrel) / tau_rel;

	Jserca = Vserca * (Cai/(Cai+Kmserca) - CaNSR/CaNSR_max );

	Jtr = (CaNSR-CaJSR) / tau_transfer;

	/* Buffering factors for rapid buffering approximation */
	BJSR = 1.0/(1 + CSQNtot * KmCSQN / ((KmCSQN + CaJSR)*(KmCSQN + CaJSR)));
	Bi = 1.0/(1 + ( CMDNtot * KmCMDN / ( (Cai+KmCMDN)*(Cai+KmCMDN) ) ) + (TRPNtot * KmTRPN / ( (Cai+KmTRPN)*(Cai+KmTRPN) ) ) );

	/* Total Current */
	Iion = INa + INab + ICaL + ICaL_Na + ICaL_K + ICab + ICaT + IpCa + IKr + IKs + IK1 + IKp + INCX + INaK - (input(0)/Cm);

	/* Derivatives for ionic concentration */
	dNai = -(INa + INab + ICaL_Na + 3*INCX + 3*INaK)*Acap/(Vmyo*F);
	dKi = -(IKr + IKs + IK1 + IKp + ICaL_K - 2*INaK - (input(0)/Cm))*Acap/(Vmyo*F);
	dCai = Bi * ( -Jserca*VNSR/Vmyo + Jrel*VJSR/Vmyo - (ICaL + ICaT + ICab + IpCa - 2*INCX) * Acap / (2*Vmyo*F) ) ;
	dCaJSR = BJSR * (Jtr - Jrel);
	dCaNSR = Jserca - Jtr * VJSR / VNSR;

	/* Derivative for voltage */
	dVdt = -(Iion);
		
	/* Update voltage and ionic concentrations */

	V += DT*dVdt;
	Nai += DT*dNai;
	Ki += DT*dKi;
	Cai += DT*dCai;
	CaJSR += DT*dCaJSR;
	CaNSR += DT*dCaNSR;
	Jrel += DT*dJreldt;

	/* Set gating variables using lookup table */
	ilow = fabs((V-V_min)/0.1);
	linext = -(-V+lkup[ilow+1][0])/0.1;
	hinf = (lkup[ilow+1][1]-lkup[ilow][1])*linext+lkup[ilow+1][1];
	tauh = (lkup[ilow+1][2]-lkup[ilow][2])*linext+lkup[ilow+1][2];
	tauj = (lkup[ilow+1][3]-lkup[ilow][3])*linext+lkup[ilow+1][3];
	jinf = (lkup[ilow+1][4]-lkup[ilow][4])*linext+lkup[ilow+1][4];
	minf = (lkup[ilow+1][5]-lkup[ilow][5])*linext+lkup[ilow+1][5];
	taum = (lkup[ilow+1][6]-lkup[ilow][6])*linext+lkup[ilow+1][6];
	dinf = (lkup[ilow+1][7]-lkup[ilow][7])*linext+lkup[ilow+1][7];
	taud = (lkup[ilow+1][8]-lkup[ilow][8])*linext+lkup[ilow+1][8];
	finf = (lkup[ilow+1][9]-lkup[ilow][9])*linext+lkup[ilow+1][9];
	tauf = (lkup[ilow+1][10]-lkup[ilow][10])*linext+lkup[ilow+1][10];
	binf = (lkup[ilow+1][11]-lkup[ilow][11])*linext+lkup[ilow+1][11];
	taub = (lkup[ilow+1][12]-lkup[ilow][12])*linext+lkup[ilow+1][12];
	ginf = (lkup[ilow+1][13]-lkup[ilow][13])*linext+lkup[ilow+1][13];
	taug = (lkup[ilow+1][14]-lkup[ilow][14])*linext+lkup[ilow+1][14];
	xKrinf = (lkup[ilow+1][15]-lkup[ilow][15])*linext+lkup[ilow+1][15];
	tauxKr = (lkup[ilow+1][16]-lkup[ilow][16])*linext+lkup[ilow+1][16];
	tauxs1 = (lkup[ilow+1][17]-lkup[ilow][17])*linext+lkup[ilow+1][17];
	tauxs2 = (lkup[ilow+1][18]-lkup[ilow][18])*linext+lkup[ilow+1][18];
	xsinf = (lkup[ilow+1][19]-lkup[ilow][19])*linext+lkup[ilow+1][19];

	/* Update gate variables */
	//h += DT*((hinf-h) / tauh);
	h = (hinf - (hinf-h)*fastEXP(-DT/tauh)); //Rush-Larsen approximation
	//j += DT*((jinf-j) / tauj);
	j = (jinf - (jinf-j)*fastEXP(-DT/tauj)); //Rush-Larsen approximation
	//m += DT*((minf-m) / taum);
	m = (minf - (minf-m)*fastEXP(-DT/taum)); //Rush-Larsen approximation
	d += DT*((dinf - d)/taud);
	f += DT*((finf - f)/tauf);
	b += DT*((binf-b) / taub);
	g += DT*((ginf-g) / taug);
	xKr += DT*((xKrinf - xKr)/tauxKr);
	xs1 += DT*((xsinf - xs1) / tauxs1);
	xs2 += DT*((xsinf - xs2) / tauxs2);
} 

/*** Execute function ***/
void LivR2009_Model::execute(void) {
	/*
	 * Because the real-time thread may run much slower than we want to
	 *   integrate we need to run multiple interations of the solver.
	 */
	start  = RT::OS::getTime()*1e-3;
	for(int i = 0;i < steps;++i){
		solve(dt);
	}
	output(0) = V*1e-3; //output in nV (to mimic amplifier)

	end  = RT::OS::getTime()*1e-3;

	/*
	 * if( (end-start) > (period*1e3) ) {
	 * 	printf("\nWARNING: Calculation time over 50ms\n");
	 * 	setActive(false);
	 *	}
	 */
}

/*** Update function ***/
void LivR2009_Model::update(DefaultGUIModel::update_flags_t flag) {
	if(flag == MODIFY) { // Modify button resets model to initial conditions
		// Initial conditions at 0 beats
		/*
			V = -84.7;
			Cai = 0.0822e-3;
			CaNSR = 1.25;
			CaJSR = 1.25;
			Nai = 9.71;
			Ki = 142.82;
			m = 2.46e-4;
			h = 0.99869;
			j = 0.99887;
			d = 1e-4;
			f = 0.983;
			b = 1e-4;
			g = 0.983;
			xKr = 0.229;
			xs1 = 1e-4;
			xs2 = 1e-4;
			Jrel = 1e-4;
		*////////////////////////////////

		V = -8.185872e+01;
		Cai = 1.797384e-04; 
		CaNSR = 2.463960e+00; 
		CaJSR = 1.524945e+00;
		Nai = 1.357382e+01; 
		Ki = 1.239044e+02; 
		m = 2.601169e-03; 
		h = 9.697101e-01; 
		j = 9.806867e-01; 
		d = 5.928788e-322; 
		f = 9.981991e-01; 
		b = 1.866354e-03;
		g = 9.650771e-01;
		xKr = 4.291283e-04; 
		xs1 = 2.954278e-02; 
		xs2 = 8.283927e-02;
		Jrel = 1.101473e-38;

		Nai = 13.7601813818;
		Ki = 135.590987673;
		Cai = 0.000195549061095;
		CaJSR = 1.75015936476;
		CaNSR = 2.53205784602;

		rate = getParameter("Rate").toDouble();

		cout << "\n\nRT Period: " << RT::System::getInstance()->getPeriod() << "\n";

		period = RT::System::getInstance()->getPeriod()*1e-6;                                                               
		steps = static_cast<int>(ceil(period*rate/1000.0));
		dt = (1/rate)*1000;
		std::cout << xK1 << " " << V << " " << EK << " " << IK1 << " " << GK1_ << std::endl;
	}
	else if(flag == PERIOD) {
		period = RT::System::getInstance()->getPeriod()*1e-6;
		steps = static_cast<int>(ceil(period*rate/1000.0));
	}
}

/*** fastEXP Function
*
* Fast and stable alternative to math.h exp(x) function
* Uses powFast library and contains overflow prevention
* Max error < 0.001% when compared to math.h exp(x)
*
***/
double LivR2009_Model::fastEXP(double x){
	if(x > 88.5 || x < -87){ // Overflow prevention due to powFast library's use of float

		// Close to 0 shortcut: If x < -746, math.h exp(x) function returns 0 due to truncation
		if( x < -746 ) 
			return 0;

		if( x > 710 ){ // Infinity shortcut: If x > 710, math.h exp(x) returns infinity
			return INFINITY;
		}

		// Standard overflow prevention
		// Due to powFast limitations from using floats,
		// if x is too large or small, problem is broken into smaller pieces
		// e.g. fastEXP(100) ==> e^180 = e^60 * e^60 * e^60
		powCount = 1;
		powAns = 1;
		powExp = x;

		while(powExp > 88.5 || powExp < -87){
			powCount++;
			powExp = x / powCount;
		}

		for(n = 1; n <= powCount; n++){
			powAns = powAns * powFast->e(powExp);
		}

		return powAns;

	} // end overflow prevention

	else // call powFast e^x function
		return powFast->e(x);
}
