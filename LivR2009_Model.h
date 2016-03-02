/*
 * Copyright (C) 2011 Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*** INTRO
 *
 * 2009 Livshitz Rudy Model of a ventricular guinea pig myocte
 * Biophysical Journal, 2009
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
 * Model equations taken from Matlab model created by Dr. Eric Sobie.
 * Plugin constructed using only DefaultGUIModel.
 *
 * math.h function exp(x) replaced with fastEXP(x) which uses PowFast.hpp due
 * to spikes in computation time of math.h function.
 *
 ***/

#include <default_gui_model.h>
#include "include/PowFast.hpp"

// Physical constants 
const static double F = 96485;          // Faraday's constant, C/mol
const static double R = 8314;           // gas constant, mJ/K
const static double T = 273+37;         // absolute temperature, K
                                        // const static double T = 273+22;
const static double RTF = R*T/F;
const static double FRT = 1/RTF;
const static double pi = 4.0*atan(1.0); // Pi

// Cell Geometry 
const static double length_cell = 0.01;             // Length of the cell (cm)
const static double radius = 0.0011;                // Radius of the cell (cm)
const static double Vcell = 1000*pi*radius*
                            radius*length_cell;     // 3.801e-5 uL Cell volume (uL)
const static double Ageo = 2*pi*radius*radius+
                           2*pi*radius*length_cell; // 7.671e-5 cm^2 Geometric membrane area (cm^2)
const static double Acap = 2*Ageo;                  // 1.534e-4 cm^2 Capacitive membrane area (cm^2)
const static double Vmyo = Vcell*0.68;              // Myoplasm volume (uL)
const static double Vmito = Vcell*0.24;             // Mitochondria volume (uL)
const static double VNSR = Vcell*0.0552;            // NSR volume (uL)
const static double VJSR = Vcell*0.0048;            // JSR volume (uL)
const static double Vss = Vcell*0.02;

// Cell Capacitance 
static double Cm = 1.0;

// Fixed ionic concentrations 
const static double Ko = 5.4 ;  // uM 
const static double Nao = 137;  // uM 
const static double Cao = 2.0 ; // uM 

// Na current constants 
const static double GNa_= 16; // mS/cm^2
const static double GNab = 0.004;
double GNaL_= 6.5e-3;

// Ca current constants	
const static double PCa = 5.4e-4;       // cm/s
const static double PCa_Na = 6.75e-7;   // cm/s
const static double PCa_K = 1.93e-7;    // cm/s
const static double PCab = 1.995084e-7; // cm/s
const static double gamma_Cao = 0.341;  // dimensionless
const static double gamma_Cai = 1;      // dimensionless
const static double gamma_Nao = 0.75;   // dimensionless
const static double gamma_Nai = 0.75;   // dimensionless
const static double gamma_Ko = 0.75;    // dimensionless
const static double gamma_Ki = 0.75;    // dimensionless
const double hLca = 1;                  // dimensionless, Hill coefficient
const static double KmCa = 6e-4;        // Half saturation constant, mM

// T-type & background currents	
const static double GCaT = 0.05;
const static double GCab = 0.003016;

// K Currents 
const static double GK1_ = 0.75;
const static double GKr_ = 0.02614;
const static double GKs_ = 0.433;
const static double pKNa = 0.01833; // relative permeability of IKs, Na to K
const static double GKp_ = 5.52e-3;

const static double INaK_ = 2.25;   // Max. current through Na-K pump (uA/uF)
const static double KmNa_NaK = 10;  // Half-saturation concentration of NaK pump (mM)
const static double KmK_NaK = 1.5;  // Half-saturation concentration of NaK pump (mM)

const static double kNCX = 0.00025;
const static double ksat = 0.0001;
const static double eta = 0.15;

const static double alpha_rel = 0.125;
const static double Krel_inf = 1;
const static double hrel = 9;
const static double beta_tau = 4.75;
const static double Krel_tau = 0.0123;

// Pumps and Transporters 
const static double IpCa_ = 1.15;     // Max. Ca current through sarcolemmal Ca pump (uA/uF)
const static double KmpCa = 5e-4;     // Half-saturation concentration of sarcolemmal Ca pump (mM)
const static double Vserca = 8.75e-3; // mM/ms
const static double Kmserca = 9.0e-4; // mM
const static double CaNSR_max = 15.0;
const static double tau_transfer = 120;

// Buffers in cytosol 
const static double TRPNtot = 70e-3;    
const static double KmTRPN = 0.5e-3;    

const static double CMDNtot = 50e-3;
const static double KmCMDN = 2.38e-3;    

// Buffers in JSR 
const static double CSQNtot = 10;
const static double KmCSQN = 0.8;


class LivR2009_Model : public DefaultGUIModel {

	public:
		LivR2009_Model(void);
		virtual ~LivR2009_Model(void);
		void execute(void);

	protected:
		void update(DefaultGUIModel::update_flags_t);

	private:
		void solve(double);

		// fastEXP function using PowFast
		double fastEXP(double);

		double V;
		double dVdt;

		// Lookup table variables
		double V_min;
		double Vx;
		int z;
		int ilow;
		double linext;
		double (*lkup)[20];

		// Integration variables 
		double dt;
		int steps;
		double rate;

		// Gates and intracellar concentrations 
		double Cai;
		double CaNSR;
		double CaJSR;
		double Nai;
		double Ki;
		double m;
		double h;
		double j;
		double d;
		double f;
		double b;
		double g;
		double xKr;
		double xs1;
		double xs2;
		double Jrel;

		// Gating variables in lookup table 
		double lambda_na;
		double ah;
		double bh;
		double aj;
		double bj;
		double am;
		double bm;
		double dinf_0;
		double dinf_1;
		double lambda_g;
		double xsinf;
		double tau_xs1;
		double hinf;
		double tauh;
		double tauj;
		double jinf;
		double minf;
		double taum;
		double dinf;
		double taud;
		double finf;
		double tauf;
		double binf;
		double taub;
		double ginf;
		double taug;
		double xKrinf;
		double tauxKr;
		double tauxs1;
		double tauxs2;

		// Current Variables 
		double ENa, EK, EKs, ECa;                  // Reversal potetentials 
		double INa, INab;                          // Na current 
		double ICa_, ICaK_, ICaNa_, fCa, 
		       ICaL, ICaL_K, ICaL_Na;              // L-type Ca current 
		double ICab;                               // Background calcium current 
		double IpCa;                               // Sarcolemmal calcium pump 
		double ICaT;                               // T-type calcium current
		double xK1, IK1, RKr, IKr, IKs, Kp, IKp;   // K currents
		double sigma_NaK, fNaK, INaK, INaKx, INCX; // Pumps and transporters
		double Iion;                               // Total current
		double Iinjected;                          // Injected current

		// Intracellular Ca flux 	
		double Jrelinf, tau_rel, Jserca, Jtr;	

		// Buffering 
		double BJSR, Bi;

		// Derivatives 
		double dNai, dKi, dCai, dCaJSR, dCaNSR, dJreldt;

		// Integration Rate 
		double period;

		// Loop variable 
		int i;

		// Debug variables 
		double start;
		double end;

		// PowFast variables 
		const PowFast *powFast;
		double powAns;
		double powExp;
		int powCount;
		int n;
};
