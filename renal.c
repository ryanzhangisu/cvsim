 /* This file contains the subroutines necessary to calculate the body-fluid homoeostasis 
 * including Kidney, Hormone & Electrolyes Balance and Baroreceptor Reflex systems.
 * 
 *
 * Yan Zhang, January 30, 2014
 */
#include "renal.h"

#define PI M_PI

void renal_function(Data_vector *pressure, Reflex_vector *r, double cp[13], double rp[8], double qw[8], 
double qs[8], double c[8], double cn[10], double ref[7], double coef[6], double eff[17], double vol[3], 
double tao[3], double vas[3], double dydt[8], double sodin, double waterin, double delta_t)
{
	static double cum_hr = 0.0;
	static double rel_time = 0.0; // relative time elapsed in one heart beat
	static double MA_P = 0.0, RA_P = 0.0, RV_P = 0.0, LA_P = 0.0, LV_P = 0.0; // mean arterial pressure and right atrial pressure
	static double RA_C = 0.0, RV_C = 0.0, LA_C = 0.0, LV_C = 0.0;
	static double RA_V = 0.0, RV_V = 0.0, LA_V = 0.0, LV_V = 0.0;
	static double kidney_p = 0.0, kidney_f = 0.0;
	
	double vb_old = 0.0;
	double dt = 0.0;
	double t = 0.0;
	
 	MA_P += pressure->x[0]*delta_t; // beat-by-beat average mean aortic pressure
	RA_P += pressure->x[15]*delta_t; // beat-by-beat average right atrial pressure
	RV_P += pressure->x[16]*delta_t; // ...right ventricular
	LA_P += pressure->x[19]*delta_t; // ...left atrial
	LV_P += pressure->x[20]*delta_t; // ...left ventricular
	RA_C += pressure->c[0]*delta_t; // beat-by-beat average right atrial compliance
	RV_C += pressure->c[1]*delta_t; // ...right ventricular
	LA_C += pressure->c[2]*delta_t; // ...left atrial
	LV_C += pressure->c[3]*delta_t; // ...left ventricular	
	RA_V += pressure->v[15]*delta_t; // beat-by-beat average right atrial volume
	RV_V += pressure->v[16]*delta_t; // ...right ventricular
	LA_V += pressure->v[19]*delta_t; // ...left atrial
	LV_V += pressure->v[20]*delta_t; // ...left ventricular	
	
	kidney_p += pressure->x[7]*delta_t; // average renal arterial pressure
	kidney_f += pressure->q[9]*delta_t; // average renal blood flow
	
	rel_time += delta_t;
	
	if ((*r).hr[1] < cum_hr)  // detect the onset of a new beat
	{
		static int n = 0;
		printf("Enter Kidney Cycle No. %d \n", n);
		
		cp[0] = MA_P = MA_P / rel_time;
		cp[1] = RA_P = RA_P / rel_time; // once a new heart beat starts, average MAP and RAP for the last beat
		cp[2] = RV_P = RV_P / rel_time;
		cp[3] = LA_P = LA_P / rel_time;
		cp[4] = LV_P = LV_P / rel_time;
		cp[5] = RA_C = RA_C / rel_time;
		cp[6] = RV_C = RV_C / rel_time;
		cp[7] = LA_C = LA_C / rel_time;
		cp[8] = LV_C = LV_C / rel_time;
		cp[9] = RA_V = RA_V / rel_time;
		cp[10] = RV_V = RV_V / rel_time;
		cp[11] = LA_V = LA_V / rel_time;
		cp[12] = LV_V = LV_V / rel_time;
		
		rp[0] = kidney_p / rel_time; 	// average renal arterial pressure
		//qw[3] = kidney_f / rel_time /1000.; 	// average renal blood flow
		
		//printf("currently working at %d min\n", n*1);	
		//if(t<= 24.*60.*5.) sodin = 0.126;
		//else if(t>24.*60.*5. && t<=24.*60.*10.) sodin = 0.126 * 2.0;
		//else sodin = 0.02;
		/*if(t > (24.*60.*5.  + 30.) && t< (24.*60.*5.+ 150.)) {
			waterin = 4./1000.;
			sodin = 0.126*5.55;
		}
		else{
			 waterin = 0.0;
			 sodin = 0.126;
		 }*/
		dt = rel_time;  		// unit conversion from sec. to min.
		
		rk4_routine(cp, rp, qw, qs, c, cn, ref, coef, eff, vol, tao, vas, dydt, sodin, waterin, dt);
		
		vb_old = vol[2];
		
		renal_eqns(cp, rp, qw, qs, c, cn, ref, coef, eff, vol, tao, vas, dydt, sodin, waterin);
		
		//qw[7] = (vol[2] - vb_old) / delta_t * 1000.;	// update flow rate changes to venous returen
		qw[7] = (waterin + qw[0] - qw[6]) * 1000.;
		
		r->resistance[1] = rp[7] * 0.06;						// update renal vascular resistance
		//r->resistance[0] = 4.9 * ref[1];
		//r->resistance[2] = 3.0 * ref[1];
		//r->resistance[3] = 4.5 * ref[1];
		
		t = pressure->time[0];
		
		if(n % 10 == 0)renal_write(1, n, t, cp, rp, eff, cn, qw, qs, c, ref, vol, vas);
		
		MA_P = 0.0;
		RA_P = 0.0;
		RV_P = 0.0;
		LA_P = 0.0;
		LV_P = 0.0;
		RA_C = 0.0;
		RV_C = 0.0;
		LA_C = 0.0;
		LV_C = 0.0;
		RA_V = 0.0;
		RV_V = 0.0;
		LA_V = 0.0;
		LV_V = 0.0;
		
		kidney_p = 0.0;
		kidney_f = 0.0;
		
		rel_time = 0.0;
		n++;		
	}	
	cum_hr = (*r).hr[1];
}

void renal_eqns(double cp[13], double rp[8], double qw[8], double qs[8], double c[8], double cn[10], 
double ref[7], double coef[6], double eff[17], double vol[3], double tao[3], double vas[3], double dydt[8], 
double sodin, double waterin)
{				
		// EQN(32) update blood volume 
		vol[2] = 4.560227 + 2.431217 / (1.0 + exp(-(vol[1] - 18.11278) * 0.47437));	
		// EQN(56) update sodium concentration 
		c[1] = vol[0] / vol[1];
		// EQN(38) update vascularity destruction rate
		vas[1] = vas[0] * coef[5];
		// EQN(40) update basic arterial pressure
		//r[2] = coef[4] / vas[0];
		// EQN(48) update ADH concentration		
		c[2] = 4.0 * cn[1];
		// EQN(62) update Ang-II conentration
		c[5] = 20.0 * c[6];
		// EQN(63) update ALD concentration
		c[3] = cn[3] * 85.0;
		if (c[3] > 500.0) c[3] = 500.0;

		/*
		// define intermediate variables for Steffensen's Iteration
		double p_ma[3] = {0.0};
		double p_ra[3] = {0.0};
		//double p_ra_old = 0.0, a_auto_old = 0.0;
		double a_auto_old = 0.0;
		//p_ra_old = p[1];
		a_auto_old = ref[3];
		do 
		{
			p_ma[0] = p[0];
			int i = 0;
			for(i = 1; i < 3; i++)
			{
				// update autonomous multiplier
				ref[3] = 3.079 * exp(-p[0] * 0.011);
				ref[5] = 0.25 * ref[3];
				ref[1] = ref[5] + ref[4];
				// update arterial resistance
				r[0] = r[2] * ref[1];
				// update mean filling pressure
				p[2] = (7.436 * vol[2] - 30.18) * ref[1];
				// update total peripheral resistance
				r[7] = r[0] + r[3];
				// update resistance to venous return
				r[1] = (8.0 * r[3] + r[0]) / 31.0;
				
				do
				{
					p_ra[0] = p[1];
					int j = 0;
					for(j = 1; j < 3; j++)
					{
						// update venous return
						qw[2] = (p[2] - p[1]) / r[1];
						// update cardiac output
						qw[1] = qw[2];
						// update right atrial pressure
						p[1] = 0.2787 * exp(qw[1] * 0.2281) - 0.89;
						//p[1] = 0.2787 * exp(qw[1] * 0.2281);
						p_ra[j] = p[1];
					}
					p[1] = p_ra[0] - (p_ra[1] - p_ra[0]) * (p_ra[1] - p_ra[0])/(p_ra[2] - 2.0*p_ra[1] + p_ra[0]);
				} while(fabs(p[1] - p_ra[0]) > 0.001);
				
				// update mean arterial pressure
				p[0] = qw[1] * r[7];
				p_ma[i] = p[0];
			}
			p[0] = p_ma[0] - (p_ma[1] - p_ma[0]) * (p_ma[1] - p_ma[0])/(p_ma[2] - 2.0*p_ma[1] + p_ma[0]);
		} while(fabs(p[0] - p_ma[0]) > 0.001);
		*/

		double a_auto_old = 0.0;
		a_auto_old = ref[3];
		// update autonomous multiplier
		ref[3] = 3.079 * exp(-cp[0] * 0.011);
		ref[5] = 0.25 * ref[3];
		ref[1] = ref[5] + ref[4];
		
		// update RSNA
		eff[0] = 0.5 + 1.1 / (1.0 + exp((cp[0] - 100.0) / 15.0)); 
		eff[1] = 1.0 - 0.008 * (cp[1]-2.3);
		ref[0] = cn[0] * eff[0] * eff[1];
		
		double qs_md[3] = {0.0};
		
		do
		{
			qs_md[0] = qs[3];
			int k = 0;
			for(k = 1; k < 3; k++)
			{
				// update tubuloglomerular feedback 
				ref[2] = 0.3408 + 3.449 / (3.88 + exp((3.859 - qs[3])/0.9617));
				//ref[2] = 0.3412 + 0.06296 / (0.07079 + exp(-2.064 * qs[3]));
				
				// update renal vascular resistance
				eff[2] = 1.5 * (ref[0] - 1.0) + 1.0;
				rp[5] = 31.67 * eff[2] * ref[2];
				rp[6] = 51.66 * (0.9432 + 0.1363 / (0.2069 + exp(3.108 - 1.785 * log10(c[5]))));
				rp[7] = rp[5] + rp[6];
				
				//r->resistance[1] = rp[6] * 0.06;
				// update renal blood flow
				qw[3] = rp[0] / rp[7];
				// update glomerular filtration rate
				rp[1] = rp[0] - qw[3] * rp[5];				// Bowman's capsule pressure
				rp[4] = rp[1] - (rp[2] + rp[3]);			// filtration pressure
				qw[4] = rp[4] * coef[0]; 					// glomerular filtration rate
				// update filtered sodium loads
				qs[1] = qw[4] * c[1];
				// update proximal tubular sodium reabsorption
				eff[5] = 0.8 + 0.3 / (1.0 + exp(qs[1] - 14.0) / 138.0); 
				eff[4] = 0.95 + 0.12 / (1.0 + exp(2.6 - 1.8 * log10(c[5])));
				eff[6] = 0.5 + 0.7 / (1.0 + exp(1.0 - ref[0]) / 2.18);
				qs[2] = qs[1] * cn[6] * eff[5] * eff[4] * eff[6];
				// update macula densa sodium flow
				qs[3] = qs[1] - qs[2];
				
				qs_md[k] = qs[3];
			}
			qs[3] = qs_md[0] - (qs_md[1] - qs_md[0]) * (qs_md[1] - qs_md[0])/(qs_md[2] - 2.0*qs_md[1] + qs_md[0]);
		} while(fabs(qs[3] - qs_md[0]) > 0.001);
		
		// update distal tubular sodium reabsorption
		eff[16] = 0.17 + 0.94 / (1.0 + exp((0.48 - 1.2 * log10(c[3]))/0.88));
		qs[4] = qs[3] * cn[7] * eff[16];
		// update distal tubular sodium flow
		qs[5] = qs[3] - qs[4];
		// update collecting duct sodium reabsorption
		eff[8] = 0.82 + 0.39 / (1.0 + exp((qs[5] - 1.6) / 2.0));
		eff[7] = -0.1 * c[4] + 1.1;   // change from 1.1199 to 1.1
		qs[6] = qs[5] * cn[8] * eff[8] * eff[7];
		// update urinary sodium flow 
		qs[7] = qs[5] - qs[6];	
		
		// update tubular water reabsorption
		eff[9] = 0.37 + 0.8 / (1.0 + exp(0.6 - 3.7 * log10(c[2])));
		eff[10] = 0.17 + 0.94 / (1.0 + exp((0.48 - 1.2 * log10(c[3]))/0.88));
		qw[5] = (0.025 - 0.001 / (eff[10] * eff[9])) + 0.8 * qw[4];
		// update unrine flow
		qw[6] = qw[4] - qw[5];
		
		// update ADH secretion rate
		if(c[1] > 141.0 || ref[1] > 1.0) cn[2] = ((c[1] - 141.0) + (ref[1] - 1.0) - eff[3]) / 3.0;
		//else if(cn[1] > 143.0 && ref[1] <= 1.0) cn[2] = ((c[1] - 140.0) - eff[3]) / 3.0;
		//else if(cn[1] <= 143.0 && ref[1] > 1.0) cn[2] = ((ref[1] - 1.0) - eff[3]) / 3.0;
		//else cn[2] = 0.0;
		//else cn[2] = 0.0;
		else cn[2] = 0.01;
		if(cn[2]<0.01) cn[2] = 0.01;
		//else cn[2] = 0.1;
		//if(cn[2] < 0.0) cn[2] = 0.1;
		
		// update renin secretion rate
		eff[11] = 0.2262 + 28.04 / (11.56 + exp((qs[3] - 1.667) / 0.6056));
		eff[12] = 1.89 - 2.056 / (1.358 + exp(ref[0] - 0.8667));
		cn[5] = eff[11] * eff[12];
		
		// update Aldosterone secretion rate
		eff[14] = c[0] / c[1] / 0.00347 - 9.0;
		//if(eff[14] < 0.0) eff[14] = 1.0;
		eff[13] = 0.4 + 2.4 / (1.0 + exp((2.82 - 1.5 * log10(c[5]))/0.8));
		if(cp[0] > 100.) eff[15] = 1.0;
		else eff[15] = 69.03 * exp(-0.0425 * cp[0]);
		cn[4] = eff[14] * eff[15] * eff[13];

		// update vascularity (independent part)
		vas[2] = 11.312 * exp(-qw[1] * 0.4799) / 100000.0;
		
		// update ANP concentration
		c[4] = 7.4052 - 6.554 / (1.0 + exp(cp[1] - 2.3 - 3.762)); // normalized ANP concentration
	
		// update sodium intake
		qs[0] = sodin;
		
		// update thirst-driven water flow
		qw[0] = 0.008 / (1.0 + 1.822 * pow(c[2], -1.607)) - 0.0053;
		//qw[0] = 0.01 * (0.37 + 0.8 /(1.0 + exp(0.6 - 3.7 * log10(c[2])))) - 0.0094;
		if (qw[0] < 0.0) qw[0] = 0.0;
		//qw[0] = 0.001;
		
		// Update Derivatives of Denpendent Varibales in Modules 14, 20, 25, 26, 30, 32, 34
		dydt[0] = waterin + qw[0] - qw[6];						// change of extracellular fluid volume
		//dydt[0] = waterin - qw[6];
		dydt[1] = vas[2] - vas[1];						// change of vascularity
		//dydt[2] = 0.75 * (ref[3] - 0.0005 * ref[4] - 1.0);				// change of baroreflex
		dydt[2] = 0.75 * (ref[3] - a_auto_old) - 0.0005 * (ref[4] - 0.75);
		dydt[3] = (cn[2] - cn[1]) / tao[0];				// change of normal ADH concentration
		dydt[4] = 0.0007 * (0.2 * (cp[1]-2.3) - eff[3]); 		// change of effect of right atrial pressure on ADH
		//dydt[4] = 0.2 * (p[1] - p_ra_old) - 0.0007 * eff[3];
		dydt[5] = qs[0] - qs[7];						// change of total sodium amount
		dydt[6] = (cn[5] - c[6]) / tao[2];				// change of renin concentration
		dydt[7] = (cn[4] - cn[3]) / tao[1]; 			// change of ALD concentration
}

void rk4_routine(double cp[13], double rp[8], double qw[8], double qs[8], double c[8], double cn[10], double ref[7], 
double coef[6], double eff[17], double vol[3], double tao[3], double vas[3], double dydt[8], double sodin, 
double waterin, double dt)
{
	double h = 0.0, hh = 0.0, h6 = 0.0;
	double y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0, y5 = 0.0, y6 = 0.0, y7 = 0.0, y8 = 0.0;
	double k1_1 = 0.0, k2_1 = 0.0, k3_1 = 0.0, k4_1 = 0.0, k5_1 = 0.0, k6_1 = 0.0, k7_1 = 0.0, k8_1 = 0.0;
	double k1_2 = 0.0, k2_2 = 0.0, k3_2 = 0.0, k4_2 = 0.0, k5_2 = 0.0, k6_2 = 0.0, k7_2 = 0.0, k8_2 = 0.0;
	double k1_3 = 0.0, k2_3 = 0.0, k3_3 = 0.0, k4_3 = 0.0, k5_3 = 0.0, k6_3 = 0.0, k7_3 = 0.0, k8_3 = 0.0;
	double k1_4 = 0.0, k2_4 = 0.0, k3_4 = 0.0, k4_4 = 0.0, k5_4 = 0.0, k6_4 = 0.0, k7_4 = 0.0, k8_4 = 0.0;
	
	/* time steps for RK4 */
	h = dt;
	hh = dt / 2.0;
	h6 = dt / 6.0;
	
	/* first step */
	y1 = vol[1];
	y2 = vas[0];
	y3 = ref[4];
	y4 = cn[1];
	y5 = eff[3];
	y6 = vol[0];
	y7 = c[6];
	y8 = cn[3];
	
	k1_1 = dydt[0];
	k2_1 = dydt[1];
	k3_1 = dydt[2];
	k4_1 = dydt[3];
	k5_1 = dydt[4];
	k6_1 = dydt[5];
	k7_1 = dydt[6];
	k8_1 = dydt[7];
	
	/* take half time step to get K2 */
	vol[1] = y1 + hh * k1_1;
	vas[0] = y2 + hh * k2_1;
	ref[4] = y3 + hh * k3_1;
	cn[1] = y4 + hh * k4_1;
	eff[3] = y5 + hh * k5_1;
	vol[0] = y6 + hh * k6_1;
	c[6] = y7 + hh * k7_1;
	cn[3] = y8 + hh * k8_1;
	
	renal_eqns(cp, rp, qw, qs, c, cn, ref, coef, eff, vol, tao, vas, dydt, sodin, waterin);

	k1_2 = dydt[0];
	k2_2 = dydt[1];
	k3_2 = dydt[2];
	k4_2 = dydt[3];
	k5_2 = dydt[4];
	k6_2 = dydt[5];
	k7_2 = dydt[6];
	k8_2 = dydt[7];	
	
	/* take another half time step for K3 */
	vol[1] = y1 + hh * k1_2;
	vas[0] = y2 + hh * k2_2;
	ref[4] = y3 + hh * k3_2;
	cn[1] = y4 + hh * k4_2;
	eff[3] = y5 + hh * k5_2;
	vol[0] = y6 + hh * k6_2;
	c[6] = y7 + hh * k7_2;
	cn[3] = y8 + hh * k8_2;
	
	
	renal_eqns(cp, rp, qw, qs, c, cn, ref, coef, eff, vol, tao, vas, dydt, sodin, waterin);
	
	k1_3 = dydt[0];
	k2_3 = dydt[1];
	k3_3 = dydt[2];
	k4_3 = dydt[3];
	k5_3 = dydt[4];
	k6_3 = dydt[5];
	k7_3 = dydt[6];
	k8_3 = dydt[7];	
	
	/* take a final full time step for K4 */
	vol[1] = y1 + h * k1_3;
	vas[0] = y2 + h * k2_3;
	ref[4] = y3 + h * k3_3;
	cn[1] = y4 + h * k4_3;
	eff[3] = y5 + h * k5_3;
	vol[0] = y6 + h * k6_3;
	c[6] = y7 + h * k7_3;
	cn[3] = y8 + h * k8_3;	
	
	renal_eqns(cp, rp, qw, qs, c, cn, ref, coef, eff, vol, tao, vas, dydt, sodin, waterin);

	k1_4 = dydt[0];
	k2_4 = dydt[1];
	k3_4 = dydt[2];
	k4_4 = dydt[3];
	k5_4 = dydt[4];
	k6_4 = dydt[5];
	k7_4 = dydt[6];
	k8_4 = dydt[7];	 
	
	/* update value for next time step */
	vol[1] = y1 + h6 * (k1_1 + 2.0 * k1_2 + 2.0 * k1_3 + k1_4);
	vas[0] = y2 + h6 * (k2_1 + 2.0 * k2_2 + 2.0 * k2_3 + k2_4);
	ref[4] = y3 + h6 * (k3_1 + 2.0 * k3_2 + 2.0 * k3_3 + k3_4);
	cn[1] = y4 + h6 * (k4_1 + 2.0 * k4_2 + 2.0 * k4_3 + k4_4);
	eff[3] = y5 + h6 * (k5_1 + 2.0 * k5_2 + 2.0 * k5_3 + k5_4);
	vol[0] = y6 + h6 * (k6_1 + 2.0 * k6_2 + 2.0 * k6_3 + k6_4);
	c[6] = y7 + h6 * (k7_1 + 2.0 * k7_2 + 2.0 * k7_3 + k7_4);
	cn[3] = y8 + h6 * (k8_1 + 2.0 * k8_2 + 2.0 * k8_3 + k8_4);	
	
}

void renal_init(Data_vector *pressure, double cp[13], double rp[8], double qw[8], double qs[8], double c[8], double cn[10], double ref[7], double coef[6], 
double eff[17], double vol[3], double tao[3], double vas[3])
{
	/* initialize pressure variables (mmHg) */
	/*p[0] = 100.;						// [P_ma], mean arterial pressure
	p[1] = 0.;							// [P_ra], right arterial pressure
	p[2] = 7.;							// [P_mf], mean filling pressure
	p[3] = 62.;							// [P_gh], glomerular filtration pressure
	p[4] = 18.;							// [P_B], Bowman hydro. pressure
	p[5] = 28.;							// [P_go], glomerular osmotic pressure
	p[6] = 16.;							// [P_f], net filtration pressure
	*/
	
	/* initialize resistance variables (mmHg*min/L) */
	/*r[0] = 16.6;						// [R_a], arterial pressure
	r[1] = 1.4;							// [R_vr], venous return resistance
	r[2] = 16.6;						// [R_ba], basic arterial pressure
	r[3] = 3.4;							// [R_bv], basic venous resistance
	r[4] = 31.67;						// [R_aa], afferent arteriolar resistance
	r[5] = 51.66;						// [R_ea], efferect arteriolar resistance
	r[6] = 83.33;						// [R_r], renal vascular resistance
	r[7] = 20.;							// [R_tp], total peripheral resistance
	*/
	
	cp[0] = pressure->x[0];
	cp[1] = pressure->x[15];
	
	rp[0] = pressure->x[7];							// renal arterial pressure
	rp[1] = pressure->x[7]-1.2*37.67;				// glomerular hydrostatic pressure
	rp[2] = 18.;						// Bowman's capsule pressure
	rp[3] = 28.;						// glomerular osmotic pressure
	rp[4] = 0.;							// filtration pressure
	rp[5] = 37.67;						// [R_aa], afferent arteriolar resistance
	rp[6] = 51.66;						// [R_ea], efferect arteriolar resistance
	rp[7] = 83.33;						// [R_r], renal vascular resistance
	
	/* initialize water/blood flow rate (L/min) */
	qw[0] = 0.001;						// [Q_win], thirst-driven water intake
	qw[1] = 5.;							// [Q_co], cardiac output
	qw[2] = 5.;							// [Q_vr], venous return
	qw[3] = 1.2;						// [Q_rb], renal blood flow
	qw[4] = 0.125;						// [Q_gfilt], glomerular filtraton rate
	qw[5] = 0.124;						// [Q_t-wreab], tubular water reabsorption rate
	qw[6] = 0.001;						// [Q_u], urine flow rate
	qw[7] = 0.0;						// change of flow rate to CVSim Code
	
	/* initialize sodium flow rate (mEq/min) */
	qs[0] = 0.126;						// [Q_sodin], sodium intake
	qs[1] = 18.;						// [Q_filsod], filtered sodium load
	qs[2] = 14.4;						// [Q_pt-sodreab], absolute prox. tubule sodium reab
	qs[3] = 3.6; 						// [Q_md-sod], Macula Densa sodium flow
	qs[4] = 1.8;						// [Q_dt-sodreab], absolute dist. tubule sodium reab
	qs[5] = 1.8;						// [Q_dt-sod], distal tubule sodium outflow
	qs[6] = 1.674;						// [Q_cd-sodreab], absolute collecting duct sodium reab
	qs[7] = 0.126;						// [Q_u-sod], urine sodium flow
	
	/* initialize electrolyte and hormone concentrations (mEq/L, ng/L) */
	c[0] = 5.;							// [C_k], potassium concentration
	c[1] = 144.;						// [C_sod], sodium concentration
	c[2] = 4.;							// [C_adh], ADH concentration
	c[3] = 85.;							// [C_al], Aldosterone concentration
	c[4] = 1.;							// [C_anp], normalized ANP concentration
	c[5] = 20.;							// [C_at], Angiotensin concentration
	c[6] = 1.;							// [C_r], normalized renin concentration
	c[7] = 0.0;							// Alcohol concentration (mmol/L)
	
	/* initialize normalized concentraions (no units) */
	cn[0] = 1.;							// [N_rsna], normalized RSNA
	cn[1] = 1.;							// [N_adh], normalized ADH concentration
	cn[2] = 1.;							// [N_adhs], normalized ADH secretion 
	cn[3] = 1.;							// [N_al], normalized ALD concentration
	cn[4] = 1.;							// [N_als], normalied ALD secretion
	cn[5] = 1.;							// [N_rs], normalized renin secretion
	cn[6] = 0.8;						// [n_eta-pt], normal value of fractional prox. sodium reab
	cn[7] = 0.5;						// [n_epsilon-dt], normal value of fractional dist. sodium reab
	cn[8] = 0.93;						// [n_et-cd], normal value of fractional collecting duct sodium reab
	cn[9] = 1.;							// fraction of ADH reduction due to alcohol intake
	
	/* initialize effects parameters (no units) */
	eff[0] = 1.;						// [alpha_map], effect of mean arterial pressure on RSNA
	eff[1] = 1.;						// [alpha_rap], effect of right atrial pressure on RSNA
	eff[2] = 1.;						// [beta_rsna], effect of RSNA on R_aa;
	eff[3] = 0.;						// [sigma_ra], effect of right atrial pressure on normalized ADH secretion
	eff[4] = 1.;						// [gamma_at], effect of C_at on fractional prox. sodium reab
	eff[5] = 1.;						// [gamma_filsod], effect of Q_filsod on fractional prox. sodium reab
	eff[6] = 1.;						// [gamma_rsna], effect of RSNA on fractional prox. sodium reab
	eff[7] = 1.;						// [lambda_anp], effect of ANP on Q_cd-sodreab
	eff[8] = 1.;						// [lambda_dt], effect of Q_dt-sod on Q_cd-sodreab
	eff[9] = 1.;						// [miu_adh], effect of C_adh on Q_t-wreab
	eff[10] = 1.;						// [miu_al], effect of C_al on Q_t-wreab
	eff[11] = 1.;						// [niu_md-sod], effect of Q_md-sod on normalized renin secretion
	eff[12] = 1.;						// [niu_rsna], effect of RSNA on normalized renin secretion
	eff[13] = 1.;						// [ksi_at], effect of ANG on normalized ALD secretion
	eff[14] = 1.;						// [ksi_k/sod], effect of K to Na ratio on normalized ALD secretion
	eff[15] = 1.;						// [ksi_map], effect of mean arterial pressure on normalized ALD secretion
	eff[16] = 1.;						// [psi_al], effect of C_ald on fractional distal sodium reab
	
	/* initialize reflex parameters (no units) */
	ref[0] = 1.; 					// [rsna], renal sympathetic nerve activity
	ref[1] = 1.;						// [epsilon_aum], autonomic multiplier effect
	ref[2] = 1.;						// [Sigma_tgf], tubulogomerular feedback signal
	ref[3] = 1.;						// [a_auto], autonomus system activity
	ref[4] = 0.75;					// [a_baro], baroreceptor activity
	ref[5] = 0.25;					// [a_chemo], chemorecepter activity
	ref[6] = 0.0;                   // [delta_baro], baro-reflex change, added by Yan Zhang (2014)
	
	/* initialize coefficients and fractions (varied units) */
	coef[0] = 0.00781;					// [C_gcf], glomerular capillary filtration rate
	coef[1] = 0.8;						// [eta_pt-sodreab], fractional proximal sodium reab
	coef[2] = 0.5;						// [eta_dt-sodreab], fractional distal sodium reab
	coef[3] = 0.93;						// [eta_cd-sodreab], fractional collecting duct sodium reab
	coef[4] = 16.6;						// [K_bar], coefficient relating basic arterial resistance to vascularity
	coef[5] = 0.00001;					// [K_vd], vascularity destruction coeffeicient
	
	/* initialize fluid volumes (L) and total sodium amount (mEq) */
	vol[0] = 2160.;						// [M_sod], total sodium amount
	vol[1] = 15.;						// [V_ecf], extracellular fluid volume
	vol[2] = 5.;						// [V_b], blood volume
	
	/* initialize time constancs (minutes) */
	tao[0] = 6.;						// [T_adh], time constant for ADH secretion
	tao[1] = 60.;						// [T_al], time constant for ALD secretion
	tao[2] = 15.;						// [T_r], time constant for renin secretion
	
	/* initialize vascularity parameters (no units) */
	vas[0] = 1.; 						// [vas], vascularity
	vas[1] = 0.00001;					// [vas_d], vascularity destruction rate
	vas[2] = 0.00001;					// [vas_f], vascularity formation rate
	
}


void renal_write(int k, int i, double t, double cp[13], double rp[8], 
double eff[17], double cn[9], double qw[8], double qs[8], double c[7], 
double ref[7], double vol[3], double vas[3])
{
  	char buf[256];
  	
	FILE *fp_cp;
  	(void) sprintf(buf, "cardiac_parameters-%d.dat", k);
	fp_cp= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_cp, "Time Pma Pra Prv Pla Plv Cra Crv Cla Clv Vra Vrv Vla Vlv\n");
		fprintf(fp_cp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", t, cp[0], cp[1], cp[2], cp[3], cp[4], 
		cp[5], cp[6], cp[7], cp[8], cp[9], cp[10], cp[11], cp[12]);
	}
	else fprintf(fp_cp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", t, cp[0], cp[1], cp[2], cp[3], cp[4], 
		cp[5], cp[6], cp[7], cp[8], cp[9], cp[10], cp[11], cp[12]);
	fclose(fp_cp);
	
  	FILE *fp_rp;
  	(void) sprintf(buf, "renal_parameters-%d.dat", k);
	fp_rp= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_rp, "Time Pr Pgh PB Pgo PB Pf Raa Rea Rr\n");
		fprintf(fp_rp, "%e %e %e %e %e %e %e %e %e\n", t, rp[0], rp[1], rp[2], rp[3], rp[4], rp[5], rp[6], rp[7]);
	}
	else fprintf(fp_rp, "%e %e %e %e %e %e %e %e %e\n", t, rp[0], rp[1], rp[2], rp[3], rp[4], rp[5], rp[6], rp[7]);
	fclose(fp_rp);
	
	/*FILE *fp_r;
  	(void) sprintf(buf, "resistance-%d.dat", k);
	fp_r= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_r, "Time Ra Rvr Rba Rbv Raa Rea Rr Rtp\n");
		fprintf(fp_r, "%e %e %e %e %e %e %e %e %e\n", t, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
	}
	else fprintf(fp_r, "%e %e %e %e %e %e %e %e %e\n", t, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
	fclose(fp_r);
	*/
	
  	FILE *fp_qw;
  	(void) sprintf(buf, "water_flow-%d.dat", k);
	fp_qw= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_qw, "Time q_win q_co q_vr q_rb q_gfilt q_twreab q_u delta_q\n");
		fprintf(fp_qw, "%e %e %e %e %e %e %e %e %e\n", t, qw[0], qw[1], qw[2], qw[3], qw[4], qw[5], qw[6], qw[7]);
	}
	else fprintf(fp_qw, "%e %e %e %e %e %e %e %e %e\n", t, qw[0], qw[1], qw[2], qw[3], qw[4], qw[5], qw[6], qw[7]);
	fclose(fp_qw);
	
	FILE *fp_qs;
  	(void) sprintf(buf, "sodium_flow-%d.dat", k);
	fp_qs= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_qs, "Time q_sodin q_filsod q_ptsodre q_mdsod q_dtsodre q_dtsod q_cdsodre q_usod\n");
		fprintf(fp_qs, "%e %e %e %e %e %e %e %e %e\n", t, qs[0], qs[1], qs[2], qs[3], qs[4], qs[5], qs[6], qs[7]);
	}
	else fprintf(fp_qs, "%e %e %e %e %e %e %e %e %e\n", t, qs[0], qs[1], qs[2], qs[3], qs[4], qs[5], qs[6], qs[7]);
	fclose(fp_qs);
	
	FILE *fp_c;
  	(void) sprintf(buf, "concentrations-%d.dat", k);
	fp_c= fopen(buf, "a+");
    if(i == 0)
    {
		fprintf(fp_c, "Time c_k c_sod c_adh c_ald c_anp c_ang c_renin, c_etoh\n");
		fprintf(fp_c, "%e %e %e %e %e %e %e %e %e\n", t, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);
	}
	else fprintf(fp_c, "%e %e %e %e %e %e %e %e %e\n", t, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);
	fclose(fp_c);
	
	FILE *fp_ref;
  	(void) sprintf(buf, "reflex-%d.dat", k);
	fp_ref= fopen(buf, "a+");
    if(i == 0)
    {
		fprintf(fp_ref, "Time rsna am tgf 3 4\n");
		fprintf(fp_ref, "%e %e %e %e %e %e\n", t, ref[0], ref[1], ref[2], ref[3], ref[4]);
	}
	else fprintf(fp_ref, "%e %e %e %e %e %e\n", t, ref[0], ref[1], ref[2], ref[3], ref[4]);
	fclose(fp_ref);
	
	FILE *fp_vol;
  	(void) sprintf(buf, "volume-%d.dat", k);
	fp_vol= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_vol, "Time Sod_Total Vecf Vb\n");
		fprintf(fp_vol, "%e %e %e %e\n", t, vol[0], vol[1], vol[2]);
	}
	else fprintf(fp_vol, "%e %e %e %e\n", t, vol[0], vol[1], vol[2]);
	fclose(fp_vol);
	
	FILE *fp_eff;
  	(void) sprintf(buf, "effects-%d.dat", k);
	fp_eff= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_eff, "Time eff0 eff1 eff2 eff3 eff4 eff5 eff6 eff7 eff8 eff9 eff10 eff11 eff12 eff13 eff14 eff15 eff16\n");
		fprintf(fp_eff, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", t, 
		eff[0], eff[1], eff[2], eff[3], eff[4], eff[5], eff[6], eff[7], eff[8], eff[9], 
		eff[10], eff[11], eff[12], eff[13], eff[14], eff[15], eff[16]);
	}
	else fprintf(fp_eff, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", t, 
		eff[0], eff[1], eff[2], eff[3], eff[4], eff[5], eff[6], eff[7], eff[8], eff[9], 
		eff[10], eff[11], eff[12], eff[13], eff[14], eff[15], eff[16]);
	fclose(fp_eff);
	
	FILE *fp_cn;
  	(void) sprintf(buf, "normal_concentrations-%d.dat", k);
	fp_cn= fopen(buf, "a+");
    if(i == 0)
    {
		fprintf(fp_cn, "Time cn_rsna cn_adh cn_adhs cn_ald cn_alds cn_renin cn_ptsodre cn_dtsodre cn_cdsodre frac_etoh\n");
		fprintf(fp_cn, "%e %e %e %e %e %e %e %e %e %e %e\n", t, cn[0], cn[1], cn[2], cn[3], cn[4], cn[5], cn[6], cn[7], cn[8], cn[9]);
	}
	else fprintf(fp_cn, "%e %e %e %e %e %e %e %e %e %e %e\n", t, cn[0], cn[1], cn[2], cn[3], cn[4], cn[5], cn[6], cn[7], cn[8], cn[9]);
	fclose(fp_cn);
	
	FILE *fp_vas;
  	(void) sprintf(buf, "vascularity-%d.dat", k);
	fp_vas= fopen(buf, "a+");
    if(i == 0) 
    {
		fprintf(fp_vas, "Time vas0 vasd vasf\n");
		fprintf(fp_vas, "%e %e %e %e\n", t, vas[0], vas[1], vas[2]);
	}
	else fprintf(fp_vas, "%e %e %e %e\n", t, vas[0], vas[1], vas[2]);
	fclose(fp_vas);
	
}

