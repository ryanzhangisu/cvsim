# cvsim
CVSim source code and modifications for simulating cardiovascular flow waveforms
[INTRODUCTION]
This code combines the original CVSim code with renal function model in the following way:
  1.	CVSim code is original without changes. The code is heavily commented by the original authors. A brief introduction of the structure of the CVSim:
    a.	Main.c is the highest-level file, with the control of on/off of the reflex system and the gravity orthostatic modules.
    b.	Simulation.c is the main file for the simulation and calls every other submodules, including the renal function we added.
    c.	Initial.c is the file defines all initial parameters
    d.	Estimate.c is the file that gives an estimation of the dependent variable before the iteration starts.
    e.	Reflex.c is the file defines baroreflex functions and heart rhythm and cardiac timing
    f.	Rkqc.c is the file that contains the RK routines
    g.	Renal.c is the renal function module we added.
    2.	Renal functions models are contained in separate files renal.c and renal.h
  3.	Renal function module contains subroutines: renal_function(), renal_eqns(), rk4_routine(), renal_init(), renal_write() that are declared in simulation.h and called in simulation.c. The subroutine functions are mostly the same with those in the original renal function model, but with different subroutine names.
  4.	Only changes in renal function module: The renal_eqns() module removes the do loop for calculating mean arterial pressure and mean right atrial pressure, because it takes these parameters directly from the CVSim code, as a part of the coupling scheme.
  Please refer to the master thesis (Heldt 2004) for the detailed parameters and equations and to (Zhang et al. 2014) for the application of the coupled model.
[USAGE]
  1.	If you only wish to run the original CVSim  function, simply comment out the renal function subroutines called by the simulator_ptr() function in simulation.c.
  2.	The time-step of CVSim is fixed at 0.001 sec. The time scale between it and the renal function model is hugely different. In the current application, the CVSim parameters were averaged and passed onto renal function on a beat-by-beat basis, which does not provide a true match of the time scales.
  3.	One could choose to turn on/off the barorefex control by setting values (1 vs 0) of ABReflexon and CPReflexon in main.c. 
  4.	For studies of sodium and water intakes, the input for intake values are defined as sodin and waterin, in the simulator_ptr() subroutine. The values should be given as an averaged value with units of mEq/min and L/min, respectively.
REFERENCE:
Heldt, T., 2004. Computational models of cardiovascular response to orthostatic stress.
Zhang, Y., Liou, W.W. & Gupta, V., 2014. Modeling of High Sodium Intake Effects on Left Ventricular Hypertrophy. Computers in Biology and Medicine, 58, pp.31–39.


