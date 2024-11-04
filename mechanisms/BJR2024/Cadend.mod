: Calcium buffer shell model with instant Calcium handling for clustered channels by Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"

NEURON {
	SUFFIX Ca_conc_dend
	USEION ca READ ica WRITE eca, cai
	GLOBAL depth
	RANGE cai
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(um) = (micrometer)
	F = (faraday) (coulomb)
        R = 8.3134      (joule/degC)
}
    
FUNCTION ktf() (mV) {
        ktf = (1000)*R*(celsius +273.15)/(2*F)
}

PARAMETER {
	depth = .2 		(um)
	Fa = 96485.3365 (coulomb)
 	cao = 2.0 		(mM)
	tau = 10 (ms)
	cai0 = 7e-5 	(mM)
}

ASSIGNED {
	diam	(um)
	VSR (um)
        eca		(mV)
	ica		(mA/cm2) : instantaneous calcium current of l-type calcium channel
	B (mM*cm2/mA)
}

STATE {
	cai (mM) <1e-5>
}

BREAKPOINT {
	SOLVE state METHOD cnexp
        eca = ktf() * log(cao/cai) 
}

DERIVATIVE state {
    cai' = -ica * B - (cai-cai0)/tau
}


INITIAL {
        cai = cai0
	if (2*depth >= diam) {
		VSR = 0.25*diam : if diameter gets less than double the depth, the surface to volume ratio (here volume to surface ratio VSR) cannot be less than 1/(0.25diam) (instead of diam/(d*(diam-d)) )
	}else{
		VSR = depth*(1-depth/diam)
	}
	B = (1e4)/(2*Fa*VSR)
}