: Na channel


NEURON {
	SUFFIX Na_Ket2001
	USEION na READ nai, nao WRITE ina
	RANGE pbar, p, sh
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
        FARADAY = (faraday) (coulomb)
        R = (k-mole) (joule/degC)

}

PARAMETER {
  pbar = 0 (cm/s)
}

ASSIGNED {
        celsius
        p               (cm/s)
	v		(mV)
	ina		(mA/cm2)
        am	        (/ms)	
	bm		(/ms)
        nai             (mM)
        nao             (mM)
}

STATE { m }

INITIAL { 
	rates(v)
	m = am/(am + bm)
}
    
BREAKPOINT {
  SOLVE states METHOD cnexp
  p = pbar * m 
  ina = p * ghk(v, nai, nao)
:  printf("v = %g m = %g ina = %g nai = %g nao = %g am = %g bm = %g\n", v, m, ina, nai, nao, am, bm)
} 

DERIVATIVE states {
   rates(v)
   m' = am*(1 - m) - bm*m  
   printf("at t %g: v = %g m = %g m rhs = %g\n", t, v, m, am*(1 - m) - bm*m)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
        LOCAL z, eci, eco
        z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
        eco = co*efun(z)
        eci = ci*efun(-z)
        ghk = (.001)*2*FARADAY*(eci - eco)
        }

FUNCTION efun(z) {
        if (fabs(z) < 1e-4) {
                efun = 1 - z/2
        }else{
                efun = z/(exp(z) - 1)
        }
}


PROCEDURE rates(v (mV)) { LOCAL q10
  q10 = 3^((celsius - 20)/10)  : the fit was performed at 20C
        
  
  am = q10*1/(1 + exp(0.3*(-44.04 - v)))
  bm = q10*1.1/(1 + exp(-0.11*(-75.42 - v)))
}


 
