NEURON {

    POINT_PROCESS Gap
    POINTER vgap
    RANGE g, i
    NONSPECIFIC_CURRENT i
}

UNITS {

  (nA) = (nanoamp)
  (mV) = (millivolt)
  (nS) = (nanosiemens)
}

PARAMETER { g = 0 (uS) }
    
ASSIGNED {

    v    (mV)
    vgap (mV)
    i    (nA)
}
 
BREAKPOINT { 

  if (g>0) {i = g * (v-vgap) }

}
