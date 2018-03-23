NEURON {
	POINT_PROCESS gap
	NONSPECIFIC_CURRENT i
	RANGE r, i
	POINTER vgap
}
PARAMETER {
	v (millivolt)
	vgap (millivolt)
	r = 30.6 (megohm)		:30.6
}
ASSIGNED {
	i (nanoamp)
}
BREAKPOINT {
	i = (v - vgap)/r
}

