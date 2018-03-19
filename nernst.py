
from math import log
T=310
F = 96.4846         # microcoulomb_per_nanomole
R = 8.3144


cli_ = 88
clo_ = 134

ki_= 120	#23
ko_= 7	#23

nai_ = 30
nao_ = 137

ecl_=(R*T/F)*log(cli_/clo_)
ena_=(R*T/F)*log(nao_/nai_)
ek_=(R*T/F)*log(ko_/ki_)
