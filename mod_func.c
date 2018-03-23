#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ERG_reg();
extern void _Kb_reg();
extern void _Na_reg();
extern void _ano1_reg();
extern void _bk_reg();
extern void _cacl_reg();
extern void _concyto_reg();
extern void _conpu_reg();
extern void _gap_reg();
extern void _kv_11_reg();
extern void _ltype_reg();
extern void _nscc_reg();
extern void _pmca_reg();
extern void _vddr_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ERG.mod");
fprintf(stderr," Kb.mod");
fprintf(stderr," Na.mod");
fprintf(stderr," ano1.mod");
fprintf(stderr," bk.mod");
fprintf(stderr," cacl.mod");
fprintf(stderr," concyto.mod");
fprintf(stderr," conpu.mod");
fprintf(stderr," gap.mod");
fprintf(stderr," kv_11.mod");
fprintf(stderr," ltype.mod");
fprintf(stderr," nscc.mod");
fprintf(stderr," pmca.mod");
fprintf(stderr," vddr.mod");
fprintf(stderr, "\n");
    }
_ERG_reg();
_Kb_reg();
_Na_reg();
_ano1_reg();
_bk_reg();
_cacl_reg();
_concyto_reg();
_conpu_reg();
_gap_reg();
_kv_11_reg();
_ltype_reg();
_nscc_reg();
_pmca_reg();
_vddr_reg();
}
