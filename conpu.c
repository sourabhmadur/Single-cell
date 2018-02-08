/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__conpu
#define _nrn_initial _nrn_initial__conpu
#define nrn_cur _nrn_cur__conpu
#define _nrn_current _nrn_current__conpu
#define nrn_jacob _nrn_jacob__conpu
#define nrn_state _nrn_state__conpu
#define _net_receive _net_receive__conpu 
#define settables settables__conpu 
#define states states__conpu 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define J_max_leak _p[0]
#define Jmax_serca _p[1]
#define Jmax_IP3 _p[2]
#define J_ERleak _p[3]
#define Jmax_NaCa _p[4]
#define Jmax_uni _p[5]
#define h _p[6]
#define cai _p[7]
#define FoRT _p[8]
#define capui _p[9]
#define Dcapui _p[10]
#define caeri _p[11]
#define Dcaeri _p[12]
#define cami _p[13]
#define Dcami _p[14]
#define Dh _p[15]
#define v _p[16]
#define _g _p[17]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_caeri	*_ppvar[1]._pval
#define _style_caer	*((int*)_ppvar[2]._pvoid)
#define _ion_capui	*_ppvar[3]._pval
#define _style_capu	*((int*)_ppvar[4]._pvoid)
#define _ion_cami	*_ppvar[5]._pval
#define _style_cam	*((int*)_ppvar[6]._pvoid)
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_J_leak(void);
 static void _hoc_J_SERCA(void);
 static void _hoc_J_ERout(void);
 static void _hoc_J_uni(void);
 static void _hoc_J_NaCa(void);
 static void _hoc_MWC(void);
 static void _hoc_settables(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_conpu", _hoc_setdata,
 "J_leak_conpu", _hoc_J_leak,
 "J_SERCA_conpu", _hoc_J_SERCA,
 "J_ERout_conpu", _hoc_J_ERout,
 "J_uni_conpu", _hoc_J_uni,
 "J_NaCa_conpu", _hoc_J_NaCa,
 "MWC_conpu", _hoc_MWC,
 "settables_conpu", _hoc_settables,
 0, 0
};
#define J_leak J_leak_conpu
#define J_SERCA J_SERCA_conpu
#define J_ERout J_ERout_conpu
#define J_uni J_uni_conpu
#define J_NaCa J_NaCa_conpu
#define MWC MWC_conpu
 extern double J_leak( _threadargsprotocomma_ double , double );
 extern double J_SERCA( _threadargsprotocomma_ double );
 extern double J_ERout( _threadargsprotocomma_ double , double , double );
 extern double J_uni( _threadargsprotocomma_ double );
 extern double J_NaCa( _threadargsprotocomma_ double );
 extern double MWC( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define F F_conpu
 double F = 96.4846;
#define IP3 IP3_conpu
 double IP3 = 0.0006;
#define K_trans K_trans_conpu
 double K_trans = 0.006;
#define K_act K_act_conpu
 double K_act = 0.00038;
#define K_Ca K_Ca_conpu
 double K_Ca = 0.003;
#define K_Na K_Na_conpu
 double K_Na = 9.4;
#define L L_conpu
 double L = 50;
#define Na_i Na_i_conpu
 double Na_i = 30;
#define Pmito Pmito_conpu
 double Pmito = 0.12871;
#define Per Per_conpu
 double Per = 0.1;
#define P_PU P_PU_conpu
 double P_PU = 0.001;
#define P_cyto P_cyto_conpu
 double P_cyto = 0.7;
#define R R_conpu
 double R = 8.3144;
#define Vol Vol_conpu
 double Vol = 1e-012;
#define b b_conpu
 double b = 0.5;
#define conc conc_conpu
 double conc = 0.001;
#define deltaPsi deltaPsi_conpu
 double deltaPsi = 164;
#define deltaPsi_star deltaPsi_star_conpu
 double deltaPsi_star = 91;
#define d_INH d_INH_conpu
 double d_INH = 0.0014;
#define d_ACT d_ACT_conpu
 double d_ACT = 0.001;
#define d_IP3 d_IP3_conpu
 double d_IP3 = 0.00025;
#define fm fm_conpu
 double fm = 0.0003;
#define fe fe_conpu
 double fe = 0.01;
#define fc fc_conpu
 double fc = 0.01;
#define k_serca k_serca_conpu
 double k_serca = 0.00042;
#define naa naa_conpu
 double naa = 2.8;
#define n n_conpu
 double n = 2;
#define tauh tauh_conpu
 double tauh = 4;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Vol_conpu", "litre",
 "k_serca_conpu", "mM",
 "IP3_conpu", "millimolar",
 "d_IP3_conpu", "millimolar",
 "d_ACT_conpu", "millimolar",
 "d_INH_conpu", "millimolar",
 "tauh_conpu", "s",
 "F_conpu", "microcoulomb/nanomole",
 "deltaPsi_star_conpu", "mV",
 "K_Na_conpu", "millimolar",
 "K_Ca_conpu", "millimolar",
 "Na_i_conpu", "millimolar",
 "deltaPsi_conpu", "mV",
 "K_act_conpu", "mM",
 "conc_conpu", "mM",
 "K_trans_conpu", "mM",
 "J_max_leak_conpu", "mM/s",
 "Jmax_serca_conpu", "mM/s",
 "Jmax_IP3_conpu", "1/s",
 "J_ERleak_conpu", "1/s",
 "Jmax_NaCa_conpu", "mM/s",
 "Jmax_uni_conpu", "mM/s",
 0,0
};
 static double cami0 = 0;
 static double caeri0 = 0;
 static double capui0 = 0;
 static double delta_t = 0.01;
 static double h0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "P_cyto_conpu", &P_cyto_conpu,
 "Vol_conpu", &Vol_conpu,
 "fc_conpu", &fc_conpu,
 "P_PU_conpu", &P_PU_conpu,
 "Per_conpu", &Per_conpu,
 "fe_conpu", &fe_conpu,
 "k_serca_conpu", &k_serca_conpu,
 "IP3_conpu", &IP3_conpu,
 "d_IP3_conpu", &d_IP3_conpu,
 "d_ACT_conpu", &d_ACT_conpu,
 "d_INH_conpu", &d_INH_conpu,
 "tauh_conpu", &tauh_conpu,
 "F_conpu", &F_conpu,
 "R_conpu", &R_conpu,
 "deltaPsi_star_conpu", &deltaPsi_star_conpu,
 "K_Na_conpu", &K_Na_conpu,
 "K_Ca_conpu", &K_Ca_conpu,
 "Na_i_conpu", &Na_i_conpu,
 "b_conpu", &b_conpu,
 "n_conpu", &n_conpu,
 "deltaPsi_conpu", &deltaPsi_conpu,
 "fm_conpu", &fm_conpu,
 "Pmito_conpu", &Pmito_conpu,
 "naa_conpu", &naa_conpu,
 "K_act_conpu", &K_act_conpu,
 "conc_conpu", &conc_conpu,
 "K_trans_conpu", &K_trans_conpu,
 "L_conpu", &L_conpu,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[7]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"conpu",
 "J_max_leak_conpu",
 "Jmax_serca_conpu",
 "Jmax_IP3_conpu",
 "J_ERleak_conpu",
 "Jmax_NaCa_conpu",
 "Jmax_uni_conpu",
 0,
 0,
 "h_conpu",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _caer_sym;
 static Symbol* _capu_sym;
 static Symbol* _cam_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	J_max_leak = 0.01;
 	Jmax_serca = 1.8333;
 	Jmax_IP3 = 50000;
 	J_ERleak = 1.66667;
 	Jmax_NaCa = 0.05;
 	Jmax_uni = 5000;
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 8, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 prop_ion = need_memb(_caer_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[1]._pval = &prop_ion->param[1]; /* caeri */
 	_ppvar[2]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for caer */
 prop_ion = need_memb(_capu_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* capui */
 	_ppvar[4]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for capu */
 prop_ion = need_memb(_cam_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[5]._pval = &prop_ion->param[1]; /* cami */
 	_ppvar[6]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for cam */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _conpu_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("caer", 2.0);
 	ion_reg("capu", 2.0);
 	ion_reg("cam", 2.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_caer_sym = hoc_lookup("caer_ion");
 	_capu_sym = hoc_lookup("capu_ion");
 	_cam_sym = hoc_lookup("cam_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 5);
  _extcall_thread = (Datum*)ecalloc(4, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 18, 8);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "caer_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "#caer_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "capu_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "#capu_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cam_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "#cam_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "cvodeieq");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 conpu C:/Users/admin/Desktop/single cell/conpu.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "conpu channel by madur";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int settables(_threadargsproto_);
 
#define _deriv1_advance _thread[0]._i
#define _dith1 1
#define _recurse _thread[2]._i
#define _newtonspace1 _thread[3]._pvoid
extern void* nrn_cons_newtonspace(int);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist2[4];
  static int _slist1[4], _dlist1[4];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   Dcapui = fc * ( ( ( J_NaCa ( _threadargscomma_ cami ) - J_uni ( _threadargscomma_ cami ) ) * ( ( Vol * Pmito ) / ( Vol * P_PU ) ) ) + ( ( J_ERout ( _threadargscomma_ capui , caeri , h ) - J_SERCA ( _threadargscomma_ capui ) ) * ( ( Vol * Per ) / ( Vol * P_PU ) ) ) - J_leak ( _threadargscomma_ capui , cai ) * ( ( Vol * P_cyto ) / ( Vol * P_PU ) ) ) ;
   Dcaeri = fe * ( J_SERCA ( _threadargscomma_ capui ) - J_ERout ( _threadargscomma_ capui , caeri , h ) ) ;
   Dcami = fm * ( J_uni ( _threadargscomma_ cami ) - J_NaCa ( _threadargscomma_ cami ) ) ;
   Dh = ( ( d_INH - ( ( capui + d_INH ) * h ) ) / tauh ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 Dcapui = Dcapui  / (1. - dt*( (( fc * ( ( ( J_NaCa ( _threadargscomma_ cami ) - J_uni ( _threadargscomma_ cami ) ) * ( ( Vol * Pmito ) / ( Vol * P_PU ) ) ) + ( ( J_ERout ( _threadargscomma_ ( capui  + .001) , caeri , h ) - J_SERCA ( _threadargscomma_ ( capui  + .001) ) ) * ( ( Vol * Per ) / ( Vol * P_PU ) ) ) - J_leak ( _threadargscomma_ ( capui  + .001) , cai ) * ( ( Vol * P_cyto ) / ( Vol * P_PU ) ) ) ) - ( fc * ( ( ( J_NaCa ( _threadargscomma_ cami ) - J_uni ( _threadargscomma_ cami ) ) * ( ( Vol * Pmito ) / ( Vol * P_PU ) ) ) + ( ( J_ERout ( _threadargscomma_ capui , caeri , h ) - J_SERCA ( _threadargscomma_ capui ) ) * ( ( Vol * Per ) / ( Vol * P_PU ) ) ) - J_leak ( _threadargscomma_ capui , cai ) * ( ( Vol * P_cyto ) / ( Vol * P_PU ) ) )  )) / .001 )) ;
 Dcaeri = Dcaeri  / (1. - dt*( (( fe * ( J_SERCA ( _threadargscomma_ capui ) - J_ERout ( _threadargscomma_ capui , ( caeri  + .001) , h ) ) ) - ( fe * ( J_SERCA ( _threadargscomma_ capui ) - J_ERout ( _threadargscomma_ capui , caeri , h ) )  )) / .001 )) ;
 Dcami = Dcami  / (1. - dt*( (( fm * ( J_uni ( _threadargscomma_ ( cami  + .001) ) - J_NaCa ( _threadargscomma_ ( cami  + .001) ) ) ) - ( fm * ( J_uni ( _threadargscomma_ cami ) - J_NaCa ( _threadargscomma_ cami ) )  )) / .001 )) ;
 Dh = Dh  / (1. - dt*( ( ( ( ( - ( ( ( capui + d_INH ) )*( 1.0 ) ) ) ) ) / tauh ) )) ;
 return 0;
}
 /*END CVODE*/
 
static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0; int error = 0;
 { double* _savstate1 = _thread[_dith1]._pval;
 double* _dlist2 = _thread[_dith1]._pval + 4;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 4; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = nrn_newton_thread(_newtonspace1, 4,_slist2, _p, states, _dlist2, _ppvar, _thread, _nt);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   Dcapui = fc * ( ( ( J_NaCa ( _threadargscomma_ cami ) - J_uni ( _threadargscomma_ cami ) ) * ( ( Vol * Pmito ) / ( Vol * P_PU ) ) ) + ( ( J_ERout ( _threadargscomma_ capui , caeri , h ) - J_SERCA ( _threadargscomma_ capui ) ) * ( ( Vol * Per ) / ( Vol * P_PU ) ) ) - J_leak ( _threadargscomma_ capui , cai ) * ( ( Vol * P_cyto ) / ( Vol * P_PU ) ) ) ;
   Dcaeri = fe * ( J_SERCA ( _threadargscomma_ capui ) - J_ERout ( _threadargscomma_ capui , caeri , h ) ) ;
   Dcami = fm * ( J_uni ( _threadargscomma_ cami ) - J_NaCa ( _threadargscomma_ cami ) ) ;
   Dh = ( ( d_INH - ( ( capui + d_INH ) * h ) ) / tauh ) ;
   {int _id; for(_id=0; _id < 4; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 
double J_leak ( _threadargsprotocomma_ double _lcapui , double _lcai ) {
   double _lJ_leak;
 _lJ_leak = J_max_leak * ( _lcapui - _lcai ) ;
   
return _lJ_leak;
 }
 
static void _hoc_J_leak(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  J_leak ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double J_SERCA ( _threadargsprotocomma_ double _lcapui ) {
   double _lJ_SERCA;
 _lJ_SERCA = Jmax_serca * ( ( pow( _lcapui , 2.0 ) ) / ( ( pow( k_serca , 2.0 ) ) + ( pow( _lcapui , 2.0 ) ) ) ) ;
   
return _lJ_SERCA;
 }
 
static void _hoc_J_SERCA(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  J_SERCA ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double J_ERout ( _threadargsprotocomma_ double _lcapui , double _lcaeri , double _lh ) {
   double _lJ_ERout;
 _lJ_ERout = ( Jmax_IP3 * ( pow( ( IP3 / ( IP3 + d_IP3 ) ) , 3.0 ) ) * ( pow( ( _lcapui / ( _lcapui + d_ACT ) ) , 3.0 ) ) * ( pow( _lh , 3.0 ) ) + J_ERleak ) * ( _lcaeri - _lcapui ) ;
   
return _lJ_ERout;
 }
 
static void _hoc_J_ERout(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  J_ERout ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double J_NaCa ( _threadargsprotocomma_ double _lcami ) {
   double _lJ_NaCa;
 _lJ_NaCa = Jmax_NaCa * ( exp ( b * FoRT * ( deltaPsi - deltaPsi_star ) ) / ( ( ( 1.0 + ( pow( ( K_Na / Na_i ) , n ) ) ) ) * ( 1.0 + ( K_Ca / _lcami ) ) ) ) ;
   
return _lJ_NaCa;
 }
 
static void _hoc_J_NaCa(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  J_NaCa ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double MWC ( _threadargsprotocomma_ double _lcapui ) {
   double _lMWC;
 _lMWC = conc * ( ( _lcapui / K_trans ) * ( pow( ( 1.0 + _lcapui / K_trans ) , 3.0 ) ) ) / ( ( pow( ( 1.0 + _lcapui / K_trans ) , 4.0 ) ) + ( L / ( pow( ( 1.0 + _lcapui / K_act ) , naa ) ) ) ) ;
   
return _lMWC;
 }
 
static void _hoc_MWC(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  MWC ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double J_uni ( _threadargsprotocomma_ double _lcami ) {
   double _lJ_uni;
 _lJ_uni = Jmax_uni * ( ( MWC ( _threadargscomma_ capui ) - ( _lcami * exp ( - 2.0 * FoRT * ( deltaPsi - deltaPsi_star ) ) ) ) * ( 2.0 * FoRT * ( deltaPsi - deltaPsi_star ) ) / ( 1.0 - exp ( - 2.0 * FoRT * ( deltaPsi - deltaPsi_star ) ) ) ) ;
   
return _lJ_uni;
 }
 
static void _hoc_J_uni(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  J_uni ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  settables ( _threadargsproto_ ) {
   FoRT = F / ( R * ( celsius + 273.0 ) ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  capui = _ion_capui;
  capui = _ion_capui;
  cami = _ion_cami;
  cami = _ion_cami;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  _ion_caeri = caeri;
  _ion_capui = capui;
  _ion_cami = cami;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 	_pv[1] = &(_ion_caeri);
 	_pv[0] = &(_ion_capui);
 	_pv[2] = &(_ion_cami);
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  capui = _ion_capui;
  capui = _ion_capui;
  cami = _ion_cami;
  cami = _ion_cami;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[_dith1]._pval = (double*)ecalloc(8, sizeof(double));
   _newtonspace1 = nrn_cons_newtonspace(4);
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[_dith1]._pval));
   nrn_destroy_newtonspace(_newtonspace1);
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_caer_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_capu_sym, _ppvar, 3, 1);
   nrn_update_ion_pointer(_cam_sym, _ppvar, 5, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
 {
   settables ( _threadargs_ ) ;
   capui = 0.0000902 ;
   caeri = 0.007299 ;
   cami = 0.000136 ;
   h = 0.9397 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  capui = _ion_capui;
  capui = _ion_capui;
  cami = _ion_cami;
  cami = _ion_cami;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_caeri = caeri;
  nrn_wrote_conc(_caer_sym, (&(_ion_caeri)) - 1, _style_caer);
  _ion_capui = capui;
  nrn_wrote_conc(_capu_sym, (&(_ion_capui)) - 1, _style_capu);
  _ion_cami = cami;
  nrn_wrote_conc(_cam_sym, (&(_ion_cami)) - 1, _style_cam);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  capui = _ion_capui;
  capui = _ion_capui;
  cami = _ion_cami;
  cami = _ion_cami;
 {  _deriv1_advance = 1;
 derivimplicit_thread(4, _slist1, _dlist1, _p, states, _ppvar, _thread, _nt);
_deriv1_advance = 0;
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   }
  _ion_caeri = caeri;
  _ion_capui = capui;
  _ion_cami = cami;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(capui) - _p;  _dlist1[0] = &(Dcapui) - _p;
 _slist1[1] = &(caeri) - _p;  _dlist1[1] = &(Dcaeri) - _p;
 _slist1[2] = &(cami) - _p;  _dlist1[2] = &(Dcami) - _p;
 _slist1[3] = &(h) - _p;  _dlist1[3] = &(Dh) - _p;
 _slist2[0] = &(cami) - _p;
 _slist2[1] = &(capui) - _p;
 _slist2[2] = &(caeri) - _p;
 _slist2[3] = &(h) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
