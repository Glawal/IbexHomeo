#include "ibex_HeuristicOptimizer.h"

#include "ibex_CtcHC4.h"
#include "ibex_CtcAcid.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcLinearRelax.h"
#include "ibex_CellDoubleHeap.h"
#include "ibex_SmearFunction.h"
#include "ibex_LSmear.h"
#include "ibex_LoupFinderDefault.h"
#include "ibex_LoupFinderCertify.h"
#include "ibex_LinearizerCompo.h"
#include "ibex_Array.h"
#include "ibex_Random.h"
#include "ibex_CellBeamSearch.h"
#include "ibex_CellHeap.h"
#include "ibex_Optimizer.h"

using namespace std;

namespace ibex {

#define NORMALIZED_SYSTEM_TAG 1
#define EXTENDED_SYSTEM_TAG 2

#define default_relax_ratio 0.2
#define default_bisect_ratio 0.5

// The two next functions are necessary because we need
// the normalized and extended system to build
// arguments of the base class constructor (ctc, bsc, loup finder, etc.)
// and we don't know which argument is evaluated first

NormalizedSystem& HeuristicOptimizer::get_norm_sys(const System& sys, double eps_h) {
	if (found(NORMALIZED_SYSTEM_TAG)) {
		return get<NormalizedSystem>(NORMALIZED_SYSTEM_TAG);
	} else {
		return rec(new NormalizedSystem(sys,eps_h), NORMALIZED_SYSTEM_TAG);
	}
}

ExtendedSystem& HeuristicOptimizer::get_ext_sys(const System& sys, double eps_h) {
	if (found(EXTENDED_SYSTEM_TAG)) {
		return get<ExtendedSystem>(EXTENDED_SYSTEM_TAG);
	} else {
		return rec(new ExtendedSystem(sys,eps_h), EXTENDED_SYSTEM_TAG);
	}
}

Bsc& HeuristicOptimizer::helperBsc(const System& sys, double eps_h, double eps_x, HeuristicOptimizer::Bsctype i_bsctype){
	switch(i_bsctype){
		case HeuristicOptimizer::Bsctype::LSMEAR:
			return rec(new LSmear(get_ext_sys(sys,eps_h),eps_x,rec(new OptimLargestFirst(get_ext_sys(sys,eps_h).goal_var(),eps_x,default_bisect_ratio))));
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_SUM:
			return rec(new SmearSum(get_ext_sys(sys,eps_h),eps_x,rec(new OptimLargestFirst(get_ext_sys(sys,eps_h).goal_var(),eps_x,default_bisect_ratio))));
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_SUM_REL:
			return rec(new SmearSumRelative(get_ext_sys(sys,eps_h),eps_x,rec(new OptimLargestFirst(get_ext_sys(sys,eps_h).goal_var(),eps_x,default_bisect_ratio))));
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_MAX:
			return rec(new SmearMax(get_ext_sys(sys,eps_h),eps_x,rec(new OptimLargestFirst(get_ext_sys(sys,eps_h).goal_var(),eps_x,default_bisect_ratio))));
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_MAX_REL:
			return rec(new SmearMaxRelative(get_ext_sys(sys,eps_h),eps_x,rec(new OptimLargestFirst(get_ext_sys(sys,eps_h).goal_var(),eps_x,default_bisect_ratio))));
			break;
		default:
			cerr << "not an implemented bisector." << endl;
			exit(1);
	}
}

CellBufferOptim& HeuristicOptimizer::helperBuffer(const System& sys, double eps_h, HeuristicOptimizer::Buffertype i_buffertype){
	switch(i_buffertype){
		case HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH:
			return (CellBufferOptim&) rec (new  CellBeamSearch (
								    (CellHeap&) rec (new CellHeap (get_ext_sys(sys,eps_h))),
								    (CellHeap&) rec (new CellHeap (get_ext_sys(sys,eps_h))),
								    get_ext_sys(sys,eps_h)));
			break;
		case HeuristicOptimizer::Buffertype::DOUBLE_HEAP:
			return (CellBufferOptim&) rec (new  CellDoubleHeap (get_ext_sys(sys,eps_h)));
			break;
		default:
			cerr << "not an implemented cellbuffer." << endl;
			exit(1);
	}
}

HeuristicOptimizer::HeuristicOptimizer(const System& sys, double rel_eps_f, double abs_eps_f, double eps_h, bool rigor, bool inHC4, double random_seed, double eps_x, 
	HeuristicOptimizer::Bsctype i_bsctype, HeuristicOptimizer::Buffertype i_buffertype, ThreadOptimizer::Obj t_obj) :
	ThreadOptimizer(sys.nb_var,
 		ctc(get_ext_sys(sys,eps_h)), // warning: we don't know which argument is evaluated first
		helperBsc(sys,eps_h,eps_x,i_bsctype),
		rec(rigor? (LoupFinder*) new LoupFinderCertify(sys,rec(new LoupFinderDefault(get_norm_sys(sys,eps_h), inHC4))) :
			(LoupFinder*) new LoupFinderDefault(get_norm_sys(sys,eps_h), inHC4)),
    	helperBuffer(sys,eps_h,i_buffertype),
		get_ext_sys(sys,eps_h).goal_var(),
		t_obj,
		eps_x,
		rel_eps_f,
		abs_eps_f)
{
	RNG::srand(random_seed);
	buffertype=i_buffertype;
	bsctype=i_bsctype;
}

Ctc&  HeuristicOptimizer::ctc(const ExtendedSystem& ext_sys) {
	Array<Ctc> ctc_list(3);

	// first contractor on ext_sys : incremental HC4 (propag ratio=0.01)
	ctc_list.set_ref(0, rec(new CtcHC4 (ext_sys,0.01,true)));
	// second contractor on ext_sys : "Acid" with incremental HC4 (propag ratio=0.1)
	ctc_list.set_ref(1, rec(new CtcAcid (ext_sys,rec(new CtcHC4 (ext_sys,0.1,true)),true)));
	// the last contractor is "XNewton"

	if (ext_sys.nb_ctr > 1) {
		ctc_list.set_ref(2,rec(new CtcFixPoint
				(rec(new CtcCompo(
						rec(new CtcLinearRelax(ext_sys)),
						rec(new CtcHC4(ext_sys,0.01)))), default_relax_ratio)));
	} else {
		ctc_list.set_ref(2,rec(new CtcLinearRelax(ext_sys)));
	}

	return rec(new CtcCompo(ctc_list));
}

}//namespace ibex
