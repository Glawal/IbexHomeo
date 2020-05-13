#ifndef __IBEX_HEURISTIC_OPTIMIZER_H__
#define __IBEX_HEURISTIC_OPTIMIZER_H__

#include "ibex_Optimizer.h"
#include "ibex_CtcCompo.h"
#include "ibex_Memory.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_ExtendedSystem.h"
//#include "ibex_inputibexoptbox.h"
#include "ibex_BoundComm.h"
#include "ibex_InputVar.h"
#include "ibex_InputIbexoptbox.h"
#include "ibex_OptimizerConfig.h"
#include <chrono>

extern ibex::BoundComm global_comm;
extern std::mutex mutex_global;
extern std::mutex mutex_event;

using namespace std;
namespace ibex {

class HeuristicOptimizer : private Memory, public ThreadOptimizer {
public:
	/**
	 * \brief Create a default optimizer.
	 *
	 * \param sys         - The system to optimize.
	 * \param rel_eps_f   - Relative precision on the objective.
	 * \param abs_eps_f   - Absolute precision on the objective.
	 * \param eps_h       - Equality thickness.
	 * \param rigor       - If true, feasibility of equalities is certified. By default:
	 *                      false.
	 * \param inHC4       - If true, feasibility is also tried with LoupFinderInHC4.
	 * \param random_seed - The sequence of random numbers is reinitialized with
	 *                      this seed before calling optimize(..) (useful for
	 *                      reproducibility). Set by default to #default_random_seed.
	 * \param eps_x       - Stopping criterion for box splitting (absolute precision).
	 *                      (**deprecated**).
	 */

	enum class Bsctype {LSMEAR,SMEAR_SUM,SMEAR_SUM_REL,SMEAR_MAX,SMEAR_MAX_REL};

	enum class Buffertype {CELL_BEAM_SEARCH,DOUBLE_HEAP};	

	HeuristicOptimizer(const System& sys, 
		double rel_eps_f=OptimizerConfig::default_rel_eps_f,
		double abs_eps_f=OptimizerConfig::default_abs_eps_f,
		double eps_h=NormalizedSystem::default_eps_h,
		bool rigor=false, bool inHC4=true,
		double random_seed=default_random_seed,
    	double eps_x=OptimizerConfig::default_eps_x,
    	Bsctype i_bsctype=Bsctype::LSMEAR,
    	Buffertype i_buffertype=Buffertype::CELL_BEAM_SEARCH,
    	ThreadOptimizer::Obj t_obj=ThreadOptimizer::OBJ_MIN);

    static constexpr double default_random_seed = 1.0;

    Bsctype bsctype;

	Buffertype buffertype;

protected:
	
	/**
     * The contractor: HC4 + acid(HC4) + X-Newton
     */
	Ctc& ctc(const ExtendedSystem& ext_sys);

	Bsc& helperBsc(const System& sys, double eps_h, double eps_x, HeuristicOptimizer::Bsctype i_bsctype);

	CellBufferOptim& helperBuffer(const System& sys, double eps_h, HeuristicOptimizer::Buffertype i_buffertype);

	NormalizedSystem& get_norm_sys(const System& sys, double eps_h);

	ExtendedSystem& get_ext_sys(const System& sys, double eps_h);
};

struct Output{
	Optimizer::Status issuccessful;
	double o_lastloup;
	double o_lastuplo;
	IntervalVector o_lastlouppoint;
	std::chrono::duration<double> timer;
	HeuristicOptimizer::Bsctype bsctype;
	HeuristicOptimizer::Buffertype buffertype;
};

inline ostream& operator<<( ostream &flux, HeuristicOptimizer::Bsctype const& i_bsctype){
	switch(i_bsctype){
		case HeuristicOptimizer::Bsctype::LSMEAR:
			flux << "LSMEAR";
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_SUM:
			flux << "SMEAR_SUM";
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_SUM_REL:
			flux << "SMEAR_SUM_REL";
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_MAX:
			flux << "SMEAR_MAX";
			break;
		case HeuristicOptimizer::Bsctype::SMEAR_MAX_REL:
			flux << "SMEAR_MAX_REL";
			break;
	}
	return flux;
}

inline ostream& operator<<( ostream &flux, HeuristicOptimizer::Buffertype const& i_buffertype){
	switch(i_buffertype){
		case HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH:
			flux << "CELL_BEAM_SEARCH";
			break;
		case HeuristicOptimizer::Buffertype::DOUBLE_HEAP:
			flux << "DOUBLE_HEAP";
			break;
	}
	return flux;
}

}

#endif // __IBEX_HEURISTIC_OPTIMIZER_H__