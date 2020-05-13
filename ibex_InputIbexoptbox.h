#ifndef __IBEX_INPUT_IBEXOPTBOX_H__
#define __IBEX_INPUT_IBEXOPTBOX_H__

#include "ibex.h"
#include "args.hxx"
//#include <chrono>
//#include <condition_variable>
//#include "ibex_ThreadOptimizer.h"
//#include "ibex_Interval.h"
//#include <string>

using namespace std;
namespace ibex {

struct InputIbexoptbox{
public:
	double hom, rel_eps_f, abs_eps_f, eps_h, eps_x, timeoutloop;
	bool timeout_b;
	double timeout_v, random_seed, lastloup;
	bool rigor, trace, quiet;

	InputIbexoptbox(const double& homeostasis, const double& rel_eps_f_v, const double& abs_eps_f_v, const double& eps_h_v, 
		const double& eps_x_v, const double& timeoutloop_v, const bool& timeout_bool, const double& timeout_value, const double& random_seed_v, 
		const double initial_loup, const bool& rigor_b, const bool& trace_b, const bool& quiet_b);

	InputIbexoptbox(const InputIbexoptbox& copy);
	
};

inline InputIbexoptbox::InputIbexoptbox(const double& homeostasis, const double& rel_eps_f_v, const double& abs_eps_f_v, const double& eps_h_v, 
		const double& eps_x_v, const double& timeoutloop_v, const bool& timeout_bool, const double& timeout_value, const double& random_seed_v, 
		const double initial_loup, const bool& rigor_b, const bool& trace_b, const bool& quiet_b) :
	hom(homeostasis), rel_eps_f(rel_eps_f_v), abs_eps_f(abs_eps_f_v), eps_h(eps_h_v), eps_x(eps_x_v), timeoutloop(timeoutloop_v), timeout_b(timeout_bool), 
	timeout_v(timeout_value), random_seed(random_seed_v), lastloup(initial_loup), rigor(rigor_b), trace(trace_b), quiet(quiet_b){
}

inline InputIbexoptbox::InputIbexoptbox(const InputIbexoptbox& copy):
	hom(copy.hom), rel_eps_f(copy.rel_eps_f), abs_eps_f(copy.abs_eps_f), eps_h(copy.eps_h), eps_x(copy.eps_x), timeoutloop(copy.timeoutloop), timeout_b(copy.timeout_b), 
	timeout_v(copy.timeout_v), random_seed(copy.random_seed), lastloup(copy.lastloup), rigor(copy.rigor), trace(copy.trace), quiet(copy.quiet){
}

}


#endif // __IBEX_INPUT_IBEXOPTBOX_H__