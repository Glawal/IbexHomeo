#ifndef __IBEX_INPUTVAR_H__
#define __IBEX_INPUTVAR_H__

#include "ibex.h"
#include "args.hxx"
#include <string>
//#include <chrono>
//#include <condition_variable>
//#include "ibex_ThreadOptimizer.h"
#include "ibex_Interval.h"

using namespace std;
namespace ibex {

struct InputVar{
public:
	string name;
	Interval range;

	InputVar();

	InputVar(const string& namevar, const Interval& rangevar);

	InputVar(const InputVar& copy);

};

inline InputVar::InputVar(): name(""), range(Interval()){

}

inline InputVar::InputVar(const string& namevar, const Interval& rangevar) : name(namevar),range(rangevar){
}

inline InputVar::InputVar( const InputVar& copy): name(copy.name), range(copy.range){
	
}

}

#endif // __IBEX_INPUTVAR_H__