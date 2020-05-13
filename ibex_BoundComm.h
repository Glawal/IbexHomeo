#ifndef __IBEX_BOUND_COMM_H__
#define __IBEX_BOUND_COMM_H__

#include "ibex.h"
#include "args.hxx"
#include <chrono>
#include <condition_variable>
#include "ibex_ThreadOptimizer.h"
#include "ibex_Interval.h"
#include <string>

using namespace std;
namespace ibex {

struct BoundComm{
public:
	vector<IntervalVector> loup_point_global;//mutex obligatoire;
	atomic<int> nbopdone_for_min;
	atomic<int> nbopdone_for_max;
	atomic<double> uploformin_global;
	atomic<double> loupformin_global;
	atomic<double> uploformax_global;
	atomic<double> loupformax_global;
	atomic<double> uploforsubmin_global;
	atomic<double> loupforsubmin_global;
	atomic<double> uploforsubmax_global;
	atomic<double> loupforsubmax_global;
	atomic<double> stop_timeout_formin;
	atomic<double> stop_timeout_formax;
	atomic<ThreadOptimizer::Status> statusformin;
	atomic<ThreadOptimizer::Status> statusformax;
	condition_variable check_who_finished_first;
	bool min_finished_first;

	BoundComm();
};

inline BoundComm::BoundComm() :
	loup_point_global{},
	nbopdone_for_min(0),
	nbopdone_for_max(1),
	uploformin_global(NEG_INFINITY),
	loupformin_global(POS_INFINITY),
	uploformax_global(NEG_INFINITY),
	loupformax_global(POS_INFINITY),
	uploforsubmin_global(NEG_INFINITY),
	loupforsubmin_global(POS_INFINITY),
	uploforsubmax_global(NEG_INFINITY),
	loupforsubmax_global(POS_INFINITY),
	stop_timeout_formin(-1),
	stop_timeout_formax(-1),
	statusformin(ThreadOptimizer::Status::TIME_OUT),
	statusformax(ThreadOptimizer::Status::TIME_OUT),
	check_who_finished_first{},
	min_finished_first(true)
	{
}

inline ostream& operator<<( ostream &flux, BoundComm const& comm){
	if(comm.loup_point_global.empty()){
		flux << "loup_point_global is empty" << endl;
	}
	else{
		flux << "loup_point_global is not empty" << endl;
	}
	flux << "uploformin_global = " << comm.uploformin_global.load() << endl;
	flux << "loupformin_global = " << comm.loupformin_global.load() << endl;
	flux << "uploformax_global = " << comm.uploformax_global.load() << endl;
	flux << "loupformax_global = " << comm.loupformax_global.load() << endl;
	return flux;
}

}

#endif // __IBEX_BOUND_COMM_H__
