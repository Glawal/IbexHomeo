#ifndef __IBEX_REAL_TIMER_H__
#define __IBEX_REAL_TIMER_H__

#include "ibex_Exception.h"
#include <chrono>

using namespace std;

namespace ibex {

class RealTimeOutException : public Exception { };

class RealTimer {
 public:
	RealTimer();
	RealTimer(RealTimer& copy);

	void start();
	void stop();
	double get_time();
	void check(double timeout);
	bool get_active();
	std::chrono::high_resolution_clock::time_point get_start_time();
	std::chrono::high_resolution_clock::time_point get_current_time();
protected:
	  std::chrono::high_resolution_clock::time_point start_time;
	  std::chrono::high_resolution_clock::time_point current_time;
	  bool active;
};

}

#endif // __IBEX_REAL_TIMER_H__