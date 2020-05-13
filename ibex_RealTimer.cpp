#include "ibex_RealTimer.h"
#include <chrono>

namespace ibex {

RealTimer::RealTimer() : active(false){
	start_time=std::chrono::high_resolution_clock::now();
	current_time=std::chrono::high_resolution_clock::now();
}

RealTimer::RealTimer(RealTimer& copy){
	active=copy.get_active();
	start_time=copy.get_start_time();
	current_time=copy.get_current_time(); 
}

void RealTimer::start(){
	if(!active){
		active=true;
	}
	start_time=std::chrono::high_resolution_clock::now();
}

void RealTimer::stop(){
	active=false;
	current_time=std::chrono::high_resolution_clock::now();
}

double RealTimer::get_time(){
	if(active){
		current_time=std::chrono::high_resolution_clock::now();
	}
	return std::chrono::duration_cast<std::chrono::duration<double>>(current_time-start_time).count();
}

void RealTimer::check(double timeout){
	if (RealTimer::get_time() >= timeout) throw RealTimeOutException();
}

bool RealTimer::get_active(){
	return active;
}

std::chrono::high_resolution_clock::time_point RealTimer::get_start_time(){
	return start_time;
}

std::chrono::high_resolution_clock::time_point RealTimer::get_current_time(){
	return current_time;
}

}