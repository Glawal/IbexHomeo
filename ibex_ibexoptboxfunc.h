#ifndef __IBEX__IBEXOPTBOXFUNC_H__
#define __IBEX__IBEXOPTBOXFUNC_H__

#include "ibex.h"
#include "args.hxx"
#include <sstream>
#include <vector>
#include <regex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <thread>
#include <chrono>
#include <future>
#include <pthread.h>
#include <assert.h>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include "ibex_HeuristicOptimizer.h"
//#include "ibex_LowLevelOptimizer.cpp"
#include "ibex_ThreadOptimizer.h"
//#include "ibex_ThreadOptimizer.cpp"
//#include "ibex_inputibexoptbox.h"
#include "ibex_BoundComm.h"
#include "ibex_InputVar.h"
#include "ibex_InputIbexoptbox.h"

/*#ifndef _IBEX_WITH_OPTIM_
#error "You need to install the IbexOpt plugin (--with-optim)."
#endif*/

extern ibex::BoundComm global_comm;
extern std::mutex mutex_global;
extern std::mutex mutex_event;

using namespace std;
using namespace ibex;

void makebeginanswer(const string& originalfilename, const string& originalfilename_no_ext, string& result, const InputIbexoptbox& input);

void makemodel(const string& originalfilename, const string& firstmodelopt, vector<InputVar>& paramvar, vector<InputVar>& objvar);

bool isAParam (const string& line);

string boundmulvar(const string& mulvar, vector<InputVar> objvar);

void modifymodelopt(const string& firstmodelopt, const string& lastmodelopt, vector<InputVar>& objvar, const int& nbopdone);

void getasubsystem(const string& model, const string& submodel, vector<InputVar>& paramvar, vector<InputVar>& objvar, const int& sub_iter=0);

void threadoptimize(const ThreadOptimizer::Obj obj, const HeuristicOptimizer::Bsctype bsc, const HeuristicOptimizer::Buffertype buffer, const InputIbexoptbox input, const string modelopt);

void subthreadoptimize(const string modelopt, const string overmodelopt, const ThreadOptimizer::Obj obj, const HeuristicOptimizer::Bsctype bsc, 
	const HeuristicOptimizer::Buffertype buffer, const InputIbexoptbox input, vector<InputVar> paramvar, vector<InputVar> objvar, const int sub_iter);

void beginprint(vector<InputVar>& objvar, const int& nbopdone, const InputIbexoptbox input, const string& modelopt);

void printandoptim(HeuristicOptimizer& o, System* sys, const string& modelopt, const InputIbexoptbox& input, const HeuristicOptimizer::Bsctype& bsc, const HeuristicOptimizer::Buffertype& buffer, const double& subtimer=1e-15);

void fillvalues(const string& result, vector<InputVar>& objvarcomparisonminmax, const chrono::duration<double>& execmin, const chrono::duration<double>& execmax, const int& minnbopdone, const int& maxnbopdone);

void testhom(bool homtestok[], bool ishom[], vector<InputVar>& objvarcomparisonminmax, const int& minnbopdone, const int& maxnbopdone, const double& hom);

void intermediate_report(const string& result, bool homtestok[], bool ishom[], vector<InputVar>& objvarcomparisonminmax, const int& nbvar);

void endanswer(const string& result, const chrono::duration<double>& exectime, vector<InputVar> objvar, bool homtestok[], bool ishom[], const double& hom);

#endif // _IBEX__IBEXOPTBOXFUNC_H__