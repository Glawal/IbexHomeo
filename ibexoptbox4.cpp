//============================================================================
//                                  I B E X
//
//                               ************
//                                  IbexOptBox
//                               ************
//
// Author      : Aurélien Desoeuvres
// Copyright   : Université de Montpellier (France)
// License     : See the LICENSE file
// Last Update : Feb 14, 2020
//============================================================================

#include "args.hxx"
//#include "ibex_ibexoptboxfunc.h"
#include "ibex_ibexoptboxfunc.cpp"
//#include "ibex_HeuristicOptimizer.h"
#include "ibex_HeuristicOptimizer.cpp"
//#include "ibex_ThreadOptimizer.h"
#include "ibex_ThreadOptimizer.cpp"
//#include "ibex_inputibexoptbox.h"
#include "ibex_BoundComm.h"
#include "ibex_InputVar.h"
#include "ibex_InputIbexoptbox.h"
//#include "ibex_inputibexoptbox.cpp"
#include <thread>
#include <chrono>
#include <future>
//#include <pthread.h>
#include <assert.h>
#include <atomic>
#include <utility>
#include <tuple>
#include <cmath>
#include "ibex.h"

extern ibex::BoundComm global_comm; 
extern std::mutex mutex_global;
extern std::mutex mutex_event;

using namespace std;
using namespace ibex;
 
int main(int argc, char** argv) 
{
	stringstream _rel_eps_f, _abs_eps_f, _eps_h, _random_seed, _eps_x, _homeostasis;
	_rel_eps_f << "Relative precision on the objective. Default value is 1e-02.";
	_abs_eps_f << "Absolute precision on the objective. Default value is 1e-14.";
	_eps_h << "Equality relaxation value. Default value is 1e-08.";
	_random_seed << "Random seed (useful for reproducibility). Default value is " << 1.0 << ".";
	_eps_x << "Precision on the variable (**Deprecated**). Default value is 0.";
	_homeostasis << "Precision for homeostasis. Default value is 2.";

	args::ArgumentParser parser("********* IbexOpt (defaultoptimizer) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::Flag version(parser, "version", "Display the version of this plugin (same as the version of Ibex).", {'v',"version"});
	args::ValueFlag<double> homeostasis(parser, "float", _homeostasis.str(), {"homeostasis"});
	args::ValueFlag<double> rel_eps_f(parser, "float", _rel_eps_f.str(), {'r', "rel-eps-f"});
	args::ValueFlag<double> abs_eps_f(parser, "float", _abs_eps_f.str(), {'a', "abs-eps-f"});
	args::ValueFlag<double> eps_h(parser, "float", _eps_h.str(), {"eps-h"});
	args::ValueFlag<double> timeoutloop(parser, "float", "a timeout for each search of the first loop. Default value is 30s.", {'l', "timeoutloop"});
	args::ValueFlag<double> timeout(parser, "float", "global timeout (time in seconds). Default value is infinity.", {'t', "timeout"});//à faire.
	args::ValueFlag<double> random_seed(parser, "float", _random_seed.str(), {"random-seed"});
	args::ValueFlag<double> eps_x(parser, "float", _eps_x.str(), {"eps-x"});
	args::ValueFlag<double> initial_loup(parser, "float", "Intial \"loup\" (a priori known upper bound).", {"initial-loup"});
	args::Flag rigor(parser, "rigor", "Activate rigor mode (certify feasibility of equalities).", {"rigor"});
	args::Flag trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});
	args::Flag withlog(parser, "withlog", "make a log file", {'l',"log"});//pour plus tard.

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	clock_t cput1,cput2;
	cput1=clock();

	cout.precision(15);

	try{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&){
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e){
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e){
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	if (version) {
		cout << "IbexOpt Release " << _IBEX_RELEASE_ << endl;
		exit(0);
	}
	if (filename.Get()=="") {
		ibex_error("no input file (try ibexopt --help)");
		exit(1);
	}

	double hom, rel_eps_f_value, abs_eps_f_value, eps_h_value, eps_x_value, timeoutloop_value, timeout_value, random_seed_value, lastloup;
	bool timeout_bool, rigor_bool, trace_bool, quiet_bool;

	hom = (homeostasis ? homeostasis.Get() : 2);
	rel_eps_f_value = (rel_eps_f ? rel_eps_f.Get() : 0.01);
	abs_eps_f_value = (abs_eps_f ? abs_eps_f.Get() : 1e-20);
	eps_h_value = (eps_h ? eps_h.Get() : 0.00000001);
	eps_x_value = (eps_x ? eps_x.Get() : 0.0);
	timeoutloop_value = (timeoutloop ? timeoutloop.Get() : 30);
	timeout_bool = (timeout ? true : false);
	timeout_value = (timeout ? timeout.Get() : -1);
	random_seed_value = (random_seed ? random_seed.Get() : 1.0);
	lastloup = (initial_loup ? initial_loup.Get() : POS_INFINITY);
	rigor_bool = (rigor ? true : false);
	trace_bool = (trace ? true : false);
	quiet_bool = (quiet ? true : false);

	InputIbexoptbox input(hom, rel_eps_f_value, abs_eps_f_value, eps_h_value, eps_x_value, timeoutloop_value, timeout_bool, timeout_value, random_seed_value, lastloup, 
		rigor_bool, trace_bool, quiet_bool);

	try {
		System *originalsys;
		string extension = filename.Get().substr(filename.Get().find_last_of('.')+1);
		if (extension == "nl") {
#ifdef _IBEX_WITH_AMPL_
			AmplInterface ampl(filename.Get());
			originalsys = new System(ampl);
#else
			cerr << "\n\033[31mCannot read \".nl\" files: AMPL plugin required \033[0m(try reconfigure with --with-ampl)\n\n";
			exit(0);
#endif
		}
		else
			originalsys = new System(filename.Get().c_str());
			// Load a system of equations

		int const nbvar = originalsys->nb_var;
		//number of variables of this system (variables to be tested for homeostasis).
		delete originalsys;
		int nbopdone = 0;
		//number of optimization done for these variables (go to 2*nbvar).
		vector<InputVar> paramtovar;
		// The constants (parameters) which are varying are stored here.
		string originalfilename=filename.Get();
		string::size_type const ext(filename.Get().find_last_of('.'));
		string originalfilename_no_ext=filename.Get().substr(0, ext);
		// filename without extension
		vector<int> varstoredbytimeout, varstoredbytimeout2;
		//variables which are not solved before timeoutloop.
		string result, firstmodelopt, nextmodeloptmin, nextmodeloptmax, nextmodeloptsubmin, nextmodeloptsubmax;
		stringstream firstfilename, nextfilename_min, nextfilename_max, nextfilename_submin, nextfilename_submax;
		firstfilename << originalfilename_no_ext << "_first_input.txt";
		nextfilename_min << originalfilename_no_ext << "_last_input_min.txt";
		nextfilename_max << originalfilename_no_ext << "_last_input_max.txt";
		nextfilename_submin << originalfilename_no_ext << "_last_input_submin.txt";
		nextfilename_submax << originalfilename_no_ext << "_last_input_submax.txt";
		firstmodelopt=firstfilename.str();
		nextmodeloptmin=nextfilename_min.str();
		nextmodeloptmax=nextfilename_max.str();
		nextmodeloptsubmin=nextfilename_submin.str();
		nextmodeloptsubmax=nextfilename_submax.str();
		vector<InputVar> objvar;
		//stores variables to be optimized and there bounds.
		vector<InputVar> objvarcomparisonminmax;
		//stores variable already optimized and there loup and uplo (used for homeostasis test);
		vector<tuple<bool,InputVar,int>> tempvarminfornextloop, tempvarmaxfornextloop, varminfornextloop, varmaxfornextloop;
		bool homtestok[nbvar], ishom[nbvar];
		for(int i=0;i<nbvar;i++){
			homtestok[i]=false;
			ishom[i]=false;
		}

		cout << global_comm << endl;

		makebeginanswer(originalfilename, originalfilename_no_ext, result, input);
		makemodel(originalfilename, firstmodelopt, paramtovar, objvar);

		std::vector<thread> threadtabmin, threadtabmax, threadtabsubmin, threadtabsubmax;

		do{
			int minnbopdone=nbopdone, maxnbopdone=nbopdone+1;
			global_comm.stop_timeout_formin.store(input.timeoutloop);
			global_comm.stop_timeout_formax.store(input.timeoutloop);
			global_comm.statusformin.store(ThreadOptimizer::Status::TIME_OUT);
			global_comm.statusformax.store(ThreadOptimizer::Status::TIME_OUT);
			global_comm.uploformin_global.store(NEG_INFINITY);
			global_comm.uploformax_global.store(NEG_INFINITY);
			global_comm.uploforsubmin_global.store(NEG_INFINITY);
			global_comm.uploforsubmax_global.store(NEG_INFINITY);
			global_comm.loupformin_global.store(POS_INFINITY);
			global_comm.loupformax_global.store(POS_INFINITY);
			global_comm.loupforsubmin_global.store(POS_INFINITY);
			global_comm.loupforsubmax_global.store(POS_INFINITY);
			{
				lock_guard<mutex> lm(mutex_global);
				global_comm.min_finished_first=false;
			}
			global_comm.nbopdone_for_min.store(minnbopdone);
			global_comm.nbopdone_for_max.store(maxnbopdone);

			std::chrono::high_resolution_clock::time_point top1,topmin,topmax;
			top1 = std::chrono::high_resolution_clock::now();
			modifymodelopt(firstmodelopt, nextmodeloptmin, objvar, minnbopdone);
			beginprint(objvar, minnbopdone, input, nextmodeloptmin);
			threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
			threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
			if(!paramtovar.empty()){
				threadtabsubmin.push_back(thread(subthreadoptimize, nextmodeloptsubmin, nextmodeloptmin, ThreadOptimizer::Obj::SUB_OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
			}
			modifymodelopt(firstmodelopt, nextmodeloptmax, objvar, maxnbopdone);
			beginprint(objvar, maxnbopdone, input, nextmodeloptmax);
			threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
			threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
			if(!paramtovar.empty()){
				threadtabsubmax.push_back(thread(subthreadoptimize, nextmodeloptsubmax, nextmodeloptmax, ThreadOptimizer::Obj::SUB_OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
			} 
			bool min_first=true;
			{
				unique_lock<mutex> lck(mutex_event);
				long int stop = floor(1000*input.timeoutloop);
				if(global_comm.check_who_finished_first.wait_for(lck, chrono::milliseconds(stop))==std::cv_status::timeout){//attention, si la condition est vérifié pour l'un, il faut quand même voir l'autre (ajouter un wait for plus loin).
					global_comm.stop_timeout_formin.store(1e-15);
					global_comm.stop_timeout_formax.store(1e-15);
				}
				if(!global_comm.min_finished_first) min_first=false;
			}
			if(min_first){
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
			}
			else{
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
			}
			while(!threadtabmin.empty()){
				threadtabmin.pop_back();
			}
			while(!threadtabsubmin.empty()){
				threadtabsubmin.pop_back();
			}
			while(!threadtabmax.empty()){
				threadtabmax.pop_back();
			}
			while(!threadtabsubmax.empty()){
				threadtabsubmax.pop_back();
			}
			Interval tempintervmin =Interval(global_comm.uploformin_global.load(),global_comm.loupformin_global.load());
			Interval tempintervmax =Interval(global_comm.uploformax_global.load(),global_comm.loupformax_global.load());
			objvar[nbopdone/2].range&=Interval(global_comm.uploformin_global.load(),-global_comm.uploformax_global.load());
			objvarcomparisonminmax.push_back(InputVar(objvar[nbopdone/2].name,tempintervmin));//objvarcomparisonminmax[minnbopdone]
			cout << "objvarcomparisonminmax[" << objvarcomparisonminmax.size()-1 << "] : " << objvarcomparisonminmax.back().name << " in [" << objvarcomparisonminmax.back().range.lb() << ", " << objvarcomparisonminmax.back().range.ub() << "]." << endl;
			objvarcomparisonminmax.push_back(InputVar(objvar[nbopdone/2].name,tempintervmax));//objvarcomparisonminmax[maxnbopdone]
			cout << "objvarcomparisonminmax[" << objvarcomparisonminmax.size()-1 << "] : " << objvarcomparisonminmax.back().name << " in [" << objvarcomparisonminmax.back().range.lb() << ", " << objvarcomparisonminmax.back().range.ub() << "]." << endl;
			std::chrono::duration<double> execmin = std::chrono::duration_cast<std::chrono::duration<double>>(topmin-top1);
			std::chrono::duration<double> execmax = std::chrono::duration_cast<std::chrono::duration<double>>(topmax-top1);
			fillvalues(result, objvarcomparisonminmax, execmin, execmax, minnbopdone, maxnbopdone);
			cout << global_comm << endl; 
			/*testhom(homtestok, ishom, objvarcomparisonminmax, minnbopdone, maxnbopdone, input.hom);
			if(!homtestok[minnbopdone/2]){
				if(global_comm.statusformin.load()==ThreadOptimizer::TIME_OUT){
					if(global_comm.statusformax.load()==ThreadOptimizer::TIME_OUT){
						tempvarminfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
					else{
						tempvarminfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(false,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
				}
				else{
					if(global_comm.statusformax.load()==ThreadOptimizer::TIME_OUT){
						tempvarminfornextloop.push_back(make_tuple(false,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
				}
			}*/
			if(nbopdone==0 && (global_comm.statusformin.load()==ThreadOptimizer::INFEASIBLE || global_comm.statusformax.load()==ThreadOptimizer::INFEASIBLE)){
				ofstream answer;
				answer.precision(15);
				answer.open(result, ios::app);
				answer << endl << endl << "This problem is infeasible, you don't have any equilibrium in this box."<< endl;
				std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> exectime = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
				cput2=clock();
				double cputime=(double)(cput2-cput1)/CLOCKS_PER_SEC;
				answer << "Total time : " << exectime.count() << "s." << endl;
				answer << "Total CPU time : " << cputime << "s." << endl;	
				answer.close();
				return 0;
			}
			nbopdone+=2;
		}while(nbopdone<2*nbvar);

		ofstream answer;
		answer.precision(15);
		answer.open(result, ios::app);
		answer << endl << endl << "Current list of loup_point found (with high equality thickness) : " << endl << endl;
		{
			lock_guard<mutex> lock(mutex_global);
			for(vector<IntervalVector>::iterator it=global_comm.loup_point_global.begin(); it!=global_comm.loup_point_global.end();++it){
				answer << *it << endl << endl;
			}
			answer << endl;
			while(!global_comm.loup_point_global.empty()){
				global_comm.loup_point_global.pop_back();
			}
		}
		answer.close();

		input.eps_h*=0.0001;
		global_comm.stop_timeout_formin.store(input.timeoutloop); 
		global_comm.stop_timeout_formax.store(input.timeoutloop);
		global_comm.statusformin.store(ThreadOptimizer::Status::TIME_OUT);
		global_comm.statusformax.store(ThreadOptimizer::Status::TIME_OUT);
		global_comm.uploformin_global.store(NEG_INFINITY);
		global_comm.uploformax_global.store(NEG_INFINITY);
		global_comm.uploforsubmin_global.store(NEG_INFINITY);
		global_comm.uploforsubmax_global.store(NEG_INFINITY);
		global_comm.loupformin_global.store(POS_INFINITY);
		global_comm.loupformax_global.store(POS_INFINITY);
		global_comm.loupforsubmin_global.store(POS_INFINITY);
		global_comm.loupforsubmax_global.store(POS_INFINITY);
		{
			lock_guard<mutex> lm(mutex_global);
			global_comm.min_finished_first=false; 
		}
		nbopdone=0;

		do{
			int minnbopdone=nbopdone, maxnbopdone=nbopdone+1;
			global_comm.stop_timeout_formin.store(input.timeoutloop);
			global_comm.stop_timeout_formax.store(input.timeoutloop);
			global_comm.statusformin.store(ThreadOptimizer::Status::TIME_OUT);
			global_comm.statusformax.store(ThreadOptimizer::Status::TIME_OUT);
			global_comm.uploformin_global.store(NEG_INFINITY);
			global_comm.uploformax_global.store(NEG_INFINITY);
			global_comm.uploforsubmin_global.store(NEG_INFINITY);
			global_comm.uploforsubmax_global.store(NEG_INFINITY);
			global_comm.loupformin_global.store(POS_INFINITY);
			global_comm.loupformax_global.store(POS_INFINITY);
			global_comm.loupforsubmin_global.store(POS_INFINITY);
			global_comm.loupforsubmax_global.store(POS_INFINITY);
			{
				lock_guard<mutex> lm(mutex_global);
				global_comm.min_finished_first=false;
			}
			global_comm.nbopdone_for_min.store(minnbopdone);
			global_comm.nbopdone_for_max.store(maxnbopdone);

			std::chrono::high_resolution_clock::time_point top1,topmin,topmax;
			top1 = std::chrono::high_resolution_clock::now();
			modifymodelopt(firstmodelopt, nextmodeloptmin, objvar, minnbopdone);
			beginprint(objvar, minnbopdone, input, nextmodeloptmin);
			threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
			threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
			if(!paramtovar.empty()){
				threadtabsubmin.push_back(thread(subthreadoptimize, nextmodeloptsubmin, nextmodeloptmin, ThreadOptimizer::Obj::SUB_OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
			}
			modifymodelopt(firstmodelopt, nextmodeloptmax, objvar, maxnbopdone);
			beginprint(objvar, maxnbopdone, input, nextmodeloptmax);
			threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
			threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
			if(!paramtovar.empty()){
				threadtabsubmax.push_back(thread(subthreadoptimize, nextmodeloptsubmax, nextmodeloptmax, ThreadOptimizer::Obj::SUB_OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
			} 
			bool min_first=true;
			{
				unique_lock<mutex> lck(mutex_event);
				long int stop = floor(1000*input.timeoutloop);
				if(global_comm.check_who_finished_first.wait_for(lck, chrono::milliseconds(stop))==std::cv_status::timeout){//attention, si la condition est vérifié pour l'un, il faut quand même voir l'autre (ajouter un wait for plus loin).
					global_comm.stop_timeout_formin.store(1e-15);
					global_comm.stop_timeout_formax.store(1e-15);
				}
				if(!global_comm.min_finished_first) min_first=false;
			}
			if(min_first){
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
			}
			else{
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
			}
			while(!threadtabmin.empty()){
				threadtabmin.pop_back();
			}
			while(!threadtabsubmin.empty()){
				threadtabsubmin.pop_back();
			}
			while(!threadtabmax.empty()){
				threadtabmax.pop_back();
			}
			while(!threadtabsubmax.empty()){
				threadtabsubmax.pop_back();
			}
			Interval tempintervmin =Interval(global_comm.uploformin_global.load(),global_comm.loupformin_global.load());
			Interval tempintervmax =Interval(global_comm.uploformax_global.load(),global_comm.loupformax_global.load());
			objvar[nbopdone/2].range&=Interval(global_comm.uploformin_global.load(),-global_comm.uploformax_global.load());
			objvarcomparisonminmax[minnbopdone].range =tempintervmin;//objvarcomparisonminmax[minnbopdone]
			cout << "objvarcomparisonminmax[" << minnbopdone << "] : " << objvarcomparisonminmax[minnbopdone].name << " in [" << objvarcomparisonminmax[minnbopdone].range.lb() << ", " << objvarcomparisonminmax[minnbopdone].range.ub() << "]." << endl;
			objvarcomparisonminmax[maxnbopdone].range =tempintervmax;//objvarcomparisonminmax[maxnbopdone]
			cout << "objvarcomparisonminmax[" << maxnbopdone << "] : " << objvarcomparisonminmax[maxnbopdone].name << " in [" << objvarcomparisonminmax[maxnbopdone].range.lb() << ", " << objvarcomparisonminmax[maxnbopdone].range.ub() << "]." << endl;
			std::chrono::duration<double> execmin = std::chrono::duration_cast<std::chrono::duration<double>>(topmin-top1);
			std::chrono::duration<double> execmax = std::chrono::duration_cast<std::chrono::duration<double>>(topmax-top1);
			fillvalues(result, objvarcomparisonminmax, execmin, execmax, minnbopdone, maxnbopdone);
			cout << global_comm << endl; 
			testhom(homtestok, ishom, objvarcomparisonminmax, minnbopdone, maxnbopdone, input.hom);
			if(!homtestok[minnbopdone/2]){
				if(global_comm.statusformin.load()==ThreadOptimizer::TIME_OUT){
					if(global_comm.statusformax.load()==ThreadOptimizer::TIME_OUT){
						tempvarminfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
					else{
						tempvarminfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(false,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
				}
				else{
					if(global_comm.statusformax.load()==ThreadOptimizer::TIME_OUT){
						tempvarminfornextloop.push_back(make_tuple(false,objvarcomparisonminmax[minnbopdone],minnbopdone));
						tempvarmaxfornextloop.push_back(make_tuple(true,objvarcomparisonminmax[maxnbopdone],maxnbopdone));
					}
				}
			}
			if(nbopdone==0 && (global_comm.statusformin.load()==ThreadOptimizer::INFEASIBLE || global_comm.statusformax.load()==ThreadOptimizer::INFEASIBLE)){
				ofstream answer;
				answer.precision(15);
				answer.open(result, ios::app);
				answer << endl << endl << "This problem is infeasible, you don't have any equilibrium in this box."<< endl;
				std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> exectime = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
				cput2=clock();
				double cputime=(double)(cput2-cput1)/CLOCKS_PER_SEC;
				answer << "Total time : " << exectime.count() << "s." << endl;
				answer << "Total CPU time : " << cputime << "s." << endl;
				answer.close();
				return 0;
			}
			nbopdone+=2;
		}while(nbopdone<2*nbvar);

		intermediate_report(result, homtestok, ishom, objvarcomparisonminmax, nbvar);

		while(!tempvarminfornextloop.empty()){
			varminfornextloop.push_back(tempvarminfornextloop.back());
			tempvarminfornextloop.pop_back();
		}
		while(!tempvarmaxfornextloop.empty()){
			varmaxfornextloop.push_back(tempvarmaxfornextloop.back());
			tempvarmaxfornextloop.pop_back();
		}

		for(int i=0; i<varminfornextloop.size();i++){ 
			cout << "varminfornextloop : " << get<0>(varminfornextloop[i]) << ", " << get<1>(varminfornextloop[i]).name << " in [" << get<1>(varminfornextloop[i]).range.lb() << ", " << get<1>(varminfornextloop[i]).range.ub() << "], minnbopdone : " << get<2>(varminfornextloop[i]) << endl;
			cout << "varmaxfornextloop : " << get<0>(varmaxfornextloop[i]) << ", " << get<1>(varmaxfornextloop[i]).name << " in [" << get<1>(varmaxfornextloop[i]).range.lb() << ", " << get<1>(varmaxfornextloop[i]).range.ub() << "], maxnbopdone : " << get<2>(varmaxfornextloop[i]) << endl;
		}

		global_comm.stop_timeout_formin.store(input.timeoutloop); 
		global_comm.stop_timeout_formax.store(input.timeoutloop);
		global_comm.statusformin.store(ThreadOptimizer::Status::TIME_OUT);
		global_comm.statusformax.store(ThreadOptimizer::Status::TIME_OUT);
		global_comm.uploformin_global.store(NEG_INFINITY);
		global_comm.uploformax_global.store(NEG_INFINITY);
		global_comm.uploforsubmin_global.store(NEG_INFINITY);
		global_comm.uploforsubmax_global.store(NEG_INFINITY);
		global_comm.loupformin_global.store(POS_INFINITY);
		global_comm.loupformax_global.store(POS_INFINITY);
		global_comm.loupforsubmin_global.store(POS_INFINITY);
		global_comm.loupforsubmax_global.store(POS_INFINITY);
		{
			lock_guard<mutex> lm(mutex_global);
			global_comm.min_finished_first=false; 
		}

		while(!varminfornextloop.empty() || !varmaxfornextloop.empty()){
			global_comm.stop_timeout_formin.store(10*input.timeoutloop);
			global_comm.stop_timeout_formax.store(10*input.timeoutloop);
			global_comm.statusformin.store(ThreadOptimizer::Status::TIME_OUT); 
			global_comm.statusformax.store(ThreadOptimizer::Status::TIME_OUT);
			{
				lock_guard<mutex> lm(mutex_global);
				global_comm.min_finished_first=false;
			}
			pair<bool,InputVar> currentobjmin, currentobjmax;
			currentobjmin=make_pair(get<0>(varminfornextloop.back()),get<1>(varminfornextloop.back()));
			currentobjmax=make_pair(get<0>(varmaxfornextloop.back()),get<1>(varmaxfornextloop.back()));
			int minnbopdone=get<2>(varminfornextloop.back()), maxnbopdone=get<2>(varmaxfornextloop.back());
			varminfornextloop.pop_back();
			varmaxfornextloop.pop_back();
			global_comm.nbopdone_for_min.store(minnbopdone);
			global_comm.nbopdone_for_max.store(maxnbopdone);
			global_comm.uploformin_global.store(get<1>(currentobjmin).range.lb());
			global_comm.uploformax_global.store(get<1>(currentobjmax).range.lb());
			global_comm.uploforsubmin_global.store(get<1>(currentobjmin).range.lb());
			global_comm.uploforsubmax_global.store(get<1>(currentobjmax).range.lb());
			global_comm.loupformin_global.store(get<1>(currentobjmin).range.ub());
			global_comm.loupformax_global.store(get<1>(currentobjmax).range.ub());
			global_comm.loupforsubmin_global.store(POS_INFINITY);
			global_comm.loupforsubmax_global.store(POS_INFINITY);
			{
				lock_guard<mutex> lock(mutex_global);
				cout << global_comm << endl;
			}
			std::chrono::high_resolution_clock::time_point top1,topmin,topmax;
			top1 = std::chrono::high_resolution_clock::now();
			if (get<0>(currentobjmin)){
				modifymodelopt(firstmodelopt, nextmodeloptmin, objvar, minnbopdone);
				beginprint(objvar, minnbopdone, input, nextmodeloptmin);
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_SUM, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_SUM, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_MAX, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_MAX, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_SUM_REL, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_SUM_REL, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_MAX_REL, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmin));
				threadtabmin.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MIN, HeuristicOptimizer::Bsctype::SMEAR_MAX_REL, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmin));
				if(!paramtovar.empty()){
					threadtabsubmin.push_back(thread(subthreadoptimize, nextmodeloptsubmin, nextmodeloptmin, ThreadOptimizer::Obj::SUB_OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
				}
			}
			if (get<0>(currentobjmax)){
				modifymodelopt(firstmodelopt, nextmodeloptmax, objvar, maxnbopdone);
				beginprint(objvar, maxnbopdone, input, nextmodeloptmax);
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_SUM, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_SUM, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_MAX, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_MAX, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_SUM_REL, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_SUM_REL, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_MAX_REL, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, nextmodeloptmax));
				threadtabmax.push_back(thread(threadoptimize, ThreadOptimizer::Obj::OBJ_MAX, HeuristicOptimizer::Bsctype::SMEAR_MAX_REL, HeuristicOptimizer::Buffertype::DOUBLE_HEAP, input, nextmodeloptmax));
				if(!paramtovar.empty()){
					threadtabsubmax.push_back(thread(subthreadoptimize, nextmodeloptsubmax, nextmodeloptmax, ThreadOptimizer::Obj::SUB_OBJ_MAX, HeuristicOptimizer::Bsctype::LSMEAR, HeuristicOptimizer::Buffertype::CELL_BEAM_SEARCH, input, paramtovar, objvar, 0));
				}
			}
			bool min_first=true;
			{
				unique_lock<mutex> lck(mutex_event);
				global_comm.check_who_finished_first.wait(lck);
				if(!global_comm.min_finished_first) min_first=false;
			}
			if(min_first){
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
			}
			else{
				for(int i=0; i<threadtabmax.size();i++){
					if(threadtabmax[i].joinable())
						threadtabmax[i].join();
				}
				topmax = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmax.size();i++){
					if(threadtabsubmax[i].joinable())
						threadtabsubmax[i].join();
				}
				for(int i=0; i<threadtabmin.size();i++){
					if(threadtabmin[i].joinable())
						threadtabmin[i].join();
				}
				topmin = std::chrono::high_resolution_clock::now();
				for(int i=0; i<threadtabsubmin.size();i++){
					if(threadtabsubmin[i].joinable())
						threadtabsubmin[i].join();
				}
			}
			while(!threadtabmin.empty()){
				threadtabmin.pop_back();
			}
			while(!threadtabsubmin.empty()){
				threadtabsubmin.pop_back();
			}
			while(!threadtabmax.empty()){
				threadtabmax.pop_back();
			}
			while(!threadtabsubmax.empty()){
				threadtabsubmax.pop_back();
			}
			objvar[minnbopdone/2].range&=Interval(global_comm.uploformin_global.load(),-global_comm.uploformax_global.load());
			objvarcomparisonminmax[minnbopdone].range=Interval(global_comm.uploformin_global.load(),global_comm.loupformin_global.load());
			objvarcomparisonminmax[maxnbopdone].range=Interval(global_comm.uploformax_global.load(),global_comm.loupformax_global.load());
			std::chrono::duration<double> execmin = std::chrono::duration_cast<std::chrono::duration<double>>(topmin-top1);
			std::chrono::duration<double> execmax = std::chrono::duration_cast<std::chrono::duration<double>>(topmax-top1);
			fillvalues(result, objvarcomparisonminmax, execmin, execmax, minnbopdone, maxnbopdone);
			testhom(homtestok, ishom, objvarcomparisonminmax, minnbopdone, maxnbopdone, input.hom);
		}
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> exectime = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		endanswer(result, exectime, objvar, homtestok, ishom, input.hom);;
		//ofstream answer;
		answer.precision(15);
		answer.open(result, ios::app);
		cput2=clock();
		double cputime=(double)(cput2-cput1)/CLOCKS_PER_SEC;
		answer << "Total CPU time : " << cputime << "s." << endl;
		answer.close();
		return 0;
	}//end try
	catch(ibex::UnknownFileException& e) {
	cerr << "Error: cannot read file '" << filename.Get() << "'" << endl;
	}
	catch(ibex::SyntaxError& e) {
		cout << e << endl;
	}
}