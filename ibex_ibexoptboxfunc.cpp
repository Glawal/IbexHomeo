#include "ibex_ibexoptboxfunc.h"
#include <chrono>
#include <cmath>

using namespace std;
using namespace ibex;

void makebeginanswer(const string& originalfilename, const string& originalfilename_no_ext, string& result, const InputIbexoptbox& input)
{
	stringstream resultss;
	resultss << originalfilename_no_ext << "_answer.txt";
	result = resultss.str();
	ofstream answer;
	ifstream testanswer;
	testanswer.open(result, ios::in);

	if (testanswer.is_open()){
		stringstream ss;
		ss << "_old_answer" << originalfilename_no_ext << ".txt";
		string resultcopy=ss.str();
		ofstream cop(resultcopy, ios::binary);

		istreambuf_iterator<char> begin_source(testanswer);
	    istreambuf_iterator<char> end_source;
	    ostreambuf_iterator<char> begin_cop(cop);
	    copy(begin_source, end_source, begin_cop);
	}
	testanswer.close();

	answer.precision(15);

	ifstream modeluser;
	modeluser.open(originalfilename);
	answer.open(result, ios::out);
	answer << "The original model is : " << endl << endl;
	istreambuf_iterator<char> begin_source(modeluser);
	istreambuf_iterator<char> end_source;
	ostreambuf_iterator<char> begin_answer(answer);
	copy(begin_source, end_source, begin_answer);
	answer << endl << endl;
	answer << "We solve it with ibexoptbox with arguments : " << endl;
	answer << "Ibex Release : " << _IBEX_RELEASE_ << endl;
	answer << "rel_eps_f : " << input.rel_eps_f << endl;
	answer << "abs_eps_f : " << input.abs_eps_f << endl;
	answer << "eps_h : " << input.eps_h << endl;
	answer << "eps_x : " << input.eps_x << endl;
	answer << "timeout : " << input.timeout_v << endl;
	answer << "rigor : " << (input.rigor? "ON" : "OFF") << endl;
	answer << endl << endl;
	answer.close();
	modeluser.close();
}

void makemodel(const string& originalfilename, const string& firstmodelopt, vector<InputVar>& paramvar, vector<InputVar>& objvar)
{
	vector<string> multiplevar, storageline, paramtovar;
	string name;
	double min,max;
	ifstream originalmodel;
	ofstream modelopt;
	originalmodel.open(originalfilename);
	modelopt.open(firstmodelopt);
	bool keyWordConstants=false, keyWordVariables=false, keyWordConstraints=false;
	string line="";
	smatch sm,sm2;
	regex Constants("Constants|constants"), Variables("Variables|variables"), Constraints("Constraints|constraints"), searchminmax("(\\[)(.+)(,)(.+)(\\])"), 
		searchnamevar("(x\\d{1,3})"), searchmulvar("((x\\d{1,3})(\\^\\d{1,2})*([[:space:]]*\\*[[:space:]]))+(x\\d{1,3})(\\^\\d{1,2})*"), searchforparam("([[:alnum:]_-]+)(\\s+in\\s+)(\\[)(.+)(,)(.+)(\\])");
	while(getline(originalmodel,line)){
		if(!keyWordConstraints){
			if (!keyWordVariables){
				if (!keyWordConstants){
					keyWordConstants=regex_search(line,Constants);
					if (keyWordConstants)
						modelopt << line << endl;
					else{
						keyWordVariables=regex_search(line,Variables);//au cas où pas de constante.
						modelopt << line << endl;//on écrit de toute façon ce qu'il y a (il s'agit peut-être d'un commentaire).
					}
				}
				//sinon on a déjà trouvé (et écrit) Constants.
				else{
					keyWordVariables=regex_search(line,Variables);
					if (!keyWordVariables){
						bool param=isAParam(line);//on cherche si ce paramètre varie, auquel cas on le traitera comme une variable.
						if (!param){
							paramtovar.push_back(line);//mise en cache pour l'écrire plus tard (quand on trouvera Variables)
							if(regex_search(line,sm,searchforparam)){
								name=sm[1];
								min=stod(sm[4]);
								max=stod(sm[6]);
								Interval minmax(min,max);
								paramvar.push_back(InputVar(name,minmax));
							}
						}
						else
							modelopt << line << endl;
					}
					else{
						modelopt << line << endl;
					}
				}
			}
			else{//on a touvé variable mais pas Constraints
				keyWordConstraints=regex_search(line,Constraints);
				if(!keyWordConstraints){
					modelopt << line << endl;
					if(regex_search(line,sm,searchnamevar)){
						name=sm[1];
						regex_search(line,sm,searchminmax);
						min=stod(sm[2]);
						max=stod(sm[4]);
						Interval minmax(min,max);
						objvar.push_back(InputVar(name, minmax));
					}
				}
				else{//on vient de trouver Constraints
					storageline.push_back(line);
				}
			}
		}
		else{
			string tempstorage=line;
			while(regex_search(line,sm,searchmulvar)){
				string matched=sm[0];
				bool ispresent=false;
				int position=multiplevar.size()+1;
				for(int i=0;i<multiplevar.size();i++){
					if (matched==multiplevar[i]){
						ispresent=true;
						position=i+1;
					}
				}
				if(!ispresent){
					multiplevar.push_back(matched);
				}
				string newvar="v"+to_string(position);
				line=regex_replace(line,searchmulvar,newvar,regex_constants::format_first_only);
			}
			storageline.push_back(line);
		}
	}
	for (int i=0;i<paramtovar.size();i++){
		modelopt << paramtovar[i] << endl;
	}
	for (int i=0;i<multiplevar.size();i++)
	{
		string newvarline="v"+to_string(i+1)+boundmulvar(multiplevar[i], objvar);
		modelopt << newvarline << endl;
	}
	modelopt << storageline[0] << endl;//Constraints
	for(int i=0;i<multiplevar.size();i++)
	{
		string newvarconstraint="v"+to_string(i+1)+"="+multiplevar[i]+";";
		modelopt << newvarconstraint << endl;
	}
	for(int i=1;i<storageline.size();i++){
		modelopt << storageline[i] << endl;
	}
	modelopt.close();
	originalmodel.close();
}

bool isAParam (const string& line)
{
	bool isParam=true;
	regex searchForParamOrVar("(\\[)(.+)(,)(.+)(\\])");
	smatch sm;
	bool havefind=regex_search(line,sm,searchForParamOrVar);
	if (havefind){
		if (sm[2]!=sm[4])
			isParam=false;
	}
	return isParam;
}

string boundmulvar(const string& mulvar, vector<InputVar> objvar){//be careful, here range are supposed positives.
	regex onevar("(x\\d{1,3})(\\^\\d{1,2})*"), nb("\\d{1,2}");
	smatch sm, sm2;
	string copy, res;
	copy=mulvar;
	double valuemin=1, valuemax=1;
	bool namefound=false;
	while(regex_search(copy, sm, onevar)){
		if(sm[2]==""){
			int i=0;
			while(!namefound){
				if(objvar[i].name==sm[1]){
					namefound=true;
				}
				else{
					i++;
				}
			}
			valuemin=valuemin*objvar[i].range.lb();
			valuemax=valuemax*objvar[i].range.ub();
		}
		else{
			string st=sm[2];
			regex_search(st,sm2,nb);
			int i=0;
			while(!namefound){
				if(objvar[i].name==sm[1]){
					namefound=true;
				}
				else{
					i++;
				}
			}
			valuemin=valuemin*pow(objvar[i].range.lb(),stoi(sm2[0]));
			valuemax=valuemax*pow(objvar[i].range.ub(),stoi(sm2[0]));
		}
		string emptystring="";
		copy=regex_replace(copy,onevar,emptystring,regex_constants::format_first_only);
	}
	stringstream ss;
	ss.precision(15);
	ss << " in [" << valuemin << "," << valuemax << "];";
	res=ss.str();
	return res;
}

void modifymodelopt(const string& firstmodelopt, const string& lastmodelopt, vector<InputVar>& objvar, const int& nbopdone)
{
	ifstream classicmodel;
	ofstream modelopt;
	classicmodel.open(firstmodelopt);
	modelopt.open(lastmodelopt);
	bool keyWordVariables=false, keyWordConstraints=false;
	regex searchgroups("(\\[)(.+)(,)(.+)(\\])"), Variables("Variables|variables"), Constraints("Constraints|constraints"), searchnamevar("(x\\d{1,3})");
	string line="";
	int countlinevar=0;
	while(getline(classicmodel,line)){
		if(!keyWordConstraints){
			if(!keyWordVariables){
				keyWordVariables=regex_search(line, Variables);
			}
			else{
				keyWordConstraints=regex_search(line,Constraints);
				if(keyWordConstraints){
					modelopt << "Minimize" << endl;
					if(nbopdone%2==0){
						modelopt << objvar[nbopdone/2].name << endl;
					}
					else{
						modelopt << "-" << objvar[nbopdone/2].name << endl;
					}
				}
				else{
					if(countlinevar<objvar.size()){
						if(regex_search(line,searchnamevar)){
							stringstream minv,maxv;
							minv.precision(15);
							maxv.precision(15);
							minv << objvar[countlinevar].range.lb();
							maxv << objvar[countlinevar].range.ub();
							line=regex_replace(line,searchgroups,"$1 "+minv.str()+"$3 "+maxv.str()+"$5");
							countlinevar++;
						}
					}
				}
			}
		}
		modelopt << line << endl;
	}
	classicmodel.close();
	modelopt.close();
}

void getasubsystem(const string& model, const string& submodel, vector<InputVar>& paramvar, vector<InputVar>& objvar, const int& sub_iter)
{
	ifstream modelopt(model);
	ofstream submodelopt(submodel);
	vector<string> storageline;
	bool keyWordVariables=false;
	regex Variables("Variables|variables");
	string line="",towrite="";
	int countlinevar=0,countlineparam=0;
	while(getline(modelopt,line)){
		if(!keyWordVariables){
			keyWordVariables=regex_search(line,Variables);
			if(keyWordVariables){
				storageline.push_back(line);
			}
			else{
				submodelopt << line << endl;
			}
		}
		else{
			if(countlinevar<objvar.size()){
				storageline.push_back(line);
				countlinevar++;
			}
			else {
				if(countlineparam<paramvar.size()){
					stringstream ss;
					if(sub_iter==0){
						ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.lb() << ";";
					}
					else {
						if(sub_iter==1){
							ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.ub() << ";";
						}
						else{
							if(sub_iter<paramvar.size()){
								if(countlineparam<sub_iter-1){
									ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.lb() << ";";
								}
								else{
									ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.ub() << ";";
								}
							}
							else{
								if(sub_iter<2*paramvar.size()){
									if(countlineparam<(sub_iter%paramvar.size())){
										ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.lb() << ";";
									}
									else{
										ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.mid() << ";";
									}
								}
								else{
									if(sub_iter<3*paramvar.size()){
										if(countlineparam<(sub_iter%paramvar.size())){
											ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.ub() << ";";
										}
										else{
											ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.mid() << ";";
										}
									}
									else{
										ss << paramvar[countlineparam].name << " = " << paramvar[countlineparam].range.lb() << ";";//à modifier avec un aléatoire reproducible;
									}
								}
							}
						}
					}
					towrite=ss.str();
					submodelopt << towrite << endl;
					countlineparam++;
					if(countlineparam==paramvar.size())
					{
						for(int i=0; i<storageline.size();i++){
							submodelopt << storageline[i] << endl;
						}
					}
				}
				else{
					submodelopt << line << endl;
				}
			}
		}
	}
	modelopt.close();
	submodelopt.close();
}

void threadoptimize(const ThreadOptimizer::Obj obj, const HeuristicOptimizer::Bsctype bsc, const HeuristicOptimizer::Buffertype buffer, const InputIbexoptbox input, const string modelopt)
{
	switch(obj){
		case ThreadOptimizer::Obj::OBJ_MIN:
		{
			System *sys;
			sys = new System(modelopt.c_str());
			if (!sys->goal)
				ibex_error(" input file has not goal (it is not an optimization problem).");
			bool inHC4=true;
			if (sys->nb_ctr>0 && sys->nb_ctr<sys->f_ctrs.image_dim())
				inHC4=false;
			HeuristicOptimizer o(*sys, input.rel_eps_f, input.abs_eps_f, input.eps_h, input.rigor, inHC4, input.random_seed, input.eps_x, bsc, buffer, obj);
			printandoptim(o, sys, modelopt, input, bsc, buffer);
			{
				unique_lock<mutex> lck(mutex_event);
				global_comm.min_finished_first=true;
				global_comm.stop_timeout_formin.store(1e-15);
				global_comm.check_who_finished_first.notify_all();
			}
			delete sys;
			break;
		}
		case ThreadOptimizer::Obj::OBJ_MAX:
		{
			System *sys;
			sys = new System(modelopt.c_str());
			if (!sys->goal)
				ibex_error(" input file has not goal (it is not an optimization problem).");
			bool inHC4=true;
			if (sys->nb_ctr>0 && sys->nb_ctr<sys->f_ctrs.image_dim())
				inHC4=false;
			HeuristicOptimizer o(*sys, input.rel_eps_f, input.abs_eps_f, input.eps_h, input.rigor, inHC4, input.random_seed, input.eps_x, bsc, buffer, obj);
			printandoptim(o, sys, modelopt, input, bsc, buffer);
			{
				unique_lock<mutex> lck(mutex_event);
				global_comm.min_finished_first=false;
				global_comm.stop_timeout_formax.store(1e-15);
				global_comm.check_who_finished_first.notify_all();
			}
			delete sys;
			break;
		}
		default:
			cerr << "not an implemented type of objective and system." << endl;
			exit(1);
	}
}

void subthreadoptimize(const string modelopt, const string overmodelopt, const ThreadOptimizer::Obj obj, const HeuristicOptimizer::Bsctype bsc, const HeuristicOptimizer::Buffertype buffer, 
	const InputIbexoptbox input, vector<InputVar> paramvar, vector<InputVar> objvar, const int sub_iter)
{
	std::chrono::high_resolution_clock::time_point n,r;
	switch(obj){
		case ThreadOptimizer::Obj::SUB_OBJ_MIN:
		{
			int iter=sub_iter;
			n = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> execsubmin(0);
			while(iter<3*paramvar.size() && ((global_comm.stop_timeout_formin.load()<=0) || (global_comm.stop_timeout_formin.load()>1)) && (execsubmin.count()<=std::max(0.0,
				global_comm.stop_timeout_formin.load()))){//attention si le glob est négatif, ça stop instant.
				getasubsystem(overmodelopt, modelopt, paramvar, objvar, iter);
				global_comm.uploforsubmin_global.store(NEG_INFINITY);
				global_comm.loupforsubmin_global.store(POS_INFINITY);
				System *sys;
				sys = new System(modelopt.c_str());
				if (!sys->goal)
					ibex_error(" input file has not goal (it is not an optimization problem).");
				bool inHC4=true;
				if (sys->nb_ctr>0 && sys->nb_ctr<sys->f_ctrs.image_dim())
					inHC4=false;
				HeuristicOptimizer o(*sys, input.rel_eps_f, input.abs_eps_f, input.eps_h, input.rigor, inHC4, input.random_seed, input.eps_x, bsc, buffer, obj);
				double t = global_comm.stop_timeout_formin.load()-execsubmin.count();
				printandoptim(o, sys, modelopt, input, bsc, buffer, t);
				delete sys;
				iter++;
				r = std::chrono::high_resolution_clock::now();
				execsubmin = std::chrono::duration_cast<std::chrono::duration<double>>(r-n);
				//cout << "SUB_OBJ_MIN : check = " << global_comm.stop_timeout_formin.load() << endl;
			}
			break;
		}
		case ThreadOptimizer::Obj::SUB_OBJ_MAX:
		{
			int iter=sub_iter;
			n = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> execsubmax(0);
			while(iter<3*paramvar.size() && ((global_comm.stop_timeout_formax.load()<=0) || (global_comm.stop_timeout_formax.load()>1)) && (execsubmax.count()<=std::max(0.0,
				global_comm.stop_timeout_formax.load()))){//attention si le glob est négatif, ça stop instant.
				getasubsystem(overmodelopt, modelopt, paramvar, objvar, iter);
				global_comm.uploforsubmax_global.store(NEG_INFINITY);
				global_comm.loupforsubmax_global.store(POS_INFINITY);
				System *sys;
				sys = new System(modelopt.c_str());
				if (!sys->goal)
					ibex_error(" input file has not goal (it is not an optimization problem).");
				bool inHC4=true;
				if (sys->nb_ctr>0 && sys->nb_ctr<sys->f_ctrs.image_dim())
					inHC4=false;
				HeuristicOptimizer o(*sys, input.rel_eps_f, input.abs_eps_f, input.eps_h, input.rigor, inHC4, input.random_seed, input.eps_x, bsc, buffer, obj);
				double t = global_comm.stop_timeout_formax.load()-execsubmax.count();
				printandoptim(o, sys, modelopt, input, bsc, buffer, t);
				delete sys;
				iter++;
				r = std::chrono::high_resolution_clock::now();
				execsubmax = std::chrono::duration_cast<std::chrono::duration<double>>(r-n);
				//cout << "SUB_OBJ_MAX : check = " << global_comm.stop_timeout_formax.load() << endl;
			}
			break;
		}
		default:
			cerr << "not an implemented type of objective and system." << endl;
			exit(1);
	}
}

void beginprint(vector<InputVar>& objvar, const int& nbopdone, const InputIbexoptbox input, const string& modelopt)
{
	if(!input.quiet){
		cout << endl << "************************ setup ************************" << endl;
		cout << "  file loaded:\t\t" << modelopt << endl;
		cout << "  homeostasis:\t\t" << input.hom << endl;
		cout << "  rel-eps-f:\t\t" << input.rel_eps_f << "\t(relative precision on objective)" << endl;
		cout << "  abs-eps-f:\t\t" << input.abs_eps_f << "\t(absolute precision on objective)" << endl;
		cout << "  eps-h:\t\t" << input.eps_h << "\t(equality thickening)" << endl;
		cout << "  eps-x:\t\t" << input.eps_x << "\t(precision on variables domain)" << endl;
		cout << "  random_seed:\t\t" << input.random_seed << endl;
		if(input.rigor)
			cout << "  rigor mode:\t\tON\t(feasibility of equalities certified)" << endl;
		if(input.lastloup!=POS_INFINITY)
			cout << "  initial loup:\t\t" << input.lastloup << " (a priori upper bound of the minimum)" << endl;
		cout << "  timeout:\t\t" << global_comm.stop_timeout_formin.load() << "s" << endl;
		if(input.trace){
			cout << "  trace:\t\tON" << endl;
		}
		if (nbopdone%2==0){
			cout << "  Objective : Minimize " << objvar[nbopdone/2].name << endl;
		}
		else{
			cout << "  Objective : Minimize -" << objvar[nbopdone/2].name << endl;
		}
		cout << "*******************************************************" << endl << endl;
		cout << "running............" << endl << endl;
	}
}

void printandoptim(HeuristicOptimizer& o, System* sys, const string& modelopt, const InputIbexoptbox& input, const HeuristicOptimizer::Bsctype& bsc, const HeuristicOptimizer::Buffertype& buffer, const double& subtimer)
{
	if(o.get_obj()==ThreadOptimizer::Obj::OBJ_MIN){
		o.timeout=global_comm.stop_timeout_formin.load();
	}
	if(o.get_obj()==ThreadOptimizer::Obj::SUB_OBJ_MIN){
		o.timeout=subtimer;
	}
	if (o.get_obj()==ThreadOptimizer::Obj::OBJ_MAX){
		o.timeout=global_comm.stop_timeout_formax.load();
	}
	if (o.get_obj()==ThreadOptimizer::Obj::SUB_OBJ_MAX){
		o.timeout=subtimer;
	}
	if(!input.quiet){
		if(o.get_obj()==ThreadOptimizer::Obj::OBJ_MIN)
			cout << "  MIN : Launch a thread with bisector: " << o.bsctype << " and buffer: " << o.buffertype << endl;
		if(o.get_obj()==ThreadOptimizer::Obj::OBJ_MAX)
			cout << "  MAX : Launch a thread with bisector: " << o.bsctype << " and buffer: " << o.buffertype << endl;
		if(o.get_obj()==ThreadOptimizer::Obj::SUB_OBJ_MIN)
			cout << "  SUB_MIN : Launch a thread." << endl;
		if(o.get_obj()==ThreadOptimizer::Obj::SUB_OBJ_MAX)
			cout << "  SUB_MAX : Launch a thread." << endl;
		if(input.trace){
			o.trace=1;
		}
		o.optimize(sys->box,input.lastloup);
		if (input.trace) cout << endl;
		o.report();
	}
	else{
		if(input.trace) o.trace=1;
		o.optimize(sys->box, input.lastloup);
		if(input.trace) cout << endl;
	}
}

void fillvalues(const string& result, vector<InputVar>& objvarcomparisonminmax, const chrono::duration<double>& execmin, const chrono::duration<double>& execmax, const int& minnbopdone, const int& maxnbopdone)
{
	ofstream answer;
	answer.precision(15);
	answer.open(result, ios::app);
	switch(global_comm.statusformin.load()){
		case ThreadOptimizer::SUCCESS:
			answer << "The min of " << objvarcomparisonminmax[minnbopdone].name << " is " << objvarcomparisonminmax[minnbopdone].range.lb() << " (result obtained in : " << execmin.count() << "s)." << endl;
			break;
		case ThreadOptimizer::TIME_OUT:
			answer << "TIME_OUT (" << execmin.count() << "s) : the min of " <<objvarcomparisonminmax[minnbopdone].name << " is between " 
			<< objvarcomparisonminmax[minnbopdone].range.lb() << " and " << objvarcomparisonminmax[minnbopdone].range.ub() << endl;
			break;
		case ThreadOptimizer::INFEASIBLE:
			answer << "Infeasible problem. (Minimize " << objvarcomparisonminmax[minnbopdone].name << ")." << endl;
			break;
		case ThreadOptimizer::NO_FEASIBLE_FOUND:
			answer << "No feasible point found (the problem (Minimize " << objvarcomparisonminmax[minnbopdone].name << ") may be infeasible)." << endl;
			break;
		case ThreadOptimizer::UNBOUNDED_OBJ:
			answer << "Possibly unbounded objective (Minimize " << objvarcomparisonminmax[minnbopdone].name << "--> -oo)" << endl;
			break;
		case ThreadOptimizer::UNREACHED_PREC:
			answer << "Unreached precision (Minimize " << objvarcomparisonminmax[minnbopdone].name << ")." << endl;
			break;

	}
	switch(global_comm.statusformax.load()){
		case ThreadOptimizer::SUCCESS:
			answer << "The max of " << objvarcomparisonminmax[maxnbopdone].name << " is " << -objvarcomparisonminmax[maxnbopdone].range.lb() << " (result obtained in : " << execmax.count() << "s)." << endl;
			break;
		case ThreadOptimizer::TIME_OUT:
			answer << "TIME_OUT (" << execmax.count() << "s) : the max of " <<objvarcomparisonminmax[maxnbopdone].name << " is between " 
			<< -objvarcomparisonminmax[maxnbopdone].range.ub() << " and " << -objvarcomparisonminmax[maxnbopdone].range.lb() << endl;
			break;
		case ThreadOptimizer::INFEASIBLE:
			answer << "Infeasible problem. (Maximize " << objvarcomparisonminmax[maxnbopdone].name << ")." << endl;
			break;
		case ThreadOptimizer::NO_FEASIBLE_FOUND:
			answer << "No feasible point found (the problem (Minimize " << objvarcomparisonminmax[maxnbopdone].name << ") may be infeasible)." << endl;
			break;
		case ThreadOptimizer::UNBOUNDED_OBJ:
			answer << "Possibly unbounded objective (Minimize " << objvarcomparisonminmax[maxnbopdone].name << "--> -oo)" << endl;
			break;
		case ThreadOptimizer::UNREACHED_PREC:
			answer << "Unreached precision (Minimize " << objvarcomparisonminmax[maxnbopdone].name << ")." << endl;
			break;
	}
	answer.close();
}

void testhom(bool homtestok[], bool ishom[], vector<InputVar>& objvarcomparisonminmax, const int& minnbopdone, const int& maxnbopdone, const double& hom)
{
	double minmin=objvarcomparisonminmax[minnbopdone].range.lb();
	double minmax=objvarcomparisonminmax[minnbopdone].range.ub();
	double maxmin=-objvarcomparisonminmax[maxnbopdone].range.ub();
	double maxmax=-objvarcomparisonminmax[maxnbopdone].range.lb();
	if((minmin!=0 && (maxmax/minmin <= hom)) || (maxmax==0 && minmin==0)){
		homtestok[minnbopdone/2]=true;
		ishom[minnbopdone/2]=true;
	}
	else{ 
		if ((minmax!=0 && (maxmin/minmax > hom)) || (minmax==0 && maxmin!=0)){
		homtestok[minnbopdone/2]=true;
		ishom[minnbopdone/2]=false;
		}
		else{
		homtestok[minnbopdone/2]=false;
		ishom[minnbopdone/2]=false;
		}
	}
}

void intermediate_report(const string& result, bool homtestok[], bool ishom[], vector<InputVar>& objvarcomparisonminmax, const int& nbvar)
{
	ofstream answer;
	answer.precision(15);
	answer.open(result, ios::app);
	answer << endl << endl << "INTERMEDIATE REPORT : " << endl << endl;
	vector<string> varhom, varnothom;
	for(int i=0;i<nbvar;i++){
		if(homtestok[i])
		{
			if(ishom[i]){
				varhom.push_back(objvarcomparisonminmax[2*i].name);
			}
			else{
				varnothom.push_back(objvarcomparisonminmax[2*i].name);
			}
		}
	}
	int j=0, k=0;
	if(varhom.size()!=0){
		answer << "The variables ";
		while(j<varhom.size()-1){
			answer << varhom[j] << ", ";
			j++;
		}
		answer << varhom[varhom.size()-1] << " are homeostatics." << endl << endl;
	}
	if(varnothom.size()!=0){
		answer << "The variables ";
		while(k<varnothom.size()-1){
			answer << varnothom[k] << ", ";
			k++;
		}
		answer << varnothom[varnothom.size()-1] << " are not homeostatics." << endl << endl;
	}
	if(varhom.empty() && varnothom.empty())
	{
		answer << "Nothing already known." << endl << endl;
	}
	answer << endl;
	answer.close();
}

void endanswer(const string& result, const chrono::duration<double>& exectime, vector<InputVar> objvar, bool homtestok[], bool ishom[], const double& hom)
{
	ofstream answer;
	answer.precision(15);
	answer.open(result, ios::app);
	answer << endl << endl << "This is the list of loup_point found : " << endl << endl;
	{
		lock_guard<mutex> lock(mutex_global);
		for(vector<IntervalVector>::iterator it=global_comm.loup_point_global.begin(); it!=global_comm.loup_point_global.end();++it){
			answer << *it << endl << endl;
		}
	}
	answer << endl << endl << "The results are : " << endl;		
	for(int i=0;i<objvar.size();i++)
	{
		if (homtestok[i]){
			if(ishom[i]){
				answer << objvar[i].name << " in [ " << objvar[i].range.lb() << " , " << objvar[i].range.ub() << " ], this variable is " << hom << "-homeostatic." << endl;
			}
			else
			{
				answer << objvar[i].name << " in [ " << objvar[i].range.lb() << " , " << objvar[i].range.ub() << " ]." << endl;
			}
		}
		else{
			answer << objvar[i].name << " in [ " << objvar[i].range.lb() << " , " << objvar[i].range.ub() << " ], Possible " << hom << "-homeostatis." << endl;
		}
	}
	
	answer << endl << "total time : " << exectime.count() << "s." << endl;	
	answer.close();
	cout << "total time : " << exectime.count() << "s." << endl;
}
