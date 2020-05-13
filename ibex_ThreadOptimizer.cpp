//============================================================================
//                                  I B E X
// File        : ibex_ThreadOptimizer.cpp
// Author      : Aurélien Desoeuvres
// Copyright   : Université de Montpellier (France)
// License     : See the LICENSE file
// Created     : Feb 19, 2020
// Last Update : Feb 19, 2020
//============================================================================

#include "ibex_ThreadOptimizer.h"
//#include "ibex_RealTimer.h"
#include "ibex_RealTimer.cpp"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_BxpOptimData.h"
#include "ibex_CovOptimData.h"

#include <float.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

mutex mutex_global;
mutex mutex_event;
ibex::BoundComm global_comm{};

namespace ibex {

/*
 * TODO: redundant with ExtendedSystem.
 */
void ThreadOptimizer::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		ext_box[i2]=box[i];
	}
}

void ThreadOptimizer::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}

ThreadOptimizer::ThreadOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
		int goal_var, ThreadOptimizer::Obj t_obj, double eps_x, double rel_eps_f, double abs_eps_f) :
                						n(n), goal_var(goal_var), obj(t_obj),
										ctc(ctc), bsc(bsc), loup_finder(finder), buffer(buffer),
										eps_x(eps_x), rel_eps_f(rel_eps_f), abs_eps_f(abs_eps_f),
										trace(0), timeout(-1), anticipated_upper_bounding(true),
										status(SUCCESS),
										uplo(NEG_INFINITY), uplo_of_epsboxes(1,POS_INFINITY), loup(POS_INFINITY),
										loup_point(IntervalVector::empty(n)), initial_loup(POS_INFINITY), loup_changed(false),
										time(0), nb_cells(0){
cout.precision(15);
}

ThreadOptimizer::~ThreadOptimizer() {
}

// compute the value ymax (decreasing the loup with the precision)
// the heap and the current box are contracted with y <= ymax
double ThreadOptimizer::compute_ymax() {
	if (anticipated_upper_bounding) {
		double ymax = loup>0 ?
				1/(1+rel_eps_f)*loup
		:
				1/(1-rel_eps_f)*loup;

		if (loup - abs_eps_f < ymax)
			ymax = loup - abs_eps_f;
		return next_float(ymax);
	} else
		return loup;
}

bool ThreadOptimizer::best_loup_found(){
	bool b=false;
	lock_guard<mutex> lock(mutex_global);
	if(global_comm.loup_point_global.empty()){
		return b;
	}
	else{
		for(int i=0; i<global_comm.loup_point_global.size(); i++){
			switch(obj){
				case OBJ_MIN:{
					if(global_comm.loup_point_global[i][(global_comm.nbopdone_for_min.load())/2].lb()<loup){
						loup=global_comm.loup_point_global[i][(global_comm.nbopdone_for_min.load())/2].lb();
						loup_point=global_comm.loup_point_global[i];
						b=true;
					}
					break;
				}
				case OBJ_MAX:{
					if(-global_comm.loup_point_global[i][(global_comm.nbopdone_for_max.load())/2].lb()<loup){
						loup=-global_comm.loup_point_global[i][(global_comm.nbopdone_for_max.load())/2].lb();
						loup_point=global_comm.loup_point_global[i];
						b=true;
					}
					break;
				}
				case SUB_OBJ_MIN:{
					break;
				}
				case SUB_OBJ_MAX:{
					break;
				}
			}
		}
		switch(obj){
			case OBJ_MIN:{
				double loup_glob=global_comm.loupformin_global.load();
				while(loup<loup_glob && !global_comm.loupformin_global.compare_exchange_weak(loup_glob,loup));
				if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
					cout << po->bsctype << "  " << po->buffertype  << "  ";
				}
				cout << "OBJ_MIN : best_loup_found : loup = " << loup << endl;
				break;
			}
			case OBJ_MAX:{
				double loup_glob=global_comm.loupformax_global.load();
				while(loup<loup_glob && !global_comm.loupformax_global.compare_exchange_weak(loup_glob,loup));
				if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
					cout << po->bsctype << "  " << po->buffertype  << "  ";
				}
				cout << "OBJ_MAX : best_loup_found : loup = " << loup << endl;
				break;
			}
			case SUB_OBJ_MIN:{
				break;
			}
			case SUB_OBJ_MAX:{
				break;
			}
		}
	}
	return b;
}

bool ThreadOptimizer::update_loup(const IntervalVector& box, BoxProperties& prop) {

	try {
		pair<IntervalVector,double> p=loup_finder.find(box,loup_point,loup,prop);
		loup_point = p.first;
		{
			lock_guard<mutex> lock(mutex_global);
			global_comm.loup_point_global.push_back(loup_point);
		}
		loup = p.second;
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
			{
				double loup_glob_min=global_comm.loupformin_global.load();
				double loup_glob_max=-global_comm.loupformax_global.load();
				double loup_glob;
				if (loup_glob_max!=NEG_INFINITY){
					loup_glob=std::min(loup_glob_max,loup_glob_min);
				}
				else{
					loup_glob=loup_glob_min;
				}
				double minofnew=std::min(loup,loup_glob);
				if (trace){
					if (minofnew<loup_glob_min) {
						if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
							cout << po->bsctype << "  " << po->buffertype  << "  ";
						}
						cout << "                    ";
						cout << "\033[32m loup for min = " << minofnew << "\033[0m" << endl;
					}
				}
				while(minofnew<loup_glob_min && !global_comm.loupformin_global.compare_exchange_weak(loup_glob_min,minofnew));
				loup=global_comm.loupformin_global.load();
				break;
			}
			case ThreadOptimizer::Obj::OBJ_MAX:
			{
				double loup_glob_max=global_comm.loupformax_global.load();
				double loup_glob_min=-global_comm.loupformin_global.load();
				double loup_glob;
				if (loup_glob_min!=NEG_INFINITY){
					loup_glob=std::min(loup_glob_max,loup_glob_min);
				}
				else{
					loup_glob=loup_glob_max;
				}
				double minofnew=std::min(loup,loup_glob);
				if (minofnew<loup_glob_max){
					if (trace) {
						if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
							cout << po->bsctype << "  " << po->buffertype  << "  ";
						}
						cout << "                    ";
						cout << "\033[32m loup for max = " << minofnew << "\033[0m" << endl;
					}
				}
				while(minofnew<loup_glob_max && !global_comm.loupformax_global.compare_exchange_weak(loup_glob_max,minofnew));
				loup=global_comm.loupformax_global.load();
				break;
			}
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
			{
				double loup_glob=global_comm.loupforsubmin_global.load();
				while(loup<loup_glob && !global_comm.loupforsubmin_global.compare_exchange_weak(loup_glob,loup));
				loup=global_comm.loupforsubmin_global.load();
				double big_loup_glob=global_comm.loupformin_global.load();
				if (loup<big_loup_glob){
					if (trace) {
						cout << "                    ";
						cout << "\033[32m loup for min (from sub) = " << loup << "\033[0m" << endl;
					}
				}
				while(loup<big_loup_glob && !global_comm.loupformin_global.compare_exchange_weak(big_loup_glob,loup));
				break;
			}
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
			{
				double loup_glob=global_comm.loupforsubmax_global.load();
				while(loup<loup_glob && !global_comm.loupforsubmax_global.compare_exchange_weak(loup_glob,loup));
				loup=global_comm.loupforsubmax_global.load();
				double big_loup_glob=global_comm.loupformax_global.load();
				if (loup<big_loup_glob){
					if (trace) {
						cout << "                    ";
						cout << "\033[32m loup for max (from sub) = " << loup << "\033[0m" << endl;
					}
				}
				while(loup<big_loup_glob && !global_comm.loupformax_global.compare_exchange_weak(big_loup_glob,loup));
				break;
			}
			default :
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
		return true;
	} catch(LoupFinder::NotFound&) {
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
			{
				double loup_glob_min=global_comm.loupformin_global.load();
				double loup_glob_max=-global_comm.loupformax_global.load();
				double loup_glob;
				if (loup_glob_max!=NEG_INFINITY){
					loup_glob=std::min(loup_glob_max,loup_glob_min);
				}
				else{
					loup_glob=loup_glob_min;
				}
				if(loup>loup_glob){
					loup=loup_glob;
					return true;
				}
				break;
			}
			case ThreadOptimizer::Obj::OBJ_MAX:
			{
				double loup_glob_max=global_comm.loupformax_global.load();
				double loup_glob_min=-global_comm.loupformin_global.load();
				double loup_glob;
				if (loup_glob_min!=NEG_INFINITY){
					loup_glob=std::min(loup_glob_max,loup_glob_min);
				}
				else{
					loup_glob=loup_glob_max;
				}
				if(loup>loup_glob){
					loup=loup_glob;
					return true;
				}
				break;
			}
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
				break;
			default :
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
		return false;
	}
}

void ThreadOptimizer::update_uplo() {
	double new_uplo=POS_INFINITY;

	if (! buffer.empty()) {
		double new_temp_uplo=buffer.minimum();
		bool uplofromglobal=true;
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
				{
					double uplo_glob=global_comm.uploformin_global.load();
					while(new_temp_uplo>uplo_glob && !global_comm.uploformin_global.compare_exchange_weak(uplo_glob,new_temp_uplo));
					if(new_temp_uplo>uplo_glob){uplofromglobal=false;}
					new_uplo=global_comm.uploformin_global.load();
					break;
				}
			case ThreadOptimizer::Obj::OBJ_MAX:
				{
					double uplo_glob=global_comm.uploformax_global.load();
					while(new_temp_uplo>uplo_glob && !global_comm.uploformax_global.compare_exchange_weak(uplo_glob,new_temp_uplo));
					if(new_temp_uplo>uplo_glob){uplofromglobal=false;}
					new_uplo=global_comm.uploformax_global.load();
					break;
				}
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
				{
					new_uplo=new_temp_uplo;
					break;
				}
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
				{
					new_uplo=new_temp_uplo;
					break;
				}
			default:
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
		if(uplofromglobal){
			while(new_uplo>uplo_of_epsboxes.back()){
				uplo_of_epsboxes.pop_back();
			}
		}
		if (new_uplo > loup && uplo_of_epsboxes.back() > loup) {
			switch(obj){
				case ThreadOptimizer::Obj::OBJ_MIN:
					cout << "OBJ_MIN :";
					break;
				case ThreadOptimizer::Obj::OBJ_MAX:
					cout << "OBJ_MAX :";
					break;
				case ThreadOptimizer::Obj::SUB_OBJ_MIN:
					cout << "SUB_OBJ_MIN :";
					break;
				case ThreadOptimizer::Obj::SUB_OBJ_MAX:
					cout << "SUB_OBJ_MAX :";
					break;
			}
			cout << " loup = " << loup << " new_uplo=" << new_uplo <<  " uplo_of_epsboxes=" << uplo_of_epsboxes.back() << endl;
			ibex_error("optimizer: new_uplo>loup (please report bug)");
		}
		if (new_uplo < uplo) {
			switch(obj){
				case ThreadOptimizer::Obj::OBJ_MIN:
					cout << "OBJ_MIN :";
					break;
				case ThreadOptimizer::Obj::OBJ_MAX:
					cout << "OBJ_MAX :";
					break;
				case ThreadOptimizer::Obj::SUB_OBJ_MIN:
					cout << "SUB_OBJ_MIN :";
					break;
				case ThreadOptimizer::Obj::SUB_OBJ_MAX:
					cout << "SUB_OBJ_MAX :";
					break;
			}
			cout << "uplo= " << uplo << " new_uplo=" << new_uplo << endl;
			ibex_error("optimizer: new_uplo<uplo (please report bug)");
		}
		// uplo <- max(uplo, min(new_uplo, uplo_of_epsboxes))
		if (new_uplo < uplo_of_epsboxes.back()) {
			if (new_uplo > uplo) {
				uplo = new_uplo;
				if(!uplofromglobal){
					if (trace){
						switch(obj){
							case ThreadOptimizer::Obj::OBJ_MIN:
							{
								if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
									cout << po->bsctype << "  " << po->buffertype  << "  ";
								}
								cout << "\033[33m uplo for min = " << uplo << "\033[0m" << endl;
								break;
							}
							case ThreadOptimizer::Obj::OBJ_MAX:
							{
								if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
									cout << po->bsctype << "  " << po->buffertype  << "  ";
								}
								cout << "\033[33m uplo for max = " << uplo << "\033[0m" << endl;
								break;
							}
							case ThreadOptimizer::Obj::SUB_OBJ_MIN:
								break;
							case ThreadOptimizer::Obj::SUB_OBJ_MAX:
								break;
							default:
								cerr << "not an implemented type of objective and system." << endl;
								exit(1);
						}
					}
				}
			}
		}
		else uplo = uplo_of_epsboxes.back();
	}
	else if (buffer.empty() && loup != POS_INFINITY) {
		// empty buffer : new uplo is set to ymax (loup - precision) if a loup has been found
		new_uplo=compute_ymax(); // not new_uplo=loup, because constraint y <= ymax was enforced
		double m = (new_uplo < uplo_of_epsboxes.back()) ? new_uplo :  uplo_of_epsboxes.back();
		if (uplo < m) uplo = m; // warning: hides the field "m" of the class
		// note: we always have uplo <= uplo_of_epsboxes but we may have uplo > new_uplo, because
		// ymax is strictly lower than the loup.
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:{
				double uplo_glob=global_comm.uploformin_global.load();
				while(uplo>uplo_glob && !global_comm.uploformin_global.compare_exchange_weak(uplo_glob,uplo));
				break;
			}
			case ThreadOptimizer::Obj::OBJ_MAX:{
				double uplo_glob=global_comm.uploformax_global.load();
				while(uplo>uplo_glob && !global_comm.uploformax_global.compare_exchange_weak(uplo_glob,uplo));
				break;
			}
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:{
				break;
			}
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:{
				break;
			}
		}
	}
}

void ThreadOptimizer::update_uplo_of_epsboxes(double ymin) {//améliorer l'affichage.

	// the current box cannot be bisected.  ymin is a lower bound of the objective on this box
	// uplo of epsboxes can only go down, but not under uplo : it is an upperbound for uplo,
	// that indicates a lowerbound for the objective in all the small boxes
	// found by the precision criterion
	assert (uplo_of_epsboxes.back() >= uplo);
	assert(ymin >= uplo);
	if (uplo_of_epsboxes.back() > ymin) {
		uplo_of_epsboxes.push_back(ymin);
		if (trace) {
			switch(obj){
				case ThreadOptimizer::Obj::OBJ_MIN:
				{
					if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
						cout << po->bsctype << "  " << po->buffertype  << "  ";
					}
					cout << " unprocessable tiny box (OBJ_MIN) : now uplo <= " << setprecision(15) <<  uplo_of_epsboxes.back() << " uplo = " << uplo << endl;
					break;
				}
				case ThreadOptimizer::Obj::OBJ_MAX:
				{
					if (HeuristicOptimizer* po = dynamic_cast<HeuristicOptimizer*>(this)){
						cout << po->bsctype << "  " << po->buffertype  << "  ";
					}
					cout << " unprocessable tiny box (OBJ_MAX) : now uplo <= " << setprecision(15) <<  uplo_of_epsboxes.back() << " uplo = " << uplo << endl;
					break;
				}
				case ThreadOptimizer::Obj::SUB_OBJ_MIN:
					break;
				case ThreadOptimizer::Obj::SUB_OBJ_MAX:
					break;
				default:
					cerr << "not an implemented type of objective and system." << endl;
					exit(1);
			}
		}
	}
}

void ThreadOptimizer::handle_cell(Cell& c) {

	contract_and_bound(c);

	if (c.box.is_empty()) {
		delete &c;
	} else {
		buffer.push(&c);
	}
}

void ThreadOptimizer::contract_and_bound(Cell& c) {

	/*======================== contract y with y<=loup and y>=uplo ========================*/
	Interval& y=c.box[goal_var];

	double ymax;
	if (loup==POS_INFINITY) ymax = POS_INFINITY;
	// ymax is slightly increased to favour subboxes of the loup
	// TODO: useful with double heap??
	else ymax = next_float(compute_ymax());//+1.e-15;

	y &= Interval(uplo,ymax);

	if (y.is_empty()) {
		c.box.set_empty();
		return;
	} else {
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	/*================ contract x with f(x)=y and g(x)<=0 ================*/
	//cout << " [contract]  x before=" << c.box << endl;
	//cout << " [contract]  y before=" << y << endl;

	ContractContext context(c.prop);
	if (c.bisected_var!=-1) {
		context.impact.clear();
		context.impact.add(c.bisected_var);
		context.impact.add(goal_var);
	}

	ctc.contract(c.box, context);
	//cout << c.prop << endl;
	if (c.box.is_empty()) return;

	//cout << " [contract]  x after=" << c.box << endl;
	//cout << " [contract]  y after=" << y << endl;
	/*====================================================================*/

	/*========================= update loup =============================*/

	IntervalVector tmp_box(n);
	read_ext_box(c.box,tmp_box);

	c.prop.update(BoxEvent(c.box,BoxEvent::CHANGE));

	bool loup_ch=update_loup(tmp_box, c.prop);

	// update of the upper bound of y in case of a new loup found
	if (loup_ch) {
		y &= Interval(NEG_INFINITY,compute_ymax());
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	//TODO: should we propagate constraints again?

	loup_changed |= loup_ch;

	if (y.is_empty()) { // fix issue #44
		c.box.set_empty();
		return;
	}

	/*====================================================================*/
	// Note: there are three different cases of "epsilon" box,
	// - NoBisectableVariableException raised by the bisector (---> see optimize(...)) which
	//   is independent from the optimizer
	// - the width of the box is less than the precision given to the optimizer ("prec" for the original variables
	//   and "goal_abs_prec" for the goal variable)
	// - the extended box has no bisectable domains (if prec=0 or <1 ulp)
	if ((tmp_box.max_diam()<=eps_x && y.diam() <=abs_eps_f) || !c.box.is_bisectable()) {
		update_uplo_of_epsboxes(y.lb());
		c.box.set_empty();
		return;
	}

	// ** important: ** must be done after upper-bounding
	//kkt.contract(tmp_box);

	if (tmp_box.is_empty()) {
		c.box.set_empty();
	} else {
		// the current extended box in the cell is updated
		write_ext_box(tmp_box,c.box);
	}
}

ThreadOptimizer::Status ThreadOptimizer::optimize(const IntervalVector& init_box, double obj_init_bound) {
	start(init_box, obj_init_bound);
	return optimize();
}

void ThreadOptimizer::start(const IntervalVector& init_box, double obj_init_bound) {

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes.back()=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	write_ext_box(init_box, root->box);

	// add data required by the bisector
	bsc.add_property(init_box, root->prop);

	// add data required by the contractor
	ctc.add_property(init_box, root->prop);

	// add data required by the buffer
	buffer.add_property(init_box, root->prop);

	// add data required by the loup finder
	loup_finder.add_property(init_box, root->prop);

	//cout << "**** Properties ****\n" << root->prop << endl;

	loup_changed=false;
	initial_loup=obj_init_bound;

	loup_point = init_box; //.set_empty();
	time=0;

	handle_cell(*root);
}

ThreadOptimizer::Status ThreadOptimizer::optimize() {
	bool yneginf=false;
	RealTimer timer;
	timer.start();
	bool isnttimeout=true;
	//update_uplo(); //lead to bugs

	loup_changed |= best_loup_found();	

	if (loup_changed) {
		// In case of a new upper bound (loup_changed == true), all the boxes
		// with a lower bound greater than (loup - goal_prec) are removed and deleted.
		// Note: if contraction was before bisection, we could have the problem
		// that the current cell is removed by contractHeap. See comments in
		// older version of the code (before revision 284).

		double ymax=compute_ymax();

		buffer.contract(ymax);
		update_uplo();
	
		//cout << " now buffer is contracted and min=" << buffer.minimum() << endl;

		// TODO: check if happens. What is the return code in this case?
		if (ymax <= NEG_INFINITY) {
			if (trace) cout << " infinite value for the minimum " << endl;
			yneginf=true;
		}
	}

	try {
	     while (!buffer.empty() && !yneginf) {
		  
			loup_changed=false;
			// for double heap , choose randomly the buffer : top  has to be called before pop
			Cell *c = buffer.top(); 
			if (trace >= 2) cout << " current box " << c->box << endl;

			try {

				pair<Cell*,Cell*> new_cells=bsc.bisect(*c);
				buffer.pop();
				delete c; // deletes the cell.

				nb_cells+=2;  // counting the cells handled ( in previous versions nb_cells was the number of cells put into the buffer after being handled)
                
				handle_cell(*new_cells.first);
				handle_cell(*new_cells.second);

				if (uplo_of_epsboxes.back() == NEG_INFINITY) {
					break;
				}
				if (loup_changed) {
					// In case of a new upper bound (loup_changed == true), all the boxes
					// with a lower bound greater than (loup - goal_prec) are removed and deleted.
					// Note: if contraction was before bisection, we could have the problem
					// that the current cell is removed by contractHeap. See comments in
					// older version of the code (before revision 284).

					double ymax=compute_ymax();

					buffer.contract(ymax);
				
					//cout << " now buffer is contracted and min=" << buffer.minimum() << endl;

					// TODO: check if happens. What is the return code in this case?
					if (ymax <= NEG_INFINITY) {
						if (trace) cout << " infinite value for the minimum " << endl;
						break;
					}
				}
				update_uplo();

				if (!anticipated_upper_bounding) // useless to check precision on objective if 'true'
					if (get_obj_rel_prec()<rel_eps_f || get_obj_abs_prec()<abs_eps_f)
						break;

				if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
				switch(obj){
					case ThreadOptimizer::Obj::OBJ_MIN:
						if (global_comm.stop_timeout_formin.load()>0) timer.check(global_comm.stop_timeout_formin.load());
						break;
					case ThreadOptimizer::Obj::SUB_OBJ_MIN:
						if (global_comm.stop_timeout_formin.load()>0) timer.check(global_comm.stop_timeout_formin.load());
						break;
					case ThreadOptimizer::Obj::OBJ_MAX:
						if (global_comm.stop_timeout_formax.load()>0) timer.check(global_comm.stop_timeout_formax.load());
						break;
					case ThreadOptimizer::Obj::SUB_OBJ_MAX:
						if (global_comm.stop_timeout_formax.load()>0) timer.check(global_comm.stop_timeout_formax.load());
						break;
					default:
						cerr << "not an implemented type of objective and system." << endl;
						exit(1);
				}
				time = timer.get_time();
			}
			catch (NoBisectableVariableException& ) {
				Interval& y=c->box[goal_var];
				y &= Interval(uplo,POS_INFINITY);
				if (y.is_empty()) {
					c->box.set_empty();
				} else {
					update_uplo_of_epsboxes((c->box)[goal_var].lb());
				}
				buffer.pop();
				delete c; // deletes the cell.
				update_uplo(); // the heap has changed -> recalculate the uplo (eg: if not in best-first search)
			}
		}

	 	timer.stop();
	 	switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
				global_comm.stop_timeout_formin.store(1e-15);
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
				break;
			case ThreadOptimizer::Obj::OBJ_MAX:
				global_comm.stop_timeout_formax.store(1e-15);
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
				break;
			default:
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
	 	time = timer.get_time();

		// No solution found and optimization stopped with empty buffer
		// before the required precision is reached => means infeasible problem
	 	if (uplo_of_epsboxes.back() == NEG_INFINITY)
	 		status = UNBOUNDED_OBJ;
	 	else if (uplo_of_epsboxes.back() == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0)))
	 		status = INFEASIBLE;
	 	else if (loup==initial_loup)
	 		status = NO_FEASIBLE_FOUND;
	 	else if (get_obj_rel_prec()>rel_eps_f && get_obj_abs_prec()>abs_eps_f)
	 		status = UNREACHED_PREC;
	 	else
	 		status = SUCCESS;
	}
	catch (RealTimeOutException& ) {
		status = TIME_OUT;
		isnttimeout=false;
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
				global_comm.stop_timeout_formin.store(1e-15);
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
				break;
			case ThreadOptimizer::Obj::OBJ_MAX:
				global_comm.stop_timeout_formax.store(1e-15);
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
				break;
			default:
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
	}
	if(isnttimeout){
		switch(obj){
			case ThreadOptimizer::Obj::OBJ_MIN:
				global_comm.statusformin.store(status);
				break;
			case ThreadOptimizer::Obj::OBJ_MAX:
				global_comm.statusformax.store(status);
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MIN:
				break;
			case ThreadOptimizer::Obj::SUB_OBJ_MAX:
				break;
			default:
				cerr << "not an implemented type of objective and system." << endl;
				exit(1);
		}
	}
	while (!buffer.empty()) {
		delete buffer.pop();
	}
	return status;
}

namespace {
const char* green() {
#ifndef _WIN32
	return "\033[32m";
#else
	return "";
#endif
}

const char* red(){
#ifndef _WIN32
	return "\033[31m";
#else
	return "";
#endif
}

const char* white() {
#ifndef _WIN32
	return "\033[0m";
#else
	return "";
#endif
}

}

void ThreadOptimizer::report() {

	if (/*!cov || */!buffer.empty()) { // not started
		cout << " not started." << endl;
		return;
	}

	switch(obj){
		case OBJ_MIN:
		{
			switch(status) {
				case SUCCESS: 
					cout << green() << "OBJ_MIN : optimization successful!" << endl;
					break;
				case INFEASIBLE: 
					cout << red() << "OBJ_MIN : infeasible problem" << endl;
					break;
				case NO_FEASIBLE_FOUND: 
					cout << red() << "OBJ_MIN : no feasible point found (the problem may be infeasible)" << endl;
					break;
				case UNBOUNDED_OBJ: 
					cout << red() << "OBJ_MIN : possibly unbounded objective (f*=-oo)" << endl;
					break;
				case TIME_OUT: 
					cout << red() << "OBJ_MIN : time limit " << timeout << "s. reached " << endl;
					break;
				case UNREACHED_PREC: 
					cout << red() << "OBJ_MIN : unreached precision" << endl;
					break;
			}
			break;
		}
		case OBJ_MAX:
		{
			switch(status) {
				case SUCCESS: 
					cout << green() << "OBJ_MAX : optimization successful!" << endl;
					break;
				case INFEASIBLE: 
					cout << red() << "OBJ_MAX : infeasible problem" << endl;
					break;
				case NO_FEASIBLE_FOUND: 
					cout << red() << "OBJ_MAX : no feasible point found (the problem may be infeasible)" << endl;
					break;
				case UNBOUNDED_OBJ: 
					cout << red() << "OBJ_MAX : possibly unbounded objective (f*=-oo)" << endl;
					break;
				case TIME_OUT: 
					cout << red() << "OBJ_MAX : time limit " << timeout << "s. reached " << endl;
					break;
				case UNREACHED_PREC: 
					cout << red() << "OBJ_MAX : unreached precision" << endl;
					break;
			}
			break;
		}
		case SUB_OBJ_MIN:
		{
			switch(status) {
				case SUCCESS: 
					cout << green() << "SUB_OBJ_MIN : optimization successful!" << endl;
					break;
				case INFEASIBLE: 
					cout << red() << "SUB_OBJ_MIN : infeasible problem" << endl;
					break;
				case NO_FEASIBLE_FOUND: 
					cout << red() << "SUB_OBJ_MIN : no feasible point found (the problem may be infeasible)" << endl;
					break;
				case UNBOUNDED_OBJ: 
					cout << red() << "SUB_OBJ_MIN : possibly unbounded objective (f*=-oo)" << endl;
					break;
				case TIME_OUT: 
					cout << red() << "SUB_OBJ_MIN : time limit " << timeout << "s. reached " << endl;
					break;
				case UNREACHED_PREC: 
					cout << red() << "SUB_OBJ_MIN : unreached precision" << endl;
					break;
			}
			break;
		}
		case SUB_OBJ_MAX:
		{
			switch(status) {
				case SUCCESS: 
					cout << green() << "SUB_OBJ_MAX : optimization successful!" << endl;
					break;
				case INFEASIBLE: 
					cout << red() << "SUB_OBJ_MAX : infeasible problem" << endl;
					break;
				case NO_FEASIBLE_FOUND: 
					cout << red() << "SUB_OBJ_MAX : no feasible point found (the problem may be infeasible)" << endl;
					break;
				case UNBOUNDED_OBJ: 
					cout << red() << "SUB_OBJ_MAX : possibly unbounded objective (f*=-oo)" << endl;
					break;
				case TIME_OUT: 
					cout << red() << "SUB_OBJ_MAX : time limit " << timeout << "s. reached " << endl;
					break;
				case UNREACHED_PREC: 
					cout << red() << "SUB_OBJ_MAX : unreached precision" << endl;
					break;
			}
			break;
		}
	}
	cout << white() <<  endl;

	// No solution found and optimization stopped with empty buffer
	// before the required precision is reached => means infeasible problem
	if (status==INFEASIBLE) {
		cout << " infeasible problem " << endl;
	} else {
		cout << " f* in\t[" << uplo << "," << loup << "]" << endl;
		cout << "\t(best bound)" << endl << endl;

		if (loup==initial_loup)
			cout << " x* =\t--\n\t(no feasible point found)" << endl;
		else {
			if (loup_finder.rigorous())
				cout << " x* in\t" << loup_point << endl;
			else
				cout << " x* =\t" << loup_point.lb() << endl;
			cout << "\t(best feasible point)" << endl;
		}
		cout << endl;
		/*double rel_prec=get_obj_rel_prec();
		double abs_prec=get_obj_abs_prec();

		cout << " relative precision on f*:\t" << rel_prec;
		if (rel_prec <= rel_eps_f)
			cout << green() << " [passed] " << white();
		cout << endl;

		cout << " absolute precision on f*:\t" << abs_prec;
		if (abs_prec <= abs_eps_f)
			cout << green() << " [passed] " << white();
		cout << endl;*/
	}

	//cout << " cpu time used:\t\t\t" << time << "s";
	/*if (cov->time()!=time)
		cout << " [total=" << cov->time() << "]";*/
	//cout << endl;
	//cout << " number of cells:\t\t" << nb_cells;
	/*if (cov->nb_cells()!=nb_cells)
		cout << " [total=" << cov->nb_cells() << "]";*/
	//cout << endl << endl;
}

} // end namespace ibex
