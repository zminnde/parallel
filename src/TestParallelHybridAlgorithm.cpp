
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <chrono> 
#include <random> 
#include <mpi.h>
#include <time.h>

#include <unistd.h>




#include "Solution.h"
#include "Function.h"
#include "OptimizedFunction.h"
#include "SurogateFunction.h"
#include "SurogateOneCriteriaFunction.h"
#include "Hypercube.h"
#include "ParallelHybridAlgorithm.h"
#include "LocalSearch.h"
#include "HookeJeeves.h"
#include "AllTestFunctions.h"
#include <iomanip>

#define ACCEPT_SUB_OPTIMAL_PARETO_SOLUTIONS true

using namespace std;


void calculate_mean_and_variance_and_standart_deviation(vector<double> &IGD_values, double &mean, double &variance, double &std_deviance){
    mean=0;
    for(int i=0;i<IGD_values.size();i++){
        mean+=IGD_values.at(i);
    }
    mean/=IGD_values.size();

    variance=0;
    for(int i=0;i<IGD_values.size();i++){
        variance+=(IGD_values.at(i)-mean)*(IGD_values.at(i)-mean);
    }
    variance/=(IGD_values.size());
    std_deviance=sqrt(variance);
}
  

double max_min(vector<Solution*>* P1,vector<Solution*>* P2){
    double max=-1;
    for (vector<Solution*>::iterator it = P1->begin(); it != P1->end(); ++it) {   
        double min_distance=-1;
        for (vector<Solution*>::iterator it1 = P2->begin(); it1 != P2->end(); ++it1) {
            if(min_distance<0)
                min_distance=(*it1)->distance_in_objective_space_not_normalized(*it);
            else{
                double distance=(*it1)->distance_in_objective_space_not_normalized(*it);
                if(min_distance>distance)
                    min_distance=distance;
            }
        }
        if(max<0)
            max=min_distance;
        else{
            if(max<min_distance)
                max=min_distance;
        }
    }
    return max;
}  



void put_to_Pareto_front( vector<Solution*>* P, Solution* s){
    if(P->empty()){
        P->push_back(s);
        double *evaluation_of_objective_functions = s->get_evaluation_of_objective_functions();
        int m=s->get_f()->get_m();
        if(s->get_f()->max_values_of_objective_functions == nullptr && s->get_f()->min_values_of_objective_functions == nullptr){
            s->get_f()->max_values_of_objective_functions = new double[m];
            s->get_f()->min_values_of_objective_functions = new double[m];
            for(int i=0;i<m;i++){
                s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
                s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
            }
        }
    }else{
        bool dominated=false;
        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
            if( s->is_dominated_by(*it) || s->equals(*it)){
                dominated=true;
                break;
            }
        }
        if(!dominated){
            double *evaluation_of_objective_functions = s->get_evaluation_of_objective_functions();
            int m=s->get_f()->get_m();
            for(int i=0;i<m;i++){
                s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
                s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions[i];
            }
            for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ) {
                if((*it)->is_dominated_by(s)){
                    it=P->erase(it);   
                }else{
                    double *evaluation_of_objective_functions1 = (*it)->get_evaluation_of_objective_functions();
                    for(int i=0;i<m;i++){
                        if(s->get_f()->max_values_of_objective_functions[i] < evaluation_of_objective_functions1[i])
                            s->get_f()->max_values_of_objective_functions[i] = evaluation_of_objective_functions1[i];
                        if(s->get_f()->min_values_of_objective_functions[i] > evaluation_of_objective_functions1[i])
                            s->get_f()->min_values_of_objective_functions[i] = evaluation_of_objective_functions1[i];
                    }

                    ++it;
                }
            }
            P->push_back(s);
        }
    }


    /*if(P->empty()){
        P->push_back(s);
    }else{
        bool dominated=false;
        for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ++it) {
            if( s->is_dominated_by(*it) || s->equals(*it)){
                dominated=true;
                delete s;
                break;
            }
        }
        if(!dominated){
            for (vector<Solution*>::iterator it = P->begin(); it != P->end(); ) {
                if((*it)->is_dominated_by(s)){
                    it=P->erase(it);   
                }else{
                    ++it;
                }
            }
            P->push_back(s);
        }
    }*/


}



/*double Euclidean_distance(Solution* p1,Solution* p2){
    double distance=0;
    double *val1=p1->get_evaluation_of_objective_functions();
    double *val2=p2->get_evaluation_of_objective_functions();
    for(int i=0; i<p1->get_f()->get_m(); i++){
            distance+=(val1[i]-val2[i])*(val1[i]-val2[i]);
            //cout<<"Not normalized";
    }
    distance=pow(distance,0.5);
    return distance;
}*/




void put_ZDT_function_to_vector(Function* m_f, vector<vector<Solution*>*> &true_Pareto_fronts, double  optimal_normalized_const){
    
    
    vector<Solution*>* true_Pareto_front = new vector<Solution*>;
    
    vector<double> normalized_loc;
    normalized_loc.assign(m_f->get_d(), optimal_normalized_const);// for ZDT4 should be 0.5
                                          
    int number_of_points=10000;//10000
    for(int i=0;i<number_of_points+1;i++){
        normalized_loc.at(0)=1.0*i/number_of_points;
        Solution* s = new Solution(m_f);
        s->evaluate_function(normalized_loc);
        put_to_Pareto_front( true_Pareto_front,  s);
    }
    //cout<<"true_Pareto_front->size() : "<<true_Pareto_front->size()<<endl; 
    /*for (vector<Solution*>::iterator it = true_Pareto_front->begin(); it != true_Pareto_front->end(); ++it) {
        //cout<<"*it->to_string()"<<(*it)->to_string()<<endl;
    } */
    //>>>>>>>>>>>  REMOVE EXTRA POINTS  >>>>>>>>>>>>>>>>>>>

    
    for(int i=0;i<true_Pareto_front->size()-1;i++)
        for(int i1=i+1;i1<true_Pareto_front->size();i1++){
            if(true_Pareto_front->at(i)->get_evaluation_of_objective_functions()[0] > true_Pareto_front->at(i1)->get_evaluation_of_objective_functions()[0]){
                Solution *p_tmps=true_Pareto_front->at(i);
                true_Pareto_front->at(i)=true_Pareto_front->at(i1);
                true_Pareto_front->at(i1)=p_tmps;
            }
        }
    //cout<<"Euclidean_distance(true_Pareto_front->at(0),true_Pareto_front->at(true_Pareto_front->size()-1)) : "<<Euclidean_distance(true_Pareto_front->at(0),true_Pareto_front->at(true_Pareto_front->size()-1))<<endl;
    double initial_distance=true_Pareto_front->at(0)->distance_in_objective_space_not_normalized(true_Pareto_front->at(true_Pareto_front->size()-1))/500;
    //cout<<"true_Pareto_front->at(0)->distance_in_objective_space_not_normalized(true_Pareto_front->at(true_Pareto_front->size()-1)) : "<<true_Pareto_front->at(0)->distance_in_objective_space_not_normalized(true_Pareto_front->at(true_Pareto_front->size()-1))<<endl;
    //cout<<"Euclidean_distance(true_Pareto_front->at(0),true_Pareto_front->at(true_Pareto_front->size()-1)) : "<<Euclidean_distance(true_Pareto_front->at(0),true_Pareto_front->at(true_Pareto_front->size()-1))<<endl;
    //cout<<"initial_distance"<<initial_distance<<endl;
    double step_size=initial_distance*0.02;

    vector<Solution*>* true_Pareto_front_500 = new vector<Solution*>;
    do {
        delete true_Pareto_front_500;
        true_Pareto_front_500 = new vector<Solution*>;
        true_Pareto_front_500->push_back(true_Pareto_front->at(0));

        for(int i=1;i<true_Pareto_front->size()-1;i++){
            if(initial_distance<true_Pareto_front_500->at(true_Pareto_front_500->size()-1)->distance_in_objective_space_not_normalized(true_Pareto_front->at(i))){
                true_Pareto_front_500->push_back(true_Pareto_front->at(i));
                //true_Pareto_front_500->push_back(true_Pareto_front->at(0));
            }
            //cout<<"true_Pareto_front->at(i)->get_val().data[0] : "<<true_Pareto_front->at(i)->get_val().data[0]<<endl;
        }
        true_Pareto_front_500->push_back(true_Pareto_front->at(true_Pareto_front->size()-1));
        if(true_Pareto_front_500->size()>505){
            initial_distance+=step_size;
        }
        if(true_Pareto_front_500->size()<495){
            initial_distance-=step_size;
        }
        //cout<<"11true_Pareto_front_500->size()"<<true_Pareto_front_500->size()<<endl;
    }
    while ( (true_Pareto_front_500->size()>505 || true_Pareto_front_500->size()<495));

    //cout<<"true_Pareto_front->size()"<<true_Pareto_front->size()<<endl;
    //cout<<"true_Pareto_front_500->size()"<<true_Pareto_front_500->size()<<endl;
    //exit(0);
    //<<<<<<<<<<<  END REMOVE EXTRA POINTS  <<<<<<<<<<<<<<<

    true_Pareto_fronts.push_back(true_Pareto_front_500);
}




void optimise_ZDT(int world_rank) { 
    int number_of_experiments=10;/*10*/
    
    
    vector<vector<Solution*>*> true_Pareto_fronts;

    
    put_ZDT_function_to_vector(new F_ZDT1(),true_Pareto_fronts,0);//1

    put_ZDT_function_to_vector(new F_ZDT2(),true_Pareto_fronts,0);//2
    put_ZDT_function_to_vector(new F_ZDT3(),true_Pareto_fronts,0);//3

    put_ZDT_function_to_vector(new F_ZDT4(),true_Pareto_fronts,0.5);//4
    put_ZDT_function_to_vector(new F_ZDT6(),true_Pareto_fronts,0);//6

    //cout<<"true_Pareto_fronts.size() : "<<true_Pareto_fronts.size()<<endl;
    for(int i=0;i<true_Pareto_fronts.size();i++){
        vector<double> GD_values;
        vector<double> actual_number_of_function_evaluations_values;
        vector<double> run_time_values;
        //Function * m_f = functions.at(i);
        //cout<<"Function number: "<<i<<endl;
        vector<Solution*> *true_Pareto_front = true_Pareto_fronts.at(i);
        //cout<<"true_Pareto_front->size()"<<true_Pareto_front->size()<<endl;

        for(int i1 = 0;i1<number_of_experiments;i1++){
            //cout<<"i1 : "<<i1<<endl;
            Function *of1;
            if(i==0){
                //cout<<"optimize F_ZDT1"<<endl;
                of1 = new F_ZDT1();
            }
            if(i==1){
                //cout<<"optimize F_ZDT2"<<endl;
                of1 = new F_ZDT2();
            }
            if(i==2){
                //cout<<"optimize F_ZDT3"<<endl;
                of1 = new F_ZDT3();
            }
            if(i==3){
                //cout<<"optimize F_ZDT4"<<endl;
                of1 = new F_ZDT4();
            }
            if(i==4){
                //cout<<"optimize F_ZDT6"<<endl;
                of1 = new F_ZDT6();
            }
            ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
            params->f=of1;
            
            params->I_max=2000;
            params->N_max=23000;//90

            params->N=100;//
            params->q=1;//
            params->p=0.8;//

            params->h0=2;//2//4
            params->update_h0_and_hn=true;
            params->hn=8;//
            params->accept_sub_optimal_Pareto_solutions=ACCEPT_SUB_OPTIMAL_PARETO_SOLUTIONS;

            ParallelHybridAlgorithm *parallel1 = new ParallelHybridAlgorithm(params); 
            parallel1->optimize();

            if(world_rank==0){
                
                double GD_avg=0;
                for (vector<Solution*>::iterator it = parallel1->get_P()->begin(); it != parallel1->get_P()->end(); ++it) {
                    //cout<<(*it)->to_string()<<endl;
                    double min=-1;
                    for (vector<Solution*>::iterator it1 = true_Pareto_front->begin(); it1 != true_Pareto_front->end(); ++it1) {
                        //(**it1)
                        double ed=(*it)->distance_in_objective_space_not_normalized(*it1);
                        if(min<0)
                            min=ed;
                        else if(min>ed)
                            min=ed;
                    }
                    GD_avg+=min/parallel1->get_P()->size();
                    
                }
                
                
                GD_values.push_back(GD_avg);
                actual_number_of_function_evaluations_values.push_back(parallel1->get_total_function_evaluation_count());
                run_time_values.push_back(parallel1->get_run_time());
            }
            delete parallel1;
            
        }
        
        double mean = 0;
        double variance = 0;
        double standart_deviation = 0;
        if(world_rank==0){
            calculate_mean_and_variance_and_standart_deviation(GD_values, mean,variance, standart_deviation);
            cout<<"--------------------------- GD mean : "<<mean<<" GD variance : "<<variance<<" GD standart deviation : "<<standart_deviation<<endl;

            calculate_mean_and_variance_and_standart_deviation(actual_number_of_function_evaluations_values, mean,variance, standart_deviation);
            cout<<"--------------------------- actual_number_of_function_evaluations mean : "<<mean<<" actual_number_of_function_evaluations variance : "<<variance<<" actual_number_of_function_evaluations standart deviation : "<<standart_deviation<<endl;
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            cout<<"actual_number_of_function_evaluations mean / number_of_CPU : "<<mean/world_size<<endl;
            calculate_mean_and_variance_and_standart_deviation(run_time_values, mean,variance, standart_deviation);
            cout<<"--------------------------- run_time_values mean : "<<mean<<" run_time_values variance : "<<variance<<" run_time_values standart deviation : "<<standart_deviation<<endl<<endl<<endl;
        }
    }

    return;
}


void optimise_ZDT1(int world_rank) { 
    int number_of_experiments=51;/*51*/
    vector<vector<Solution*>*> true_Pareto_fronts;


    put_ZDT_function_to_vector(new F_ZDT1(),true_Pareto_fronts,0);//1

    put_ZDT_function_to_vector(new F_ZDT2(),true_Pareto_fronts,0);//2
    put_ZDT_function_to_vector(new F_ZDT3(),true_Pareto_fronts,0);//3

    put_ZDT_function_to_vector(new F_ZDT4(),true_Pareto_fronts,0.5);//4
    put_ZDT_function_to_vector(new F_ZDT6(),true_Pareto_fronts,0);//6

    //cout<<"true_Pareto_fronts.size() : "<<true_Pareto_fronts.size()<<endl;
    for(int i=0;i<true_Pareto_fronts.size();i++){
        vector<double> IGD_values;
        vector<double> actual_number_of_function_evaluations_values;
        vector<double> run_time_values;
        //Function * m_f = functions.at(i);
        //cout<<"Function number: "<<i<<endl;
        vector<Solution*> *true_Pareto_front = true_Pareto_fronts.at(i);
        //cout<<"true_Pareto_front->size()"<<true_Pareto_front->size()<<endl;
        
        for(int i1 = 0;i1<number_of_experiments;i1++){
            //cout<<"i1 : "<<i1<<endl;
            Function *of1;
            ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
            if(i==0){
                //cout<<"optimize F_ZDT1"<<endl;
                of1 = new F_ZDT1();
                params->N_max=17098.0*0.9;
            }
            if(i==1){
                //cout<<"optimize F_ZDT2"<<endl;
                of1 = new F_ZDT2();
                params->N_max=17657.0*0.9;
            }
            if(i==2){
                //cout<<"optimize F_ZDT3"<<endl;
                of1 = new F_ZDT3();
                params->N_max=16559.0*0.9;
            }
            if(i==3){
                //cout<<"optimize F_ZDT4"<<endl;
                of1 = new F_ZDT4();
                params->N_max=24451.0*0.9;
            }
            if(i==4){
                //cout<<"optimize F_ZDT6"<<endl;
                of1 = new F_ZDT6();
                params->N_max=24833.0*0.9;
            }
            
            params->f=of1;
            
            params->I_max=2000;
            //params->N_max=230000;//90

            params->N=100;//
            params->q=1;//
            params->p=0.8;//
            params->update_h0_and_hn=true;
            params->h0=2;//2//4
            params->hn=8;//
            params->accept_sub_optimal_Pareto_solutions=ACCEPT_SUB_OPTIMAL_PARETO_SOLUTIONS;

            ParallelHybridAlgorithm *parallel1 = new ParallelHybridAlgorithm(params); 
            parallel1->optimize();


            if(world_rank==0){
                double IGD_avg=0;
                for (vector<Solution*>::iterator it = true_Pareto_front->begin(); it != true_Pareto_front->end(); ++it) {
                    //cout<<(*it)->to_string()<<endl;
                    double min=-1;
                    for (vector<Solution*>::iterator it1 = parallel1->get_P()->begin(); it1 != parallel1->get_P()->end(); ++it1) {
                        //(**it1)
                        double ed=(*it)->distance_in_objective_space(*it1);
                        if(min<0)
                            min=ed;
                        else if(min>ed)
                            min=ed;
                    }
                    IGD_avg+=min/true_Pareto_front->size();
                    
                }
                
                IGD_values.push_back(IGD_avg);
                actual_number_of_function_evaluations_values.push_back(parallel1->get_total_function_evaluation_count());
                run_time_values.push_back(parallel1->get_run_time());
            }

            delete parallel1;
        }
        double mean = 0;
        double variance = 0;
        double standart_deviation = 0;
        if(world_rank==0){
            calculate_mean_and_variance_and_standart_deviation(IGD_values, mean,variance, standart_deviation);
            cout<<"--------------------------- IGD mean : "<<mean<<" IGD variance : "<<variance<<" IGD standart deviation : "<<standart_deviation<<endl;

            calculate_mean_and_variance_and_standart_deviation(actual_number_of_function_evaluations_values, mean,variance, standart_deviation);
            cout<<"--------------------------- actual_number_of_function_evaluations mean : "<<mean<<" actual_number_of_function_evaluations variance : "<<variance<<" actual_number_of_function_evaluations standart deviation : "<<standart_deviation<<endl;
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            cout<<"actual_number_of_function_evaluations mean / number_of_CPU : "<<mean/world_size<<endl;
            calculate_mean_and_variance_and_standart_deviation(run_time_values, mean,variance, standart_deviation);
            cout<<"--------------------------- run_time_values mean : "<<mean<<" run_time_values variance : "<<variance<<" run_time_values standart deviation : "<<standart_deviation<<endl<<endl<<endl;
        }
    }

    return;
}







void conference_function1(int world_rank) {

    int number_of_experiments=100;
    double array_of_GD[number_of_experiments];
    double array_of_EI[number_of_experiments];
    double array_of_U_size[number_of_experiments];
    double array_of_run_times[number_of_experiments];

    Function *of1 = new F_conference1();

    ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
    params->f=of1;
    params->I_max=1;
    params->N_max=1;

    if(world_rank==0)
        params->N=1000000;
    else
        params->N=100;


    ParallelHybridAlgorithm *parallel = new ParallelHybridAlgorithm(params); 
    parallel->optimize();
    vector<Solution*>* P_full = parallel->get_P();


    for(int i=0;i<number_of_experiments;i++){
        //cout<<"i : "<<i<<endl;
        Function *of1 = new F_conference1();

        ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
        params->f=of1;
        params->I_max=10000;
        params->N_max=90;//90

        params->N=20;//
        params->q=10000;//
        params->p=0.8;//

        params->h0=4;//2//4
        params->hn=4;//
        params->update_h0_and_hn=true;
        params->accept_sub_optimal_Pareto_solutions=ACCEPT_SUB_OPTIMAL_PARETO_SOLUTIONS;

        ParallelHybridAlgorithm *parallel1 = new ParallelHybridAlgorithm(params); 
        parallel1->optimize();
        if(world_rank==0){
            //cout<<"GD "<<max_min(parallel1->get_P(), P_full)<<endl;
            array_of_GD[i]=max_min(parallel1->get_P(), P_full);

            //cout<<"EI "<<max_min( P_full, parallel1->get_P())<<endl;
            array_of_EI[i]=max_min( P_full, parallel1->get_P());
            array_of_U_size[i]=parallel1->get_total_function_evaluation_count();
            array_of_run_times[i]=parallel1->get_run_time();
        }
        if(IS_TESTING){
            cout<<"FULL PARETO OPTIMAL POINTS:"<<endl;
            for (vector<Solution*>::iterator it = parallel->get_P()->begin(); it != parallel->get_P()->end(); ++it) {
                //break;
                //double *upper_bound=(*it)->get_f()->upper_bound;
                //double *lower_bound=(*it)->get_f()->lower_bound;
                cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
                //cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
            }
            cout<<"PARETO OPTIMAL POINTS:"<<endl;
            cout<<endl<<endl;
            for (vector<Solution*>::iterator it = parallel1->get_P()->begin(); it != parallel1->get_P()->end(); ++it) {
                //break;
                //double *upper_bound=(*it)->get_f()->upper_bound;
                //double *lower_bound=(*it)->get_f()->lower_bound;
                cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
                //cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
            }
        }
        /*cout<<endl;
        cout<<endl;
        for (vector<Solution*>::iterator it = parallel->get_U()->begin(); it != parallel->get_U()->end(); ++it) {
            break;
            double *upper_bound=(*it)->get_f()->upper_bound;
            double *lower_bound=(*it)->get_f()->lower_bound;
            //cout<<"(*it)->get_foud_by_local_search()"<<(*it)->get_foud_by_local_search()<<endl;
            //cout<<"(*it)->to_string() : "<< (*it)->to_string()<<endl;
            //cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
            cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
        }
        cout<<endl;
        cout<<endl;
        for (vector<Solution*>::iterator it = parallel->get_P()->begin(); it != parallel->get_P()->end(); ++it) {
            break;
            double *upper_bound=(*it)->get_f()->upper_bound;
            double *lower_bound=(*it)->get_f()->lower_bound;
            //cout<<"(*it)->get_foud_by_local_search()"<<(*it)->get_foud_by_local_search()<<endl;
            //cout<<"(*it)->to_string() : "<< (*it)->to_string()<<endl;
            //cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
            cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
        }*/
        delete parallel1;
    }
    if(world_rank==0){
        double mean_GD=0;
        double mean_EI=0;
        double mean_U_size=0;
        double mean_run_time=0;
        for(int i=0;i<number_of_experiments;i++){
            mean_GD+=array_of_GD[i];
            mean_EI+=array_of_EI[i];
            mean_U_size+=array_of_U_size[i];
            mean_run_time+=array_of_run_times[i];
            //cout<<"array_of_U_size[i] : "<<array_of_U_size[i];
        }
        mean_GD/=number_of_experiments;
        mean_EI/=number_of_experiments;
        mean_U_size/=number_of_experiments;
        mean_run_time/=number_of_experiments;

        cout<<"F_conference1 mean_GD : "<<mean_GD<<endl;
        cout<<"mean_EI : "<<mean_EI<<endl;
        cout<<"mean_U_size : "<<mean_U_size<<endl;
        cout<<"mean_run_time (miliseconds): "<<mean_run_time<<endl;
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        //cout<<"world_size : "<<world_size<<endl;
        cout<<"mean_U_size / number_of_CPU : "<<mean_U_size/world_size<<endl;
    }

    delete parallel;

}

void conference_function(int world_rank) {

    int number_of_experiments=100;//must be 100
    double array_of_GD[number_of_experiments];
    double array_of_EI[number_of_experiments];
    double array_of_U_size[number_of_experiments];
    double array_of_run_times[number_of_experiments];

    Function *of1 = new F_conference();

    ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
    params->f=of1;
    params->I_max=1;
    params->N_max=1;
    if(world_rank==0)
        params->N=1000000;
    else
        params->N=100;

    ParallelHybridAlgorithm *parallel = new ParallelHybridAlgorithm(params); 
    
    parallel->optimize();
    vector<Solution*>* P_full = parallel->get_P();
    

    for(int i=0;i<number_of_experiments;i++){
        //cout<<"i : "<<i<<endl;
        Function *of1 = new F_conference();

        ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
        params->f=of1;
        params->I_max=1000000;
        params->N_max=90;//90

        params->N=20;//
        params->q=10000;//
        params->p=0.8;//

        params->h0=2;//2//4
        params->hn=4;//
        params->update_h0_and_hn=true;
        params->accept_sub_optimal_Pareto_solutions=ACCEPT_SUB_OPTIMAL_PARETO_SOLUTIONS;


        ParallelHybridAlgorithm *parallel1 = new ParallelHybridAlgorithm(params); 
        parallel1->optimize();
        //cout<<"parallel1->get_P()->size(): "<<parallel1->get_P()->size()<<" world_rank: "<<world_rank<<endl;
        if(world_rank==0){
            //cout<<"GD "<<max_min(parallel1->get_P(), P_full)<<endl;
            array_of_GD[i]=max_min(parallel1->get_P(), P_full);

            //cout<<"EI "<<max_min( P_full, parallel1->get_P())<<endl;
            array_of_EI[i]=max_min( P_full, parallel1->get_P());
            array_of_U_size[i]=parallel1->get_total_function_evaluation_count();
            array_of_run_times[i]=parallel1->get_run_time();
        }
        if(IS_TESTING){
            cout<<"FULL PARETO OPTIMAL POINTS:"<<endl;
            for (vector<Solution*>::iterator it = parallel->get_P()->begin(); it != parallel->get_P()->end(); ++it) {
                //break;
                //double *upper_bound=(*it)->get_f()->upper_bound;
                //double *lower_bound=(*it)->get_f()->lower_bound;
                cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
                //cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
            }
            cout<<"PARETO OPTIMAL POINTS:"<<endl;
            cout<<endl<<endl;
            for (vector<Solution*>::iterator it = parallel1->get_P()->begin(); it != parallel1->get_P()->end(); ++it) {
                //break;
                //double *upper_bound=(*it)->get_f()->upper_bound;
                //double *lower_bound=(*it)->get_f()->lower_bound;
                cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
                //cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
            }
        }
        /*cout<<endl;
        cout<<endl;
        for (vector<Solution*>::iterator it = parallel->get_U()->begin(); it != parallel->get_U()->end(); ++it) {
            break;
            double *upper_bound=(*it)->get_f()->upper_bound;
            double *lower_bound=(*it)->get_f()->lower_bound;
            //cout<<"(*it)->get_foud_by_local_search()"<<(*it)->get_foud_by_local_search()<<endl;
            //cout<<"(*it)->to_string() : "<< (*it)->to_string()<<endl;
            //cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
            cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
        }
        cout<<endl;
        cout<<endl;
        for (vector<Solution*>::iterator it = parallel->get_P()->begin(); it != parallel->get_P()->end(); ++it) {
            break;
            double *upper_bound=(*it)->get_f()->upper_bound;
            double *lower_bound=(*it)->get_f()->lower_bound;
            //cout<<"(*it)->get_foud_by_local_search()"<<(*it)->get_foud_by_local_search()<<endl;
            //cout<<"(*it)->to_string() : "<< (*it)->to_string()<<endl;
            //cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;//   get_decision_point()[0]
            cout<<lower_bound[0]+(upper_bound[0]-lower_bound[0])* (*it)->get_decision_point()[0]<<" "<<lower_bound[1]+(upper_bound[1]-lower_bound[1])* (*it)->get_decision_point()[1]<<endl;
        }*/
        delete parallel1;
    }
    if(world_rank==0){
        double mean_GD=0;
        double mean_EI=0;
        double mean_U_size=0;
        double mean_run_time=0;
        for(int i=0;i<number_of_experiments;i++){
            mean_GD+=array_of_GD[i];
            mean_EI+=array_of_EI[i];
            mean_U_size+=array_of_U_size[i];
            mean_run_time+=array_of_run_times[i];
            //cout<<"array_of_U_size[i] : "<<array_of_U_size[i];
        }
        mean_GD/=number_of_experiments;
        mean_EI/=number_of_experiments;
        mean_U_size/=number_of_experiments;
        mean_run_time/=number_of_experiments;

        cout<<"F_conference mean_GD : "<<mean_GD<<endl;
        cout<<"mean_EI : "<<mean_EI<<endl;
        cout<<"mean_U_size : "<<mean_U_size<<endl;
        cout<<"mean_run_time (miliseconds): "<<mean_run_time<<endl;
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        cout<<"mean_U_size / number_of_CPU : "<<mean_U_size/world_size<<endl;
    }


    delete parallel;

}


// ---------------------------------------------------------------------------
// Modeled on conference_java(); objective values come from scripts/worker.py.
//
//   f1(x) = (x1 - 1)^2 + (x2 - 2)^2
//   f2(x) = (x1 + 1)^2 + (x2 + 2)^2
//
// Analytic Pareto set: the straight segment between (1,2) and (-1,-2).
// ---------------------------------------------------------------------------

// SHORT TEST: 10000.  Set to 10000000 for the full run.
// Override without editing: mpicxx -DPYTHON_N_MAX=...
#ifndef PYTHON_N_MAX
#define PYTHON_N_MAX 10000
#endif

void conference_python(int world_rank) {

    int number_of_experiments=1;//must be 1

    for(int i=0;i<number_of_experiments;i++){

        Function *of1 = new PythonFunction();

        ParallelHybridAlgorithmParams *params=new ParallelHybridAlgorithmParams();
        params->f=of1;
        params->I_max=1;
        params->N_max=PYTHON_N_MAX;

        params->N=1000;//
        params->q=10;//
        params->p=0.1;//

        params->h0=4;//2//4
        params->hn=6;//
        params->update_h0_and_hn=true;
        params->accept_sub_optimal_Pareto_solutions=false;


        ParallelHybridAlgorithm *parallel1 = new ParallelHybridAlgorithm(params);
        parallel1->optimize();
        if(1||world_rank==0)
            cout<<"parallel1->get_P()->size(): "<<parallel1->get_P()->size()<<" world_rank: "<<world_rank<<endl;

        if(world_rank==0){

            cout<<"PARETO OPTIMAL POINTS (f1 f2):"<<endl;
            cout<<endl<<endl;
            for (vector<Solution*>::iterator it = parallel1->get_P()->begin(); it != parallel1->get_P()->end(); ++it) {
                cout<<(*it)->get_evaluation_of_objective_functions()[0]<<" "<<(*it)->get_evaluation_of_objective_functions()[1]<<endl;
            }
            cout<<endl<<endl;

            // denormalized decision points: these are what must lie on the
            // segment between (1,2) and (-1,-2)
            cout<<"PARETO OPTIMAL DECISION POINTS (x1 x2):"<<endl;
            for (vector<Solution*>::iterator it = parallel1->get_P()->begin(); it != parallel1->get_P()->end(); ++it) {
                double *upper_bound=(*it)->get_f()->upper_bound;
                double *lower_bound=(*it)->get_f()->lower_bound;
                for(int j=0;j<2;j++){
                    cout<<lower_bound[j]+(upper_bound[j]-lower_bound[j])*(*it)->get_decision_point()[j]<<" ";
                }
                cout<<endl;
            }
            cout<<endl<<endl;

            cout<<"Total function evaluations: "<<parallel1->get_total_function_evaluation_count()<<endl;
            cout<<"Run time (miliseconds)"<<parallel1->get_run_time()<<endl;
            cout<<"Run time (seconds)"<<parallel1->get_run_time()/1000.0<<endl;
            cout<<"Run time (minutes)"<<parallel1->get_run_time()/60000.0<<endl;
        }
        delete parallel1;
    }

}


int main() {
    clock_t start = clock();
    MPI_Init(NULL, NULL);
    cout<<"--------------------------------------------------------------------------------"<<endl;
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //conference_function( world_rank);
    //conference_function1(world_rank);
    

    //optimise_ZDT(world_rank);

    
    //optimise_ZDT1(world_rank);

    conference_python( world_rank);
    

    
    
    MPI_Finalize();

    clock_t end = clock();
    double run_time = double(end - start)/(CLOCKS_PER_SEC);
    
    if(world_rank==0){

	   cout<<endl << "END 19  !!!" << endl;
        cout<<"Is testing: "<<IS_TESTING<<endl;
        cout<<"run_time (seconds): "<<run_time<<endl;
    }
    
    
}



























