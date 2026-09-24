#ifndef SRC_ALL_TEST_FUNCTIONS_H_
#define SRC_ALL_TEST_FUNCTIONS_H_

#include <iostream>
#include <cassert>

#include "Function.h"
#include "PythonFunction.h"
#include <cmath>

using namespace std;



#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif





class OptimizedFunction11: public Function { 
protected:
	
	double* evaluate_func_values(double* decision_point){
		//static int i=0;
		//i++;
		//cout<<"i : "<<i<<endl;
		double * func_values= new double[m];
		func_values[0]=(decision_point[0]-0.5)*(decision_point[0]-0.5)+(decision_point[1]-1)*(decision_point[1]-1);
		func_values[1]=(decision_point[0]+1)*(decision_point[0]+1)+(decision_point[1]+1)*(decision_point[1]+1);
		//func_values[0]=100*(decision_point[0]*decision_point[0]+decision_point[1]*decision_point[1]);
		//func_values[0]=decision_point[0]+decision_point[1];
		return func_values;
	}	
  	
public:

	OptimizedFunction11()  {
		d=2;//  dimension of decision space
		m=2;//  number of objective functions
		double upper_bound[]={2,3}; // decision  space bounds
	  	double lower_bound[]={-1,-4}; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
	}

	~OptimizedFunction11() {
		
	}
};








class F_conference: public Function { 
protected:
	
	double* evaluate_func_values(double* decision_point){
		double * func_values= new double[m];
		func_values[0]=1  -exp(-(  (decision_point[0]-1/sqrt(2))*(decision_point[0]-1/sqrt(2)) + (decision_point[1]-1/sqrt(2))*(decision_point[1]-1/sqrt(2)) ));      
		func_values[1]=1  -exp(-(  (decision_point[0]+1/sqrt(2))*(decision_point[0]+1/sqrt(2)) + (decision_point[1]+1/sqrt(2))*(decision_point[1]+1/sqrt(2)) ));

		return func_values;
	}	
  	
public:

	F_conference()  {
		d=2;//  dimension of decision space
		m=2;//  number of objective functions
		double upper_bound[]={ 4, 4 }; // decision  space bounds
	  double lower_bound[]={ -4, -4 }; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
	}

	~F_conference() {
		
	}
};






class F_conference1: public Function { 
protected:
	
	double* evaluate_func_values(double* decision_point){
		double * func_values= new double[m];
		func_values[0]=-0.1/(0.1+(decision_point[0]-0.1)*(decision_point[0]-0.1) +2*(decision_point[1]-0.1)*(decision_point[1]-0.1)) 
			-0.1/(0.14+20*(decision_point[0]-0.45)*(decision_point[0]-0.45) +20*(decision_point[1]-0.55)*(decision_point[1]-0.55))  ;      
		func_values[1]=-0.1/(0.15+40*(decision_point[0]-0.55)*(decision_point[0]-0.55) +40*(decision_point[1]-0.45)*(decision_point[1]-0.45)) 
			-0.1/(0.1+(decision_point[0]-0.3)*(decision_point[0]-0.3) +(decision_point[1]-0.95)*(decision_point[1]-0.95));

		return func_values;
	}	
  	
public:

	F_conference1()  {
		d=2;//  dimension of decision space
		m=2;//  number of objective functions
		double upper_bound[]={ 1, 1 }; // decision  space bounds
	    double lower_bound[]={ 0, 0 }; // decision  space bounds
		intailize(d, upper_bound, lower_bound);
	}

	~F_conference1() {
		
	}
};


class F_ZDT: public Function {
    
public:
	virtual double f1(double x1)=0;
    virtual double g(double *x)=0;
    virtual double h(double f1, double g)=0;

protected:

    double* evaluate_func_values(double* x1){
      
        double* ats= new double[m];
        ats[0]=f1(x1[0]);      
        ats[1]=g(x1)*h(f1(x1[0]),g(x1));
        //cout<<"OKI"<<M_PI<<endl;
        return ats;
    }
};



class F_ZDT1: public F_ZDT {  
public:
    
    F_ZDT1() {
		m=2;//  number of objective functions
        d = 30;
        double lower_bound[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        double upper_bound[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

        intailize(d, upper_bound, lower_bound);

        //_xmin.push_back(new Vec_double());
        //init_vector(*_xmin.at(0), xmin0, _dim);
    }
    virtual ~F_ZDT1() {
    }
    virtual double f1(double x1){
        return x1;
    };
    virtual double g(double *x){
        double sum;
        sum=0;
        for(int i =1;i<d;i++){
            sum+=x[i]/(d-1);
        }
        return 1+9*sum;
    };
    virtual double h(double f1, double g){
        assert(f1/g>=0);
        return 1-sqrt(f1/g);
    };

};



class F_ZDT2: public F_ZDT {  
public:

    F_ZDT2() {
		m=2;//  number of objective functions
        d = 30;
        double lower_bound[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        double upper_bound[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

        intailize(d, upper_bound, lower_bound);

        //_xmin.push_back(new Vec_double());
        //init_vector(*_xmin.at(0), xmin0, _dim);
    }




    virtual ~F_ZDT2() {
    }
    virtual double f1(double x1){
        return x1;
    };
    virtual double g(double *x){
        double sum;
        sum=0;
        for(int i =1;i<d;i++){
            sum+=x[i]/(d-1);
        }
        return 1+9*sum;
    };
    virtual double h(double f1, double g){
        assert(f1/g>=0);
        return 1-pow(f1/g,2);
    };
    
};



class F_ZDT3: public F_ZDT {  
public:
    
    F_ZDT3() {
		m=2;//  number of objective functions
        d = 30;
        double lower_bound[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
        double upper_bound[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

        intailize(d, upper_bound, lower_bound);

        //_xmin.push_back(new Vec_double());
        //init_vector(*_xmin.at(0), xmin0, _dim);
    }
    virtual ~F_ZDT3() {
    }

    virtual double f1(double x1){
        return x1;
    };
    virtual double g(double *x){
        double sum;
        sum=0;
        for(int i =1;i<d;i++){
            sum+=x[i]/(d-1);
        }
        return 1+9*sum;
    };
    virtual double h(double f1, double g){
        
        assert(f1/g>=0);
        return 1-sqrt(f1/g)-sin(10*M_PI*f1)*f1/g;
    };



};









class F_ZDT4: public F_ZDT {  
public:
    F_ZDT4() {
        m=2;
        d = 10;
        //
        double lower_bound[] = {0,-5,-5,-5,-5,-5,-5,-5,-5,-5 };
        double upper_bound[] = {1,5,5,5,5,5,5,5,5,5 };

        intailize(d, upper_bound, lower_bound);

        //_xmin.push_back(new Vec_double());
        //init_vector(*_xmin.at(0), xmin0, _dim);
    }
    

    virtual ~F_ZDT4() {
    }

    virtual double f1(double x1){
        return x1;
    };
    virtual double g(double *x){
        double sum;
        sum=1+10*(d-1);
        for(int i =1;i<d;i++){
            sum+=x[i]*x[i]-10*cos ( 4*M_PI*x[i] );;
        }
        return sum;
    };
    virtual double h(double f1, double g){
        
        assert(f1/g>=0);
        return 1-sqrt(f1/g);
    };

};




class F_ZDT6: public F_ZDT {  
public:

    F_ZDT6() {
        m = 2;
        d = 10;
        //
        double lower_bound[] = {0,0,0,0,0,0,0,0,0,0 };
        double upper_bound[] = {1,1,1,1,1,1,1,1,1,1 };

        intailize(d, upper_bound, lower_bound);

        //_xmin.push_back(new Vec_double());
        //init_vector(*_xmin.at(0), xmin0, _dim);
    }

    virtual ~F_ZDT6() {
    }

    virtual double f1(double x1){
        return 1-exp(-4*x1)*pow (sin(6*M_PI*x1),6);
    };
    virtual double g(double *x){
        double sum;
        sum=0;
        for(int i =1;i<d;i++){
            sum+=x[i];
        }
        sum/=(d-1);
        sum=1+9*pow(sum,0.25);
        return sum;
    };
    virtual double h(double f1, double g){
        
        assert(f1/g>=0);
        return 1-pow(f1/g,2);
    };

};






#endif /* SRC_FUNCTION_H_ */
