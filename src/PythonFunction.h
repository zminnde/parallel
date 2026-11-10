#ifndef SRC_PYTHON_FUNCTION_H_
#define SRC_PYTHON_FUNCTION_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <sys/wait.h>
#include <unistd.h>

#include "Function.h"

// Override without editing: mpicxx -DPYTHON_INTERPRETER='"..."' -DPYTHON_WORKER_SCRIPT='"..."'
#ifndef PYTHON_INTERPRETER
#define PYTHON_INTERPRETER "/usr/bin/python3"
#endif
#ifndef PYTHON_WORKER_SCRIPT
#define PYTHON_WORKER_SCRIPT "/opt/multiobjective-optimisation/LondonAirQuality/2026-11-10/optimized/worker.py"
#endif

class PythonFunction: public Function {

private:
	int pipe_cpp_to_python[2];
	int pipe_python_to_cpp[2];

	FILE* write_to_python;
	FILE* read_from_python;

	pid_t python_pid;

	// one child process + one pipe pair per instance: copying would
	// hand two objects the same pid and the same FILE*s
	PythonFunction(const PythonFunction&) = delete;
	PythonFunction& operator=(const PythonFunction&) = delete;

protected:

	void init_python_function_call() {
		if (pipe(pipe_cpp_to_python) == -1 || pipe(pipe_python_to_cpp) == -1) {
			std::cerr << "Failed to create pipes.\n";
			exit(1);
		}

		// drain our own buffers first: anything still pending would be
		// duplicated in the child, and the child's stdout is the protocol pipe
		std::cout.flush();
		std::cerr.flush();
		fflush(nullptr);

		python_pid = fork();

		if (python_pid == -1) {
			std::cerr << "Failed to fork process.\n";
			exit(1);
		}

		if (python_pid == 0) {

			dup2(pipe_cpp_to_python[0], STDIN_FILENO);
			dup2(pipe_python_to_cpp[1], STDOUT_FILENO);

			close(pipe_cpp_to_python[0]); close(pipe_cpp_to_python[1]);
			close(pipe_python_to_cpp[0]); close(pipe_python_to_cpp[1]);

			execlp(PYTHON_INTERPRETER,
			       PYTHON_INTERPRETER,
			       "-u",
			       PYTHON_WORKER_SCRIPT,
			       nullptr);

			// exec failed; _exit so no inherited stdio buffer is flushed into the pipe
			std::cerr << "Failed to launch Python interpreter.\n";
			_exit(127);

		} else {

			close(pipe_cpp_to_python[0]);
			close(pipe_python_to_cpp[1]);

			write_to_python = fdopen(pipe_cpp_to_python[1], "w");
			read_from_python = fdopen(pipe_python_to_cpp[0], "r");

			if (write_to_python == nullptr || read_from_python == nullptr) {
				std::cerr << "Failed to open pipe streams.\n";
				exit(1);
			}

			return;
		}
	}

	// signature must match Function::evaluate_func_values(double*) exactly
	double* evaluate_func_values(double* decision_point) override {

		double * func_values= new double[m];

		// one line per evaluation: the d decision values separated by spaces
		for (int i = 0; i < d; i++) {
			fprintf(write_to_python, i == 0 ? "%.6f" : " %.6f", decision_point[i]);
		}
		fprintf(write_to_python, "\n");
		fflush(write_to_python);

		double f1, f2;
		if (fscanf(read_from_python, "%lf %lf", &f1, &f2) == 2) {

			func_values[0]=f1;
			func_values[1]=f2;

		} else {
			// no per-evaluation output; this path is fatal anyway
			std::cerr<<"ERROR while getting values from Python function"<<std::endl;
			exit(1);
		}

		return func_values;
	}

public:

	PythonFunction()  {
		d=6;//  dimension of decision space: 5 distance weights + threshold S
		m=2;//  number of objective functions: mean NMAE12 and refused prediction count

		// S upper bound 14: with every weight at 1 the largest combined distance
		// in the data is 13.81, so above it every candidate is always a neighbour.
		double upper_bound[]={1, 1, 1, 1, 1, 14}; // decision  space bounds
		double lower_bound[]={0, 0, 0, 0, 0, 0};   // decision  space bounds

		intailize(d, upper_bound, lower_bound);

		init_python_function_call();
	}

	~PythonFunction() {
		// closing the write end sends EOF, which ends the worker's read loop
		fclose(write_to_python);
		fclose(read_from_python);

		waitpid(python_pid, nullptr, 0);
	}

	virtual string to_string(){
		static string name;
		name = "PythonFunction";
		return name;
	}
};

#endif /* SRC_PYTHON_FUNCTION_H_ */
