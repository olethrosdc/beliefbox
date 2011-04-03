/* -*- Mode: C++; -*- */
// copyright (c) 2010-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "ContextTreeKDTree.h"
#include "ConditionalKDContextTree.h"
#include "Random.h"
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "ReadFile.h"
#include "KFold.h"
#include <vector>
#include <cstring>
#include <getopt.h>


static const char* const help_text = "Usage: conditional_density_estimation [options]\n\
\nOptions:\n\
    --T:              number of examples to load\n\
    --max_depth:      maximum depth of the local density estimator tree\n\
    --max_depth_cond: maximum depth of the conditional tree\n\
    --joint:          perform simple density estimation\n\
    --data:           filename\n\
    --n_inputs:       number of columns to condition on\n\
    --grid_size:      grid size for plot\n\
    --pdf_test:       test against the actual pdf at given locations\n\
    --test:           test log loss on additional data\n\
    --kfold:         test on a k-fold rather than a separate set.\n\
\n";



struct LocalOptions{
	int T;
	int max_depth;
	int max_depth_cond;
	bool joint;
	char* filename;
	char* test_filename;
	char* pdf_test_filename;
	int n_inputs;
	int grid_size;
	int kfold;
};

void train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options);

int main (int argc, char** argv)
{
	LocalOptions options;
    options.T = 0;
    options.max_depth = 8;
    options.max_depth_cond = 8;
    options.joint = false;
    options.filename = NULL;
    options.test_filename = NULL;
    options.pdf_test_filename = NULL;
    options.n_inputs = 1;
    options.grid_size = 0;
	options.kfold = 0;

    {
        // options
        int c;
        int digit_optind = 0;
        while (1) {
            int this_option_optind = optind ? optind : 1;
            int option_index = 0;
            static struct option long_options[] = {
                {"T", required_argument, 0, 0}, //0
                {"max_depth", required_argument, 0, 0}, //1
                {"max_depth_cond", required_argument, 0, 0}, //2
                {"joint", no_argument, 0, 0}, //3
                {"data", required_argument, 0, 0}, // 4
                {"n_inputs", required_argument, 0, 0}, // 5
                {"grid_size", required_argument, 0, 0}, // 6
                {"test", required_argument, 0, 0}, // 7
                {"pdf_test", required_argument, 0, 0}, // 8
                {"kfold", required_argument, 0, 0}, // 9
                {0, 0, 0, 0}
            };
            c = getopt_long (argc, argv, "",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch (c) {
            case 0:
#if 0
                printf ("option %s (%d)", long_options[option_index].name, option_index);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
#endif
                switch (option_index) {
                case 0: options.T = atoi(optarg); break;
                case 1: options.max_depth = atoi(optarg); break;
                case 2: options.max_depth_cond = atoi(optarg); break;
                case 3: options.joint = true; break;
                case 4: options.filename = optarg; break;
                case 5: options.n_inputs = atoi(optarg); break;
                case 6: options.grid_size = atoi(optarg); assert(options.grid_size >= 0); break;
                case 7: options.test_filename = optarg; break;
                case 8: options.pdf_test_filename = optarg; break;
                case 9: options.kfold = atoi(optarg); break;
                default:
                    fprintf (stderr, "%s", help_text);
                    exit(0);
                    break;
                }
                break;
            case '0':
            case '1':
            case '2':
                if (digit_optind != 0 && digit_optind != this_option_optind)
                    printf ("digits occur in two different argv-elements.\n");
                digit_optind = this_option_optind;
                printf ("option %c\n", c);
                break;
            default:
                std::cout << help_text;
                exit (-1);
            }
        }
	
        if (optind < argc) {
            printf ("non-option ARGV-elements: ");
            while (optind < argc) {
                printf ("%s ", argv[optind++]);
                
            }
            printf ("\n");
        }
    }

    Matrix data;
    int n_records = ReadFloatDataASCII(data, options.filename, options.T);
    assert(n_records == data.Rows());

    if (options.max_depth <= 0) {
        Serror("max_depth should be >= 0\n");
        exit(-1);
    }


    if (options.max_depth_cond <= 0) {
        Serror("max_depth_cond should be >= 0\n");
        exit(-1);
    }
	
	if (options.test_filename && options.kfold > 0) {
		Matrix test_data;
        int n_test = ReadFloatDataASCII(test_data, options.test_filename);
		if (n_test > 0) {
			train_and_test(data, test_data, options);
		} else {
			fprintf(stderr, "Failed to load test data\n");
			exit(-1);
		}
		KFold k_fold(data, options.kfold);
		for (int i=0; i<options.kfold; ++i) {
			Matrix train_data = k_fold.getTrainFold(i);
			train_and_test(train_data, test_data, options);
		}
	} else if (options.test_filename && !options.kfold) {
		Matrix test_data;
        int n_test = ReadFloatDataASCII(test_data, options.test_filename);
		if (n_test > 0) {
			train_and_test(data, test_data, options);
		} else {
			fprintf(stderr, "Failed to load test data\n");
			exit(-1);
		}
	} else if (options.kfold > 0) {
		KFold k_fold(data, options.kfold);
		for (int i=0; i<options.kfold; ++i) {
			Matrix train_data = k_fold.getTrainFold(i);
			Matrix test_data = k_fold.getTestFold(i);
			train_and_test(train_data, test_data, options);
		}
	}

}

void train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options)
{ 
	int n_inputs = options.n_inputs;
    int data_dimension = data.Columns();
	Vector lower_bound(data_dimension);
	Vector upper_bound(data_dimension);
	Vector lower_bound_x(n_inputs);
	Vector upper_bound_x(n_inputs);
    int n_outputs = data_dimension - n_inputs;
	Vector lower_bound_y(n_outputs);
	Vector upper_bound_y(n_outputs);
	for (int t=0; t<data.Rows(); ++t) {
        for (int i=0; i<data_dimension; ++i) {
            real x = data(t, i);
            //printf ("%f ", x);
            lower_bound(i) = std::min<real>(x, lower_bound(i));
            upper_bound(i) = std::max<real>(x, upper_bound(i));
        }
        //printf ("\n");
    }
    lower_bound -=1;
    upper_bound +=1;
    
    for (int i=0; i<n_inputs; ++i) {
        lower_bound_x(i) = lower_bound(i);
        upper_bound_x(i) = upper_bound(i);
	}
    for (int i=0; i<n_outputs; ++i) {
        lower_bound_y(i) = lower_bound(i + n_inputs);
        upper_bound_y(i) = upper_bound(i + n_inputs);
	}
	
    printf ("# T: %d\n", options.T);
    printf ("# Joint: %d\n", options.joint);
    printf ("# L: "); lower_bound.print(stdout);
    printf ("# U: "); upper_bound.print(stdout);
    printf ("# LX: "); lower_bound_x.print(stdout);
    printf ("# LY: "); lower_bound_y.print(stdout);
    printf ("# UX: "); upper_bound_x.print(stdout);
    printf ("# UY: "); upper_bound_y.print(stdout);


    ContextTreeKDTree* pdf = NULL;
    ConditionalKDContextTree* cpdf = NULL;
    
    if (options.joint) {
        pdf = new ContextTreeKDTree (2, options.max_depth, lower_bound, upper_bound);
    } else {
        cpdf = new ConditionalKDContextTree(2,
                                            options.max_depth, options.max_depth_cond,
                                            lower_bound_x, upper_bound_x,
                                            lower_bound_y, upper_bound_y);
    }
	Vector z(data_dimension);
	real log_loss = 0;
    for (int t=0; t<data.Rows(); ++t) {
        z = data.getRow(t);
#if 0
        for (int i=0; i<z.Size(); ++i) {
            z(i) += urandom()*0.01;
        }
#endif
		real p = 0;
        if (pdf) {
			p = pdf->Observe(z);
        } 
        if (cpdf) {
            Vector x(n_inputs);
            for (int i=0; i<n_inputs; ++i) {
                x(i) = z(i);
            }
            Vector y(n_outputs);
            for (int i=0; i<n_outputs; ++i) {
                y(i) = z(i + n_inputs);
            }
            p = cpdf->Observe(x, y);
        }

        real log_p = log(p);
        if (log_p < -40) {
            log_p = -40;
        }
		log_loss -= log_p;
        //printf ("%f # p_t\n", p);
    }
    

	printf ("%f # AVERAGE LOG LOSS\n", log_loss / (real) data.Rows());
    if (options.grid_size) {
        if (options.joint) {
            Vector v(data_dimension);
            if (data_dimension==1) {
                real min_axis = Min(lower_bound);
                real max_axis = Max(upper_bound);
                real step = (max_axis - min_axis) / (real) options.grid_size;
                printf ("# MIN AXIS: %f\n", min_axis);
                printf ("# MAX AXIS: %f\n", max_axis);
                printf ("# STEP: %f\n", step);

                for (real z=min_axis; z<max_axis; z+=step) {
                    v(0) = z;
                    printf ("%f %f # P_XY\n", z, pdf->pdf(v));
                }
            } else {
                Vector step = (upper_bound - lower_bound) / (real) options.grid_size;

                printf ("# MIN AXIS:"); lower_bound.print(stdout);
                printf ("# MAX AXIS:"); upper_bound.print(stdout);
                printf ("# STEP"); step.print(stdout);

                Vector v = lower_bound;
                bool running = true;
                while (running) {
                    if (data_dimension != 2) {
                        for (int i=0; i<data_dimension; ++i) {
                            printf ("%f ", v(i));
                        }
                        printf ("%f # P_XY\n", pdf->pdf(v));
                    } else {
                        printf ("%f ", pdf->pdf(v));
                    }

                    int i = 0;
                    bool carry = true;
                    bool end_of_line = false;
                    while (carry) {
                        v(i) += step(i);
                        if (v(i) > upper_bound(i)) {
                            v(i) = lower_bound(i);
                            carry = true;
                            end_of_line = true;
                            ++i;
                        } else {
                            carry = false;
                        }
                        if (i == data_dimension) {
                            carry = false;
                            running = false;
                        }
                    }
				
                    if (data_dimension == 2 && end_of_line) {
                        printf ("# P_XY\n");
                    }
                }

            }
        } else {
            real min_axis = Min(lower_bound);
            real max_axis = Max(upper_bound);
            real step = (max_axis - min_axis) / (real) options.grid_size;

            printf ("# MIN AXIS: %f\n", min_axis);
            printf ("# MAX AXIS: %f\n", max_axis);
            printf ("# STEP: %f\n", step);

            for (real y=min_axis; y<max_axis; y+=step) {
                for (real x=min_axis; x<max_axis; x+=step) {
                    Vector X(1);
                    Vector Y(1);
                    X(0) = x;
                    Y(0) = y;
                    printf (" %f ", cpdf->pdf(X, Y));// distribution.pdf(x)*distribution2.pdf(y));
                }
                printf(" # P_Y_X\n");
            }
        }
    } // options.grid_size

    if (pdf) {
        printf ("PDF model\n");
        pdf->Show();
    }
    if (cpdf) {
        printf ("CPDF model\n");
        cpdf->Show();
    }

    // ------------------ tests --------------------- //
    
    // Test against a PDF
    if (test_data.Rows() > 0) {
		int n_test = test_data.Rows();
        real mse = 0;
        real abs = 0;
        if (pdf) {
            Vector x(data_dimension);
            for (int t=0; t<n_test; ++t) {
                Vector z = test_data.getRow(t);
                for (int i=0; i<data_dimension; ++i) {
                    x(i) = z(i);
                }
                real p = pdf->pdf(x);
                real p_test = z(data_dimension);
                mse += (p - p_test)*(p - p_test);
                abs += fabs(p - p_test);
            }
        }
        real Z = 1.0 / (real) n_test;
        printf ("%f %f # mismatch (MSE, L1)\n", Z * mse, Z * abs);
    }
    
    // Test against prediction.
    if (test_data.Rows()) {
        int n_test = test_data.Rows();
        real log_loss = 0;
        for (int t=0; t<n_test; ++t) {
            Vector z = test_data.getRow(t);
            if (pdf) {
                log_loss -= log(pdf->pdf(z));
            }
            
            if (cpdf) {
                Vector x(n_inputs);
                for (int i=0; i<n_inputs; ++i) {
                    x(i) = z(i);
                }
                Vector y(n_outputs);
                for (int i=0; i<n_outputs; ++i) {
                    y(i) = z(i + n_inputs);
                }
                log_loss -= log(cpdf->pdf(x, y));
            }
        }
        printf ("%f # log loss test\n", log_loss / (real) n_test);
    }

    delete cpdf;
    delete pdf;
}

#endif
