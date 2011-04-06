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
#include "KernelDensityEstimator.h"
#include "KernelConditionalDensityEstimator.h"
#include "DoubleKernelCDE.h"
#include "ContextTreeKDTree.h"
#include "ConditionalKDContextTree.h"
#include "ReadFile.h"
#include "KFold.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "Random.h"
#include "EasyClock.h"
#include "MersenneTwister.h"
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
    --kernel:       use kernel estimator\n\
    --double_kernel:       use double kernel estimator\n\
    --bandwidth:      bandwidth for kernel estimator\n\
    --tune_bandwidth: tune the bandwidth using a random hold out sample\n\
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
	bool kernel;
	real bandwidth;
	bool tune_bandwidth;
	bool double_kernel;
    real noise_level;
};

void train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options);
void tree_train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options);
void kernel_train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options);

void Noisify(Matrix& data, real noise_level)
{
    if (noise_level == 0) {
        printf("# Adding no noise\n");
        return;
    } 
    printf("# Adding uniform noise %f\n", noise_level);

    MersenneTwisterRNG rng;
    rng.manualSeed(1849435425); // to ensure consistent noise is added
    int T = data.Rows();
    int n = data.Columns();
    for (int i=0; i<T; ++i) {
        for (int j=0; j<n; ++j) {
            data(i,j) += noise_level*(rng.uniform() - 0.5);
        }
    }
}
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
	options.kernel = false;
	options.double_kernel = false;
	options.bandwidth = 1.0;
	options.tune_bandwidth = false;
    options.noise_level = 0;
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
                {"kernel", no_argument, 0, 0}, // 10
                {"bandwidth", required_argument, 0, 0}, // 11
                {"tune_bandwidth", no_argument, 0, 0}, // 12
                {"double_kernel", no_argument, 0, 0}, // 13
                {"noise", required_argument, 0, 0}, // 14
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
                case 10: options.kernel = true; break;
                case 11: options.bandwidth = atof(optarg); break;
                case 12: options.tune_bandwidth = true; break;
                case 13: options.kernel = true; options.double_kernel = true; printf("double\n"); break;
                case 14: options.noise_level = atof(optarg); break;
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
	printf("# OPTIONS\n");
    printf("#T: %d \n", options.T);
    printf("# depth: %d \n", options.max_depth);
    printf("# c-depth: %d \n", options.max_depth_cond);
    printf("# joint: %d\n", options.joint);
    printf("# filename: %s\n", options.filename);
    printf("# test: %s \n", options.test_filename);
    printf("# pdf test: %s \n", options.pdf_test_filename);
    printf("# inputs: %d \n", options.n_inputs);
    printf("# grid: %d \n", options.grid_size);
	printf("# folds: %d \n", options.kfold);
	printf("# kernel: %d \n", options.kernel);
	printf("# double: %d \n", options.double_kernel);
	printf("# bw: %f \n", options.bandwidth);
	printf("# tune: %d \n", options.tune_bandwidth);



    if (options.max_depth <= 0) {
        Serror("max_depth should be >= 0\n");
        exit(-1);
    }


    if (options.max_depth_cond <= 0) {
        Serror("max_depth_cond should be >= 0\n");
        exit(-1);
    }
	
	if (options.test_filename && options.kfold > 0) {
		Matrix data;
		int n_records = ReadFloatDataASCII(data, options.filename, options.T);
		if (n_records <= 0) {
			Serror("Failed to read train data\n");
			exit(-1);
		}
		Matrix test_data;
        int n_test = ReadFloatDataASCII(test_data, options.test_filename);
		if (n_test <= 0) {
			Serror("Failed to read test data\n");
			exit(-1);
		}
		KFold k_fold(data, options.kfold);
		for (int i=0; i<options.kfold; ++i) {
			printf("# Fold %d\n", i);
			Matrix train_data = k_fold.getTrainFold(i);
            Noisify(train_data, options.noise_level);
			train_and_test(train_data, test_data, options);
		}
	} else if (options.test_filename && !options.kfold) {
		Matrix data;
		int n_records = ReadFloatDataASCII(data, options.filename, options.T);
        Noisify(data, options.noise_level);
		if (n_records <= 0) {
			Serror("Failed to read train data\n");
			exit(-1);
		}
		Matrix test_data;
        int n_test = ReadFloatDataASCII(test_data, options.test_filename);
		if (n_test > 0) {
			train_and_test(data, test_data, options);
		} else {
			fprintf(stderr, "Failed to load test data\n");
			exit(-1);
		}
	} else if (options.kfold > 0) {
		Matrix data;
		int n_records = ReadFloatDataASCII(data, options.filename);
		if (n_records <= 0) {
			Serror("Failed to read train data\n");
			exit(-1);
		}
		KFold k_fold(data, options.kfold);
		for (int i=0; i<options.kfold; ++i) {
			printf("# Fold %d\n", i);
			Matrix train_data = k_fold.getTrainFold(i, options.T);
            Noisify(train_data, options.noise_level);
			Matrix test_data = k_fold.getTestFold(i);
			train_and_test(train_data, test_data, options);
		}
	}

}


void train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options)
{
	if (options.kernel) {
		kernel_train_and_test(data, test_data, options);
	} else {
		tree_train_and_test(data, test_data, options);
	}
}

void tree_train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options)
{ 
	printf ("# TREE\n");
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
	
    printf ("# T: %d\n", data.Rows());
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
		fflush(stdout);
		fflush(stderr);
    }

    delete cpdf;
    delete pdf;
}


void kernel_train_and_test(Matrix& data, Matrix& test_data, LocalOptions& options)
{ 
	printf ("# KERNEL\n");
	int n_inputs = options.n_inputs;
	real bandwidth = options.bandwidth;
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
	
    printf ("# T: %d\n", data.Rows());
    printf ("# Joint: %d\n", options.joint);
    printf ("# L: "); lower_bound.print(stdout);
    printf ("# U: "); upper_bound.print(stdout);
    printf ("# LX: "); lower_bound_x.print(stdout);
    printf ("# LY: "); lower_bound_y.print(stdout);
    printf ("# UX: "); upper_bound_x.print(stdout);
    printf ("# UY: "); upper_bound_y.print(stdout);


    KernelDensityEstimator* pdf = NULL;
    KernelConditionalDensityEstimator* cpdf = NULL;
    DoubleKernelCDE* dcpdf = NULL;
    int knn = 0;
    if (options.joint) {
 			printf ("# JOINT KERNEL\n");
        pdf = new KernelDensityEstimator(lower_bound.Size(), bandwidth, knn);
    } else {
		if (options.double_kernel) {
 			printf ("# DOUBLE_KERNEL\n");
			dcpdf = new DoubleKernelCDE(lower_bound_x.Size(),
									   lower_bound_y.Size(),
									   bandwidth);
		} else {
 			printf ("# SINGLE_KERNEL\n");
			cpdf = new KernelConditionalDensityEstimator(lower_bound_x.Size(),
														 lower_bound_y.Size(),
														 bandwidth,
														 knn);
		}
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
        if (cpdf || dcpdf) {
            Vector x(n_inputs);
            for (int i=0; i<n_inputs; ++i) {
                x(i) = z(i);
            }
            Vector y(n_outputs);
            for (int i=0; i<n_outputs; ++i) {
                y(i) = z(i + n_inputs);
            }
			if (cpdf) {
				p = cpdf->Observe(x, y);
			}
			if (dcpdf) {
				p = dcpdf->Observe(x, y);
			}
        }

        real log_p = log(p);
        if (log_p < -40) {
            log_p = -40;
        }
		log_loss -= log_p;
        //printf ("%f # p_t\n", p);
    }
    

	printf ("%f # AVERAGE LOG LOSS\n", log_loss / (real) data.Rows());

	if (options.tune_bandwidth) {
        if (pdf) {
            pdf->BootstrapBandwidth();
        }
        if (cpdf) {
            cpdf->BootstrapBandwidth();
        }
        if (dcpdf) {
            dcpdf->BootstrapBandwidth();
        }
	}
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
					if (cpdf) {
						printf (" %f ", cpdf->pdf(X, Y));// distribution.pdf(x)*distribution2.pdf(y));
					} 
					if (dcpdf) {
						printf (" %f ", dcpdf->pdf(X, Y));
					} 
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
    if (dcpdf) {
        printf ("DCPDF model\n");
        dcpdf->Show();
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
            
            if (cpdf || dcpdf) {
                Vector x(n_inputs);
                for (int i=0; i<n_inputs; ++i) {
                    x(i) = z(i);
                }
                Vector y(n_outputs);
                for (int i=0; i<n_outputs; ++i) {
                    y(i) = z(i + n_inputs);
                }
                if (cpdf) {
					log_loss -= log(cpdf->pdf(x, y));
				}
				if (dcpdf) {
					log_loss -= log(dcpdf->pdf(x, y));
				}
            }
        }
        printf ("%f # log loss test\n", log_loss / (real) n_test);
    }

    delete cpdf;
    delete pdf;
}

#endif
