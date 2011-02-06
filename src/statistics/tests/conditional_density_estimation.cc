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
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "ReadFile.h"
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
    --test:           test file for pdf\n\
\n";


int main (int argc, char** argv)
{
    int T = 0;
    int max_depth = 8;
    int max_depth_cond = 8;
    bool joint = false;
    char* filename = NULL;
    char* test_filename = NULL;
    int n_inputs = 1;
    int grid_size = 256;

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
                case 0: T = atoi(optarg); break;
                case 1: max_depth = atoi(optarg); break;
                case 2: max_depth_cond = atoi(optarg); break;
                case 3: joint = true; break;
                case 4: filename = optarg; break;
                case 5: n_inputs = atoi(optarg); break;
                case 6: grid_size = atoi(optarg); assert(grid_size > 0); break;
                case 7: test_filename = optarg; break;
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
    int n_records = ReadFloatDataASCII(data, filename);
    assert(n_records == data.Rows());
    if (T > 0) {
        T = std::min<int>(T, n_records);
    } else {
        T = n_records;
    }
    int data_dimension = data.Columns();

    if (T < 0) {
        Serror("T should be >= 0\n");
        exit(-1);
    }


    if (max_depth <= 0) {
        Serror("max_depth should be >= 0\n");
        exit(-1);
    }


    if (max_depth_cond <= 0) {
        Serror("max_depth_cond should be >= 0\n");
        exit(-1);
    }

    
	Vector lower_bound(data_dimension);
	Vector upper_bound(data_dimension);
	Vector lower_bound_x(n_inputs);
	Vector upper_bound_x(n_inputs);
    int n_outputs = data_dimension - n_inputs;
	Vector lower_bound_y(n_outputs);
	Vector upper_bound_y(n_outputs);
	for (int t=0; t<T; ++t) {
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

    printf ("# T: %d\n", T);
    printf ("# Joint: %d\n", joint);
    printf ("# L: "); lower_bound.print(stdout);
    printf ("# U: "); upper_bound.print(stdout);
    printf ("# LX: "); lower_bound_x.print(stdout);
    printf ("# LY: "); lower_bound_y.print(stdout);
    printf ("# UX: "); upper_bound_x.print(stdout);
    printf ("# UY: "); upper_bound_y.print(stdout);


    ContextTreeKDTree* pdf = NULL;
    ConditionalKDContextTree* cpdf = NULL;
    
    if (joint) {
        pdf = new ContextTreeKDTree (2, max_depth, lower_bound, upper_bound);
    } else {
        cpdf = new ConditionalKDContextTree(2,
                                            max_depth, max_depth_cond,
                                            lower_bound_x, upper_bound_x,
                                            lower_bound_y, upper_bound_y);
    }
	Vector z(data_dimension);
	real log_loss = 0;
    for (int t=0; t<T; ++t) {
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
    

	printf ("%f # AVERAGE LOG LOSS\n", log_loss / (real) T);
    if (joint) {
        Vector v(data_dimension);
        if (data_dimension==1) {
            real min_axis = Min(lower_bound);
            real max_axis = Max(upper_bound);
            real step = (max_axis - min_axis) / (real) grid_size;
            printf ("# MIN AXIS: %f\n", min_axis);
            printf ("# MAX AXIS: %f\n", max_axis);
            printf ("# STEP: %f\n", step);

            for (real z=min_axis; z<max_axis; z+=step) {
                v(0) = z;
                printf ("%f %f # P_XY\n", z, pdf->pdf(v));
            }
        } else {
            Vector step = (upper_bound - lower_bound) / (real) grid_size;

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
        real step = (max_axis - min_axis) / (real) grid_size;

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
    if (pdf) {
        printf ("PDF model\n");
        pdf->Show();
    }
    if (cpdf) {
        printf ("CPDF model\n");
        cpdf->Show();
    }

    if (test_filename) {
        Matrix test_data;
        int n_test = ReadFloatDataASCII(test_data, test_filename);
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
        printf ("%f %f # mismatch\n", Z * mse, Z * abs);
    }
    

    delete cpdf;
    delete pdf;
    return 0;
}

#endif
