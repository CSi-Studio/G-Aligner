// STDLIB:
#include <cstdio>
#include <iostream>
#include "string.h"
#include "math.h"

// 3RD PARTY 
#include <argtable2.h>

// MINE
#include "vec.h"
#include "mat.h"

#define DEBUG (0)

using namespace VEC;

int mymain(double *deviations, int deviations_cnt, const char **infiles, int infile_cnt );

int main (int argc, char **argv) {
    struct arg_lit  *help = arg_lit0("h", "help", "prints this help and exits");
    struct arg_dbl *deviations = arg_dbl0("d", "dev", "<deviations>", "deviations cutoff (default 4.0)");
    struct arg_file *infiles = arg_filen(NULL, NULL, NULL,1,1024, "files to align (first is template)");
    struct arg_end  *end     = arg_end(20);
    void* argtable[] = {help, deviations, infiles, end};
    const char* progname = "outliers";
    int nerrors;
    int exitcode=0;

    /* verify the argtable[] entries were allocated sucessfully */
    if (arg_nullcheck(argtable) != 0)
    {
        /* NULL entries were detected, some allocations must have failed */
        printf("%s: insufficient memory\n",progname);
        exitcode=1;
        goto exit;
    }

    /* set any command line default values prior to parsing */

    /* Parse the command line as defined by argtable[] */
    nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0 || argc < 2)
    {
        printf("*************************************************************************\n");
        printf("Usage: %s ", progname);
        arg_print_syntax(stdout,argtable,"\n");
        printf("\n");
        printf("tosses out outliers from regression line beyond a certain deviation\n");
        printf("\n");
        arg_print_glossary(stdout,argtable,"  %-26s %s\n");
        printf("\n");
        printf("*************************************************************************\n");

        exitcode=0;
        goto exit;
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
    {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        exitcode=1;
        goto exit;
    }

    /* normal case: take the command line options at face value */
    exitcode = mymain(deviations->dval, deviations->count,  infiles->filename, infiles->count);

    exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return exitcode;
}

int mymain(double *deviations, int deviations_cnt, const char **infiles, int infile_cnt ) {
    int i;
  
    for (int cnt = 0; cnt < infile_cnt; ++cnt) {
        double dev;
        if (deviations_cnt == 1) {
            dev = deviations[0];  // if they give one deviation, use it for all
        }
        else if (deviations_cnt == 0) {
            // Set deviations default:
            dev = 4.0;
        }
        else {
            dev = deviations[cnt]; // multiple deviatsion, one for each file
        }
        char file[1024];
        strcpy(file, infiles[cnt]);

        if (DEBUG) {
            std::cerr << "**********************************************\n";
            std::cerr << "file: " << file << "\n";
            std::cerr << "**********************************************\n";
        }

        MatD mat;
        mat.set_from_ascii(file, 1);
        VecD vecs[2];
        int cnt;
        MatD as_rows;
        mat.transpose(as_rows);
        as_rows.row_vecs(cnt, vecs);
        //vecs[0].print();
        //vecs[1].print();

        double rsq, slope, y_intercept;
        VecD::rsq_slope_intercept(vecs[0], vecs[1], rsq, slope, y_intercept);
        //printf("rsq %f slope %f intercept %f\n", rsq, slope, y_intercept);
        //mat.print();
        // mx + b
        // y = slope(x) + intercept

        // Get the differences from the regression line
        // expected_y = slope(x) + intercept
        // ydiff = actual_y - expected_y
        // run = ydiff/slope
        ////////////////// ydiff = abs(actual_y - expected_y)
        // ydiff / run = tan a
        // sin a = run / x
        // final = run/( sin(arctan(ydiff/run)) )
        /////////////////// if (actual_y - expected_y) < 0 the diff should be (-)
        VecD residuals(vecs[0].length());
        for (i = 0; i < vecs[0].length(); ++i) {
            double expected_y = (slope*vecs[0][i]) + y_intercept;
            double ydiff = vecs[1][i] - expected_y;
            double run = ydiff/slope;
            residuals[i] = run/( sin(atan(ydiff/run)) );
        }
        //puts("RESIDUALS: ");
        //residuals.print();
        //puts("END RESIDUALS: ");

        // get the mean and standard deviation
        double mean, stddev;
        residuals.sample_stats(mean, stddev);
        //printf("m: %f std: %f\n", mean, stddev);

        // for each difference calculate standard deviations
        MatD acceptable_tmp(vecs[0].length(), 2);

        int num_accept = 0;
        int not_accept = 0;
        for (i = 0; i < residuals.length(); ++i) {
            // #stddevsaway = abs(residuals[i] - mean)/stddev );
            double point_devs = (residuals[i] - mean)/stddev;
            if (point_devs < 0.0) { point_devs = -1.0*point_devs; } // abs val

            if (point_devs <= dev) {  // acceptable
                //printf("acceptable dev: %f\n", point_devs);
                acceptable_tmp(num_accept,0) = vecs[0][i];
                acceptable_tmp(num_accept,1) = vecs[1][i];
                ++num_accept;
            }
            else {  // not acceptable! toss out
                //printf("NOT ACCEPTABLE: %f, %f\n", vecs[0][i], vecs[1][i]);
                ++not_accept;
            }
        }

        printf("TOSSED %d points > %.1f deviations from regression line (of %d total) reading file: %s\n", not_accept, dev, residuals.length(), file);
        MatD accept(num_accept,2,(double*)acceptable_tmp,1);


        // print to file named file+.4.0out
        char devs_str[10];
        sprintf(devs_str, "%.1f",dev);
        strcat(file, ".");
        strcat(file, devs_str); 
        strcat(file, "out"); 
        //puts("Accept");
        //accept.print();
        accept.print(file,1);
    }
    return 0;
}

