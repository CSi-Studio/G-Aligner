// STDLIB:
#include <cstdio>
#include <iostream>
#include "string.h"

// MINE
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "dynprog.h"
#include "pngio.h"


char file1[1024];
char file2[1024];
int mi_bins = 5;
float init_penalty = 0.f;
char toError[300];

int AXES = 0;
int BINARY = 0;
int LOCAL = 0;
int IMAGES = 0;
int LOGYESNO = 0;
char SCORE[10];


int main (int argc, char *argv[]) {
    /************************************************************
     * GET ARGUMENTS
     ************************************************************/ 
    strcpy(SCORE, "covariance");
    int file1_found_already = 0;
    if (argc == 1) { 
        std::cerr << "usage: smat_dist [-a] [-b] [-l] [-g] [-s <scoretype>] file1 file2\n" <<
        "FORMAT: \n" <<
        "Data should be in an m(rows)x n(cols) matrix (space delimited)\n" <<
        "where each line contains one row of data.  Should be same # cols.\n" << 
        "Data will be aligned along the m axis.\n" <<
        "ARGUMENTS (default marked by asterik*): \n" <<
        "b|binary = file is binary [precision?, etc] rather than *ascii\n" <<
        "a|axes = 1st line in file contains x coordinates, 2nd the y\n" << 
        "s|score = scoring function: *covariance, product (dot product)\n" << 
        "          pearsons_r, pearsons_r2, mutual_info\n" <<
        "l|local = local rather than *global alignment\n" <<
        "i|images = creates png images of the alignment process\n" <<
        "g|log = takes the log (base 2) of smat\n" <<
        "[space between argument and value, please.]\n";
        exit(1);
    }
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-a")) {
            AXES = 1;
        }
        else if (!strcmp(argv[i],"-b")) {
            BINARY = 1;
        }
        else if (!strcmp(argv[i],"-g")) {
            LOGYESNO = 1;
        }
        else if (!strcmp(argv[i],"-i")) {
            IMAGES = 1;
        }
        else if (!strcmp(argv[i],"-l")) {
            LOCAL = 1;
        }
        else if (!strcmp(argv[i],"-s")) {
            i++;
            strcpy(SCORE, argv[i]);
        }
        // if it doesn't match an option then it is our file!
        else {
            if (file1_found_already) {
                strcpy(file2, argv[i]);
            }
            else {
                strcpy(file1, argv[i]);
                file1_found_already = 1;
            }
        }
    }
    std::cerr << "**********************************************\n";
    std::cerr << "SCORE: " << SCORE << "\n";
    std::cerr << "file1: " << file1 << "\n";
    std::cerr << "file2: " << file2 << "\n";
    std::cerr << "LOCAL: " << LOCAL << "\n";
    std::cerr << "BINARY: " << BINARY << "\n";
    std::cerr << "IMAGES: " << IMAGES << "\n";
    std::cerr << "AXES: " << AXES << "\n";
    std::cerr << "LOG: " << LOGYESNO << "\n";
    std::cerr << "**********************************************\n";

    /************************************************************
     * READ IN FILES TO GET MAT 
     ************************************************************/ 
    LMat lmat1;
    LMat lmat2;
    MatF mat1;
    MatF mat2;
    MatF smat;
    DynProg dyn;
    
    if (AXES) {
        lmat1.set_from_ascii(file1);
        lmat2.set_from_ascii(file2);
        lmat1._mat->copy(mat1, 1);
        lmat2._mat->copy(mat2, 1);
    }
    else {
        //mat1.set_from_ascii(file1);  @TODO: write this guy
        //mat2.set_from_ascii(file2); 
    }

    /************************************************************
     * SCORE THE MATRICES
     ************************************************************/ 
    std::cerr << "Scoring the mats!\n";
    dyn.score(mat1, mat2, smat, SCORE, mi_bins);

    /************************************************************
     * PREPARE GAP PENALTY ARRAY
     ************************************************************/ 
    VecF gp_array;  // use default for now

    /************************************************************
     * DYNAMIC PROGRAM
     ************************************************************/ 

    int minimize = 0;
    std::cerr << "Dynamic Time Warping Score Matrix!\n";
    dyn.find_path_with_gaps(smat, gp_array, minimize, LOCAL, init_penalty);
    printf("DYNPROG SCORE: %f\n", dyn._bestScore);

    // Run through various distances:
    int reply;
    int steps;
    char steps_st[3];
    char basefilename[255];

    // strip the lmata:
    char *pointer;
    pointer = strstr(file1, ".lmata");
    *pointer = '\0';
    pointer = strstr(file2, ".lmata");
    *pointer = '\0';

    strcpy(basefilename, file1);
    strcat(basefilename, "_");
    strcat(basefilename, file2);
    strcat(basefilename, "_");
    strcat(basefilename, SCORE);
    strcat(basefilename, "_");
    strcat(basefilename, "steps");
    strcat(basefilename, "_");
    char finalfn[255];
    if (LOGYESNO) {
        strcat(basefilename, "logbase2");
        strcat(basefilename, "_");
        smat.logarithm(2);
    }

    for (int steps = 0; steps < 50; steps += steps + 1) {
        MatI tbpathe;
        dyn._tbpath.expand(tbpathe,1,steps,steps,steps,steps,0,0,0,0); 
    
        strcpy(finalfn,basefilename);
        sprintf(steps_st, "%d", steps);
        strcat(finalfn,steps_st);

        VecF result;
        smat.mask_as_vec(1, tbpathe, result);

        VecD _bins;
        VecI _freqs;
        result.hist(100, _bins, _freqs); 
        char *hist_fn = "hist.txt";
        std::ofstream fh(hist_fn);
        // print filename and title:
        fh << finalfn << "\n"; // filename
        fh << finalfn << "\n"; // title
        // print the x and y axis labels:
        fh << "score" << "\n";
        fh << "frequency" << "\n";
        // print the data:
        _bins.print(fh);
        _freqs.print(fh);
        fh.close();
        reply = system("plot_xy.rb hist.txt -b"); 
        if (reply == -1) { puts("Error!"); }
        else { puts("success"); }
    }

    if (IMAGES) {
        PngIO wrt(1);
        //char tb_fn[100];
        //strcpy(tb_fn, "tb.png");
        char *tb_fn = "tb.png";
        wrt.write(tb_fn, dyn._tb);
        char *tbpath_fn = "tbpath.png";
        wrt.write(tbpath_fn, dyn._tbpath);
        char *asmat_fn = "asmat.png";
        wrt.write(asmat_fn, dyn._asmat);
        char *smat_fn = "smat.png";
        wrt.write(smat_fn, *dyn._smat);
    }

    
    /*
    char silly[100];
    strcpy(silly, "png_");
    char tmpp[5];
    sprintf(tmpp, "%d", i);
    strcat(silly, tmpp); 
    strcat(silly, ".png");

    PngIO wrt(0);
    //wrt.write(silly, dyn._tbpath);
    wrt.write(silly, _scorepath);
    */
}

