// STDLIB:
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "string.h"

// MINE
#include "lmat.h"


using namespace VEC;
int main (int argc, char *argv[]) {

    if (argc == 1) { 
        std::cerr << 
        "**********************************************************************\n" <<
        "usage: lmat2chrms file.lmat m/z ...\n" <<
        "**********************************************************************\n";
        exit(1);
    }
    /************************************************************
     * GET ARGUMENTS
     ************************************************************/ 
    int i;
    char file[1024];

    LMat lmat;
    strcpy(file, argv[1]);
    lmat.set_from_binary(file);
    MatF trans;
    lmat.mat()->transpose(trans);
    //printf("rows %d cols %d\n", lmat.mat()->rows(), lmat.mat()->cols());
    //printf("rows %d cols %d\n", trans.rows(), trans.cols());
    VecF *vecs = new VecF[trans.rows()];
    int num_vecs;
    trans.row_vecs(num_vecs, vecs);

    char fn[1024];
    char toplotfn[1024];

    strcpy(fn, file);
    char *pch;
    pch = strstr(fn, ".lmat");
    *pch = '\0';

    char fnbase[1024];
    strcpy(fnbase, fn);
    char *start = strrchr(fnbase, '/');
    if (start != NULL) {
        strcpy(fnbase, ++start);
    }

    strcpy(toplotfn, fn);
    strcat(toplotfn, ".toplot");
    std::ofstream fh(toplotfn);
    printf("WRITING TO: %s\n", toplotfn);

    fh << "XYData" << "\n";
    fh << fnbase << "\n";
    fh << fnbase << " chromatograms" << "\n";
    fh << "time (sec)" << "\n";
    fh << "ion counts" << "\n";

    for (int i = 2; i < argc; ++i) {
        fh << "m/z " << argv[i] << "\n";
        //std::cout << "m/z " << argv[i] << "\n";
        int ind = lmat.mz()->index(atof(argv[i]));
        //printf("IND: %d\n", ind);
        lmat.tm()->print(fh, 1);
        vecs[ind].print(fh, 1);
    }
        
    delete[] vecs;
}

