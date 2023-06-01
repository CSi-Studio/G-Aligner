// STDLIB:
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "string.h"

// MINE
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "dynprog.h"
#include "pngio.h"


char file[1024];

int BINARY = 0;

int main (int argc, char *argv[]) {
    /************************************************************
     * GET ARGUMENTS
     ************************************************************/ 
    if (argc == 1) { 
        std::cerr << 
        "*****************************************************************\n" <<
        "usage: lmat2png [-b] file1 ...\n" <<
        "*****************************************************************\n" <<
        "requires png2matrix callable\n" <<
        "*****************************************************************\n";
        exit(1);
    }
    int i;
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-b")) {
            BINARY = 1;
        }
    }
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-b")) {  // if this is a file (not an arg)
            char outfile[1024];
            strcpy(file, argv[i]);
            strcpy(outfile, file);
            char *ptr;
            ptr = strstr(outfile, ".lmat");  //works for lmat and lmata
            *ptr = '\0';
            strcat(outfile, ".png");
            LMat lmat;
            if (BINARY) {
                lmat.set_from_binary(file);
            }
            else {
                lmat.set_from_ascii(file);
            }
            PngIO wrt(1);
            wrt.write(outfile, *lmat.mat());
        }
    }
}

