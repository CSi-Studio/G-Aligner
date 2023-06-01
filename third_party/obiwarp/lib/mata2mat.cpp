// STDLIB:
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "string.h"

// MINE
#include "mat.h"


using namespace VEC;
int main (int argc, char *argv[]) {

    if (argc == 1) { 
        std::cerr << 
        "**********************************************************************\n" <<
        "usage: mata2mat file.mata ... \n" <<
        "**********************************************************************\n" <<
        "outputs file.mat\n" <<
        "converts ascii mata file into binary file (does NOT delete original)\n" <<
        "\n" <<
        "FORMATS: \n" <<
        "mat = binary (all vals 4 bytes):\n" << 
        "        (int)#rows,(int)#cols,(floats)matrix_data_values...\n" <<
        "mata = ascii format (space delimited, with newlines as shown below):\n" <<
        "         #rows, #cols\n" <<
        "         matrix data row1\n" <<
        "         matrix data row2\n" <<
        "         matrix data row3 ...\n" <<
        "**********************************************************************\n";
        exit(1);
    }
    /************************************************************
     * GET ARGUMENTS
     ************************************************************/ 
    int i;
    char file[1024];
    char outfile[1024];

    MatF mat;
    for (i = 1; i < argc; i++) {
        strcpy(file, argv[i]);
        strcpy(outfile, file); 
        outfile[strlen(outfile)-1] = '\0';
        //std::cerr << "creating: " << outfile << "\n";
        mat.set_from_ascii(file);
        mat.write(outfile);
    }
}

