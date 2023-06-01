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
        "usage: mat2mata file.mat ... \n" <<
        "**********************************************************************\n" <<
        "outputs file.mata\n" <<
        "converts binary file into ascii mata file(does NOT delete original)\n" <<
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
        int outfile_strlen = strlen(outfile);
        outfile[outfile_strlen] = 'a';
        outfile[outfile_strlen+1] = '\0';
        //std::cerr << "creating: " << outfile << "\n";
        mat.set_from_binary(file);
        mat.print(outfile);
    }
}

