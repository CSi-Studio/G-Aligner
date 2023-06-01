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
        "usage: lmat2lmata file.lmat ... \n" <<
        "**********************************************************************\n" <<
        "converts binary lmat file into ascii file (does NOT delete original)\n" <<
        "outputs file.lmata\n" <<
        "\n" <<
        "FORMATS: \n" <<
        "lmat = binary (all vals 4 bytes; all on one line):\n" << 
        "        (int) # rows, (floats) m axis values (i.e. time vals),\n" <<
        "        (int) # cols, (floats) n axis values (i.e. m/z vals),\n" <<
        "        (floats) matrix data values row1, row2, row3 ...\n" <<
        "lmata = ascii format (space delimited, with newlines as shown below):\n" <<
        "         # rows\n" <<
        "         m axis values (i.e. time vals)\n" <<
        "         # cols\n" <<
        "         n axis values (i.e. m/z vals)\n" <<
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

    LMat lmat;
    for (i = 1; i < argc; i++) {
        strcpy(file, argv[i]);
        strcpy(outfile, file); 
        int outfile_strlen = strlen(outfile);
        outfile[outfile_strlen] = 'a';
        outfile[outfile_strlen+1] = '\0';
        //std::cerr << "creating: " << outfile << "\n";
        lmat.set_from_binary(file);
        lmat.print(outfile);
    }
}

