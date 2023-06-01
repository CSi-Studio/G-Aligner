// STDLIB:
#include <cstdio>
#include <iostream>
#include <fstream>
#include "string.h"

// MINE
#include "vec.h"
#include "mat.h"
#include "lmat.h"
#include "dynprog.h"
#include "pngio.h"
#include "cmdparser.h"

#define DEBUG (0)

int main (int argc, char **argv) {

    CmdParser opts(argc, argv);

    int NUM_INTERNAL_ANCHORS = 10000;   // Max this out!
    char file1[1024];
    char file2[1024];
    strcpy(file1, opts.infiles[0]);
    strcpy(file2, opts.infiles[1]); 
    char outfilename[1024 + 7];

    int mi_bins = 5;
    //char toError[300];
    
    if (DEBUG) {
        std::cerr << "**********************************************\n";
        std::cerr << "opts.local: " << opts.local << "\n";
        std::cerr << "opts.images: " << opts.images << "\n";
        std::cerr << "opts.score: " << opts.score << "\n";
        std::cerr << "opts.outfile: "<< opts.outfile << "\n";
        std::cerr << "opts.timefile: " << opts.timefile << "\n";
        std::cerr << "file1: " << file1 << "\n";
        std::cerr << "file2: " << file2 << "\n";
        std::cerr << "**********************************************\n";
    }


    
    if (opts.outfile == NULL) {
        strcpy(outfilename, file2);
        strcat(outfilename, ".warped");
    }
    else {
        strcpy(outfilename, opts.outfile);
    }

    // ************************************************************
    // * READ IN FILES TO GET MAT 
    // ************************************************************
    LMat lmat1;
    LMat lmat2;
    MatF mat1;
    MatF mat2;
    MatF smat;
    DynProg dyn;

//    if (opts.axes) {
//        if (opts.binary) {
//            lmat1.set_from_binary(file1);
//            lmat2.set_from_binary(file2);
//            lmat1._mat->copy(mat1, 1);
//            lmat2._mat->copy(mat2, 1);
//        }
//        else {
//            lmat1.set_from_ascii(file1);
//            lmat2.set_from_ascii(file2);
//            lmat1._mat->copy(mat1, 1);
//            lmat2._mat->copy(mat2, 1);
//        }
//    }
//    else {
//        //mat1.set_from_ascii(file1);  @TODO: write this guy
//        //mat2.set_from_ascii(file2); 
//    }
//
//    // ************************************************************
//    // * SCORE THE MATRICES
//    // ************************************************************
//    if (DEBUG) {
//        std::cerr << "Scoring the mats!\n";
//    }
//    if (opts.smat_in != NULL) {
//        smat.set_from_binary(opts.smat_in);
//        dyn._smat = &smat;
//    }
//    else {
//        dyn.score(mat1, mat2, smat, opts.score, mi_bins);
//        // SETTING THE SMAT TO BE std normal
//        smat -= smat.avg();
//        double mean, stdev;
//        smat._dat.sample_stats(mean, stdev);
//        smat /= stdev;
//        if (!strcmp(opts.score,"euclidean")) {
//            smat *= -1; // inverting euclidean
//        }
//    }
//    if (opts.smat_out != NULL) {
//        printf("Writing binary smat to '%s'\n", opts.smat_out);
//        smat.write(opts.smat_out);
//        //smat.print(smat_out_files[0]);
//        exit(0);
//    }
//
//    // ************************************************************
//    // * PREPARE GAP PENALTY ARRAY
//    // ************************************************************
//   
//    MatF tester;
//    MatF tester_trans;
//    VecF mpt;
//    VecF npt;
//    VecF mOut_tm;
//    VecF nOut_tm;
//
//    double average = smat.avg();
//    int gp_length = smat.rows() + smat.cols();
//
//    VecF gp_array;
//    dyn.linear_less_before(opts.gap_extend,opts.gap_init,gp_length,gp_array);
//
//    // ************************************************************
//    // * DYNAMIC PROGRAM
//    // ************************************************************ 
//    int minimize = 0;
//    if (DEBUG) {
//        std::cerr << "Dynamic Time Warping Score Matrix!\n";
//    }
//    dyn.find_path(smat, gp_array, minimize, opts.factor_diag, opts.factor_gap, opts.local, opts.init_penalty);
//
//    VecI mOut;
//    VecI nOut;
//    dyn.warp_map(mOut, nOut, minimize, NUM_INTERNAL_ANCHORS);
//
//    if (opts.timefile != NULL) {
//        tester.set_from_ascii(opts.timefile, 1);  // no headers on the files
//        tester.transpose(tester_trans);
//        mpt.set(tester_trans.cols(), tester_trans.pointer(0));
//        npt.set(tester_trans.cols(), tester_trans.pointer(1));
//        float ssr, asr, sad, aad;
//        dyn.path_accuracy((*lmat1._tm), (*lmat2._tm), mOut, nOut, mpt, npt, ssr, asr, sad, aad);
//        //printf("average residual^2 (sec): %f\n", asr);
//        //printf("average abs time diff (sec): %f\n", aad);
//        printf("%f %f %f %f\n", ssr, asr, sad, aad);
//    }
//
//    // Warp the second lmat run!
//    if (opts.axes) {
//        VecF nOutF;
//        VecF mOutF;
//        lmat1.tm_axis_vals(mOut, mOutF);
//        lmat2.tm_axis_vals(nOut, nOutF); //
//        lmat2.warp_tm(nOutF, mOutF); 
//    }
//    else {
//        // or warp the mat itself!
//        // @TODO: write the warping of mat itself!
//    }
//
//    if (opts.binary) {
//        //lmat2.print(outfilename);
//        lmat2.write(outfilename);
//    }
//    else {
//        lmat2.print(outfilename);
//    }
//
//
//
//    if (opts.images) {
//        PngIO wrt(1);
//        char base_fn[1024];
//        strcpy(base_fn, "obi-warp_");
//        char tb_fn[1024];
//        strcpy(tb_fn, base_fn);
//        strcat(tb_fn, "tb.png");
//        //char *tb_fn = "tb.png";
//        wrt.write(tb_fn, dyn._tb);
//        char tbpath_fn[1024];
//        strcpy(tbpath_fn, base_fn);
//        strcat(tbpath_fn, "tbpath.png");
//        wrt.write(tbpath_fn, dyn._tbpath);
//
//        char asmat_fn[1024];
//        strcpy(asmat_fn, base_fn);
//        strcat(asmat_fn, "asmat.png");
//        //wrt.write(asmat_fn, dyn._asmat);
//
//        //strcpy(base_fn, "tb.png");
//        //char *tbpath_fn = "tbpath.png";
//        //char *tbscores_fn = "tbscores.png";
//        //wrt.write(tbscores_fn, dyn._tbscores);
//        //char *asmat_fn = "asmat.png";
//        //wrt.write(asmat_fn, dyn._asmat);
//        char *smat_fn = "smat.png";
//        //wrt.write(smat_fn, *dyn._smat);
//    }
//
///*
//   char silly[100];
//   strcpy(silly, "png_");
//   char tmpp[5];
//   sprintf(tmpp, "%d", i);
//   strcat(silly, tmpp); 
//   strcat(silly, ".png");
//
//   PngIO wrt(0);
////wrt.write(silly, dyn._tbpath);
//wrt.write(silly, _scorepath);
//*/
//

return 0;
}

