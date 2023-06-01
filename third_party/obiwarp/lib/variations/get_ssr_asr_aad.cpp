// STDLIB:
#include <cstdio>
#include "string.h"

// MINE
#include "vec.h"
#include "mat.h"

#define DEBUG (0)

using namespace VEC;

int main (int argc, char **argv) {
    MatF tester;
    MatF tester_trans;
    VecF mpt;
    VecF npt;
    for (int c = 1; c < argc; ++c) {
        tester.set_from_ascii(argv[c], 1);  // no headers on the files
        tester.transpose(tester_trans);
        mpt.set(tester_trans.cols(), tester_trans.pointer(0));
        npt.set(tester_trans.cols(), tester_trans.pointer(1));
        double ssr = VecF::sum_sq_res_yeqx(mpt, npt);
        double asr = VecF::avg_sq_res_yeqx(mpt, npt);
        double aad = VecF::avg_abs_diff(mpt, npt);
        printf("%s %f %f %f\n", argv[c], ssr, asr, aad);
    }
    return 1;
}

