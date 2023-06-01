#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include "mat.h"
#include "vec.h"


using namespace VEC;
class AlignTestSuite : public CxxTest::TestSuite 
{
    public:
        // @TODO: either get these smat tests working or move them over to ruby
        //
        // Assures that the same data is representated before and after
        // conversions
        void test_smat_io( void ) {
#ifdef WIN32
            //system("../bin/obiwarp.exe -a -s product -x smat_product.tmp tfiles/tmp1.lmata tfiles/tmp2.lmata");
#else
            //system("../bin/obiwarp -a -s product -x smat_product.tmp tfiles/tmp1.lmata tfiles/tmp2.lmata");
#endif
            TS_ASSERT_EQUALS(0, 0);
//            MatF smat_f;
//            char *tmpfile = "smat_product.tmp";
//
//            smat_f.set_from_binary(tmpfile);
//            TS_ASSERT_DELTA(smat_f(0,0), 1.63204e+14, 1.0e11);
//            TS_ASSERT_DELTA(smat_f(0,1), 1.46614e+14, 1.0e11);
//            TS_ASSERT_DELTA(smat_f(39,45), 2.02115e+14, 1.0e11);
//            TS_ASSERT_EQUALS(smat_f.rows(), 40);
//            TS_ASSERT_EQUALS(smat_f.cols(), 46);
            /*
            if (WIN32) {
                system("obiwarp -a --smat_in smat_product.tmp tfiles/tmp1.lmata tfiles/tmp2.lmata");
            }
            else {
                system("./obiwarp -a --smat_in smat_product.tmp tfiles/tmp1.lmata tfiles/tmp2.lmata");
            }
            */
//            remove(tmpfile);
        }
        
};

