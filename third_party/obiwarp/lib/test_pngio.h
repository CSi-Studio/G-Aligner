#include <cxxtest/TestSuite.h>
#include <stdlib.h>

#include "mat.h"
#include "vec.h"
#include "pngio.h"


using namespace VEC;
class PngIOTestSuite : public CxxTest::TestSuite 
{
    public:
        void test_simple( void ) {
            MatI silly(4,4);
            //silly = 8;
            for (int m = 0; m < silly.rows(); ++m) {
                for (int n = 0; n < silly.cols(); ++n) {
                    silly(m,n) = 0;
                }
            }
            silly(0,2) = 1;
            silly(1,3) = 1;
            silly(2,4) = 1;

            PngIO ioguy(1);
            ioguy.write("trial.png", silly);
            
            remove("trial.png");
        }
};

