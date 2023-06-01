#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include "mat.h"
#include "vec.h"


using namespace VEC;
class OutliersTestSuite : public CxxTest::TestSuite 
{
    public:
        // Assures that the same data is representated before and after
        // conversions
        void test_outliers( void ) {
            system("./outliers -d 1.2 tfiles/tmptimes.txt");
            //TS_ASSERT_EQUALS(fromascii.mzlen(), ch_mz_vals);
        }


};

