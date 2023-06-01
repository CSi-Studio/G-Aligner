#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include "lmat.h"
#include "mat.h"
#include "vec.h"


using namespace VEC;
class LMatTestSuite : public CxxTest::TestSuite 
{
    public:
        // Assures that the same data is representated before and after
        // conversions
        void test_conversions( void ) {
            // Depends on three file conversions performed elsewhere:
            // sh "./lmata2lmat tfiles/tmp1.lmata"
            // File.copy('tfiles/tmp1.lmat', 'tfiles/tmp1B.lmat')
            // sh "./lmat2lmata tfiles/tmp1B.lmat"

            int ch_mz_vals = 30;
            int ch_tm_vals = 40;

            float _mz[30] = { 400,401,402,403,404,405,406,407,408,409,
                              410,411,412,413,414,415,416,417,418,419,
                              420,421,422,423,424,425,426,427,428,429
            };
            float _tm[40] = { 1200.34, 1212.34, 1224.34, 1236.34, 1248.34, 
                1260.34, 1272.34, 1284.34, 1296.34, 1308.34, 1320.34, 1332.34, 
                1344.34, 1356.34, 1368.34, 1380.34, 1392.34, 1404.34, 1416.34, 
                1428.34, 1440.34, 1452.34, 1464.34, 1476.34, 1488.34, 1500.34, 
                1512.34, 1524.34, 1536.34, 1548.34, 1560.34, 1572.34, 1584.34, 
                1596.34, 1608.34, 1620.34, 1632.34, 1644.34, 1656.34, 1668.34 
            };
            VecF mzv(ch_mz_vals,_mz,1);
            VecF tmv(ch_tm_vals,_tm,1);

            // ************************************************** 
            // Set from ascii
            // ************************************************** 
            LMat fromascii;
            fromascii.set_from_ascii("tfiles/tmp1.lmata");

            // Assert that this guy is like we expect
            TS_ASSERT_EQUALS(fromascii.mzlen(), ch_mz_vals);
            TS_ASSERT_EQUALS(fromascii.tmlen(), ch_tm_vals);
            TS_ASSERT_EQUALS( mzv, *(fromascii.mz()) );
            TS_ASSERT_EQUALS( tmv, *(fromascii.tm()) );
            TS_ASSERT_DELTA( (*fromascii.mat())(0,0), 6139950.06794636, 0.1 );
            TS_ASSERT_DELTA( (*fromascii.mat())(7,8),1397963.17842461, 0.1 );
            TS_ASSERT_DELTA( (*fromascii.mat())(39,29), 2292810.65100822, 0.1 );

            // ************************************************** 
            // Read from binary 
            // ************************************************** 
            LMat readnew("tfiles/tmp1.lmat");

            // Assert that it is identical to 'fromascii'
            TS_ASSERT_EQUALS(fromascii.mzlen(), readnew.mzlen());
            TS_ASSERT_EQUALS(fromascii.tmlen(), readnew.tmlen());
            // Problems in WINDOWS HERE::
            TS_ASSERT_SAME_DATA((float*)(*fromascii.mat()),(float*)(*readnew.mat()),ch_mz_vals*ch_tm_vals);
            TS_ASSERT_SAME_DATA((float*)(*fromascii.mz()),(float*)(*readnew.mz()),ch_mz_vals);
            TS_ASSERT_SAME_DATA((float*)(*fromascii.tm()),(float*)(*readnew.tm()),ch_tm_vals);

            // ************************************************** 
            // read from ascii
            // ************************************************** 
            LMat fromascii2;
            fromascii2.set_from_ascii("tfiles/tmp1B.lmata");
      
            TS_ASSERT_EQUALS(fromascii.mzlen(), fromascii2.mzlen());
            TS_ASSERT_EQUALS(fromascii.tmlen(), fromascii2.tmlen());
//            TS_ASSERT_SAME_DATA((float*)(*fromascii.mat()),(float*)(*fromascii2.mat()),ch_mz_vals*ch_tm_vals);
            TS_ASSERT_SAME_DATA((float*)(*fromascii.mz()),(float*)(*fromascii2.mz()),ch_mz_vals);
            TS_ASSERT_SAME_DATA((float*)(*fromascii.tm()),(float*)(*fromascii2.tm()),ch_tm_vals);

        }


};

