#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include "lmat.h"
#include "mat.h"
#include "vec.h"


using namespace VEC;
class LMatTestSuite : public CxxTest::TestSuite 
{
    public:
        void test_creation( void ) {
            LMat testing;
            TS_ASSERT_EQUALS(testing._mz_vals, 0);
            TS_ASSERT_EQUALS(testing._tm_vals, 0);
        }

        void test_warp_tm( void ) {
            float self_arr[5] = {0,3,5,6,12};
            float other_arr[5] = {1,3,7,9,20};
            VecF selfy(5,self_arr,1);
            VecF other(5,other_arr,1);
            float time_arr[16] = {0,1,2,3,4,4.1,5,6,7,8,9,10,10.001,11,12.1, 13.1};
            float answ_arr[16] = {1, 1.28863, 2.02171, 3, 4.75862, 4.98493, 7, 9, 10.9224, 12.8188, 14.6819, 16.5046, 16.5064, 18.2797, 20.1687, 21.8204};
            VecF answ(16, answ_arr,1);
            LMat obj;
            obj._tm->set(16, time_arr);
            obj.warp_tm(selfy, other); 
            for (int i = 0; i < 16; i++) {
                TS_ASSERT_DELTA(answ[1], obj.tm()->at(1), 0.001);
            }
        }

        void test_ascii_read_write( void ) {
            LMat readin;
            readin.set_from_ascii("tfiles/tmp1.lmata");
            TS_ASSERT_EQUALS(readin.mzlen(), 30);
            TS_ASSERT_EQUALS(readin.tmlen(), 40);
            float _mz[30] = { 400,401,402,403,404,405,406,407,408,409,
                              410,411,412,413,414,415,416,417,418,419,
                              420,421,422,423,424,425,426,427,428,429
            };
            VecF mzv(30,_mz,1);
            TS_ASSERT_EQUALS( mzv, *(readin.mz()) );
            TS_ASSERT_DELTA( (*readin.mat())(0,0), 6139950.06794636, 0.1 );
            TS_ASSERT_DELTA( (*readin.mat())(39,29), 2292810.65100822, 0.1 );
            TS_ASSERT_DELTA((*readin.mat())(7,8), 1397963.17842461, 0.1 );

            // *******************************************
            // TEST COORDS
            // *******************************************
            VecI obj1(4);
            obj1[0] = 1;
            obj1[1] = 3;
            obj1[2] = 4;
            obj1[3] = 8;
            VecF out;
            readin.mz_axis_vals(obj1,out);
            TS_ASSERT_DELTA(out[0], 401, 0.001);
            TS_ASSERT_DELTA(out[1], 403, 0.001);
            TS_ASSERT_DELTA(out[2], 404, 0.001);
            TS_ASSERT_DELTA(out[3], 408, 0.001);
            readin.tm_axis_vals(obj1,out);
            TS_ASSERT_DELTA(out[0], 1212.34, 0.001);
            TS_ASSERT_DELTA(out[1], 1236.34, 0.001);
            TS_ASSERT_DELTA(out[2], 1248.34, 0.001);
            TS_ASSERT_DELTA(out[3], 1296.34, 0.001);
            // *******************************************

            char *tmpfile = (char *)"tmp.tmp.tmp";
            readin.print(tmpfile);
            
            LMat readnew;
            readnew.set_from_ascii(tmpfile);
            TS_ASSERT_EQUALS(readnew.mzlen(), 30);
            TS_ASSERT_EQUALS(readnew.tmlen(), 40);
            TS_ASSERT_EQUALS(readnew.tmlen(), 40);
            TS_ASSERT_EQUALS( mzv, *(readnew.mz()) );
            TS_ASSERT_DELTA( (*readnew.mat())(0,0), 6139950.06794636, 0.1 );
            TS_ASSERT_DELTA( (*readnew.mat())(39,29), 2292810.65100822, 0.1 );
            TS_ASSERT_DELTA((*readnew.mat())(7,8), 1397963.17842461, 0.1 );
            remove(tmpfile); 
                
            // Test printing to stdout
            //readin.set_from_ascii("tfiles/tmp1.lmata");
            //readin.print();
        }

        void test_creation_from_mat( void ) {
            LMat obj;
            obj.set_from_binary_mat("tfiles/file1.mat");
            TS_ASSERT_EQUALS(obj.tmlen(), 4);
            TS_ASSERT_EQUALS(obj.mzlen(), 3);
            TS_ASSERT_EQUALS((*obj.mat())(0,0), 1.0);
            TS_ASSERT_EQUALS((*obj.mat())(3,2), 12.0);
        }

        void test_creation_from_mata( void ) {
            LMat obj;
            obj.set_from_ascii_mat("tfiles/file1.mata");
            TS_ASSERT_EQUALS(obj.tmlen(), 4);
            TS_ASSERT_EQUALS(obj.mzlen(), 3);
            TS_ASSERT_EQUALS((*obj.mat())(0,0), 1.0);
            TS_ASSERT_EQUALS((*obj.mat())(3,2), 12.0);
        }

        void test_binary_read_write( void ) {
            int ch_mz_vals = 30;
            int ch_tm_vals = 40;
            LMat readin;
            readin.set_from_ascii("tfiles/tmp1.lmata");
            TS_ASSERT_EQUALS(readin.mzlen(), ch_mz_vals);
            TS_ASSERT_EQUALS(readin.tmlen(), ch_tm_vals);
            float *mptr = (float*)(*readin.mat());
            TS_ASSERT_DELTA(mptr[0], 6139950.06794636, 0.1 );
            TS_ASSERT_DELTA(mptr[(ch_mz_vals*ch_tm_vals)-1], 2292810.65100822, 0.1 );
            TS_ASSERT_DELTA((*readin.mat())(7,8), 1397963.17842461, 0.1 );

            char *tmpfile = (char *)"tmp2.tmp.tmp";
            readin.write(tmpfile);
            LMat readnew(tmpfile);

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
            TS_ASSERT_EQUALS(readnew.mzlen(), ch_mz_vals);
            TS_ASSERT_EQUALS(readnew.tmlen(), ch_tm_vals);
            TS_ASSERT_EQUALS( mzv, *(readnew.mz()) );
            TS_ASSERT_EQUALS( tmv, *(readnew.tm()) );
            TS_ASSERT_DELTA( (*readnew.mat())(0,0), 6139950.06794636, 0.1 );
            TS_ASSERT_DELTA( (*readnew.mat())(7,8),1397963.17842461, 0.1 );
            TS_ASSERT_DELTA( (*readnew.mat())(39,29), 2292810.65100822, 0.1 );
            remove(tmpfile);

            // Test writing binary file to stdout
            //readin.write();
        }
      

};

