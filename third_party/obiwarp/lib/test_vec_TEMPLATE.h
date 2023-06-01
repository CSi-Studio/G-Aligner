#include <cxxtest/TestSuite.h>
#include "vec.h"
#include <cstdlib>
#include <fstream>
//#include "vec_utils.h"
#include "tnt_stopwatch.h"

using namespace VEC;


class VecTestSuite : public CxxTest::TestSuite 
{
    public:

        // answers derived from chim on PDL::Slatec
        void test_VecF_linear_interpolation( void ) {
            float seq[10] = {  0,1, 2, 3,4,5,6,7,8, 9 };
            VecF x(10, seq, 1);
            float seq2[10] = { 0,10,12,4,5,2,7,9,10,4 };
            VecF y(10, seq2, 1);
            float seq3[5] = { 1.0f, 2.5f, 4.0f, 6.2f, 9.0f };
            VecF new_x(5,seq3,1);
            VecF out_y;
            // derivs           ?,-8,?,2?
            float ans_arr[5] = {10.0f,8.0f,5.0f,7.4f,4.0f};
            VecF ans(5,ans_arr,1);
            VecF::linear_interp(x,y,new_x,out_y,1);
            //TS_ASSERT_EQUALS(ans, out_y);
            TS_ASSERT_EQUALS(ans[0], out_y[0]);
            TS_ASSERT_EQUALS(ans[1], out_y[1]);
            TS_ASSERT_EQUALS(ans[2], out_y[2]);
            TS_ASSERT_DELTA(ans[3], out_y[3], 0.0001);
            TS_ASSERT_EQUALS(ans[4], out_y[4]);

            VecF out_y2;
            float seq4[5] = { 6.2, 9, 1, 2.5, 4 };
            VecF new_x2(5,seq4,1);
            float ans_arr2[5] = {7.4,4,10,8,5};
            VecF ans2(5,ans_arr2,1);
            VecF::linear_interp(x,y,new_x2,out_y2);
            TS_ASSERT_DELTA(ans2[0], out_y2[0], 0.0001);
            TS_ASSERT_EQUALS(ans2[1], out_y2[1]);
            TS_ASSERT_EQUALS(ans2[2], out_y2[2]);
            TS_ASSERT_EQUALS(ans2[3], out_y2[3]);
            TS_ASSERT_EQUALS(ans2[4], out_y2[4]);
        }
 

        void test_VecF_and_VecD_abs_val( void ) {
            VecF obj1(4,-2.5f);
            float arr1[4] = {2.5f,2.5f,2.5f,2.5f};
            VecF answ(4,arr1,1);
            TS_ASSERT_DIFFERS(obj1, answ);
            obj1.abs_val();
            TS_ASSERT_EQUALS(obj1, answ);

            VecD obj2(4,-2.5);
            obj2[0] = 2.2;
            obj2[3] = 4.2;
            double arr2[4] = {2.2,2.5,2.5,4.2};
            VecD answ2(4,arr2,1);
            TS_ASSERT_DIFFERS(obj2, answ2);
            obj2.abs_val();
            TS_ASSERT_EQUALS(obj2, answ2);
        }

        void test_to_f( void ) {
            VecI test1(3,2);
            VecF out;
            test1.to_f(out);
            out[0] = 0.1f;
            out[1] = 3.6f;
            out[2] += 0.5f;
            VecF answ(3);
            answ[0] = 0.1f;
            answ[1] = 3.6f;
            answ[2] = 2.5f;
            TS_ASSERT_EQUALS(answ, out);
        }

        // Tests for specific Templates:
        void test_VecF_chfe( void ) {
            float seq[10] = { 0,1,2,3,4,5,6,7,8,9 };
            VecF x(10, seq, 1);
            float seq2[10] = { 0,10,12,4,5,2,7,9,10,4 };
            VecF y(10, seq2, 1);

            float seq3[7] = {-4, 1, 2.5, 4, 6.5, 9, 13 };
            VecF new_x(7,seq3,1);
            VecF out_y;
            //x.print(); y.print(); new_x.print();
            VecF::chfe(x,y,new_x,out_y,1);
            //out_y.print();
            TS_ASSERT_DELTA(out_y[0],93.3333333, 0.0001); 
            TS_ASSERT_EQUALS(out_y[1],10); 
            TS_ASSERT_EQUALS(out_y[2],8); 
            TS_ASSERT_EQUALS(out_y[3],5); 
            TS_ASSERT_DELTA(out_y[4],8.1904762, 0.0001); 
            TS_ASSERT_EQUALS(out_y[5],4); 
            TS_ASSERT_EQUALS(out_y[6],110); 
        }

        void test_VecD_chfe( void ) {
            double seq[10] = { 0,1,2,3,4,5,6,7,8,9 };
            VecD x(10, seq, 1);
            double seq2[10] = { 0,10,12,4,5,2,7,9,10,4 };
            VecD y(10, seq2, 1);

            double seq3[7] = {-4, 1, 2.5, 4, 6.5, 9, 13 };
            VecD new_x(7,seq3,1);
            VecD out_y;
            //x.print(); y.print(); new_x.print();
            VecD::chfe(x,y,new_x,out_y);
            //out_y.print();
            TS_ASSERT_DELTA(out_y[0],93.3333333, 0.000001); 
            TS_ASSERT_EQUALS(out_y[1],10); 
            TS_ASSERT_EQUALS(out_y[2],8); 
            TS_ASSERT_EQUALS(out_y[3],5); 
            TS_ASSERT_DELTA(out_y[4],8.1904762, 0.000001); 
            TS_ASSERT_EQUALS(out_y[5],4); 
            TS_ASSERT_EQUALS(out_y[6],110); 
        }



        // answers derived from chim on PDL::Slatec
        void test_VecF_interpolation( void ) {
            float seq[10] = { 0,1,2,3,4,5,6,7,8,9 };
            VecF x(10, seq, 1);
            float seq2[10] = { 0,10,12,4,5,2,7,9,10,4 };
            VecF y(10, seq2, 1);
            VecF d;
            VecF::chim(x,y,d);
            //derivs.print();
            TS_ASSERT_EQUALS(d[0],14); 
            TS_ASSERT_DELTA(d[1], 3.3333333, 0.0001); 
            TS_ASSERT_EQUALS(d[2],0); 
            TS_ASSERT_EQUALS(d[3],0); 
            TS_ASSERT_EQUALS(d[4],0); 
            TS_ASSERT_EQUALS(d[5],0); 
            TS_ASSERT_DELTA(d[6], 2.8571429, 0.0001); 
            TS_ASSERT_DELTA(d[7], 1.3333333, 0.0001); 
            TS_ASSERT_EQUALS(d[8],0); 
            TS_ASSERT_EQUALS(d[9],-9.5); 

            float seq3[5] = { 1, 2.5, 4, 6.5, 9 };
            VecF new_x(5,seq3,1);
            VecF out_y;
            VecF::chfe(x,y,new_x,out_y);
            TS_ASSERT_EQUALS(out_y[0],10); 
            TS_ASSERT_EQUALS(out_y[1],8); 
            TS_ASSERT_EQUALS(out_y[2],5); 
            TS_ASSERT_DELTA(out_y[3],8.1904762, 0.0001); 
            TS_ASSERT_EQUALS(out_y[4],4); 

            VecF out_y_alloced(5);
             VecF::chfe(x,y,new_x,out_y_alloced);
            TS_ASSERT_EQUALS(out_y_alloced[0],10); 
            TS_ASSERT_EQUALS(out_y_alloced[1],8); 
            TS_ASSERT_EQUALS(out_y_alloced[2],5); 
            TS_ASSERT_DELTA(out_y_alloced[3],8.1904762, 0.0001); 
            TS_ASSERT_EQUALS(out_y_alloced[4],4); 

            // The test set above is not adequate for this function:
            VecF::chfe_xy(x,y,new_x,out_y);
            TS_ASSERT_EQUALS(out_y[0],10); 
            TS_ASSERT_EQUALS(out_y[1],8); 
            TS_ASSERT_EQUALS(out_y[2],5); 
            TS_ASSERT_DELTA(out_y[3],8.1904762, 0.01); 
            TS_ASSERT_EQUALS(out_y[4],4); 

        }
        
        void Xtest_VecF_interpolation_speed( void ) {
            // A timed test to see how much our method for interpolation
            // speeds things up:
            // THE ANSWER: A LOT!
            puts("");
            int sorted = 1;
            printf("SORTED: %d\n", sorted);
            for (int factor = 10; factor < 500000; factor *= 2) {
                int i;
                int length = factor;
                float *x_arr = new float[length];
                float *y_arr = new float[length];
                float *newx_arr = new float[length];
                float float_i;
                for (i = 0; i < length; ++i) {
                    float_i = (float)i;
                    newx_arr[i] = float_i+0.5;
                    x_arr[i] = float_i;
                    y_arr[i] = float_i*2;
                }
                VecF in_x(length,x_arr);
                VecF in_y(length,y_arr);
                VecF new_x3(length,newx_arr);
                VecF new_out_y;

                TNT::Stopwatch st;
                st.start();
                for (i = 0; i < 1; ++i) {
                    VecF::chfe(in_x, in_y, new_x3, new_out_y, sorted);
                }
                st.stop();
                printf("%d %f\n", factor, st.read());
                //new_out_y.print();
          }


        }  // End test interpolation
        

        void test_VecD_interpolation( void ) {
            double seq[10] = { 0,1,2,3,4,5,6,7,8,9 };
            VecD x(10, seq, 1);
            double seq2[10] = { 0,10,12,4,5,2,7,9,10,4 };
            VecD y(10, seq2, 1);
            VecD d;
            VecD::chim(x,y,d);
            //derivs.print();
            TS_ASSERT_EQUALS(d[0],14); 
            TS_ASSERT_DELTA(d[1], 3.3333333, 0.0001); 
            TS_ASSERT_EQUALS(d[2],0); 
            TS_ASSERT_EQUALS(d[3],0); 
            TS_ASSERT_EQUALS(d[4],0); 
            TS_ASSERT_EQUALS(d[5],0); 
            TS_ASSERT_DELTA(d[6], 2.8571429, 0.0001); 
            TS_ASSERT_DELTA(d[7], 1.3333333, 0.0001); 
            TS_ASSERT_EQUALS(d[8],0); 
            TS_ASSERT_EQUALS(d[9],-9.5); 

            double seq3[5] = { 1, 2.5, 4, 6.5, 9 };
            VecD new_x(5,seq3,1);
            VecD out_y;
            VecD::chfe(x,y,new_x,out_y);
            TS_ASSERT_EQUALS(out_y[0],10); 
            TS_ASSERT_EQUALS(out_y[1],8); 
            TS_ASSERT_EQUALS(out_y[2],5); 
            TS_ASSERT_DELTA(out_y[3],8.1904762, 0.0001); 
            TS_ASSERT_EQUALS(out_y[4],4); 
            // The test set above is not adequate for this function:
            VecD::chfe_xy(x,y,new_x,out_y);
            TS_ASSERT_EQUALS(out_y[0],10); 
            TS_ASSERT_EQUALS(out_y[1],8); 
            TS_ASSERT_EQUALS(out_y[2],5); 
            TS_ASSERT_DELTA(out_y[3],8.1904762, 0.01); 
            TS_ASSERT_EQUALS(out_y[4],4); 
        }







// BEGIN TEMPLATE

        /************************************************************
         * TESTING VecABR (type: FLOAT)
         ************************************************************/ 
        void test_VecABR_creation( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(VecABR test1(-1));

            // This will abort the tests if bounds checking is on!
            //VecABR test2(-1);

            VecABR test3(10);
            TS_ASSERT_EQUALS(test3.length(), 10);

            VecABR test4;
            TS_ASSERT_EQUALS(test4.len(), 0);

            VecABR *test5 = new VecABR(0);
            VecABR *test99 = new VecABR(0);
            TS_ASSERT_EQUALS(test5->len(), 0);
            delete test5;
            delete test99;

            VecABR test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);

            FLOAT *arr = new FLOAT[10];
            // ****** Don't delete arr (its ownership has been taken!)
            for (int i = 0; i < 10; ++i) {
                arr[i] = 8;
            }
            VecABR test7(10, arr);
            TS_ASSERT_EQUALS(test7[9],arr[9]);
            test7[9] = 4;
            TS_ASSERT_EQUALS(test7[9],arr[9]);  // Still the same!
            // ****** Don't delete arr (its ownership has been taken!)

            // Let newly created object take ownership over deleting the array
            FLOAT *arr2 = new FLOAT[10];
            for (int i = 0; i < 10; ++i) {
                arr2[i] = 8;
            }
            VecABR test8(10, arr2, 1);  // 8 is shallow
            TS_ASSERT_EQUALS(test8[9],arr2[9]);
            test8[9] = 4;
            TS_ASSERT_EQUALS(test8[9],arr2[9]);  // Still the same!

            // Creation from another object
            // DEEP
            VecABR test9(test8, 0);
            ref_diff_dat(test9, test8);
        
            // Creation from another object
            // SHALLOW!
            VecABR test10(test8, 1);
            ref_same_dat(test10, test8);
            delete[]arr2;
        }


        void test_VecABR_vec_comparisons( void ) {
            FLOAT arr1[4] = {2,5,8,1};
            FLOAT arr2[4] = {1,6,9,1};
            VecABR x(4, arr1, 1);
            VecABR y(4, arr2, 1);

            TS_ASSERT_EQUALS(105, VecABR::dot_product(x,y));
            TS_ASSERT_DELTA(0.987985, VecABR::pearsons_r(x,y), 0.0001);
            TS_ASSERT_DELTA(1.732051, VecABR::euclidean(x,y), 0.0001);
            TS_ASSERT_EQUALS(9.25, VecABR::covariance(x,y));
        }

        void test_VecABR_all_equal( void ) {
            VecABR obj1(10,16);
            TS_ASSERT(obj1.all_equal());
            obj1[3] = 5;
            TS_ASSERT(!obj1.all_equal());
        }

        void test_VecABR_math( void ) {
            VecABR obj1(10,16);
            obj1.square_root();
            VecABR obj2(10,4);
            TS_ASSERT_EQUALS(obj1,obj2);
        }

        void test_VecABR_remove( void ) {
            FLOAT data1[] = {0,2,6,5,1,9,3,8,-1,1};
            VecABR obj1(10,data1,1);
            FLOAT ans_arr[] = {0,2,5,1,9,3,8,-1,1};
            VecABR answ(9,ans_arr,1);
            obj1.remove(2);
            TS_ASSERT_EQUALS(obj1,answ);
        }

        void test_VecABR_sort( void ) {
            FLOAT ans_arr[] = {-1,0,1,1,2,3,5,5,8,9};
            VecABR ans(10, ans_arr,1);
            FLOAT data1[] = {0,2,5,5,1,9,3,8,-1,1};
            VecABR obj1(10,data1,1);
            obj1.sort();
            TS_ASSERT_EQUALS(ans, obj1);
        }

        void test_VecABR_set( void ) {
            VecABR obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecABR obj2;
            obj2.set(obj1);
            TS_ASSERT(!obj1.shallow());
            TS_ASSERT(obj2.shallow());
            ref_same_dat(obj1, obj2);

            FLOAT *arr = new FLOAT[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            VecABR obj3;
            obj3.set(10,arr);
            TS_ASSERT_SAME_DATA(arr, (FLOAT*)obj3, 1);
            delete[] arr;
        }

        void test_VecABR_take( void ) {
            VecABR obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecABR obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            VecABR obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2[9], 5);
            TS_ASSERT_EQUALS(obj3[9], 5);
            TS_ASSERT_EQUALS(obj2.len(), 10);
            TS_ASSERT_EQUALS(obj3.len(), 10);

            // If the other one is deep too!
            VecABR obj4(10,5);
            VecABR obj5(4,2);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.len(), 10);
            TS_ASSERT_EQUALS(obj5[9], 5);

            // Test take(n,FLOAT *arr)
            VecABR obj6;
            FLOAT *arr = new FLOAT[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(10,arr);
            TS_ASSERT_EQUALS(obj6[0], 0);
            TS_ASSERT_EQUALS(obj6[9], 9);

            VecABR obj7(10,5);  // what if we are not shallow
            FLOAT *arr2 = new FLOAT[10];
            for (int i = 0; i < 10; ++i) { arr2[i] = i; }
            obj7.take(10,arr2);
            TS_ASSERT_EQUALS(obj7[0], 0);
            TS_ASSERT_EQUALS(obj7[9], 9);
        }

        void test_VecABR_abs_val( void ) {
            VecABR obj1(4,-2);
            obj1[0] = 2;
            obj1[3] = 4;
            FLOAT arr1[4] = {2,2,2,4};
            VecABR answ(4,arr1,1);
            TS_ASSERT_DIFFERS(obj1, answ);
            obj1.abs_val();
            TS_ASSERT_EQUALS(obj1, answ);
        }

        void test_VecABR_FUNCTIONS( void ) {
            // test avg:
            VecABR obj1(4,2);
            obj1[0] = 3;
            obj1[3] = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            FLOAT sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // test hist
            VecABR obj2(11);
            for (int i = 0; i < 11; ++i) {
                obj2[i] = i;
            }
            VecD b;
            VecI f;
            obj2.hist(5,b,f);
            // I'm not sure that this is the ONLY acceptable behavior here...
            TS_ASSERT_EQUALS(b.dim(), 5);
            TS_ASSERT_EQUALS(f.dim(), 5);
            TS_ASSERT_EQUALS(b[0], 1.0); 
            TS_ASSERT_EQUALS(b[1], 3.0);
            TS_ASSERT_EQUALS(b[2], 5.0);
            TS_ASSERT_EQUALS(b[3], 7.0);
            TS_ASSERT_EQUALS(b[4], 9.0);
            TS_ASSERT_EQUALS(f[0], 2);
            TS_ASSERT_EQUALS(f[1], 2);
            TS_ASSERT_EQUALS(f[2], 2);
            TS_ASSERT_EQUALS(f[3], 2);
            TS_ASSERT_EQUALS(f[4], 3);

            // test mask_as_vec
            VecABR obj3(10);
            VecI obj4(10);
            for (int i = 0; i < 10; ++i) {
                obj3[i] = i;
                obj4[i] = i;
            }
            VecABR result;
            obj3.mask_as_vec(3,obj4,result); 
            TS_ASSERT_EQUALS(result.dim(),1); 
            TS_ASSERT_EQUALS(result[0],3); 
            obj4 = 1;
            obj3.mask_as_vec(0,obj4,result); 
            TS_ASSERT_EQUALS(result.len(), 0); 
            obj3.mask_as_vec(1,obj4,result); 
            TS_ASSERT_EQUALS(result.len(), 10); 
            TS_ASSERT_EQUALS(result[0], 0); 
            TS_ASSERT_EQUALS(result[9], 9); 

            // test log
            VecABR obj5(3,2);
            obj5.logarithm(2);
            TS_ASSERT_EQUALS(obj5[0],1); 
            TS_ASSERT_EQUALS(obj5[2],1); 
            
        }

        void test_VecABR_min_max( void ) {
            VecABR obj1(10,5);
            FLOAT min; FLOAT max;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 5);
            TS_ASSERT_EQUALS(max, 5);
            for (int i = 0; i < 10; ++i) {
                obj1[i] = i;
            }
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 0);
            TS_ASSERT_EQUALS(max, 9);
            for (int i = 0; i < 10; ++i) {
                obj1[i] = 9-i;
            }
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 0);
            TS_ASSERT_EQUALS(max, 9);
        }

        void test_VecABR_indexing( void ) {
            VecABR test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);
            // Indexing a pointer to an object!
            VecABR *test7 = new VecABR(10,5);
            //test7[]
            FLOAT val = (*test7)[3];
            TS_ASSERT_EQUALS(val, 5);
            delete test7;
        }

        void test_VecABR_operator_thing_splat( void ) {
            //TS_WARN("NEED TO TEST IF OUR GUY HAS NO LENGTH (NEWLY CONSTRUCTED)");
            VecABR test7(10,5);
            FLOAT *silly = (FLOAT*)test7;
            TS_ASSERT_EQUALS(silly[2], 5);
            VecABR test8;
            FLOAT *silly2 = (FLOAT*)test8;
            TS_ASSERT(silly2 == 0);
        }

        void test_VecABR_assignment_operator( void ) {
            // First type, equals FLOAT
            VecABR test8(10,5);
            test8 = 3;
            TS_ASSERT_EQUALS(test8[2], 3);
            // Second type, equals another VecABR object:
            // A lot of this will just be memory testing (valgrind)
            VecABR test9 = test8;
            TS_ASSERT_EQUALS(test9[2], 3);
            ref_diff_dat(test9, test8);
            VecABR test10(10,2);
            test10 = test9;
            TS_ASSERT_EQUALS(test10[2], 3);
            VecABR test11;
            test11 = test9;
            TS_ASSERT_EQUALS(test11[2], 3);
        }

        void test_VecABR_copy_constructor( void ) {
            VecABR test10(10,5);
            // Two ways to do copy constructor:
            VecABR test11 = test10;
            TS_ASSERT_EQUALS(test10[0], test11[0]);
            ref_diff_dat(test10, test11);
            VecABR test12(test10);  // Will use assignmen!?
            TS_ASSERT_EQUALS(test10[0], test12[0]);
            ref_diff_dat(test10, test12);
        }

        void Xtest_VecABR_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            VecABR testA(1000000, 10);
//            st.stop();
//            //printf("\nINIT: %fs\n", st.read());
//            st.start();
//            VecABR testB;
//            testA.copy(testB);
//            st.stop();
            //printf("COPY: %fs\n", st.read());
            int length_vec = 100000000;
            VecABR testC(length_vec, 10);
            //FLOAT *arr = new FLOAT[length_vec];
            st.start();
            for (int i = 0; i < testC.len(); ++i) {
                testC[i] = 4;   // Averages 3.0s
                //arr[i] = 4;   // Averages 2.8s
            }
            st.stop();
            printf("ACCESS: %fs\n", st.read());

        }

        void test_VecABR_copy( void ) {
            VecABR test10(10,5);
            VecABR test11;
            test10.copy(test11);
            TS_ASSERT_EQUALS(test10, test11);
            TS_ASSERT_EQUALS(test11[0], 5);
            test10[2] = 3;
            TS_ASSERT_DIFFERS(test10, test11);
            VecABR test13;
            {
                VecABR test12(10,3);
                test12.copy(test13);
            }
            TS_ASSERT_EQUALS(test13[0], 3);
            // The following is a compiletime error if uncommented
            //test12[3] = 7;
        }

        void test_VecABR_equals_operator( void ) {
            VecABR test10(10,5);
            VecABR test11 = test10; 
            TS_ASSERT_EQUALS(test10, test11);
            VecABR test12(10,5); 
            TS_ASSERT_EQUALS(test10, test12);
            test12[2] = 3;
            TS_ASSERT_DIFFERS(test10, test12);
            TS_ASSERT_EQUALS(test10, test10);
            VecABR test13(11,5); 
            TS_ASSERT_DIFFERS(test10, test13);
        }

        void ref_same_dat(VecABR &one, VecABR &two) {
            FLOAT tmp = one[0];
            one[0] = 99;
            TS_ASSERT_EQUALS(one[0], two[0]);
            one[0] = tmp;
        }

        void ref_diff_dat(VecABR &one, VecABR &two) {
            FLOAT tmp = one[0];
            one[0] = 9999;
            TS_ASSERT_DIFFERS(one[0], two[0]);
            one[0] = tmp;
        }

        void test_VecABR_equal_math_ops( void ) {
            VecABR test10(10,5);
            VecABR test11(10,5);
            test10 += test11;
            TS_ASSERT_EQUALS(test10[0], 10);
            TS_ASSERT_EQUALS(test10[9], 10);
            test10 -= test11;
            TS_ASSERT_EQUALS(test10[0], 5);
            TS_ASSERT_EQUALS(test10[9], 5);
            test10 *= test11;
            TS_ASSERT_EQUALS(test10[0], 25);
            TS_ASSERT_EQUALS(test10[9], 25);
            test10 /= test11;
            TS_ASSERT_EQUALS(test10[0], 5);
            TS_ASSERT_EQUALS(test10[9], 5);
            //std::cout << test10;
            test10 += 7;
            TS_ASSERT_EQUALS(test10[0], 12);
            TS_ASSERT_EQUALS(test10[9], 12);
            test10 -= 7;
            TS_ASSERT_EQUALS(test10[0], 5);
            TS_ASSERT_EQUALS(test10[9], 5);
            test10 *= 5;
            TS_ASSERT_EQUALS(test10[0], 25);
            TS_ASSERT_EQUALS(test10[9], 25);
            test10 /= 5;
            TS_ASSERT_EQUALS(test10[0], 5);
            TS_ASSERT_EQUALS(test10[9], 5);
        }

        void test_VecABR_math_ops( void ) {
            VecABR test10(10,5);
            VecABR test11(10,5);
            VecABR test12;
            test10.add(test11, test12);
            TS_ASSERT_EQUALS(test12[0], 10);
            TS_ASSERT_EQUALS(test12[9], 10);
            // Test if the receiver already has data!
            VecABR test13(10,6);
            test10.add(test11, test13);
            TS_ASSERT_EQUALS(test13[0], 10);
            TS_ASSERT_EQUALS(test13[9], 10);
        }

        void test_VecABR_print( void ) {
            VecABR test10(10,5);
            //test10.print(); // prints 5 5 5 5 5 5 5 5 5 5
            std::ofstream fh("tmp.tmp");
            test10.print(fh);
            fh.close();
            test10.print("tmp.tmp.tmp");

             // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        void test_VecABR_sum_sq_res_yeqx_and_avg_abs_diff( void ) {
            FLOAT data1[] = {0,2,5,5};
            FLOAT data2[] = {0,1,2,3};
            VecABR d1(4,data1,1);
            VecABR d2(4,data2,1);
            double answ = VecABR::sum_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ, 7.0);
            double answ2 = VecABR::sum_sq_res_yeqx(d2,d1);
            TS_ASSERT_EQUALS(answ2, 7.0);
            double answ3 = VecABR::avg_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ3, (7.0/4.0));
            double answ4 = VecABR::avg_abs_diff(d1,d2);
            TS_ASSERT_EQUALS(answ4, (6.0/4.0));
        }

        void test_VecABR_stats( void ) {
            FLOAT data1[] = {0,2,5,5,6,7};
            VecABR d1(6,data1,1);
            double mean, std_dev;
            d1.sample_stats(mean, std_dev);
            TS_ASSERT_DELTA(mean, 4.166666, 0.00001);
            TS_ASSERT_DELTA(std_dev, 2.639444, 0.00001);

            // THIS will have to wait until public domain _normalCDF
            //double prob;
            //prob = d1.prob_one_side_right(4.166666);
            //TS_ASSERT_DELTA(prob, 0.5, 0.00001);
            //prob = d1.prob_one_side_right(1.527222);
            //TS_ASSERT_DELTA(prob, 0.158655, 0.00001);
            //prob = d1.prob_one_side_right(6.80611);
            //TS_ASSERT_DELTA(prob, 0.841345, 0.00001);
        }
        
        void test_VecABR_rsq_slope_intercept( void ) {
            FLOAT data1[] = {0,2,5,5};
            FLOAT data2[] = {1,3,4,7};
            VecABR d1(4,data1,1);
            VecABR d2(4,data2,1);
            double rsq, slope, y_intercept;
            VecABR::rsq_slope_intercept(d1,d2,rsq,slope,y_intercept);
            TS_ASSERT_DELTA(rsq, 0.758519, 0.000001);
            TS_ASSERT_DELTA(slope, 0.888888888, 0.000001);
            TS_ASSERT_DELTA(y_intercept, 1.083333333, 0.000001);
        }


        //*****************************************************************
        //*****************************************************************

// END TEMPLATE


};


