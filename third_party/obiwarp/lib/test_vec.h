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
         * TESTING VecF (type: float)
         ************************************************************/ 
        void test_VecF_creation( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(VecF test1(-1));

            // This will abort the tests if bounds checking is on!
            //VecF test2(-1);

            VecF test3(10);
            TS_ASSERT_EQUALS(test3.length(), 10);

            VecF test4;
            TS_ASSERT_EQUALS(test4.len(), 0);

            VecF *test5 = new VecF(0);
            VecF *test99 = new VecF(0);
            TS_ASSERT_EQUALS(test5->len(), 0);
            delete test5;
            delete test99;

            VecF test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);

            float *arr = new float[10];
            // ****** Don't delete arr (its ownership has been taken!)
            for (int i = 0; i < 10; ++i) {
                arr[i] = 8;
            }
            VecF test7(10, arr);
            TS_ASSERT_EQUALS(test7[9],arr[9]);
            test7[9] = 4;
            TS_ASSERT_EQUALS(test7[9],arr[9]);  // Still the same!
            // ****** Don't delete arr (its ownership has been taken!)

            // Let newly created object take ownership over deleting the array
            float *arr2 = new float[10];
            for (int i = 0; i < 10; ++i) {
                arr2[i] = 8;
            }
            VecF test8(10, arr2, 1);  // 8 is shallow
            TS_ASSERT_EQUALS(test8[9],arr2[9]);
            test8[9] = 4;
            TS_ASSERT_EQUALS(test8[9],arr2[9]);  // Still the same!

            // Creation from another object
            // DEEP
            VecF test9(test8, 0);
            ref_diff_dat(test9, test8);
        
            // Creation from another object
            // SHALLOW!
            VecF test10(test8, 1);
            ref_same_dat(test10, test8);
            delete[]arr2;
        }


        void test_VecF_vec_comparisons( void ) {
            float arr1[4] = {2,5,8,1};
            float arr2[4] = {1,6,9,1};
            VecF x(4, arr1, 1);
            VecF y(4, arr2, 1);

            TS_ASSERT_EQUALS(105, VecF::dot_product(x,y));
            TS_ASSERT_DELTA(0.987985, VecF::pearsons_r(x,y), 0.0001);
            TS_ASSERT_DELTA(1.732051, VecF::euclidean(x,y), 0.0001);
            TS_ASSERT_EQUALS(9.25, VecF::covariance(x,y));
        }

        void test_VecF_all_equal( void ) {
            VecF obj1(10,16);
            TS_ASSERT(obj1.all_equal());
            obj1[3] = 5;
            TS_ASSERT(!obj1.all_equal());
        }

        void test_VecF_math( void ) {
            VecF obj1(10,16);
            obj1.square_root();
            VecF obj2(10,4);
            TS_ASSERT_EQUALS(obj1,obj2);
        }

        void test_VecF_remove( void ) {
            float data1[] = {0,2,6,5,1,9,3,8,-1,1};
            VecF obj1(10,data1,1);
            float ans_arr[] = {0,2,5,1,9,3,8,-1,1};
            VecF answ(9,ans_arr,1);
            obj1.remove(2);
            TS_ASSERT_EQUALS(obj1,answ);
        }

        void test_VecF_sort( void ) {
            float ans_arr[] = {-1,0,1,1,2,3,5,5,8,9};
            VecF ans(10, ans_arr,1);
            float data1[] = {0,2,5,5,1,9,3,8,-1,1};
            VecF obj1(10,data1,1);
            obj1.sort();
            TS_ASSERT_EQUALS(ans, obj1);
        }

        void test_VecF_set( void ) {
            VecF obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecF obj2;
            obj2.set(obj1);
            TS_ASSERT(!obj1.shallow());
            TS_ASSERT(obj2.shallow());
            ref_same_dat(obj1, obj2);

            float *arr = new float[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            VecF obj3;
            obj3.set(10,arr);
            TS_ASSERT_SAME_DATA(arr, (float*)obj3, 1);
            delete[] arr;
        }

        void test_VecF_take( void ) {
            VecF obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecF obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            VecF obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2[9], 5);
            TS_ASSERT_EQUALS(obj3[9], 5);
            TS_ASSERT_EQUALS(obj2.len(), 10);
            TS_ASSERT_EQUALS(obj3.len(), 10);

            // If the other one is deep too!
            VecF obj4(10,5);
            VecF obj5(4,2);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.len(), 10);
            TS_ASSERT_EQUALS(obj5[9], 5);

            // Test take(n,float *arr)
            VecF obj6;
            float *arr = new float[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(10,arr);
            TS_ASSERT_EQUALS(obj6[0], 0);
            TS_ASSERT_EQUALS(obj6[9], 9);

            VecF obj7(10,5);  // what if we are not shallow
            float *arr2 = new float[10];
            for (int i = 0; i < 10; ++i) { arr2[i] = i; }
            obj7.take(10,arr2);
            TS_ASSERT_EQUALS(obj7[0], 0);
            TS_ASSERT_EQUALS(obj7[9], 9);
        }

        void test_VecF_abs_val( void ) {
            VecF obj1(4,-2);
            obj1[0] = 2;
            obj1[3] = 4;
            float arr1[4] = {2,2,2,4};
            VecF answ(4,arr1,1);
            TS_ASSERT_DIFFERS(obj1, answ);
            obj1.abs_val();
            TS_ASSERT_EQUALS(obj1, answ);
        }

        void test_VecF_FUNCTIONS( void ) {
            // test avg:
            VecF obj1(4,2);
            obj1[0] = 3;
            obj1[3] = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            float sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // test hist
            VecF obj2(11);
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
            VecF obj3(10);
            VecI obj4(10);
            for (int i = 0; i < 10; ++i) {
                obj3[i] = i;
                obj4[i] = i;
            }
            VecF result;
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
            VecF obj5(3,2);
            obj5.logarithm(2);
            TS_ASSERT_EQUALS(obj5[0],1); 
            TS_ASSERT_EQUALS(obj5[2],1); 
            
        }

        void test_VecF_min_max( void ) {
            VecF obj1(10,5);
            float min; float max;
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

        void test_VecF_indexing( void ) {
            VecF test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);
            // Indexing a pointer to an object!
            VecF *test7 = new VecF(10,5);
            //test7[]
            float val = (*test7)[3];
            TS_ASSERT_EQUALS(val, 5);
            delete test7;
        }

        void test_VecF_operator_thing_splat( void ) {
            //TS_WARN("NEED TO TEST IF OUR GUY HAS NO LENGTH (NEWLY CONSTRUCTED)");
            VecF test7(10,5);
            float *silly = (float*)test7;
            TS_ASSERT_EQUALS(silly[2], 5);
            VecF test8;
            float *silly2 = (float*)test8;
            TS_ASSERT(silly2 == 0);
        }

        void test_VecF_assignment_operator( void ) {
            // First type, equals float
            VecF test8(10,5);
            test8 = 3;
            TS_ASSERT_EQUALS(test8[2], 3);
            // Second type, equals another VecF object:
            // A lot of this will just be memory testing (valgrind)
            VecF test9 = test8;
            TS_ASSERT_EQUALS(test9[2], 3);
            ref_diff_dat(test9, test8);
            VecF test10(10,2);
            test10 = test9;
            TS_ASSERT_EQUALS(test10[2], 3);
            VecF test11;
            test11 = test9;
            TS_ASSERT_EQUALS(test11[2], 3);
        }

        void test_VecF_copy_constructor( void ) {
            VecF test10(10,5);
            // Two ways to do copy constructor:
            VecF test11 = test10;
            TS_ASSERT_EQUALS(test10[0], test11[0]);
            ref_diff_dat(test10, test11);
            VecF test12(test10);  // Will use assignmen!?
            TS_ASSERT_EQUALS(test10[0], test12[0]);
            ref_diff_dat(test10, test12);
        }

        void Xtest_VecF_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            VecF testA(1000000, 10);
//            st.stop();
//            //printf("\nINIT: %fs\n", st.read());
//            st.start();
//            VecF testB;
//            testA.copy(testB);
//            st.stop();
            //printf("COPY: %fs\n", st.read());
            int length_vec = 100000000;
            VecF testC(length_vec, 10);
            //float *arr = new float[length_vec];
            st.start();
            for (int i = 0; i < testC.len(); ++i) {
                testC[i] = 4;   // Averages 3.0s
                //arr[i] = 4;   // Averages 2.8s
            }
            st.stop();
            printf("ACCESS: %fs\n", st.read());

        }

        void test_VecF_copy( void ) {
            VecF test10(10,5);
            VecF test11;
            test10.copy(test11);
            TS_ASSERT_EQUALS(test10, test11);
            TS_ASSERT_EQUALS(test11[0], 5);
            test10[2] = 3;
            TS_ASSERT_DIFFERS(test10, test11);
            VecF test13;
            {
                VecF test12(10,3);
                test12.copy(test13);
            }
            TS_ASSERT_EQUALS(test13[0], 3);
            // The following is a compiletime error if uncommented
            //test12[3] = 7;
        }

        void test_VecF_equals_operator( void ) {
            VecF test10(10,5);
            VecF test11 = test10; 
            TS_ASSERT_EQUALS(test10, test11);
            VecF test12(10,5); 
            TS_ASSERT_EQUALS(test10, test12);
            test12[2] = 3;
            TS_ASSERT_DIFFERS(test10, test12);
            TS_ASSERT_EQUALS(test10, test10);
            VecF test13(11,5); 
            TS_ASSERT_DIFFERS(test10, test13);
        }

        void ref_same_dat(VecF &one, VecF &two) {
            float tmp = one[0];
            one[0] = 99;
            TS_ASSERT_EQUALS(one[0], two[0]);
            one[0] = tmp;
        }

        void ref_diff_dat(VecF &one, VecF &two) {
            float tmp = one[0];
            one[0] = 9999;
            TS_ASSERT_DIFFERS(one[0], two[0]);
            one[0] = tmp;
        }

        void test_VecF_equal_math_ops( void ) {
            VecF test10(10,5);
            VecF test11(10,5);
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

        void test_VecF_math_ops( void ) {
            VecF test10(10,5);
            VecF test11(10,5);
            VecF test12;
            test10.add(test11, test12);
            TS_ASSERT_EQUALS(test12[0], 10);
            TS_ASSERT_EQUALS(test12[9], 10);
            // Test if the receiver already has data!
            VecF test13(10,6);
            test10.add(test11, test13);
            TS_ASSERT_EQUALS(test13[0], 10);
            TS_ASSERT_EQUALS(test13[9], 10);
        }

        void test_VecF_print( void ) {
            VecF test10(10,5);
            //test10.print(); // prints 5 5 5 5 5 5 5 5 5 5
            std::ofstream fh("tmp.tmp");
            test10.print(fh);
            fh.close();
            test10.print("tmp.tmp.tmp");

             // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        void test_VecF_sum_sq_res_yeqx_and_avg_abs_diff( void ) {
            float data1[] = {0,2,5,5};
            float data2[] = {0,1,2,3};
            VecF d1(4,data1,1);
            VecF d2(4,data2,1);
            double answ = VecF::sum_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ, 7.0);
            double answ2 = VecF::sum_sq_res_yeqx(d2,d1);
            TS_ASSERT_EQUALS(answ2, 7.0);
            double answ3 = VecF::avg_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ3, (7.0/4.0));
            double answ4 = VecF::avg_abs_diff(d1,d2);
            TS_ASSERT_EQUALS(answ4, (6.0/4.0));
        }

        void test_VecF_stats( void ) {
            float data1[] = {0,2,5,5,6,7};
            VecF d1(6,data1,1);
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
        
        void test_VecF_rsq_slope_intercept( void ) {
            float data1[] = {0,2,5,5};
            float data2[] = {1,3,4,7};
            VecF d1(4,data1,1);
            VecF d2(4,data2,1);
            double rsq, slope, y_intercept;
            VecF::rsq_slope_intercept(d1,d2,rsq,slope,y_intercept);
            TS_ASSERT_DELTA(rsq, 0.758519, 0.000001);
            TS_ASSERT_DELTA(slope, 0.888888888, 0.000001);
            TS_ASSERT_DELTA(y_intercept, 1.083333333, 0.000001);
        }


        //*****************************************************************
        //*****************************************************************


        /************************************************************
         * TESTING VecD (type: double)
         ************************************************************/ 
        void test_VecD_creation( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(VecD test1(-1));

            // This will abort the tests if bounds checking is on!
            //VecD test2(-1);

            VecD test3(10);
            TS_ASSERT_EQUALS(test3.length(), 10);

            VecD test4;
            TS_ASSERT_EQUALS(test4.len(), 0);

            VecD *test5 = new VecD(0);
            VecD *test99 = new VecD(0);
            TS_ASSERT_EQUALS(test5->len(), 0);
            delete test5;
            delete test99;

            VecD test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);

            double *arr = new double[10];
            // ****** Don't delete arr (its ownership has been taken!)
            for (int i = 0; i < 10; ++i) {
                arr[i] = 8;
            }
            VecD test7(10, arr);
            TS_ASSERT_EQUALS(test7[9],arr[9]);
            test7[9] = 4;
            TS_ASSERT_EQUALS(test7[9],arr[9]);  // Still the same!
            // ****** Don't delete arr (its ownership has been taken!)

            // Let newly created object take ownership over deleting the array
            double *arr2 = new double[10];
            for (int i = 0; i < 10; ++i) {
                arr2[i] = 8;
            }
            VecD test8(10, arr2, 1);  // 8 is shallow
            TS_ASSERT_EQUALS(test8[9],arr2[9]);
            test8[9] = 4;
            TS_ASSERT_EQUALS(test8[9],arr2[9]);  // Still the same!

            // Creation from another object
            // DEEP
            VecD test9(test8, 0);
            ref_diff_dat(test9, test8);
        
            // Creation from another object
            // SHALLOW!
            VecD test10(test8, 1);
            ref_same_dat(test10, test8);
            delete[]arr2;
        }


        void test_VecD_vec_comparisons( void ) {
            double arr1[4] = {2,5,8,1};
            double arr2[4] = {1,6,9,1};
            VecD x(4, arr1, 1);
            VecD y(4, arr2, 1);

            TS_ASSERT_EQUALS(105, VecD::dot_product(x,y));
            TS_ASSERT_DELTA(0.987985, VecD::pearsons_r(x,y), 0.0001);
            TS_ASSERT_DELTA(1.732051, VecD::euclidean(x,y), 0.0001);
            TS_ASSERT_EQUALS(9.25, VecD::covariance(x,y));
        }

        void test_VecD_all_equal( void ) {
            VecD obj1(10,16);
            TS_ASSERT(obj1.all_equal());
            obj1[3] = 5;
            TS_ASSERT(!obj1.all_equal());
        }

        void test_VecD_math( void ) {
            VecD obj1(10,16);
            obj1.square_root();
            VecD obj2(10,4);
            TS_ASSERT_EQUALS(obj1,obj2);
        }

        void test_VecD_remove( void ) {
            double data1[] = {0,2,6,5,1,9,3,8,-1,1};
            VecD obj1(10,data1,1);
            double ans_arr[] = {0,2,5,1,9,3,8,-1,1};
            VecD answ(9,ans_arr,1);
            obj1.remove(2);
            TS_ASSERT_EQUALS(obj1,answ);
        }

        void test_VecD_sort( void ) {
            double ans_arr[] = {-1,0,1,1,2,3,5,5,8,9};
            VecD ans(10, ans_arr,1);
            double data1[] = {0,2,5,5,1,9,3,8,-1,1};
            VecD obj1(10,data1,1);
            obj1.sort();
            TS_ASSERT_EQUALS(ans, obj1);
        }

        void test_VecD_set( void ) {
            VecD obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecD obj2;
            obj2.set(obj1);
            TS_ASSERT(!obj1.shallow());
            TS_ASSERT(obj2.shallow());
            ref_same_dat(obj1, obj2);

            double *arr = new double[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            VecD obj3;
            obj3.set(10,arr);
            TS_ASSERT_SAME_DATA(arr, (double*)obj3, 1);
            delete[] arr;
        }

        void test_VecD_take( void ) {
            VecD obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecD obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            VecD obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2[9], 5);
            TS_ASSERT_EQUALS(obj3[9], 5);
            TS_ASSERT_EQUALS(obj2.len(), 10);
            TS_ASSERT_EQUALS(obj3.len(), 10);

            // If the other one is deep too!
            VecD obj4(10,5);
            VecD obj5(4,2);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.len(), 10);
            TS_ASSERT_EQUALS(obj5[9], 5);

            // Test take(n,double *arr)
            VecD obj6;
            double *arr = new double[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(10,arr);
            TS_ASSERT_EQUALS(obj6[0], 0);
            TS_ASSERT_EQUALS(obj6[9], 9);

            VecD obj7(10,5);  // what if we are not shallow
            double *arr2 = new double[10];
            for (int i = 0; i < 10; ++i) { arr2[i] = i; }
            obj7.take(10,arr2);
            TS_ASSERT_EQUALS(obj7[0], 0);
            TS_ASSERT_EQUALS(obj7[9], 9);
        }

        void test_VecD_abs_val( void ) {
            VecD obj1(4,-2);
            obj1[0] = 2;
            obj1[3] = 4;
            double arr1[4] = {2,2,2,4};
            VecD answ(4,arr1,1);
            TS_ASSERT_DIFFERS(obj1, answ);
            obj1.abs_val();
            TS_ASSERT_EQUALS(obj1, answ);
        }

        void test_VecD_FUNCTIONS( void ) {
            // test avg:
            VecD obj1(4,2);
            obj1[0] = 3;
            obj1[3] = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            double sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // test hist
            VecD obj2(11);
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
            VecD obj3(10);
            VecI obj4(10);
            for (int i = 0; i < 10; ++i) {
                obj3[i] = i;
                obj4[i] = i;
            }
            VecD result;
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
            VecD obj5(3,2);
            obj5.logarithm(2);
            TS_ASSERT_EQUALS(obj5[0],1); 
            TS_ASSERT_EQUALS(obj5[2],1); 
            
        }

        void test_VecD_min_max( void ) {
            VecD obj1(10,5);
            double min; double max;
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

        void test_VecD_indexing( void ) {
            VecD test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);
            // Indexing a pointer to an object!
            VecD *test7 = new VecD(10,5);
            //test7[]
            double val = (*test7)[3];
            TS_ASSERT_EQUALS(val, 5);
            delete test7;
        }

        void test_VecD_operator_thing_splat( void ) {
            //TS_WARN("NEED TO TEST IF OUR GUY HAS NO LENGTH (NEWLY CONSTRUCTED)");
            VecD test7(10,5);
            double *silly = (double*)test7;
            TS_ASSERT_EQUALS(silly[2], 5);
            VecD test8;
            double *silly2 = (double*)test8;
            TS_ASSERT(silly2 == 0);
        }

        void test_VecD_assignment_operator( void ) {
            // First type, equals double
            VecD test8(10,5);
            test8 = 3;
            TS_ASSERT_EQUALS(test8[2], 3);
            // Second type, equals another VecD object:
            // A lot of this will just be memory testing (valgrind)
            VecD test9 = test8;
            TS_ASSERT_EQUALS(test9[2], 3);
            ref_diff_dat(test9, test8);
            VecD test10(10,2);
            test10 = test9;
            TS_ASSERT_EQUALS(test10[2], 3);
            VecD test11;
            test11 = test9;
            TS_ASSERT_EQUALS(test11[2], 3);
        }

        void test_VecD_copy_constructor( void ) {
            VecD test10(10,5);
            // Two ways to do copy constructor:
            VecD test11 = test10;
            TS_ASSERT_EQUALS(test10[0], test11[0]);
            ref_diff_dat(test10, test11);
            VecD test12(test10);  // Will use assignmen!?
            TS_ASSERT_EQUALS(test10[0], test12[0]);
            ref_diff_dat(test10, test12);
        }

        void Xtest_VecD_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            VecD testA(1000000, 10);
//            st.stop();
//            //printf("\nINIT: %fs\n", st.read());
//            st.start();
//            VecD testB;
//            testA.copy(testB);
//            st.stop();
            //printf("COPY: %fs\n", st.read());
            int length_vec = 100000000;
            VecD testC(length_vec, 10);
            //double *arr = new double[length_vec];
            st.start();
            for (int i = 0; i < testC.len(); ++i) {
                testC[i] = 4;   // Averages 3.0s
                //arr[i] = 4;   // Averages 2.8s
            }
            st.stop();
            printf("ACCESS: %fs\n", st.read());

        }

        void test_VecD_copy( void ) {
            VecD test10(10,5);
            VecD test11;
            test10.copy(test11);
            TS_ASSERT_EQUALS(test10, test11);
            TS_ASSERT_EQUALS(test11[0], 5);
            test10[2] = 3;
            TS_ASSERT_DIFFERS(test10, test11);
            VecD test13;
            {
                VecD test12(10,3);
                test12.copy(test13);
            }
            TS_ASSERT_EQUALS(test13[0], 3);
            // The following is a compiletime error if uncommented
            //test12[3] = 7;
        }

        void test_VecD_equals_operator( void ) {
            VecD test10(10,5);
            VecD test11 = test10; 
            TS_ASSERT_EQUALS(test10, test11);
            VecD test12(10,5); 
            TS_ASSERT_EQUALS(test10, test12);
            test12[2] = 3;
            TS_ASSERT_DIFFERS(test10, test12);
            TS_ASSERT_EQUALS(test10, test10);
            VecD test13(11,5); 
            TS_ASSERT_DIFFERS(test10, test13);
        }

        void ref_same_dat(VecD &one, VecD &two) {
            double tmp = one[0];
            one[0] = 99;
            TS_ASSERT_EQUALS(one[0], two[0]);
            one[0] = tmp;
        }

        void ref_diff_dat(VecD &one, VecD &two) {
            double tmp = one[0];
            one[0] = 9999;
            TS_ASSERT_DIFFERS(one[0], two[0]);
            one[0] = tmp;
        }

        void test_VecD_equal_math_ops( void ) {
            VecD test10(10,5);
            VecD test11(10,5);
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

        void test_VecD_math_ops( void ) {
            VecD test10(10,5);
            VecD test11(10,5);
            VecD test12;
            test10.add(test11, test12);
            TS_ASSERT_EQUALS(test12[0], 10);
            TS_ASSERT_EQUALS(test12[9], 10);
            // Test if the receiver already has data!
            VecD test13(10,6);
            test10.add(test11, test13);
            TS_ASSERT_EQUALS(test13[0], 10);
            TS_ASSERT_EQUALS(test13[9], 10);
        }

        void test_VecD_print( void ) {
            VecD test10(10,5);
            //test10.print(); // prints 5 5 5 5 5 5 5 5 5 5
            std::ofstream fh("tmp.tmp");
            test10.print(fh);
            fh.close();
            test10.print("tmp.tmp.tmp");

             // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        void test_VecD_sum_sq_res_yeqx_and_avg_abs_diff( void ) {
            double data1[] = {0,2,5,5};
            double data2[] = {0,1,2,3};
            VecD d1(4,data1,1);
            VecD d2(4,data2,1);
            double answ = VecD::sum_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ, 7.0);
            double answ2 = VecD::sum_sq_res_yeqx(d2,d1);
            TS_ASSERT_EQUALS(answ2, 7.0);
            double answ3 = VecD::avg_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ3, (7.0/4.0));
            double answ4 = VecD::avg_abs_diff(d1,d2);
            TS_ASSERT_EQUALS(answ4, (6.0/4.0));
        }

        void test_VecD_stats( void ) {
            double data1[] = {0,2,5,5,6,7};
            VecD d1(6,data1,1);
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
        
        void test_VecD_rsq_slope_intercept( void ) {
            double data1[] = {0,2,5,5};
            double data2[] = {1,3,4,7};
            VecD d1(4,data1,1);
            VecD d2(4,data2,1);
            double rsq, slope, y_intercept;
            VecD::rsq_slope_intercept(d1,d2,rsq,slope,y_intercept);
            TS_ASSERT_DELTA(rsq, 0.758519, 0.000001);
            TS_ASSERT_DELTA(slope, 0.888888888, 0.000001);
            TS_ASSERT_DELTA(y_intercept, 1.083333333, 0.000001);
        }


        //*****************************************************************
        //*****************************************************************


        /************************************************************
         * TESTING VecI (type: int)
         ************************************************************/ 
        void test_VecI_creation( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(VecI test1(-1));

            // This will abort the tests if bounds checking is on!
            //VecI test2(-1);

            VecI test3(10);
            TS_ASSERT_EQUALS(test3.length(), 10);

            VecI test4;
            TS_ASSERT_EQUALS(test4.len(), 0);

            VecI *test5 = new VecI(0);
            VecI *test99 = new VecI(0);
            TS_ASSERT_EQUALS(test5->len(), 0);
            delete test5;
            delete test99;

            VecI test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);

            int *arr = new int[10];
            // ****** Don't delete arr (its ownership has been taken!)
            for (int i = 0; i < 10; ++i) {
                arr[i] = 8;
            }
            VecI test7(10, arr);
            TS_ASSERT_EQUALS(test7[9],arr[9]);
            test7[9] = 4;
            TS_ASSERT_EQUALS(test7[9],arr[9]);  // Still the same!
            // ****** Don't delete arr (its ownership has been taken!)

            // Let newly created object take ownership over deleting the array
            int *arr2 = new int[10];
            for (int i = 0; i < 10; ++i) {
                arr2[i] = 8;
            }
            VecI test8(10, arr2, 1);  // 8 is shallow
            TS_ASSERT_EQUALS(test8[9],arr2[9]);
            test8[9] = 4;
            TS_ASSERT_EQUALS(test8[9],arr2[9]);  // Still the same!

            // Creation from another object
            // DEEP
            VecI test9(test8, 0);
            ref_diff_dat(test9, test8);
        
            // Creation from another object
            // SHALLOW!
            VecI test10(test8, 1);
            ref_same_dat(test10, test8);
            delete[]arr2;
        }


        void test_VecI_vec_comparisons( void ) {
            int arr1[4] = {2,5,8,1};
            int arr2[4] = {1,6,9,1};
            VecI x(4, arr1, 1);
            VecI y(4, arr2, 1);

            TS_ASSERT_EQUALS(105, VecI::dot_product(x,y));
            TS_ASSERT_DELTA(0.987985, VecI::pearsons_r(x,y), 0.0001);
            TS_ASSERT_DELTA(1.732051, VecI::euclidean(x,y), 0.0001);
            TS_ASSERT_EQUALS(9.25, VecI::covariance(x,y));
        }

        void test_VecI_all_equal( void ) {
            VecI obj1(10,16);
            TS_ASSERT(obj1.all_equal());
            obj1[3] = 5;
            TS_ASSERT(!obj1.all_equal());
        }

        void test_VecI_math( void ) {
            VecI obj1(10,16);
            obj1.square_root();
            VecI obj2(10,4);
            TS_ASSERT_EQUALS(obj1,obj2);
        }

        void test_VecI_remove( void ) {
            int data1[] = {0,2,6,5,1,9,3,8,-1,1};
            VecI obj1(10,data1,1);
            int ans_arr[] = {0,2,5,1,9,3,8,-1,1};
            VecI answ(9,ans_arr,1);
            obj1.remove(2);
            TS_ASSERT_EQUALS(obj1,answ);
        }

        void test_VecI_sort( void ) {
            int ans_arr[] = {-1,0,1,1,2,3,5,5,8,9};
            VecI ans(10, ans_arr,1);
            int data1[] = {0,2,5,5,1,9,3,8,-1,1};
            VecI obj1(10,data1,1);
            obj1.sort();
            TS_ASSERT_EQUALS(ans, obj1);
        }

        void test_VecI_set( void ) {
            VecI obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecI obj2;
            obj2.set(obj1);
            TS_ASSERT(!obj1.shallow());
            TS_ASSERT(obj2.shallow());
            ref_same_dat(obj1, obj2);

            int *arr = new int[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            VecI obj3;
            obj3.set(10,arr);
            TS_ASSERT_SAME_DATA(arr, (int*)obj3, 1);
            delete[] arr;
        }

        void test_VecI_take( void ) {
            VecI obj1(10,5); 
            TS_ASSERT(!obj1.shallow());
            VecI obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            VecI obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2[9], 5);
            TS_ASSERT_EQUALS(obj3[9], 5);
            TS_ASSERT_EQUALS(obj2.len(), 10);
            TS_ASSERT_EQUALS(obj3.len(), 10);

            // If the other one is deep too!
            VecI obj4(10,5);
            VecI obj5(4,2);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.len(), 10);
            TS_ASSERT_EQUALS(obj5[9], 5);

            // Test take(n,int *arr)
            VecI obj6;
            int *arr = new int[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(10,arr);
            TS_ASSERT_EQUALS(obj6[0], 0);
            TS_ASSERT_EQUALS(obj6[9], 9);

            VecI obj7(10,5);  // what if we are not shallow
            int *arr2 = new int[10];
            for (int i = 0; i < 10; ++i) { arr2[i] = i; }
            obj7.take(10,arr2);
            TS_ASSERT_EQUALS(obj7[0], 0);
            TS_ASSERT_EQUALS(obj7[9], 9);
        }

        void test_VecI_abs_val( void ) {
            VecI obj1(4,-2);
            obj1[0] = 2;
            obj1[3] = 4;
            int arr1[4] = {2,2,2,4};
            VecI answ(4,arr1,1);
            TS_ASSERT_DIFFERS(obj1, answ);
            obj1.abs_val();
            TS_ASSERT_EQUALS(obj1, answ);
        }

        void test_VecI_FUNCTIONS( void ) {
            // test avg:
            VecI obj1(4,2);
            obj1[0] = 3;
            obj1[3] = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            int sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // test hist
            VecI obj2(11);
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
            VecI obj3(10);
            VecI obj4(10);
            for (int i = 0; i < 10; ++i) {
                obj3[i] = i;
                obj4[i] = i;
            }
            VecI result;
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
            VecI obj5(3,2);
            obj5.logarithm(2);
            TS_ASSERT_EQUALS(obj5[0],1); 
            TS_ASSERT_EQUALS(obj5[2],1); 
            
        }

        void test_VecI_min_max( void ) {
            VecI obj1(10,5);
            int min; int max;
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

        void test_VecI_indexing( void ) {
            VecI test6(10,5);
            TS_ASSERT_EQUALS(test6.len(), 10);
            TS_ASSERT_EQUALS(test6[4], 5);
            // Indexing a pointer to an object!
            VecI *test7 = new VecI(10,5);
            //test7[]
            int val = (*test7)[3];
            TS_ASSERT_EQUALS(val, 5);
            delete test7;
        }

        void test_VecI_operator_thing_splat( void ) {
            //TS_WARN("NEED TO TEST IF OUR GUY HAS NO LENGTH (NEWLY CONSTRUCTED)");
            VecI test7(10,5);
            int *silly = (int*)test7;
            TS_ASSERT_EQUALS(silly[2], 5);
            VecI test8;
            int *silly2 = (int*)test8;
            TS_ASSERT(silly2 == 0);
        }

        void test_VecI_assignment_operator( void ) {
            // First type, equals int
            VecI test8(10,5);
            test8 = 3;
            TS_ASSERT_EQUALS(test8[2], 3);
            // Second type, equals another VecI object:
            // A lot of this will just be memory testing (valgrind)
            VecI test9 = test8;
            TS_ASSERT_EQUALS(test9[2], 3);
            ref_diff_dat(test9, test8);
            VecI test10(10,2);
            test10 = test9;
            TS_ASSERT_EQUALS(test10[2], 3);
            VecI test11;
            test11 = test9;
            TS_ASSERT_EQUALS(test11[2], 3);
        }

        void test_VecI_copy_constructor( void ) {
            VecI test10(10,5);
            // Two ways to do copy constructor:
            VecI test11 = test10;
            TS_ASSERT_EQUALS(test10[0], test11[0]);
            ref_diff_dat(test10, test11);
            VecI test12(test10);  // Will use assignmen!?
            TS_ASSERT_EQUALS(test10[0], test12[0]);
            ref_diff_dat(test10, test12);
        }

        void Xtest_VecI_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            VecI testA(1000000, 10);
//            st.stop();
//            //printf("\nINIT: %fs\n", st.read());
//            st.start();
//            VecI testB;
//            testA.copy(testB);
//            st.stop();
            //printf("COPY: %fs\n", st.read());
            int length_vec = 100000000;
            VecI testC(length_vec, 10);
            //int *arr = new int[length_vec];
            st.start();
            for (int i = 0; i < testC.len(); ++i) {
                testC[i] = 4;   // Averages 3.0s
                //arr[i] = 4;   // Averages 2.8s
            }
            st.stop();
            printf("ACCESS: %fs\n", st.read());

        }

        void test_VecI_copy( void ) {
            VecI test10(10,5);
            VecI test11;
            test10.copy(test11);
            TS_ASSERT_EQUALS(test10, test11);
            TS_ASSERT_EQUALS(test11[0], 5);
            test10[2] = 3;
            TS_ASSERT_DIFFERS(test10, test11);
            VecI test13;
            {
                VecI test12(10,3);
                test12.copy(test13);
            }
            TS_ASSERT_EQUALS(test13[0], 3);
            // The following is a compiletime error if uncommented
            //test12[3] = 7;
        }

        void test_VecI_equals_operator( void ) {
            VecI test10(10,5);
            VecI test11 = test10; 
            TS_ASSERT_EQUALS(test10, test11);
            VecI test12(10,5); 
            TS_ASSERT_EQUALS(test10, test12);
            test12[2] = 3;
            TS_ASSERT_DIFFERS(test10, test12);
            TS_ASSERT_EQUALS(test10, test10);
            VecI test13(11,5); 
            TS_ASSERT_DIFFERS(test10, test13);
        }

        void ref_same_dat(VecI &one, VecI &two) {
            int tmp = one[0];
            one[0] = 99;
            TS_ASSERT_EQUALS(one[0], two[0]);
            one[0] = tmp;
        }

        void ref_diff_dat(VecI &one, VecI &two) {
            int tmp = one[0];
            one[0] = 9999;
            TS_ASSERT_DIFFERS(one[0], two[0]);
            one[0] = tmp;
        }

        void test_VecI_equal_math_ops( void ) {
            VecI test10(10,5);
            VecI test11(10,5);
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

        void test_VecI_math_ops( void ) {
            VecI test10(10,5);
            VecI test11(10,5);
            VecI test12;
            test10.add(test11, test12);
            TS_ASSERT_EQUALS(test12[0], 10);
            TS_ASSERT_EQUALS(test12[9], 10);
            // Test if the receiver already has data!
            VecI test13(10,6);
            test10.add(test11, test13);
            TS_ASSERT_EQUALS(test13[0], 10);
            TS_ASSERT_EQUALS(test13[9], 10);
        }

        void test_VecI_print( void ) {
            VecI test10(10,5);
            //test10.print(); // prints 5 5 5 5 5 5 5 5 5 5
            std::ofstream fh("tmp.tmp");
            test10.print(fh);
            fh.close();
            test10.print("tmp.tmp.tmp");

             // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        void test_VecI_sum_sq_res_yeqx_and_avg_abs_diff( void ) {
            int data1[] = {0,2,5,5};
            int data2[] = {0,1,2,3};
            VecI d1(4,data1,1);
            VecI d2(4,data2,1);
            double answ = VecI::sum_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ, 7.0);
            double answ2 = VecI::sum_sq_res_yeqx(d2,d1);
            TS_ASSERT_EQUALS(answ2, 7.0);
            double answ3 = VecI::avg_sq_res_yeqx(d1,d2);
            TS_ASSERT_EQUALS(answ3, (7.0/4.0));
            double answ4 = VecI::avg_abs_diff(d1,d2);
            TS_ASSERT_EQUALS(answ4, (6.0/4.0));
        }

        void test_VecI_stats( void ) {
            int data1[] = {0,2,5,5,6,7};
            VecI d1(6,data1,1);
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
        
        void test_VecI_rsq_slope_intercept( void ) {
            int data1[] = {0,2,5,5};
            int data2[] = {1,3,4,7};
            VecI d1(4,data1,1);
            VecI d2(4,data2,1);
            double rsq, slope, y_intercept;
            VecI::rsq_slope_intercept(d1,d2,rsq,slope,y_intercept);
            TS_ASSERT_DELTA(rsq, 0.758519, 0.000001);
            TS_ASSERT_DELTA(slope, 0.888888888, 0.000001);
            TS_ASSERT_DELTA(y_intercept, 1.083333333, 0.000001);
        }


        //*****************************************************************
        //*****************************************************************

// END TEMPLATE


};


