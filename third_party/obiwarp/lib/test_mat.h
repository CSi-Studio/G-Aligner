#include <cxxtest/TestSuite.h>
#include "mat.h"
#include <cstdlib>
#include <fstream>
#include "stdio.h"
//#include "vec_utils.h"
#include "tnt_stopwatch.h"

using namespace VEC;


class MatTestSuite : public CxxTest::TestSuite 
{
    public:

// BEGIN TEMPLATE

        /************************************************************
         * TESTING MatF (type: float)
         ************************************************************/ 
        void test_MatF_creation_and_indexing( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(MatF test1(-1));

            // This will abort the tests if bounds checking is on!
            //MatF test2(-1);

            MatF test3(10,2);
            TS_ASSERT_EQUALS(test3.dim1(), 10);
            TS_ASSERT_EQUALS(test3.dim2(), 2);

            MatF test4;
            TS_ASSERT_EQUALS(test4.dim1(), 0);
            TS_ASSERT_EQUALS(test4.dim2(), 0);

            MatF *test5 = new MatF(0,0);
            TS_ASSERT_EQUALS(test5->dim1(), 0);
            TS_ASSERT_EQUALS(test5->dim2(), 0);
            delete test5;

            MatF test6(10,5,4);
            TS_ASSERT_EQUALS(test6.dim1(), 10);
            TS_ASSERT_EQUALS(test6.dim2(), 5);
            TS_ASSERT_EQUALS(test6(0,0), 4);
            TS_ASSERT_EQUALS(test6(9,4), 4);

            // set with float array!
            float *arr = new float[50];
            for (int i = 0; i < 50; ++i) {
                arr[i] = 9;
            }
            MatF test7(10,5,arr);
            TS_ASSERT_EQUALS(test7(0,0), 9);
            test7(0,0) = 8;
            TS_ASSERT_EQUALS(test7(0,0), 8);
            TS_ASSERT_EQUALS(test7(0,0), arr[0]);

            float *arr2 = new float[50];
            for (int i = 0; i < 50; ++i) {
                arr2[i] = 9;
            }
            MatF test8(10,5,arr2);
            TS_ASSERT_EQUALS(test8(0,0), 9);
            test8(0,0) = 8;
            TS_ASSERT_EQUALS(test8(0,0), 8);
            TS_ASSERT_EQUALS(test8(0,0), arr2[0]);
            // Notice we don't have to delete this guy, because we gave him
            // away! (have to do alloc/free test with valgrind!

            MatF test9(10,5,2);
            MatF test10(test9);
            ref_diff_dat(test9, test10);

            MatF test12(test9, 1);
            ref_same_dat(test9, test12);
        }

        void test_MatF_to_vec( void ) {
            MatF obj1(4,3,3);
            VecF out;
            obj1.to_vec(out);  // deep
            TS_ASSERT_EQUALS(out, obj1._dat);
            out[0] = 6;
            TS_ASSERT_DIFFERS(out, obj1._dat);
            VecF outshallow;
            obj1.to_vec(outshallow, 1);  // shallow
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
            outshallow[0] = 6;
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
        }

        void test_MatF_row_vecs( void ) {
            MatF obj1(10,5,3);

            // Row vecs:
            VecF rvecs[10];
            int cnt;
            obj1.row_vecs(cnt, rvecs);
            TS_ASSERT_EQUALS(cnt, 10);
            VecF rvec_answ(5,3);
            for (int i = 0; i < cnt; ++i) {
                TS_ASSERT_EQUALS(rvecs[i], rvec_answ);
            }
        }

        void test_MatF_set_from_ascii( void ) {
            MatF obj1;
            obj1.set_from_ascii("tfiles/tmp1.mata");
            TS_ASSERT_EQUALS(obj1.rows(), 6);
            TS_ASSERT_EQUALS(obj1.cols(), 10);
            TS_ASSERT_EQUALS(obj1(0,0), 230);
            TS_ASSERT_EQUALS(obj1(5,9), 657);

            MatF obj2;
            // This file is perfect
            obj2.set_from_ascii("tfiles/tmp1_no_header.mata", 1);
            TS_ASSERT_EQUALS(obj2.rows(), 6);
            TS_ASSERT_EQUALS(obj2.cols(), 10);
            TS_ASSERT_EQUALS(obj2(0,0), 230);
            TS_ASSERT_EQUALS(obj2(5,9), 657);

            MatF obj3;
            // This file has an extra space at the end of the first line
            // extra spaces at the end of some other lines
            // and an extra line at the end
            obj3.set_from_ascii("tfiles/tmp1_no_header_messy.mata", 1);
            TS_ASSERT_EQUALS(obj3.rows(), 6);
            TS_ASSERT_EQUALS(obj3.cols(), 10);
            TS_ASSERT_EQUALS(obj3(0,0), 230);
            TS_ASSERT_EQUALS(obj3(5,9), 657);
        }
        
        void test_MatF_transpose( void ) {
            MatF obj1(4,2,5);
            obj1(3,0) = 10;
            MatF trans;
            obj1.transpose(trans);
            TS_ASSERT_EQUALS(obj1.rows(), trans.cols());
            TS_ASSERT_EQUALS(obj1.cols(), trans.rows());
            TS_ASSERT_EQUALS(trans(0,0), 5);
            TS_ASSERT_EQUALS(trans(0,3), obj1(3,0));
            TS_ASSERT_EQUALS(trans(0,3), 10);
            TS_ASSERT_EQUALS(trans(1,0), 5);
        }

        void test_MatF_FUNCTIONS( void ) {
            // test avg:
            MatF obj1(2,2,2);
            obj1(0,1) = 3;
            obj1(1,0) = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            // test sum:
            float sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // sum of a row:
            // 2 3 2 2 2  =  11
            // 2 2 2 2 2  =  10
            MatF obj2(2,5,2);
            obj2(0,2) = 3;
            TS_ASSERT_EQUALS(obj2.sum(0),11);
            TS_ASSERT_EQUALS(obj2.sum(1),10);

            // test mask_as_vec:
            MatF obj3(5,2);
            MatI obj4(5,2);
            int cnt = 0;
            for (int m = 0; m < 5; ++m) {
                for (int n = 0; n < 2; ++n) {
                    obj3(m,n) = cnt;
                    obj4(m,n) = cnt;
                    cnt++;
                }
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

            // test expand:
            MatF result2;
            float tmparr[15] = { 
                0,0,0,0,0,
                1,0,0,1,0,
                0,0,1,0,0
            };
            MatF obj5(3,5,tmparr,1);
            MatF obj6 = obj5; // deep copy

            obj5.expand(result2, 1, 1,0,0,0,0,0,0,0);
            obj6(1,2) = 1; obj6(2,1) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);
            TS_ASSERT_DIFFERS(obj6, obj5);
            
            obj5.expand(result2, 1, 0,0,0,1,0,0,0,0);
            obj6 = obj5; // deep copy
            obj6(2,0) = 1; obj6(2,3) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            obj5.expand(result2, 1, 1,1,1,1,1,0,1,1);
            obj6 = 1;
            obj6(0,4) = 0; obj6(0,1) = 0;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            // test log
            MatF obj7(3,3,2);
            obj7.logarithm(2);
            TS_ASSERT_EQUALS(obj7(0,0),1); 
            TS_ASSERT_EQUALS(obj7(2,2),1); 

        }

        void test_MatF_min_max( void ) {
            MatF obj1(10,5,4);
            float min; float max;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 4);
            TS_ASSERT_EQUALS(max, 4);

            obj1(0,0) = 11;
            obj1(9,4) = 2;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);

            obj1(0,0) = 2;
            obj1(9,4) = 11;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);
        }

        void test_MatF_set( void ) {
            MatF obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatF obj2;
            obj2.set(obj1);
            ref_same_dat(obj1, obj2);

            float *arr = new float[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            MatF obj3;
            obj3.set(5,2,arr);
            TS_ASSERT_SAME_DATA(arr, (float*)obj3, 1);
            delete[] arr;
        }

        void test_MatF_take( void ) {
            MatF obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatF obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            MatF obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2(9,4), 2);
            TS_ASSERT_EQUALS(obj3(9,4), 2);
            TS_ASSERT_EQUALS(obj2.dim1(), 10);
            TS_ASSERT_EQUALS(obj3.dim1(), 10);
            TS_ASSERT_EQUALS(obj2.dim2(), 5);
            TS_ASSERT_EQUALS(obj3.dim2(), 5);

            // If the other one is deep too!
            MatF obj4(10,5,2);
            MatF obj5(4,2,1);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.dim1(), 10);
            TS_ASSERT_EQUALS(obj5.dim2(), 5);
            TS_ASSERT_EQUALS(obj4.dim1(), 10);
            TS_ASSERT_EQUALS(obj4.dim2(), 5);
            TS_ASSERT_EQUALS(obj5(9,4), 2);

            // Take from a shallow?
            // UNCOMMENT THIS CODE AND THE TESTING WILL ABORT ON THIS:
            // MatF obj6;
            // obj6.take(obj4);

                  // Test take(n,float *arr)
            MatF obj6;
            float *arr = new float[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(5,2,arr);
            TS_ASSERT_EQUALS(obj6(4,1), 9);
            TS_ASSERT_EQUALS(obj6(0,0), 0);

            MatF obj7(10,5,4);  // what if we are not shallow
            float *arr2 = new float[50];
            for (int i = 0; i < 50; ++i) { arr2[i] = i; }
            obj7.take(10,5,arr2);
            TS_ASSERT_EQUALS(obj7(9,4), 49);
            TS_ASSERT_EQUALS(obj7(0,0), 0);
        }

        

        void test_MatF_operator_thing_splat( void ) {
            MatF test7(10,5,4);
            float *silly = (float*)test7;
            TS_ASSERT_EQUALS(silly[0], 4);
            TS_ASSERT_EQUALS(silly[49], 4);
        }

        void test_MatF_assignment_operator( void ) {
            // First type, equals float
            MatF test8(10,5,2);
            TS_ASSERT_EQUALS(test8(0,0), 2);
            TS_ASSERT_EQUALS(test8(9,4), 2);
            test8 = 3;
            TS_ASSERT_EQUALS(test8(0,0), 3);
            TS_ASSERT_EQUALS(test8(9,4), 3);
            // Second type, equals another Mat object:
            MatF test9 = test8;
            TS_ASSERT_EQUALS(test9(0,0), 3);
            TS_ASSERT_EQUALS(test9(9,4), 3);
            TS_ASSERT_EQUALS(test9.dim1(), test8.dim1());
            ref_diff_dat(test9, test8);
        }

        void test_MatF_copy_constructor( void ) {
            MatF obj1(10,5,2);
            // Two ways to do copy constructor:
            MatF obj2 = obj1;
            TS_ASSERT_EQUALS(obj2(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj2);
            MatF obj3(obj1);  // Will use assignmen!?
            TS_ASSERT_EQUALS(obj3(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj3);
        }

        void Xtest_MatF_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            MatF testA(10000,1000, 10);
//            st.stop();
//            printf("\nTIME INIT: %fs\n", st.read());
//            st.start();
//            MatF testB;
//            testA.copy(testB);
//            st.stop();
//            printf("TIME COPY: %fs\n", st.read());

            
            // My data access vs. raw data access:
            int ml = 10000;
            int nl = 1000;
            MatF testC(ml,nl, 10);
            //float *arr = new float[ml*nl];
            float *arr = (float*)testC;

            int m;
            int n;
            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    //testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("BRUTEFORCE INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING (same as before): %fs\n", st.read());

        }

        void test_MatF_copy( void ) {
            MatF obj1(10,5,2);
            MatF obj2;
            obj1.copy(obj2);
            TS_ASSERT_EQUALS(obj1, obj2);
            TS_ASSERT_EQUALS(obj2(0,0), 2);
            obj2(0,0) = 3;
            TS_ASSERT_DIFFERS(obj1, obj2);
            MatF obj3;
            {
                MatF obj4(10,5,4);
                obj4.copy(obj3);
            }
            TS_ASSERT_EQUALS(obj3(9,4), 4);
            // The following is a compiletime error if uncommented
            //obj3[3] = 7;
        }

        void test_MatF_equals_operator( void ) {
            MatF obj1(10,5,3);
            MatF obj2 = obj1; 
            TS_ASSERT_EQUALS(obj1, obj2);
            MatF obj3(10,5,3); 
            TS_ASSERT_EQUALS(obj1, obj3);
            obj3(2,2) = 2;
            TS_ASSERT_DIFFERS(obj1, obj3);
            TS_ASSERT_EQUALS(obj1, obj1);
            MatF obj4(11,5,3); 
            TS_ASSERT_DIFFERS(obj1, obj4);
        }

        void ref_same_dat(MatF &one, MatF &two) {
            float tmp = one(0,0);
            one(0,0) = 99;
            TS_ASSERT_EQUALS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void ref_diff_dat(MatF &one, MatF &two) {
            float tmp = one(0,0);
            one(0,0) = 9999;
            TS_ASSERT_DIFFERS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void test_MatF_equal_math_ops( void ) {
            MatF obj1(10,5,3);
            MatF obj2(10,5,3);
            obj1 += obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 6);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 -= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 *= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 9);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 /= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            //std::cout << obj1;
            obj1 += 5;
            TS_ASSERT_EQUALS(obj1(0,0), 8);
            TS_ASSERT_EQUALS(obj1(9,4), 8);
            obj1 -= 5;
            TS_ASSERT_EQUALS(obj1(0,0), 3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            obj1 *= 4;
            TS_ASSERT_EQUALS(obj1(0,0), 12);
            TS_ASSERT_EQUALS(obj1(9,4), 12);
            obj1 /= 3;
            TS_ASSERT_EQUALS(obj1(0,0), 4);
            TS_ASSERT_EQUALS(obj1(9,4), 4);
        }

        void test_MatF_math_ops( void ) {
            MatF obj1(10,5,3);
            MatF obj2(10,5,3);
            MatF obj3;
            obj1.add(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            // Test if the receiver already has data!
            MatF obj4(10,5,8);
            obj1.add(obj2, obj4);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            obj1.sub(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 0);

            obj1.mul(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 9);

            obj1.div(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 1);

        }

        void test_MatF_write_and_set_from_binary( void ) {
            MatF obj1(7,3,10);
            obj1(0,1) = 5;
            obj1(1,0) = 4;
            obj1(6,2) = 3;
            obj1.write("tmp.tmp.tmp");
            
            MatF obj2;
            obj2.set_from_binary("tmp.tmp.tmp");
            TS_ASSERT_EQUALS(obj1, obj2);
            
            // test write to stdout
            //obj1.write();
        }

        void test_MatF_print( void ) {
            float *arr = new float[10];
            for (int i = 0; i < 10; ++i) {
                arr[i] = i;
            }
            MatF obj1(2,5,arr,1);
            //obj1.print();  // prints the array 0 1 2 3 4\n5 6 7 8 9
            std::ofstream fh("tmp.tmp");
            obj1.print(fh, 1);
            fh.close();
            FILE *fptr = fopen("tmp.tmp", "r");
            char line[300];
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            obj1.print("tmp.tmp.tmp",1);
            fptr = fopen("tmp.tmp.tmp", "r"); 
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            delete[] arr;

            // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        //*****************************************************************
        //*****************************************************************


        /************************************************************
         * TESTING MatD (type: double)
         ************************************************************/ 
        void test_MatD_creation_and_indexing( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(MatD test1(-1));

            // This will abort the tests if bounds checking is on!
            //MatD test2(-1);

            MatD test3(10,2);
            TS_ASSERT_EQUALS(test3.dim1(), 10);
            TS_ASSERT_EQUALS(test3.dim2(), 2);

            MatD test4;
            TS_ASSERT_EQUALS(test4.dim1(), 0);
            TS_ASSERT_EQUALS(test4.dim2(), 0);

            MatD *test5 = new MatD(0,0);
            TS_ASSERT_EQUALS(test5->dim1(), 0);
            TS_ASSERT_EQUALS(test5->dim2(), 0);
            delete test5;

            MatD test6(10,5,4);
            TS_ASSERT_EQUALS(test6.dim1(), 10);
            TS_ASSERT_EQUALS(test6.dim2(), 5);
            TS_ASSERT_EQUALS(test6(0,0), 4);
            TS_ASSERT_EQUALS(test6(9,4), 4);

            // set with double array!
            double *arr = new double[50];
            for (int i = 0; i < 50; ++i) {
                arr[i] = 9;
            }
            MatD test7(10,5,arr);
            TS_ASSERT_EQUALS(test7(0,0), 9);
            test7(0,0) = 8;
            TS_ASSERT_EQUALS(test7(0,0), 8);
            TS_ASSERT_EQUALS(test7(0,0), arr[0]);

            double *arr2 = new double[50];
            for (int i = 0; i < 50; ++i) {
                arr2[i] = 9;
            }
            MatD test8(10,5,arr2);
            TS_ASSERT_EQUALS(test8(0,0), 9);
            test8(0,0) = 8;
            TS_ASSERT_EQUALS(test8(0,0), 8);
            TS_ASSERT_EQUALS(test8(0,0), arr2[0]);
            // Notice we don't have to delete this guy, because we gave him
            // away! (have to do alloc/free test with valgrind!

            MatD test9(10,5,2);
            MatD test10(test9);
            ref_diff_dat(test9, test10);

            MatD test12(test9, 1);
            ref_same_dat(test9, test12);
        }

        void test_MatD_to_vec( void ) {
            MatD obj1(4,3,3);
            VecD out;
            obj1.to_vec(out);  // deep
            TS_ASSERT_EQUALS(out, obj1._dat);
            out[0] = 6;
            TS_ASSERT_DIFFERS(out, obj1._dat);
            VecD outshallow;
            obj1.to_vec(outshallow, 1);  // shallow
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
            outshallow[0] = 6;
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
        }

        void test_MatD_row_vecs( void ) {
            MatD obj1(10,5,3);

            // Row vecs:
            VecD rvecs[10];
            int cnt;
            obj1.row_vecs(cnt, rvecs);
            TS_ASSERT_EQUALS(cnt, 10);
            VecD rvec_answ(5,3);
            for (int i = 0; i < cnt; ++i) {
                TS_ASSERT_EQUALS(rvecs[i], rvec_answ);
            }
        }

        void test_MatD_set_from_ascii( void ) {
            MatD obj1;
            obj1.set_from_ascii("tfiles/tmp1.mata");
            TS_ASSERT_EQUALS(obj1.rows(), 6);
            TS_ASSERT_EQUALS(obj1.cols(), 10);
            TS_ASSERT_EQUALS(obj1(0,0), 230);
            TS_ASSERT_EQUALS(obj1(5,9), 657);

            MatD obj2;
            // This file is perfect
            obj2.set_from_ascii("tfiles/tmp1_no_header.mata", 1);
            TS_ASSERT_EQUALS(obj2.rows(), 6);
            TS_ASSERT_EQUALS(obj2.cols(), 10);
            TS_ASSERT_EQUALS(obj2(0,0), 230);
            TS_ASSERT_EQUALS(obj2(5,9), 657);

            MatD obj3;
            // This file has an extra space at the end of the first line
            // extra spaces at the end of some other lines
            // and an extra line at the end
            obj3.set_from_ascii("tfiles/tmp1_no_header_messy.mata", 1);
            TS_ASSERT_EQUALS(obj3.rows(), 6);
            TS_ASSERT_EQUALS(obj3.cols(), 10);
            TS_ASSERT_EQUALS(obj3(0,0), 230);
            TS_ASSERT_EQUALS(obj3(5,9), 657);
        }
        
        void test_MatD_transpose( void ) {
            MatD obj1(4,2,5);
            obj1(3,0) = 10;
            MatD trans;
            obj1.transpose(trans);
            TS_ASSERT_EQUALS(obj1.rows(), trans.cols());
            TS_ASSERT_EQUALS(obj1.cols(), trans.rows());
            TS_ASSERT_EQUALS(trans(0,0), 5);
            TS_ASSERT_EQUALS(trans(0,3), obj1(3,0));
            TS_ASSERT_EQUALS(trans(0,3), 10);
            TS_ASSERT_EQUALS(trans(1,0), 5);
        }

        void test_MatD_FUNCTIONS( void ) {
            // test avg:
            MatD obj1(2,2,2);
            obj1(0,1) = 3;
            obj1(1,0) = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            // test sum:
            double sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // sum of a row:
            // 2 3 2 2 2  =  11
            // 2 2 2 2 2  =  10
            MatD obj2(2,5,2);
            obj2(0,2) = 3;
            TS_ASSERT_EQUALS(obj2.sum(0),11);
            TS_ASSERT_EQUALS(obj2.sum(1),10);

            // test mask_as_vec:
            MatD obj3(5,2);
            MatI obj4(5,2);
            int cnt = 0;
            for (int m = 0; m < 5; ++m) {
                for (int n = 0; n < 2; ++n) {
                    obj3(m,n) = cnt;
                    obj4(m,n) = cnt;
                    cnt++;
                }
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

            // test expand:
            MatD result2;
            double tmparr[15] = { 
                0,0,0,0,0,
                1,0,0,1,0,
                0,0,1,0,0
            };
            MatD obj5(3,5,tmparr,1);
            MatD obj6 = obj5; // deep copy

            obj5.expand(result2, 1, 1,0,0,0,0,0,0,0);
            obj6(1,2) = 1; obj6(2,1) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);
            TS_ASSERT_DIFFERS(obj6, obj5);
            
            obj5.expand(result2, 1, 0,0,0,1,0,0,0,0);
            obj6 = obj5; // deep copy
            obj6(2,0) = 1; obj6(2,3) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            obj5.expand(result2, 1, 1,1,1,1,1,0,1,1);
            obj6 = 1;
            obj6(0,4) = 0; obj6(0,1) = 0;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            // test log
            MatD obj7(3,3,2);
            obj7.logarithm(2);
            TS_ASSERT_EQUALS(obj7(0,0),1); 
            TS_ASSERT_EQUALS(obj7(2,2),1); 

        }

        void test_MatD_min_max( void ) {
            MatD obj1(10,5,4);
            double min; double max;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 4);
            TS_ASSERT_EQUALS(max, 4);

            obj1(0,0) = 11;
            obj1(9,4) = 2;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);

            obj1(0,0) = 2;
            obj1(9,4) = 11;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);
        }

        void test_MatD_set( void ) {
            MatD obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatD obj2;
            obj2.set(obj1);
            ref_same_dat(obj1, obj2);

            double *arr = new double[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            MatD obj3;
            obj3.set(5,2,arr);
            TS_ASSERT_SAME_DATA(arr, (double*)obj3, 1);
            delete[] arr;
        }

        void test_MatD_take( void ) {
            MatD obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatD obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            MatD obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2(9,4), 2);
            TS_ASSERT_EQUALS(obj3(9,4), 2);
            TS_ASSERT_EQUALS(obj2.dim1(), 10);
            TS_ASSERT_EQUALS(obj3.dim1(), 10);
            TS_ASSERT_EQUALS(obj2.dim2(), 5);
            TS_ASSERT_EQUALS(obj3.dim2(), 5);

            // If the other one is deep too!
            MatD obj4(10,5,2);
            MatD obj5(4,2,1);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.dim1(), 10);
            TS_ASSERT_EQUALS(obj5.dim2(), 5);
            TS_ASSERT_EQUALS(obj4.dim1(), 10);
            TS_ASSERT_EQUALS(obj4.dim2(), 5);
            TS_ASSERT_EQUALS(obj5(9,4), 2);

            // Take from a shallow?
            // UNCOMMENT THIS CODE AND THE TESTING WILL ABORT ON THIS:
            // MatD obj6;
            // obj6.take(obj4);

                  // Test take(n,double *arr)
            MatD obj6;
            double *arr = new double[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(5,2,arr);
            TS_ASSERT_EQUALS(obj6(4,1), 9);
            TS_ASSERT_EQUALS(obj6(0,0), 0);

            MatD obj7(10,5,4);  // what if we are not shallow
            double *arr2 = new double[50];
            for (int i = 0; i < 50; ++i) { arr2[i] = i; }
            obj7.take(10,5,arr2);
            TS_ASSERT_EQUALS(obj7(9,4), 49);
            TS_ASSERT_EQUALS(obj7(0,0), 0);
        }

        

        void test_MatD_operator_thing_splat( void ) {
            MatD test7(10,5,4);
            double *silly = (double*)test7;
            TS_ASSERT_EQUALS(silly[0], 4);
            TS_ASSERT_EQUALS(silly[49], 4);
        }

        void test_MatD_assignment_operator( void ) {
            // First type, equals double
            MatD test8(10,5,2);
            TS_ASSERT_EQUALS(test8(0,0), 2);
            TS_ASSERT_EQUALS(test8(9,4), 2);
            test8 = 3;
            TS_ASSERT_EQUALS(test8(0,0), 3);
            TS_ASSERT_EQUALS(test8(9,4), 3);
            // Second type, equals another Mat object:
            MatD test9 = test8;
            TS_ASSERT_EQUALS(test9(0,0), 3);
            TS_ASSERT_EQUALS(test9(9,4), 3);
            TS_ASSERT_EQUALS(test9.dim1(), test8.dim1());
            ref_diff_dat(test9, test8);
        }

        void test_MatD_copy_constructor( void ) {
            MatD obj1(10,5,2);
            // Two ways to do copy constructor:
            MatD obj2 = obj1;
            TS_ASSERT_EQUALS(obj2(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj2);
            MatD obj3(obj1);  // Will use assignmen!?
            TS_ASSERT_EQUALS(obj3(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj3);
        }

        void Xtest_MatD_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            MatD testA(10000,1000, 10);
//            st.stop();
//            printf("\nTIME INIT: %fs\n", st.read());
//            st.start();
//            MatD testB;
//            testA.copy(testB);
//            st.stop();
//            printf("TIME COPY: %fs\n", st.read());

            
            // My data access vs. raw data access:
            int ml = 10000;
            int nl = 1000;
            MatD testC(ml,nl, 10);
            //double *arr = new double[ml*nl];
            double *arr = (double*)testC;

            int m;
            int n;
            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    //testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("BRUTEFORCE INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING (same as before): %fs\n", st.read());

        }

        void test_MatD_copy( void ) {
            MatD obj1(10,5,2);
            MatD obj2;
            obj1.copy(obj2);
            TS_ASSERT_EQUALS(obj1, obj2);
            TS_ASSERT_EQUALS(obj2(0,0), 2);
            obj2(0,0) = 3;
            TS_ASSERT_DIFFERS(obj1, obj2);
            MatD obj3;
            {
                MatD obj4(10,5,4);
                obj4.copy(obj3);
            }
            TS_ASSERT_EQUALS(obj3(9,4), 4);
            // The following is a compiletime error if uncommented
            //obj3[3] = 7;
        }

        void test_MatD_equals_operator( void ) {
            MatD obj1(10,5,3);
            MatD obj2 = obj1; 
            TS_ASSERT_EQUALS(obj1, obj2);
            MatD obj3(10,5,3); 
            TS_ASSERT_EQUALS(obj1, obj3);
            obj3(2,2) = 2;
            TS_ASSERT_DIFFERS(obj1, obj3);
            TS_ASSERT_EQUALS(obj1, obj1);
            MatD obj4(11,5,3); 
            TS_ASSERT_DIFFERS(obj1, obj4);
        }

        void ref_same_dat(MatD &one, MatD &two) {
            double tmp = one(0,0);
            one(0,0) = 99;
            TS_ASSERT_EQUALS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void ref_diff_dat(MatD &one, MatD &two) {
            double tmp = one(0,0);
            one(0,0) = 9999;
            TS_ASSERT_DIFFERS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void test_MatD_equal_math_ops( void ) {
            MatD obj1(10,5,3);
            MatD obj2(10,5,3);
            obj1 += obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 6);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 -= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 *= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 9);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 /= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            //std::cout << obj1;
            obj1 += 5;
            TS_ASSERT_EQUALS(obj1(0,0), 8);
            TS_ASSERT_EQUALS(obj1(9,4), 8);
            obj1 -= 5;
            TS_ASSERT_EQUALS(obj1(0,0), 3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            obj1 *= 4;
            TS_ASSERT_EQUALS(obj1(0,0), 12);
            TS_ASSERT_EQUALS(obj1(9,4), 12);
            obj1 /= 3;
            TS_ASSERT_EQUALS(obj1(0,0), 4);
            TS_ASSERT_EQUALS(obj1(9,4), 4);
        }

        void test_MatD_math_ops( void ) {
            MatD obj1(10,5,3);
            MatD obj2(10,5,3);
            MatD obj3;
            obj1.add(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            // Test if the receiver already has data!
            MatD obj4(10,5,8);
            obj1.add(obj2, obj4);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            obj1.sub(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 0);

            obj1.mul(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 9);

            obj1.div(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 1);

        }

        void test_MatD_write_and_set_from_binary( void ) {
            MatD obj1(7,3,10);
            obj1(0,1) = 5;
            obj1(1,0) = 4;
            obj1(6,2) = 3;
            obj1.write("tmp.tmp.tmp");
            
            MatD obj2;
            obj2.set_from_binary("tmp.tmp.tmp");
            TS_ASSERT_EQUALS(obj1, obj2);
            
            // test write to stdout
            //obj1.write();
        }

        void test_MatD_print( void ) {
            double *arr = new double[10];
            for (int i = 0; i < 10; ++i) {
                arr[i] = i;
            }
            MatD obj1(2,5,arr,1);
            //obj1.print();  // prints the array 0 1 2 3 4\n5 6 7 8 9
            std::ofstream fh("tmp.tmp");
            obj1.print(fh, 1);
            fh.close();
            FILE *fptr = fopen("tmp.tmp", "r");
            char line[300];
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            obj1.print("tmp.tmp.tmp",1);
            fptr = fopen("tmp.tmp.tmp", "r"); 
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            delete[] arr;

            // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        //*****************************************************************
        //*****************************************************************


        /************************************************************
         * TESTING MatI (type: int)
         ************************************************************/ 
        void test_MatI_creation_and_indexing( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(MatI test1(-1));

            // This will abort the tests if bounds checking is on!
            //MatI test2(-1);

            MatI test3(10,2);
            TS_ASSERT_EQUALS(test3.dim1(), 10);
            TS_ASSERT_EQUALS(test3.dim2(), 2);

            MatI test4;
            TS_ASSERT_EQUALS(test4.dim1(), 0);
            TS_ASSERT_EQUALS(test4.dim2(), 0);

            MatI *test5 = new MatI(0,0);
            TS_ASSERT_EQUALS(test5->dim1(), 0);
            TS_ASSERT_EQUALS(test5->dim2(), 0);
            delete test5;

            MatI test6(10,5,4);
            TS_ASSERT_EQUALS(test6.dim1(), 10);
            TS_ASSERT_EQUALS(test6.dim2(), 5);
            TS_ASSERT_EQUALS(test6(0,0), 4);
            TS_ASSERT_EQUALS(test6(9,4), 4);

            // set with int array!
            int *arr = new int[50];
            for (int i = 0; i < 50; ++i) {
                arr[i] = 9;
            }
            MatI test7(10,5,arr);
            TS_ASSERT_EQUALS(test7(0,0), 9);
            test7(0,0) = 8;
            TS_ASSERT_EQUALS(test7(0,0), 8);
            TS_ASSERT_EQUALS(test7(0,0), arr[0]);

            int *arr2 = new int[50];
            for (int i = 0; i < 50; ++i) {
                arr2[i] = 9;
            }
            MatI test8(10,5,arr2);
            TS_ASSERT_EQUALS(test8(0,0), 9);
            test8(0,0) = 8;
            TS_ASSERT_EQUALS(test8(0,0), 8);
            TS_ASSERT_EQUALS(test8(0,0), arr2[0]);
            // Notice we don't have to delete this guy, because we gave him
            // away! (have to do alloc/free test with valgrind!

            MatI test9(10,5,2);
            MatI test10(test9);
            ref_diff_dat(test9, test10);

            MatI test12(test9, 1);
            ref_same_dat(test9, test12);
        }

        void test_MatI_to_vec( void ) {
            MatI obj1(4,3,3);
            VecI out;
            obj1.to_vec(out);  // deep
            TS_ASSERT_EQUALS(out, obj1._dat);
            out[0] = 6;
            TS_ASSERT_DIFFERS(out, obj1._dat);
            VecI outshallow;
            obj1.to_vec(outshallow, 1);  // shallow
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
            outshallow[0] = 6;
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
        }

        void test_MatI_row_vecs( void ) {
            MatI obj1(10,5,3);

            // Row vecs:
            VecI rvecs[10];
            int cnt;
            obj1.row_vecs(cnt, rvecs);
            TS_ASSERT_EQUALS(cnt, 10);
            VecI rvec_answ(5,3);
            for (int i = 0; i < cnt; ++i) {
                TS_ASSERT_EQUALS(rvecs[i], rvec_answ);
            }
        }

        void test_MatI_set_from_ascii( void ) {
            MatI obj1;
            obj1.set_from_ascii("tfiles/tmp1.mata");
            TS_ASSERT_EQUALS(obj1.rows(), 6);
            TS_ASSERT_EQUALS(obj1.cols(), 10);
            TS_ASSERT_EQUALS(obj1(0,0), 230);
            TS_ASSERT_EQUALS(obj1(5,9), 657);

            MatI obj2;
            // This file is perfect
            obj2.set_from_ascii("tfiles/tmp1_no_header.mata", 1);
            TS_ASSERT_EQUALS(obj2.rows(), 6);
            TS_ASSERT_EQUALS(obj2.cols(), 10);
            TS_ASSERT_EQUALS(obj2(0,0), 230);
            TS_ASSERT_EQUALS(obj2(5,9), 657);

            MatI obj3;
            // This file has an extra space at the end of the first line
            // extra spaces at the end of some other lines
            // and an extra line at the end
            obj3.set_from_ascii("tfiles/tmp1_no_header_messy.mata", 1);
            TS_ASSERT_EQUALS(obj3.rows(), 6);
            TS_ASSERT_EQUALS(obj3.cols(), 10);
            TS_ASSERT_EQUALS(obj3(0,0), 230);
            TS_ASSERT_EQUALS(obj3(5,9), 657);
        }
        
        void test_MatI_transpose( void ) {
            MatI obj1(4,2,5);
            obj1(3,0) = 10;
            MatI trans;
            obj1.transpose(trans);
            TS_ASSERT_EQUALS(obj1.rows(), trans.cols());
            TS_ASSERT_EQUALS(obj1.cols(), trans.rows());
            TS_ASSERT_EQUALS(trans(0,0), 5);
            TS_ASSERT_EQUALS(trans(0,3), obj1(3,0));
            TS_ASSERT_EQUALS(trans(0,3), 10);
            TS_ASSERT_EQUALS(trans(1,0), 5);
        }

        void test_MatI_FUNCTIONS( void ) {
            // test avg:
            MatI obj1(2,2,2);
            obj1(0,1) = 3;
            obj1(1,0) = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            // test sum:
            int sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // sum of a row:
            // 2 3 2 2 2  =  11
            // 2 2 2 2 2  =  10
            MatI obj2(2,5,2);
            obj2(0,2) = 3;
            TS_ASSERT_EQUALS(obj2.sum(0),11);
            TS_ASSERT_EQUALS(obj2.sum(1),10);

            // test mask_as_vec:
            MatI obj3(5,2);
            MatI obj4(5,2);
            int cnt = 0;
            for (int m = 0; m < 5; ++m) {
                for (int n = 0; n < 2; ++n) {
                    obj3(m,n) = cnt;
                    obj4(m,n) = cnt;
                    cnt++;
                }
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

            // test expand:
            MatI result2;
            int tmparr[15] = { 
                0,0,0,0,0,
                1,0,0,1,0,
                0,0,1,0,0
            };
            MatI obj5(3,5,tmparr,1);
            MatI obj6 = obj5; // deep copy

            obj5.expand(result2, 1, 1,0,0,0,0,0,0,0);
            obj6(1,2) = 1; obj6(2,1) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);
            TS_ASSERT_DIFFERS(obj6, obj5);
            
            obj5.expand(result2, 1, 0,0,0,1,0,0,0,0);
            obj6 = obj5; // deep copy
            obj6(2,0) = 1; obj6(2,3) = 1;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            obj5.expand(result2, 1, 1,1,1,1,1,0,1,1);
            obj6 = 1;
            obj6(0,4) = 0; obj6(0,1) = 0;
            TS_ASSERT_EQUALS(result2, obj6);
            TS_ASSERT_DIFFERS(result2, obj5);

            // test log
            MatI obj7(3,3,2);
            obj7.logarithm(2);
            TS_ASSERT_EQUALS(obj7(0,0),1); 
            TS_ASSERT_EQUALS(obj7(2,2),1); 

        }

        void test_MatI_min_max( void ) {
            MatI obj1(10,5,4);
            int min; int max;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 4);
            TS_ASSERT_EQUALS(max, 4);

            obj1(0,0) = 11;
            obj1(9,4) = 2;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);

            obj1(0,0) = 2;
            obj1(9,4) = 11;
            obj1.min_max(min,max);
            TS_ASSERT_EQUALS(min, 2);
            TS_ASSERT_EQUALS(max, 11);
        }

        void test_MatI_set( void ) {
            MatI obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatI obj2;
            obj2.set(obj1);
            ref_same_dat(obj1, obj2);

            int *arr = new int[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            MatI obj3;
            obj3.set(5,2,arr);
            TS_ASSERT_SAME_DATA(arr, (int*)obj3, 1);
            delete[] arr;
        }

        void test_MatI_take( void ) {
            MatI obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatI obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            MatI obj3;
            obj3.take(obj2);
            TS_ASSERT(obj2.shallow());
            TS_ASSERT(!obj3.shallow());
            TS_ASSERT_EQUALS(obj2(9,4), 2);
            TS_ASSERT_EQUALS(obj3(9,4), 2);
            TS_ASSERT_EQUALS(obj2.dim1(), 10);
            TS_ASSERT_EQUALS(obj3.dim1(), 10);
            TS_ASSERT_EQUALS(obj2.dim2(), 5);
            TS_ASSERT_EQUALS(obj3.dim2(), 5);

            // If the other one is deep too!
            MatI obj4(10,5,2);
            MatI obj5(4,2,1);
            obj5.take(obj4);
            TS_ASSERT(obj4.shallow());
            TS_ASSERT(!obj5.shallow());
            TS_ASSERT_EQUALS(obj5.dim1(), 10);
            TS_ASSERT_EQUALS(obj5.dim2(), 5);
            TS_ASSERT_EQUALS(obj4.dim1(), 10);
            TS_ASSERT_EQUALS(obj4.dim2(), 5);
            TS_ASSERT_EQUALS(obj5(9,4), 2);

            // Take from a shallow?
            // UNCOMMENT THIS CODE AND THE TESTING WILL ABORT ON THIS:
            // MatI obj6;
            // obj6.take(obj4);

                  // Test take(n,int *arr)
            MatI obj6;
            int *arr = new int[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(5,2,arr);
            TS_ASSERT_EQUALS(obj6(4,1), 9);
            TS_ASSERT_EQUALS(obj6(0,0), 0);

            MatI obj7(10,5,4);  // what if we are not shallow
            int *arr2 = new int[50];
            for (int i = 0; i < 50; ++i) { arr2[i] = i; }
            obj7.take(10,5,arr2);
            TS_ASSERT_EQUALS(obj7(9,4), 49);
            TS_ASSERT_EQUALS(obj7(0,0), 0);
        }

        

        void test_MatI_operator_thing_splat( void ) {
            MatI test7(10,5,4);
            int *silly = (int*)test7;
            TS_ASSERT_EQUALS(silly[0], 4);
            TS_ASSERT_EQUALS(silly[49], 4);
        }

        void test_MatI_assignment_operator( void ) {
            // First type, equals int
            MatI test8(10,5,2);
            TS_ASSERT_EQUALS(test8(0,0), 2);
            TS_ASSERT_EQUALS(test8(9,4), 2);
            test8 = 3;
            TS_ASSERT_EQUALS(test8(0,0), 3);
            TS_ASSERT_EQUALS(test8(9,4), 3);
            // Second type, equals another Mat object:
            MatI test9 = test8;
            TS_ASSERT_EQUALS(test9(0,0), 3);
            TS_ASSERT_EQUALS(test9(9,4), 3);
            TS_ASSERT_EQUALS(test9.dim1(), test8.dim1());
            ref_diff_dat(test9, test8);
        }

        void test_MatI_copy_constructor( void ) {
            MatI obj1(10,5,2);
            // Two ways to do copy constructor:
            MatI obj2 = obj1;
            TS_ASSERT_EQUALS(obj2(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj2);
            MatI obj3(obj1);  // Will use assignmen!?
            TS_ASSERT_EQUALS(obj3(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj3);
        }

        void Xtest_MatI_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            MatI testA(10000,1000, 10);
//            st.stop();
//            printf("\nTIME INIT: %fs\n", st.read());
//            st.start();
//            MatI testB;
//            testA.copy(testB);
//            st.stop();
//            printf("TIME COPY: %fs\n", st.read());

            
            // My data access vs. raw data access:
            int ml = 10000;
            int nl = 1000;
            MatI testC(ml,nl, 10);
            //int *arr = new int[ml*nl];
            int *arr = (int*)testC;

            int m;
            int n;
            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    //testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("BRUTEFORCE INDEXING: %fs\n", st.read());

            st.start();
            for (m = 0; m < testC.mlen(); ++m) {
                for (n = 0; n < testC.nlen(); ++n) {
                    // Uses less memory here, but...
                    testC(m,n) = 5;  //Avg: 4.6
                    // More memory appears to be used here, but...
                    //arr[(m*nl) + n] = 5;  //Avg: 3.2
                }
            }
            st.stop();
            printf("MY INDEXING (same as before): %fs\n", st.read());

        }

        void test_MatI_copy( void ) {
            MatI obj1(10,5,2);
            MatI obj2;
            obj1.copy(obj2);
            TS_ASSERT_EQUALS(obj1, obj2);
            TS_ASSERT_EQUALS(obj2(0,0), 2);
            obj2(0,0) = 3;
            TS_ASSERT_DIFFERS(obj1, obj2);
            MatI obj3;
            {
                MatI obj4(10,5,4);
                obj4.copy(obj3);
            }
            TS_ASSERT_EQUALS(obj3(9,4), 4);
            // The following is a compiletime error if uncommented
            //obj3[3] = 7;
        }

        void test_MatI_equals_operator( void ) {
            MatI obj1(10,5,3);
            MatI obj2 = obj1; 
            TS_ASSERT_EQUALS(obj1, obj2);
            MatI obj3(10,5,3); 
            TS_ASSERT_EQUALS(obj1, obj3);
            obj3(2,2) = 2;
            TS_ASSERT_DIFFERS(obj1, obj3);
            TS_ASSERT_EQUALS(obj1, obj1);
            MatI obj4(11,5,3); 
            TS_ASSERT_DIFFERS(obj1, obj4);
        }

        void ref_same_dat(MatI &one, MatI &two) {
            int tmp = one(0,0);
            one(0,0) = 99;
            TS_ASSERT_EQUALS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void ref_diff_dat(MatI &one, MatI &two) {
            int tmp = one(0,0);
            one(0,0) = 9999;
            TS_ASSERT_DIFFERS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void test_MatI_equal_math_ops( void ) {
            MatI obj1(10,5,3);
            MatI obj2(10,5,3);
            obj1 += obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 6);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 -= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 *= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 9);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            obj1 /= obj2;
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            //std::cout << obj1;
            obj1 += 5;
            TS_ASSERT_EQUALS(obj1(0,0), 8);
            TS_ASSERT_EQUALS(obj1(9,4), 8);
            obj1 -= 5;
            TS_ASSERT_EQUALS(obj1(0,0), 3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            obj1 *= 4;
            TS_ASSERT_EQUALS(obj1(0,0), 12);
            TS_ASSERT_EQUALS(obj1(9,4), 12);
            obj1 /= 3;
            TS_ASSERT_EQUALS(obj1(0,0), 4);
            TS_ASSERT_EQUALS(obj1(9,4), 4);
        }

        void test_MatI_math_ops( void ) {
            MatI obj1(10,5,3);
            MatI obj2(10,5,3);
            MatI obj3;
            obj1.add(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            // Test if the receiver already has data!
            MatI obj4(10,5,8);
            obj1.add(obj2, obj4);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            obj1.sub(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 0);

            obj1.mul(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 9);

            obj1.div(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 1);

        }

        void test_MatI_write_and_set_from_binary( void ) {
            MatI obj1(7,3,10);
            obj1(0,1) = 5;
            obj1(1,0) = 4;
            obj1(6,2) = 3;
            obj1.write("tmp.tmp.tmp");
            
            MatI obj2;
            obj2.set_from_binary("tmp.tmp.tmp");
            TS_ASSERT_EQUALS(obj1, obj2);
            
            // test write to stdout
            //obj1.write();
        }

        void test_MatI_print( void ) {
            int *arr = new int[10];
            for (int i = 0; i < 10; ++i) {
                arr[i] = i;
            }
            MatI obj1(2,5,arr,1);
            //obj1.print();  // prints the array 0 1 2 3 4\n5 6 7 8 9
            std::ofstream fh("tmp.tmp");
            obj1.print(fh, 1);
            fh.close();
            FILE *fptr = fopen("tmp.tmp", "r");
            char line[300];
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            obj1.print("tmp.tmp.tmp",1);
            fptr = fopen("tmp.tmp.tmp", "r"); 
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("0 1 2 3 4",line, 9);
            fgets(line,300,fptr);
            //printf("line: %s\n", line);
            TS_ASSERT_SAME_DATA("5 6 7 8 9",line, 9);
            fclose(fptr);

            delete[] arr;

            // undo these to verify printing is OK
            remove("tmp.tmp");
            remove("tmp.tmp.tmp");
        }

        //*****************************************************************
        //*****************************************************************

// END TEMPLATE

};


