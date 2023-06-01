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
         * TESTING MatABR (type: FLOAT)
         ************************************************************/ 
        void test_MatABR_creation_and_indexing( void ) {
            // This will throw if bounds check is off:
            // Otherwise it will stop the tests!
            //TS_ASSERT_THROWS_ANYTHING(MatABR test1(-1));

            // This will abort the tests if bounds checking is on!
            //MatABR test2(-1);

            MatABR test3(10,2);
            TS_ASSERT_EQUALS(test3.dim1(), 10);
            TS_ASSERT_EQUALS(test3.dim2(), 2);

            MatABR test4;
            TS_ASSERT_EQUALS(test4.dim1(), 0);
            TS_ASSERT_EQUALS(test4.dim2(), 0);

            MatABR *test5 = new MatABR(0,0);
            TS_ASSERT_EQUALS(test5->dim1(), 0);
            TS_ASSERT_EQUALS(test5->dim2(), 0);
            delete test5;

            MatABR test6(10,5,4);
            TS_ASSERT_EQUALS(test6.dim1(), 10);
            TS_ASSERT_EQUALS(test6.dim2(), 5);
            TS_ASSERT_EQUALS(test6(0,0), 4);
            TS_ASSERT_EQUALS(test6(9,4), 4);

            // set with FLOAT array!
            FLOAT *arr = new FLOAT[50];
            for (int i = 0; i < 50; ++i) {
                arr[i] = 9;
            }
            MatABR test7(10,5,arr);
            TS_ASSERT_EQUALS(test7(0,0), 9);
            test7(0,0) = 8;
            TS_ASSERT_EQUALS(test7(0,0), 8);
            TS_ASSERT_EQUALS(test7(0,0), arr[0]);

            FLOAT *arr2 = new FLOAT[50];
            for (int i = 0; i < 50; ++i) {
                arr2[i] = 9;
            }
            MatABR test8(10,5,arr2);
            TS_ASSERT_EQUALS(test8(0,0), 9);
            test8(0,0) = 8;
            TS_ASSERT_EQUALS(test8(0,0), 8);
            TS_ASSERT_EQUALS(test8(0,0), arr2[0]);
            // Notice we don't have to delete this guy, because we gave him
            // away! (have to do alloc/free test with valgrind!

            MatABR test9(10,5,2);
            MatABR test10(test9);
            ref_diff_dat(test9, test10);

            MatABR test12(test9, 1);
            ref_same_dat(test9, test12);
        }

        void test_MatABR_to_vec( void ) {
            MatABR obj1(4,3,3);
            VecABR out;
            obj1.to_vec(out);  // deep
            TS_ASSERT_EQUALS(out, obj1._dat);
            out[0] = 6;
            TS_ASSERT_DIFFERS(out, obj1._dat);
            VecABR outshallow;
            obj1.to_vec(outshallow, 1);  // shallow
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
            outshallow[0] = 6;
            TS_ASSERT_EQUALS(outshallow, obj1._dat);
        }

        void test_MatABR_row_vecs( void ) {
            MatABR obj1(10,5,3);

            // Row vecs:
            VecABR rvecs[10];
            int cnt;
            obj1.row_vecs(cnt, rvecs);
            TS_ASSERT_EQUALS(cnt, 10);
            VecABR rvec_answ(5,3);
            for (int i = 0; i < cnt; ++i) {
                TS_ASSERT_EQUALS(rvecs[i], rvec_answ);
            }
        }

        void test_MatABR_set_from_ascii( void ) {
            MatABR obj1;
            obj1.set_from_ascii("tfiles/tmp1.mata");
            TS_ASSERT_EQUALS(obj1.rows(), 6);
            TS_ASSERT_EQUALS(obj1.cols(), 10);
            TS_ASSERT_EQUALS(obj1(0,0), 230);
            TS_ASSERT_EQUALS(obj1(5,9), 657);

            MatABR obj2;
            // This file is perfect
            obj2.set_from_ascii("tfiles/tmp1_no_header.mata", 1);
            TS_ASSERT_EQUALS(obj2.rows(), 6);
            TS_ASSERT_EQUALS(obj2.cols(), 10);
            TS_ASSERT_EQUALS(obj2(0,0), 230);
            TS_ASSERT_EQUALS(obj2(5,9), 657);

            MatABR obj3;
            // This file has an extra space at the end of the first line
            // extra spaces at the end of some other lines
            // and an extra line at the end
            obj3.set_from_ascii("tfiles/tmp1_no_header_messy.mata", 1);
            TS_ASSERT_EQUALS(obj3.rows(), 6);
            TS_ASSERT_EQUALS(obj3.cols(), 10);
            TS_ASSERT_EQUALS(obj3(0,0), 230);
            TS_ASSERT_EQUALS(obj3(5,9), 657);
        }
        
        void test_MatABR_transpose( void ) {
            MatABR obj1(4,2,5);
            obj1(3,0) = 10;
            MatABR trans;
            obj1.transpose(trans);
            TS_ASSERT_EQUALS(obj1.rows(), trans.cols());
            TS_ASSERT_EQUALS(obj1.cols(), trans.rows());
            TS_ASSERT_EQUALS(trans(0,0), 5);
            TS_ASSERT_EQUALS(trans(0,3), obj1(3,0));
            TS_ASSERT_EQUALS(trans(0,3), 10);
            TS_ASSERT_EQUALS(trans(1,0), 5);
        }

        void test_MatABR_FUNCTIONS( void ) {
            // test avg:
            MatABR obj1(2,2,2);
            obj1(0,1) = 3;
            obj1(1,0) = 3;
            double average = obj1.avg();
            TS_ASSERT_EQUALS(average, 2.5);

            // test sum:
            FLOAT sm = obj1.sum();
            TS_ASSERT_EQUALS(sm, 10);

            // sum of a row:
            // 2 3 2 2 2  =  11
            // 2 2 2 2 2  =  10
            MatABR obj2(2,5,2);
            obj2(0,2) = 3;
            TS_ASSERT_EQUALS(obj2.sum(0),11);
            TS_ASSERT_EQUALS(obj2.sum(1),10);

            // test mask_as_vec:
            MatABR obj3(5,2);
            MatI obj4(5,2);
            int cnt = 0;
            for (int m = 0; m < 5; ++m) {
                for (int n = 0; n < 2; ++n) {
                    obj3(m,n) = cnt;
                    obj4(m,n) = cnt;
                    cnt++;
                }
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

            // test expand:
            MatABR result2;
            FLOAT tmparr[15] = { 
                0,0,0,0,0,
                1,0,0,1,0,
                0,0,1,0,0
            };
            MatABR obj5(3,5,tmparr,1);
            MatABR obj6 = obj5; // deep copy

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
            MatABR obj7(3,3,2);
            obj7.logarithm(2);
            TS_ASSERT_EQUALS(obj7(0,0),1); 
            TS_ASSERT_EQUALS(obj7(2,2),1); 

        }

        void test_MatABR_min_max( void ) {
            MatABR obj1(10,5,4);
            FLOAT min; FLOAT max;
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

        void test_MatABR_set( void ) {
            MatABR obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatABR obj2;
            obj2.set(obj1);
            ref_same_dat(obj1, obj2);

            FLOAT *arr = new FLOAT[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            MatABR obj3;
            obj3.set(5,2,arr);
            TS_ASSERT_SAME_DATA(arr, (FLOAT*)obj3, 1);
            delete[] arr;
        }

        void test_MatABR_take( void ) {
            MatABR obj1(10,5,2); 
            TS_ASSERT(!obj1.shallow());
            MatABR obj2;
            obj2.take(obj1);
            TS_ASSERT(obj1.shallow());
            TS_ASSERT(!obj2.shallow());
            MatABR obj3;
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
            MatABR obj4(10,5,2);
            MatABR obj5(4,2,1);
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
            // MatABR obj6;
            // obj6.take(obj4);

                  // Test take(n,FLOAT *arr)
            MatABR obj6;
            FLOAT *arr = new FLOAT[10];
            for (int i = 0; i < 10; ++i) { arr[i] = i; }
            obj6.take(5,2,arr);
            TS_ASSERT_EQUALS(obj6(4,1), 9);
            TS_ASSERT_EQUALS(obj6(0,0), 0);

            MatABR obj7(10,5,4);  // what if we are not shallow
            FLOAT *arr2 = new FLOAT[50];
            for (int i = 0; i < 50; ++i) { arr2[i] = i; }
            obj7.take(10,5,arr2);
            TS_ASSERT_EQUALS(obj7(9,4), 49);
            TS_ASSERT_EQUALS(obj7(0,0), 0);
        }

        

        void test_MatABR_operator_thing_splat( void ) {
            MatABR test7(10,5,4);
            FLOAT *silly = (FLOAT*)test7;
            TS_ASSERT_EQUALS(silly[0], 4);
            TS_ASSERT_EQUALS(silly[49], 4);
        }

        void test_MatABR_assignment_operator( void ) {
            // First type, equals FLOAT
            MatABR test8(10,5,2);
            TS_ASSERT_EQUALS(test8(0,0), 2);
            TS_ASSERT_EQUALS(test8(9,4), 2);
            test8 = 3;
            TS_ASSERT_EQUALS(test8(0,0), 3);
            TS_ASSERT_EQUALS(test8(9,4), 3);
            // Second type, equals another Mat object:
            MatABR test9 = test8;
            TS_ASSERT_EQUALS(test9(0,0), 3);
            TS_ASSERT_EQUALS(test9(9,4), 3);
            TS_ASSERT_EQUALS(test9.dim1(), test8.dim1());
            ref_diff_dat(test9, test8);
        }

        void test_MatABR_copy_constructor( void ) {
            MatABR obj1(10,5,2);
            // Two ways to do copy constructor:
            MatABR obj2 = obj1;
            TS_ASSERT_EQUALS(obj2(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj2);
            MatABR obj3(obj1);  // Will use assignmen!?
            TS_ASSERT_EQUALS(obj3(0,0), obj1(0,0));
            ref_diff_dat(obj1, obj3);
        }

        void Xtest_MatABR_speed ( void ) {
            TNT::Stopwatch st;
//            st.start();
//            MatABR testA(10000,1000, 10);
//            st.stop();
//            printf("\nTIME INIT: %fs\n", st.read());
//            st.start();
//            MatABR testB;
//            testA.copy(testB);
//            st.stop();
//            printf("TIME COPY: %fs\n", st.read());

            
            // My data access vs. raw data access:
            int ml = 10000;
            int nl = 1000;
            MatABR testC(ml,nl, 10);
            //FLOAT *arr = new FLOAT[ml*nl];
            FLOAT *arr = (FLOAT*)testC;

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

        void test_MatABR_copy( void ) {
            MatABR obj1(10,5,2);
            MatABR obj2;
            obj1.copy(obj2);
            TS_ASSERT_EQUALS(obj1, obj2);
            TS_ASSERT_EQUALS(obj2(0,0), 2);
            obj2(0,0) = 3;
            TS_ASSERT_DIFFERS(obj1, obj2);
            MatABR obj3;
            {
                MatABR obj4(10,5,4);
                obj4.copy(obj3);
            }
            TS_ASSERT_EQUALS(obj3(9,4), 4);
            // The following is a compiletime error if uncommented
            //obj3[3] = 7;
        }

        void test_MatABR_equals_operator( void ) {
            MatABR obj1(10,5,3);
            MatABR obj2 = obj1; 
            TS_ASSERT_EQUALS(obj1, obj2);
            MatABR obj3(10,5,3); 
            TS_ASSERT_EQUALS(obj1, obj3);
            obj3(2,2) = 2;
            TS_ASSERT_DIFFERS(obj1, obj3);
            TS_ASSERT_EQUALS(obj1, obj1);
            MatABR obj4(11,5,3); 
            TS_ASSERT_DIFFERS(obj1, obj4);
        }

        void ref_same_dat(MatABR &one, MatABR &two) {
            FLOAT tmp = one(0,0);
            one(0,0) = 99;
            TS_ASSERT_EQUALS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void ref_diff_dat(MatABR &one, MatABR &two) {
            FLOAT tmp = one(0,0);
            one(0,0) = 9999;
            TS_ASSERT_DIFFERS(one(0,0), two(0,0));
            one(0,0) = tmp;
        }

        void test_MatABR_equal_math_ops( void ) {
            MatABR obj1(10,5,3);
            MatABR obj2(10,5,3);
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

        void test_MatABR_math_ops( void ) {
            MatABR obj1(10,5,3);
            MatABR obj2(10,5,3);
            MatABR obj3;
            obj1.add(obj2, obj3);
            TS_ASSERT_EQUALS(obj1(9,4), 3);
            TS_ASSERT_EQUALS(obj2(9,4), 3);
            TS_ASSERT_EQUALS(obj3(9,4), 6);

            // Test if the receiver already has data!
            MatABR obj4(10,5,8);
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

        void test_MatABR_write_and_set_from_binary( void ) {
            MatABR obj1(7,3,10);
            obj1(0,1) = 5;
            obj1(1,0) = 4;
            obj1(6,2) = 3;
            obj1.write("tmp.tmp.tmp");
            
            MatABR obj2;
            obj2.set_from_binary("tmp.tmp.tmp");
            TS_ASSERT_EQUALS(obj1, obj2);
            
            // test write to stdout
            //obj1.write();
        }

        void test_MatABR_print( void ) {
            FLOAT *arr = new FLOAT[10];
            for (int i = 0; i < 10; ++i) {
                arr[i] = i;
            }
            MatABR obj1(2,5,arr,1);
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


