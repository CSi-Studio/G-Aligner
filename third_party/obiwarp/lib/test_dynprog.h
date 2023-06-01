#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include "dynprog.h"

class DynProgTestSuite : public CxxTest::TestSuite 
{
    public:
        
        void test_creation( void ) {
            DynProg dyn;
            TS_ASSERT(3 == 3);
        }

        void test_score( void ) {

            float data1[] = { 
                10.0, 12.3, 14.3, 12.0, 10.1,    // mean 11.74
                12.3, 60.3, 20.0, 12.0, 23.0,    // mean 25.52 
                12.3, 20.0, 70.2, 12.0, 23.0,    // mean 27.5
            };

            float data2[] = { 
                23.2, 34.2, 70.0, 34.0, 80.0,    // mean 48.28
                12.3, 60.3, 20.0, 12.0, 27.0,    // mean 26.32
            };

            // product answer = 
            // 2869.66 5995.62 etc...
            // 1567.39 etc..

            MatF scores;

            MatF trialm(3,5,data1,1);  
            MatF trialn(2,5,data2,1); 
            DynProg dyn;
            dyn.score_product(trialm, trialn, scores); 
            TS_ASSERT_DELTA(scores(0,0), 2869.66, 1.00001);
            TS_ASSERT_DELTA(scores(0,1), 1567.39, 1.00001);
            TS_ASSERT_DELTA(scores(1,0), 5995.62, 1.00001);

            dyn.score_covariance(trialm, trialn, scores); 
            TS_ASSERT_DELTA(scores(0,0), 7.1248, 1.00001);
            TS_ASSERT_DELTA(scores(0,1), 4.4812, 1.00001);
            TS_ASSERT_DELTA(scores(1,0), -32.9816, 1.00001);
            //TS_ASSERT_DELTA(scores(0,0), 35.62, 1.00001);
            //TS_ASSERT_DELTA(scores(0,1), 22.4, 1.00001);
            //TS_ASSERT_DELTA(scores(1,0), -164.9, 1.00001);

            dyn.score_pearsons_r(trialm, trialn, scores); 
            TS_ASSERT_DELTA(scores(0,0), 0.19994, 0.00001);
            TS_ASSERT_DELTA(scores(0,1), 0.15764, 0.00001);
            TS_ASSERT_DELTA(scores(1,0), -0.0822, 0.00001);

            dyn.score_pearsons_r2(trialm, trialn, scores); 
            TS_ASSERT_DELTA(scores(0,0), 0.03997, 0.0001);
            TS_ASSERT_DELTA(scores(0,1), 0.02485, 0.001);
            TS_ASSERT_DELTA(scores(1,0), 0.006757, 0.001);

            dyn.score_euclidean(trialm, trialn, scores); 
            TS_ASSERT_DELTA(scores(0,0), 95.5319, 0.0001);
            TS_ASSERT_DELTA(scores(0,1), 51.2580725, 0.00001);
            TS_ASSERT_DELTA(scores(1,0), 83.86310, 0.0001);
            // TEST MUTUAL INFO:
            // Don't have any exhaustive test for this, but 
            // I think it is right...
            // Some of it can be tested...
            float data3[] = { 
                0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
                1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
                0.f, 1.f, 0.f, 1.f, 0.f, 1.f,
                1.f, 1.f, 0.f, 1.f, 0.f, 1.f,
                1.f, 1.f, 0.f, 1.f, 0.f, 1.f
            };
            float data4[] = { 
                0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
                1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
                0.f, 1.f, 0.f, 1.f, 0.f, 1.f,
                0.f, 1.f, 0.f, 1.f, 0.f, 0.f
            }; 
            MatF trialMut1(5,6,data3,1);  // MAT1
            MatF trialMut2(4,6,data4,1);  // MAT1
            MatF scoresMut;  // MAT1

            dyn.score_mutual_info(trialMut1, trialMut2, scoresMut, 2);
            TS_ASSERT_EQUALS(scoresMut(1,0),0.f);
            TS_ASSERT_DELTA(scoresMut(3,0),0.f, 0.0000001);
            TS_ASSERT_DELTA(scoresMut(0,3),0.f, 0.0000001);
            TS_ASSERT_EQUALS(scoresMut(1,1),0.f);
            TS_ASSERT_EQUALS(scoresMut(2,2),1);
            TS_ASSERT_DELTA(scoresMut(3,3),0.2516,0.001);
            TS_ASSERT_DELTA(scoresMut(4,3),0.2516,0.001);
            TS_ASSERT_DELTA(scoresMut(3,2),0.4591,0.001);
            TS_ASSERT_DELTA(scoresMut(2,3),0.4591,0.001);
        }

        void test_linear_less_before( void ) {

            // Test less before
            float data1[] = { 10.0f, 12.0f, 15.0f, 16.0f, 22.0f };
            VecF arr1(5,data1,1);
            DynProg dyn;
            dyn.less_before(arr1);
            TS_ASSERT_EQUALS(arr1[0], 10.0);
            TS_ASSERT_EQUALS(arr1[1], 2.0);
            TS_ASSERT_EQUALS(arr1[2], 3.0);
            TS_ASSERT_EQUALS(arr1[3], 1.0);
            TS_ASSERT_EQUALS(arr1[4], 6.0);
            
            // Test linear
            VecF result;
            dyn.linear(2.0, 30.0, 10, result);
            TS_ASSERT_EQUALS(result[0], 30.0);
            TS_ASSERT_EQUALS(result[9], 48.0);

            dyn.less_before(result);
            TS_ASSERT_EQUALS(result[0], 30.0);
            TS_ASSERT_EQUALS(result[1], 2.0);
            TS_ASSERT_EQUALS(result[2], 2.0);
            TS_ASSERT_EQUALS(result[9], 2.0);

            // Test linear less before:
            VecF lessbefore;
            dyn.linear_less_before(2.0, 30.0, 10, lessbefore);
            TS_ASSERT_EQUALS(lessbefore[0], 30.0);
            TS_ASSERT_EQUALS(lessbefore[1], 2.0);
            TS_ASSERT_EQUALS(lessbefore[2], 2.0);
            TS_ASSERT_EQUALS(lessbefore[9], 2.0);
        }

        void test_default_gap_penalty( void ) {
            float data1[] = { 
                4.0, 4.0, 4.0, 4.0,
                4.0, 4.0, 4.0, 4.0,
                6.0, 6.0, 6.0, 6.0,
                6.0, 6.0, 6.0, 6.0,
                6.0, 6.0, 4.0, 4.0
            };
            MatF smat(5,4,data1,1);
            VecF penalty;
            DynProg dyn;
            dyn.default_gap_penalty(smat, penalty);
            TS_ASSERT_EQUALS(penalty.len(), 9);
            TS_ASSERT_EQUALS(penalty[0], 5);
            TS_ASSERT_EQUALS(penalty[1], 2);
            TS_ASSERT_EQUALS(penalty[2], 2);
        }

        void test_find_path( void ) {
            // test global alignment:
            float data1[] = { 
                23.2, 200.0, 70.0, 34.0, 43.0, 23.0,
                12.3, 60.3, 20.0, 12.0, 23.0, 12.2,
                12.3, 20.0, 70.2, 12.0, 23.0, 12.2,
                12.3, 20.0, 70.2, 12.0, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 68.9, 12.2,
                12.3, 20.0, 12.0, 68.8, 68.9, 100.3,
                12.3, 20.0, 12.0, 68.8, 86.2, 70.3 
            };
            MatF smat(9,6,data1,1); 
            DynProg dyn;
            VecF gap_penalty;
            dyn.find_path(smat, gap_penalty);
            //puts("find path using score + gaps");
            //dyn._asmat.print();
            //puts("TBPATH");
            //dyn._tbpath.print();
            TS_ASSERT_EQUALS(dyn._tbpath(8,5), 1);  // end
            TS_ASSERT_EQUALS(dyn._tbpath(8,4), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);
            TS_ASSERT_EQUALS(dyn._tbpath(3,1), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(3,2), 1);
            //TS_ASSERT_EQUALS(dyn._tbpath(3,3), 1);  //examine!
            TS_ASSERT_EQUALS(dyn._tbpath(3,4), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(0,0), 1);  // start
            TS_ASSERT_EQUALS(dyn._tbpath(0,1), 1);
            //TS_ASSERT_EQUALS(dyn._tbpath(0,2), 1);  //examine!
            TS_ASSERT_EQUALS(dyn._tbpath(0,3), 0);
        }

        void test_warp_map( void ) {
            int data1[] = {0, 0, 0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 7, 8,9,9,10,10,10,11};
            int data2[] = {0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5,6,7,8,8,8,9};
            float data3[] = {23.2, 200.0, 70.0, 20.0, 70.2, 70.2, 12.0, 68.8, 68.8, 68.8, 68.9, 68.9, 100.3, 70.3, 60.0, 20.0, 10.0, 12.0, 13.0, 14.0};
            float percent_anchs = 100.f;
            VecI mCoords(20, data1, 1);
            VecI nCoords(20, data2, 1);
            VecF scores(20, data3, 1);
            // Test warp_map:
            VecI mout;
            VecI nout;
            DynProg dyn;

            dyn._mCoords = mCoords;
            dyn._nCoords = nCoords;
            dyn._sCoords = scores;
            dyn.warp_map(mout, nout,percent_anchs,0);
            int mans[7] = {0, 3, 6, 7, 8, 9, 11};
            int nans[7] = {0, 2, 3, 4, 5, 6, 9};
            VecI manswer(7,mans,1);
            VecI nanswer(7,nans,1);
            //puts("MEXPEXT: "); manswer.print(); puts("NEXPECT: "); nanswer.print();
            //puts("MOUT: "); mout.print(); puts("NOUT: "); nout.print();
            TS_ASSERT_EQUALS(mout, manswer);
            TS_ASSERT_EQUALS(nout, nanswer);

            // Test simpler function call
            VecI mout2;
            VecI nout2;
            dyn.warp_map(mout2, nout2, percent_anchs, 0);
            TS_ASSERT_EQUALS(mout2, mout);
            TS_ASSERT_EQUALS(nout2, nout);

            // Test minimization:
            VecI mout3;
            VecI nout3;
            scores *= -1.f;  // minimize the negative should give same coords
            dyn._sCoords = scores;
            dyn.warp_map(mout3, nout3, percent_anchs, 1);
            TS_ASSERT_EQUALS(mout3, mout);
            TS_ASSERT_EQUALS(nout3, nout);
        }

        void test_find_path_with_gaps( void ) {
            // Using find_path function now with diag and gap factors!
            // test global alignment:
            float data1[] = { 
                23.2, 200.0, 70.0, 34.0, 43.0, 23.0,
                12.3, 60.3, 20.0, 12.0, 23.0, 12.2,
                12.3, 20.0, 70.2, 12.0, 23.0, 12.2,
                12.3, 20.0, 70.2, 12.0, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
                12.3, 20.0, 12.0, 68.8, 68.9, 12.2,
                12.3, 20.0, 12.0, 68.8, 68.9, 100.3,
                12.3, 20.0, 12.0, 68.8, 86.2, 70.3 
            };
            MatF smat(9,6,data1,1); 
            DynProg dyn;
            VecF gap_penalty;
            dyn.find_path(smat, gap_penalty, 0, 2.f, 0.f);   // turn off gaps

            TS_ASSERT_EQUALS(dyn._tbpath(8,5), 1);  // end
            TS_ASSERT_EQUALS(dyn._tbpath(8,4), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);
            TS_ASSERT_EQUALS(dyn._tbpath(0,0), 1);  // start
            TS_ASSERT_EQUALS(dyn._tbpath(3,2), 1);

            // Test local alignment:
            dyn.find_path(smat, gap_penalty, 0, 2.f, 0.f, 1, 5.0);
            TS_ASSERT_EQUALS(dyn._tbpath(0,1), 1);  // start at best score
            TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);  // end at best score
            TS_ASSERT_EQUALS(dyn._tbpath(8,4), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);
            TS_ASSERT_EQUALS(dyn._tbpath(3,2), 1);

            smat *= -1.f;
            smat += 200.f;
            dyn.find_path(smat, gap_penalty, 1, 1.f, 1.f, 0);
            //puts("SMAT"); smat.print(); puts("TBpath"); dyn._tbpath.print(); puts("TB"); dyn._tb.print();
            TS_ASSERT_EQUALS(dyn._tbpath(8,5), 1);  // end
            TS_ASSERT_EQUALS(dyn._tbpath(8,4), 0);
            TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);
            TS_ASSERT_EQUALS(dyn._tbpath(0,0), 1);  // start
            TS_ASSERT_EQUALS(dyn._tbpath(0,1), 1);
            TS_ASSERT_EQUALS(dyn._tbpath(3,2), 1);

            // Test local minimizing alignment:
            // This DOES not work well, because the minimal will simply be 
            //dyn.find_path(smat, gap_penalty, 1, 10.f, 0.5f, 1, -5.0);
            //puts("SMAT"); smat.print(); puts("TBpath"); dyn._tbpath.print(); puts("TB"); dyn._tb.print();
            //TS_ASSERT_EQUALS(dyn._tbpath(0,1), 1);  // start at best score
            //TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);  // end at best score
            //TS_ASSERT_EQUALS(dyn._tbpath(8,4), 0);
            //TS_ASSERT_EQUALS(dyn._tbpath(7,5), 1);
            //TS_ASSERT_EQUALS(dyn._tbpath(3,2), 1);
        }

        void test_path_accuracy( void ) {
            float marr[5] = {1.f, 2.f, 3.f, 4.f, 5.f};
            float narr[5] = {1.f, 2.f, 3.f, 4.f, 5.f};
            VecF m_tm(5,marr,1);
            VecF n_tm(5,narr,1);
            int mwarr[2] = {0,4};
            int nwarr[2] = {0,4};
            VecI mWarpMap(2, mwarr,1);
            VecI nWarpMap(2, nwarr,1);
            float mvarr[5] = {1.f, 4.f, 7.f, 4.f, 5.f};
            float nvarr[5] = {1.f, 2.f, 3.f, 4.f, 5.f};
            VecF mVals(5,mvarr,1);
            VecF nVals(5,nvarr,1);
            DynProg dyn;
            float sum = dyn.sum_sq_res_yeqx(m_tm, n_tm, mWarpMap, nWarpMap, mVals, nVals);
            TS_ASSERT_EQUALS(sum, 10);
            float ssr, asr, sad, aad;
            dyn.path_accuracy(m_tm, n_tm, mWarpMap, nWarpMap, mVals, nVals, ssr,asr,sad,aad);
            TS_ASSERT_EQUALS(ssr, 10);
            TS_ASSERT_EQUALS(asr, 2);
            TS_ASSERT_EQUALS(sad, 6);
            TS_ASSERT_DELTA(aad, 1, 0.21); // 1.2 is actual

            // Test the other invocation:
            float ssr2, asr2, sad2, aad2;
            VecF mWarpMapFt(mWarpMap.length());
            VecF nWarpMapFt(nWarpMap.length());
            for (int i = 0; i < mWarpMap.length(); ++i) {
                mWarpMapFt[i] = m_tm[mWarpMap[i]];
                nWarpMapFt[i] = n_tm[nWarpMap[i]];
            }
            dyn.path_accuracy(mWarpMapFt, nWarpMapFt, mVals, nVals, ssr2,asr2,sad2,aad2);
            
            TS_ASSERT_EQUALS(ssr, ssr2);
            TS_ASSERT_EQUALS(asr, asr2);
            TS_ASSERT_EQUALS(sad, sad2);
            TS_ASSERT_EQUALS(aad, aad2);

            // DETAILS:
            VecF sr;
            VecF ad;
            dyn.path_accuracy_details(mWarpMapFt, nWarpMapFt, mVals, nVals, sr, ad);

            float ans_sr_arr[5] = {0, 2, 8, 0, 0};
            float ans_ad_arr[5] = {0, 2, 4, 0, 0};
            VecF ans_sr(5, ans_sr_arr, 1);
            VecF ans_ad(5, ans_ad_arr, 1);
            TS_ASSERT_EQUALS(ans_sr, sr)
            TS_ASSERT_EQUALS(ans_ad, ad)
        }

        // I think I have clear the different kinds of warps I'm interested in
        // this is the most complicated one,  need to touch it up....

//        void test_warp( void ) {
//            int data1[] = {0,2,4,6,8,10};
//            int data2[] = {0,1,2,3,4,5};
//            VecI mCoords(6,data1,1);
//            VecI nCoords(6,data2,1);
//            float dataf[] = {
//                10.f,20.f,30.f,15.f,0.f,5.f,10.f,15.f,10.f,5.f,0.f,
//                1.f, 2.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 8.f, 5.f, 4.f
//            };
//            float dataansw[] = { 
//                10.f, 30.f, 0.f, 10.f, 10.f, 0.f,
//                1.f, 5.f, 7.f, 9.f, 8.f, 4.f                 
//            };
//            MatF answtr(2,6,dataansw,1); 
//            MatF answ;
//            answtr.transpose(answ);
//            MatF totrans(2,11,dataf,1);
//            MatF mguy;
//            totrans.transpose(mguy);
//            DynProg dyn;
//            MatF out;
//            dyn.warp(mCoords, nCoords, mguy, out);
//            TS_ASSERT_EQUALS(out, answ); 
//            dyn.warp(mCoords, nCoords, mguy, out,1);
//            TS_ASSERT_EQUALS(out, mguy); 
//
//            float dataf2[] = {
//                10.3f,20.f,30.f,15.3f,0.5f,5.5f,10.2f,15.3f,10.f,5.f,0.f,
//                1.f, 2.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 8.f, 5.f, 4.f
//            };
//            float dataansw2[] = { 
//                23.1031f, 0.5f, 5.5f, 10.2f, 10.f, 0.f,
//                2.3703f, 7.f, 8.f, 9.f, 8.f, 4.f
//            };
//            MatF answtr2(2,6,dataansw2,1); 
//            MatF answ2;
//            answtr.transpose(answ2);
//            MatF totrans2(2,11,dataf2,1);
//            MatF mguy2;
//            totrans2.transpose(mguy2);
//            int data3[] = {2,4,5,6,8,10};
//            int data4[] = {0,1,2,3,4,5};
//            VecI mCoords2(6,data3,1);
//            VecI nCoords2(6,data4,1);
//            MatF out2;
//            dyn.warp(mCoords2, nCoords2, mguy2, out2);
//            TS_ASSERT_DELTA(out2(0,0), 23.1031, 0.0001);
//            TS_ASSERT_EQUALS(out2(1,0), 0.5f);
//            TS_ASSERT_EQUALS(out2(2,0), 5.5f);
//            TS_ASSERT_EQUALS(out2(5,0), 0.f);
//            TS_ASSERT_DELTA(out2(0,1), 2.3703, 0.0001);
//            TS_ASSERT_EQUALS(out2(1,1), 7.f);
//            TS_ASSERT_EQUALS(out2(2,1), 8.f);
//            TS_ASSERT_EQUALS(out2(5,1), 4.f);
//
//
//            float dataf3[] = {
//                10.3f,20.f,30.f,15.3f,0.5f,5.5f,10.2f,15.3f,10.f,5.f,0.f
//            };
//            VecF towarp(11, dataf3, 1);
//            VecF out3;
//            dyn.warp(mCoords2, nCoords2, towarp, out3);
//            TS_ASSERT_DELTA(out3[0], 23.1031, 0.0001);
//            TS_ASSERT_EQUALS(out3[1], 0.5f);
//            TS_ASSERT_EQUALS(out3[2], 5.5f);
//            TS_ASSERT_EQUALS(out3[5], 0.f);
//
//            float data5[] = {2.f,4.f,5.f,6.f,8.f,10.f};
//            float data6[] = {0.f,1.f,2.f,3.f,4.f,5.f};
//            VecF out4;
//            VecF mCoords3(6,data5,1);
//            VecF nCoords3(6,data6,1);
//            dyn.warp(mCoords3, nCoords3, towarp, out4);
//            TS_ASSERT_DELTA(out4[0], 23.1031, 0.0001);
//            TS_ASSERT_EQUALS(out4[1], 0.5f);
//            TS_ASSERT_EQUALS(out4[2], 5.5f);
//            TS_ASSERT_EQUALS(out4[5], 0.f);
//        }

//            float data2[] = { 
//                23.2, 34.2, 70.0, 34.0, 80.0, 23.0,
//                12.3, 60.3, 20.0, 12.0, 27.0, 100.0,
//                12.3, 20.0, 70.2, 18.0, 23.0, 12.2,
//                12.3, 20.0, 70.2, 12.0, 23.0, 12.2,
//                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
//                12.3, 20.0, 12.0, 68.8, 23.0, 12.2,
//                12.3, 20.0, 12.0, 68.8, 68.9, 12.2,
//                12.3, 20.0, 12.0, 68.8, 68.9, 100.0,
//                12.3, 20.0, 12.0, 68.8, 86.2, 70.3 
//            };
//
//            MatF scores;
//
//            MatF trial1(9,6,data1);  
//            MatF trial2(9,6,data2); 


            //float gp [] = {10.0, 10.0, 10.0,10.0, 10.0, 10.0,10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
            //float ip = 100000.f;
            //char *type = "global";
            //int mini = 0;
            //printf("HELLO YOU!\n");
            // MatI&, float*, MatI&, int, int, char*, float)
            // float init_penalty, char *type, bool mini
            // ip, type, min);
            //DynProg dyn(trial, gp, ip, type, mini);

            //          0      1     2     3     4     5     6     7     8     9
            //          0      1     3     6     7     12    13    14    15    17
            // scores:  23.4 * 34.5  37.8  32.1  80.1  43.2  76.4  98.4  100.3 11.3
            // True eq:      1 2 2 3 3 3 4 5 5 5 5 5|6  7  8  9  9  10
            // True eq:      1 2 2 3 3 3 4 5 5 5 5 5 6  7  8  9  9  10
            // indices:      0 1 2 3 4 5 6 7 8 9 101112 13 14 15  
            /*
               int test1[18] = {0,1,2,3,3,3,4,5,6,7,8,9,10,11,12,13,13,14};
               int test2[18] = {0,1,1,2,3,4,5,6,6,6,6,6,7 , 8, 9,10,11,12};
               VecI equiv1(test1, 18);
               VecI equiv2(test2, 18);
               float scores[18] = {23.4, 34.5, 35.6, 37.8, 40.9, 23.1, 32.1, 80.1, 100.2, 300.2, 10.1, 23.1, 43.2, 76.4, 98.4, 100.3,12.1, 11.3};

               VecF scoress(18, scores);

               DynProg dyn;
               dyn._equiv = equiv;
               dyn._score_path = scoress;
               VecI wm1;
               VecI wm2;
               dyn.warpMap(5, wm1, wm2);
               wm.print();
               */

//        void Xtest_expandFlag( void ) {
//            int base[36] = {
//                1, 0, 0, 0, 0, 0,
//                1, 0, 0, 0, 0, 0,
//                1, 0, 0, 0, 0, 0,
//                1, 0, 0, 0, 0, 0,
//                1, 0, 0, 0, 0, 0,
//                1, 0, 0, 0, 0, 0
//            };
//            float floats[36];
//            for (int i = 0; i < 36; ++i) {
//                floats[i] = (float)(i+20);
//            }
//            MatI mask(6,6, base);
//            MatF vals(6,6, floats);
//
//            //DynProg::expandFlag(trial, 2, 1, newer);
//
//            DynProg dyn;
//            vals.print();
//            dyn.replaceAlignmentPathRandom(vals, mask);
//            mask.print();
//            vals.print();
//        }


};

