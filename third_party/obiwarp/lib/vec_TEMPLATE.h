#ifndef _VEC_H
#define _VEC_H

#include <fstream>

/*************************************************************
 * Creation from existing object/array is always deep!.  
 * If you want shallow, use a pointer!
 ************************************************************/ 


namespace VEC {

class VecF;
class VecD;
class VecI;

// BEGIN TEMPLATE

class VecABR {

    protected:
        // length
        int _n;
        FLOAT *_dat;
        bool _shallow;

    public:
        // Constructors:
        VecABR();
        // Data values are NOT set by default
        explicit VecABR(int n);
        VecABR(int n, const FLOAT &val);

        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        VecABR(int n, FLOAT *arr, bool shallow=0);

        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work 
        VecABR(const VecABR &A, bool shallow=0);

        operator FLOAT*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }
        operator const FLOAT*() {
            if (_n > 0) {
                return &(_dat[0]);
            }
            else {
                return 0;
            }
        }

        FLOAT first() { return _dat[0]; }
        FLOAT last() { return _dat[_n-1]; }

        FLOAT* pointer() { return &(_dat[0]); }

        // Returns the name of the class
        // del needs to be called
        char * class_name();
        // shallow ownership
        void set(int n, FLOAT *arr);
        // Deletes the object's previous memory (if not shallow) and takes
        // ownership of the array (destructor call delete[])
        // shallow in this context only refers to calling delete
        // no data is copied
        // deep ownership (no copy is performed)
        void take(int n, FLOAT *arr);

        void to_f(VecF &out);
        void to_i(VecI &out);

        // returns the first index at the value, else -1
        int index(FLOAT val);

        // shallow ownership
        void set(VecABR &A);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(VecABR &A);
        VecABR & operator=(const FLOAT &val);
        VecABR & operator=(VecABR &A);
        ~VecABR();
        // A deep copy unless shallow is set
        void copy(VecABR &receiver, bool shallow=0) const;

        bool operator==(const VecABR &A);
        
        int length() const { return _n; }
        int len() const { return _n; }
        int size() const { return _n; }
        int dim() const { return _n; }
        int dim1() const { return _n; }
        // Returns in a vector all the values matching mask value
        void mask_as_vec(FLOAT return_val, VecI &mask, VecABR &vec);

        // Returns true if all values are the same, false otherwise
        bool all_equal() { 
            FLOAT _min, _max; min_max(_min, _max);
            if (_min == _max) { return 1; } 
            else { return 0; }
        }

        FLOAT& operator[](int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { puts("index < 0 !"); exit(1); }
            if (i >= _n) { puts("i >= _n !"); exit(1); }
#endif
            return _dat[i]; 
        }
        const FLOAT& operator[](int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { puts("index < 0 !"); exit(1); }
            if (i >= _n) { puts("i >= _n !"); exit(1); }
#endif
            return _dat[i]; 
        }

        FLOAT& at(int i) {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { puts("index < 0 !"); exit(1); }
            if (i >= _n) { puts("i >= _n !"); exit(1); }
#endif
            return _dat[i]; 
        }
        const FLOAT& at(int i) const {
#ifdef JTP_BOUNDS_CHECK
            if (i < 0) { puts("index < 0 !"); exit(1); }
            if (i >= _n) { puts("i >= _n !"); exit(1); }
#endif
            return _dat[i]; 
        }

        bool shallow() {
            return _shallow;
        }

        // NOTE: All operators act on the caller!
        // Operators
        void operator+=(const VecABR &A);
        void operator-=(const VecABR &A);
        void operator*=(const VecABR &A);
        void operator/=(const VecABR &A);
        void operator+=(FLOAT val);
        void operator-=(FLOAT val);
        void operator*=(FLOAT val);
        void operator/=(FLOAT val);
    
        void add(const VecABR &toadd, VecABR &out);
        void sub(const VecABR &tosub, VecABR &out);
        void mul(const VecABR &tomul, VecABR &out);
        void div(const VecABR &todiv, VecABR &out);
        // This may be slow because we cast every value to double regardless
        void square_root();

        void logarithm(double base);
        void min_max(FLOAT &mn, FLOAT &mx);
        // alias for min_max
        void mn_mx(FLOAT &mn, FLOAT &mx) {min_max(mn,mx);}
        double avg() const;
        void hist(int num_bins, VecD &bins, VecI &freqs);
        void sample_stats(double &mean, double &std_dev);
        double prob_one_side_right(double x);
        FLOAT sum();
        FLOAT sum_of_sq();

        void abs_val();
        // converts the distribution of values into standard normal
        // Only for floating points right now!
        void std_normal();

        // uses quicksort to sort the values
        void sort();
        static int FLOATCompare( const void *a, const void *b );

        // Removes the value at index and shortens the array by one
        // not shallow anymore regardless of previous state 
        void remove(int index);
        
        //VecABR operator+(const VecABR &A);
        //void operator++();
        //void operator--();
        
        // prints the vector (space delimited on one line without any length
        // header)
        void print(bool without_length=0);
        // prints the vector to the file with the length written on the line
        // before, unless without_length == true
        void print(const char *, bool without_length=0);
        // prints the vector to the filehandle with the length written on the
        // line before, unless without_length == true
        void print(std::ostream &fout, bool without_length=0);

        // CLASS FUNCTIONS:
        static int pchst(FLOAT arg1, FLOAT arg2) {
            if      (arg1*arg2 > 0) { return  1; }
            else if (arg1*arg2 < 0) { return -1; } 
            else                    { return  0; }
        }

        static double pearsons_r(VecABR &x, VecABR &y);
        static double covariance(VecABR &x, VecABR &y);
        static double euclidean(VecABR &x, VecABR &y);
        static FLOAT dot_product(VecABR &x, VecABR &y);

        static void xy_to_x(VecABR &x, VecABR &y);
        static void x_to_xy(VecABR &x, VecABR &y);
        static void chim(VecABR &x, VecABR &y, VecABR &out_derivs);

        static void calc_cubic_coeff(VecABR &x, VecABR &y, VecABR &derivs, VecABR &c2, VecABR &c3);
        static void chfe(VecABR &xin, VecABR &yin, VecABR &xe, VecABR &out_ye, int sorted=0);
        //static void pchfe(VecABR &xin, VecABR &yin, VecABR &XE, VecABR &out_newy);
        // interpolates so that linearity is encouraged along x axis
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static inline void chfev(FLOAT X1, FLOAT F1, FLOAT D1, FLOAT C2, FLOAT C3, FLOAT XE, FLOAT &FE);
        static inline void chfev_all(FLOAT X1, FLOAT X2, FLOAT F1, FLOAT F2, FLOAT D1, FLOAT D2, FLOAT XE, FLOAT &FE);
       // static void chfev(FLOAT X1, FLOAT X2, FLOAT F1, FLOAT F2, FLOAT D1, FLOAT D2, int NE, FLOAT *XE, FLOAT *FE, int *nlr, int &ierr);

        // interpolates so that linearity is encouraged along the xy line
        // if out_new_y.length() == 0 then new memory is allocated
        // otherwise, uses whatever memory is allocated in out_new_y
        static void chfe_xy(VecABR &x, VecABR &y, VecABR &new_x, VecABR &out_new_y, int sorted=0);

        static void linear_derivs(VecABR &x, VecABR &y, VecABR &out_derivs);
        static void linear_interp(VecABR &xin, VecABR &yin, VecABR &xe, VecABR &out_ye, int sorted=0);
        //##### FOR ANY FUNCTION:
        //                B 
        //               /|
        //              / |
        //             /  |
        //          c /   |a
        //           /    |
        //          /     |
        //       A --------C
        //             b
        //  
        //  cPerp = Perpendicular from C to c
        //  (sin A)/a = (sin C)/c
        //  sin C = 1
        //  c = sqrt(a^2 + b^2)
        //  sin A = a/c
        //  cPerp = (a/c)*b      note:   = sin A * b
        //  cPerp^2 = (ab)^2/(a^2 + b^2)

        // g(x):  y = x
        // h(y):  x = y
        //  a = y - g(x)

        // b = actual x - (x = actual y)
        // a = actual y - (y = actual x)
        // avg residual = SUM(cPerp^2)/N
        // 

        //###### FOR X=Y
        // residual^2 = 1/2((Y-X)^2)        note: [or X-Y]
        static double sum_sq_res_yeqx(VecABR &x, VecABR &y);
        // divides the sum of square of residuals by the length of the vector
        static double avg_sq_res_yeqx(VecABR &x, VecABR &y);

        // returns the average of the absolute values of the differences at
        // each index
        static double avg_abs_diff(VecABR &x, VecABR &y);
        static void rsq_slope_intercept(VecABR &x, VecABR &y, double &rsq, double &slope, double &y_intercept);

    private:
        void _copy(FLOAT *p1, const FLOAT *p2, int len) const;
        static void outliers_from_regression_line(VecABR &x, VecABR &y, VecI &indices_out);
        double _zScore(double mean, double sigma, double x);

}; // End class VecABR

// END TEMPLATE

} // End namespace

#endif

