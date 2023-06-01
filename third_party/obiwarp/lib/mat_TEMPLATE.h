#ifndef _MAT_H
#define _MAT_H

#include "vec.h"

/*************************************************************
 * Creation from existing object/array is always shallow!.  
 * Will delete any memory allocated.
 * Will NOT delete any memory not allocated.
 * If you want deep then use copy function!
 ************************************************************/ 


namespace VEC {

class MatI;
class MatF;
class MatD;

// BEGIN TEMPLATE

class MatABR {

    public:
        // length
        int _m;
        int _n;
        VecABR _dat;
        // Constructors:
        MatABR();
        MatABR(int m, int n);
        MatABR(int m, int n, const FLOAT &val);
       

        // (copied from vec.h)
        // if (shallow == 1 (true)) then no memory is deleted upon destruction
        // if (shallow == 0 (false)) then delete[] is called
        // FOR THIS CONSTRUCTOR ONLY, there is no DEEP copying, EVER!
        MatABR(int m, int n, FLOAT *arr, bool shallow=0);

        // (copied from vec.h)
        // if (shallow == 0 (false)) a DEEP copy is made of the data
        // if (shallow == 1 (true)) a copy of the pointer is made
        // if (shallow) then no memory is released upon destruction
        // shallow is used for a quick copy with which to work 
        MatABR(const MatABR &A, bool shallow=0);

        operator FLOAT*() { return (FLOAT*)_dat; }
        operator const FLOAT*() { return (FLOAT*)_dat; }
        FLOAT* pointer() { return (FLOAT*)_dat; }
        FLOAT* pointer(int m) { return &_dat[m*_n]; }
        // creates vec objects 
        // caller must have allocated the array for the vec objects
        // the data is a shallow copy!
        // transpose and call row_vecs for col_vecs!
        void row_vecs(int &cnt, VecABR *vecs);

        MatABR & operator=(const FLOAT &val);
        // DEEP
        MatABR & operator=(MatABR &A);
        ~MatABR();
        // Deep copy unless shallow == true
        void copy(MatABR &receiver, bool shallow=0) const;

        void set_from_ascii(std::ifstream &stream, int m, int n, MatABR &out);
        void set_from_ascii(std::ifstream &stream, MatABR &out);
        void set_from_ascii(const char *file, bool without_axes=0);
        void set_from_binary(const char *file);
        void file_rows_cols(std::ifstream &stream, int &rows, int &cols);
        // tnt_array2d_utils.h has a good example (use ifstream)

        // shallow copy and no ownership of memory
        void set(int m, int n, FLOAT *arr);
        // shallow copy and no ownership of memory
        void set(MatABR &A);

        bool all_equal() {
            return _dat.all_equal();
        }

        // Deletes the object's memory (if not shallow) and takes ownership
        // of the array memory (we will call delete[])
        void take(int m, int n, FLOAT *arr);
        // Deletes previous memory (if not shallow) and takes ownership
        // of the other's memory.
        void take(MatABR &A);

        // flattens the matrix and returns a vector
        void to_vec(VecABR &outvec, bool shallow=0);

        bool operator==(const MatABR &A);
        
        bool shallow() { return _dat.shallow(); }
        int dim1() const { return _m; }
        int dim2() const { return _n; }
        int mlen() const { return _m; }
        int nlen() const { return _n; }
        int rows() const { return _m; }
        int cols() const { return _n; }

        FLOAT& operator()(int m, int n) {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { puts("n < 0"); exit(1); }
            if (n >= _n) { puts("n >= _n"); exit(1); }
            if (m < 0) { puts("m < 0"); exit(1); }
            if (m >= _m) { puts("m >= _m"); exit(1); }
#endif
            return _dat[m*_n + n]; 
        }
        const FLOAT& operator()(int m, int n) const {
#ifdef JTP_BOUNDS_CHECK
            if (n < 0) { puts("n < 0"); exit(1); }
            if (n >= _n) { puts("n >= _n"); exit(1); }
            if (m < 0) { puts("m < 0"); exit(1); }
            if (m >= _m) { puts("m >= _m"); exit(1); }
#endif
            return _dat[m*_n + n]; 
        }

        // NOTE: All assignment operators act on the caller!
        void operator+=(const MatABR &A);
        void operator-=(const MatABR &A);
        void operator*=(const MatABR &A);
        void operator/=(const MatABR &A);
        void operator+=(const FLOAT val) { _dat += val; }
        void operator-=(const FLOAT val) { _dat -= val; }
        void operator*=(const FLOAT val) { _dat *= val; }
        void operator/=(const FLOAT val) { _dat /= val; }

    
        void add(const MatABR &toadd, MatABR &out);
        void sub(const MatABR &tosub, MatABR &out);
        void mul(const MatABR &tomul, MatABR &out);
        void div(const MatABR &todiv, MatABR &out);

        // returns the transpose in out
        void transpose(MatABR &out);

        void std_normal() { _dat.std_normal(); }
        void logarithm(double base) { _dat.logarithm(base); }
        void expand(MatABR &result, FLOAT match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn );

        void min_max(FLOAT &_min, FLOAT &_max) { _dat.min_max(_min,_max); }
        double avg() { return _dat.avg(); }
        //void operator++();
        //void operator--();
        
        FLOAT sum() { return _dat.sum(); } // return the sum of the entire matrix
        FLOAT sum(int m);  // return the sum of a given row
        // Returns in a vector all the values matching mask value
        void mask_as_vec(FLOAT return_val, MatI &mask, VecABR &out);

        // prints the bare matrix as ascii
        void print(bool without_axes=0);
        void print(const char *file, bool without_axes=0);
        void print(std::ostream &fout, bool without_axes=0);

        // writes the matrix as binary (includes # rows first and # cols as
        // ints)
        void write(const char *file=NULL);

        
        // @TODO need to write these guys:
        // prints the matrix in binary format:
        // (int) num cols (int) num rows (FLOAT) data
//        void write(const char *file);
//        void write(std::ofstream &fout);

    private:
        void _copy(FLOAT *p1, const FLOAT *p2, int len) const;

}; // End class MatABR

// END TEMPLATE

} // End namespace

#endif

