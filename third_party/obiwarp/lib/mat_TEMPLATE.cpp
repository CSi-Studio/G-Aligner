
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "mat.h"
#include "vec.h"

//#define JTP_BOUNDS_CHECK

//#define JTP_DEBUG

namespace VEC {

// BEGIN TEMPLATE

/****************************************************************
 * MatABR
 ***************************************************************/

// Constructors:
MatABR::MatABR() : _m(0), _n(0), _dat(0) {
#ifdef JTP_DEBUG
    puts("CONSTRUCTOR MatABR()!");
#endif
}

MatABR::MatABR(int m, int n) : _m(m), _n(n), _dat(m*n) {
#ifdef JTP_BOUNDS_CHECK
    if (m < 0 || n < 0) { puts("m or n < 0"); exit(1); }
#endif
#ifdef JTP_DEBUG
    puts("CONSTRUCTOR MatABR(m,n)!");
#endif
}

MatABR::MatABR(int m, int n, const FLOAT &val) : _m(m), _n(n), _dat(m*n, val) {
#ifdef JTP_DEBUG
    puts("CONSTRUCTOR MatABR(m,n,val)!");
#endif
}

MatABR::MatABR(int m, int n, FLOAT *arr, bool shallow) : _m(m), _n(n), _dat(m*n,arr,shallow) {
#ifdef JTP_DEBUG
    printf("CONSTRUCTOR MatABR(m,n,*arr,shallow) shallow=%d!\n", this->shallow());
#endif
}

MatABR::MatABR(const MatABR &A, bool shallow) : _m(A._m), _n(A._n), _dat(A._dat, shallow) { 
#ifdef JTP_DEBUG
    printf("CONSTRUCTOR MatABR(MatABR &A,shallow) shallow=%d!\n", this->shallow());
#endif
}

void MatABR::to_vec(VecABR &outvec, bool shallow) {
    if (shallow) {
        outvec.set(_dat);
    }
    else {
        _dat.copy(outvec);
    }
}

void MatABR::set(int m, int n, FLOAT *arr) {
    _dat.set(m*n,arr);
    _m = m;
    _n = n;
}

void MatABR::set(MatABR &A) {
    _dat.set(A._dat);  
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    puts("set called!");
#endif
}


void MatABR::take(int m, int n, FLOAT *arr) {
    _dat.take(m*n,arr);
    _m = m;
    _n = n;
}

void MatABR::take(MatABR &A) {
    // Checking is done in Vec to ensure we're not taking a shallow!
    _dat.take(A._dat);  
    _m = A._m;
    _n = A._n;
#ifdef JTP_DEBUG
    puts("take called!");
#endif
}

void MatABR::row_vecs(int &cnt, VecABR *vecs) {
    cnt = rows();
    int _cols = cols();
    for (int i = 0; i < cnt; ++i) {
        FLOAT *ptr = this->pointer(i);
        vecs[i].set(_cols, ptr);  // shallow allocation
    }
}



bool MatABR::operator==(const MatABR &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n && A._m == _m) { // Same size
        return _dat == A._dat;
    }
    else {
        return false;
    }
}

void MatABR::copy(MatABR &receiver, bool shallow) const {
    receiver._m = _m;
    receiver._n = _n;
    _dat.copy(receiver._dat, shallow);
#ifdef JTP_DEBUG
    puts("copy called!");
#endif
}

void MatABR::file_rows_cols(std::ifstream &stream, int &rows, int &cols) {
    rows = 0;
    cols = 0;
    int BIGGEST_LINE = 1000000;
    char line[1000000];  // windows doesn't like that variable there
    stream.getline(line, BIGGEST_LINE);
    ++rows;
    char *ptr = line;
    int linelength = 0;
    while(*ptr != '\0') {
        if (*ptr == ' ') {
            *ptr = '\0';  // keep track of spaces
            ++cols;    
        }
        ++ptr;
        ++linelength;
    }
    ++cols; // for the last one
    // Check for spaces on the end:
    while(1) {
        --ptr;
        // char is !\n && \r && \0
        if (*ptr == '\n' || *ptr == '\r') {
            continue;
        }
        else if (*ptr == '\0') {
            --cols;  // decrement for each space at the end
        }
        else { break; }
    }
    // Count the rows:
    while( stream.getline(line, BIGGEST_LINE) ) {
        // if the line starts with a real char:
        if (line[0] != ' ' && line[0] != '\n' && line[0] != '\r' && line[0] != '\0') {
            ++rows;
        }
    }
}

void MatABR::set_from_ascii(std::ifstream &stream, int m, int n, MatABR &out) {
    MatABR tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatABR::set_from_ascii(std::ifstream &stream, MatABR &out) {
    int m,n;
    stream >> m >> n;
    MatABR tmp(m,n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            stream >> tmp(i,j);
        }
    }
    out.take(tmp);
}

void MatABR::set_from_ascii(const char *file, bool without_axes) {
    std::ifstream fh(file);
    if (fh.is_open()) {
        if (without_axes) {
            int m,n;
            file_rows_cols(fh,m,n);
            // Rewind the stream to beginning
            fh.clear(); // forget we saw the eof
            fh.seekg(0,std::ios::beg);
            set_from_ascii(fh,m,n,(*this));
        }
        else {
            set_from_ascii(fh,(*this));
        }
        fh.close(); 
    }
    else {
        printf("Couldn't open %s\n", file);
        exit(1);
    }
}

void MatABR::set_from_binary(const char *file) {
    int bytes_read;
    FILE *fh = fopen(file, "rb"); 
    if (fh == NULL) {
        printf("Could not open %s for reading\n", file);
        exit(1);
    }
    bytes_read = fread(&_m, sizeof(int), 1, fh);
    bytes_read = fread(&_n, sizeof(int), 1, fh);
    // Read the matrix:
    int rows_by_cols = _m*_n;
    //printf("rbycools: %d\n", rows_by_cols);
    FLOAT *dat_tmp = new FLOAT[rows_by_cols];
    bytes_read = fread(dat_tmp, sizeof(FLOAT), rows_by_cols, fh);
    //puts("**********************************************************");
    //puts("THIS IS THE BINARY MAT READ in:");
    //printf("First: %d:%.0f\n", 0, dat_tmp[0]);
    //printf("Last: %d:%.0f\n", rows_by_cols, dat_tmp[rows_by_cols-1]);
    //puts("**********************************************************");

    _dat.take(rows_by_cols, dat_tmp);
    fclose(fh);
}


MatABR & MatABR::operator=(const FLOAT &val) {
    _dat = val;
    return *this;
}

MatABR & MatABR::operator=(MatABR &A) {
#ifdef JTP_DEBUG
    puts("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        puts("IN ASSIGNMENT OP MID");
#endif
        _m = A._m;
        _n = A._n;
        _dat = A._dat;
    }
    return *this;
}

MatABR::~MatABR( ) {
#ifdef JTP_DEBUG
    puts("DESTRUCTOR");
#endif
}

/*************************
 * MATH OPERATORS
 ************************/
void MatABR::operator+=(const MatABR &A) {
    if (A._n == _n && A._m == _m) {
        _dat += A._dat;
    }
}

void MatABR::operator-=(const MatABR &A) {
    if (A._n == _n && A._m == _m) {
        _dat -= A._dat;
    }
}

void MatABR::operator*=(const MatABR &A) {
    if (A._n == _n && A._m == _m) {
        _dat *= A._dat;
    }
}

void MatABR::operator/=(const MatABR &A) {
    if (A._n == _n && A._m == _m) {
        _dat /= A._dat;
    }
}

void MatABR::add(const MatABR &toadd, MatABR &out) {
    if (_n == toadd._n && _m == toadd._m) {
        _dat.add(toadd._dat, out._dat);
    }
}

void MatABR::sub(const MatABR &tosub, MatABR &out) {
    if (_n == tosub._n && _m == tosub._m) {
        _dat.sub(tosub._dat, out._dat);
    }
}

void MatABR::mul(const MatABR &tomul, MatABR &out) {
    if (_n == tomul._n && _m == tomul._m) {
        _dat.mul(tomul._dat, out._dat);
    }
}

void MatABR::div(const MatABR &todiv, MatABR &out) {
    if (_n == todiv._n && _m == todiv._m) {
        _dat.div(todiv._dat, out._dat);
    }
}

void MatABR::transpose(MatABR &out) {
    MatABR me(*this, 1);
    MatABR tmp(me.nlen(), me.mlen());  // reverse m,n
    for (int m = 0; m < mlen(); ++m) {
        for (int n = 0; n < nlen(); ++n) {
            tmp(n,m) = me(m,n);
        }
    }
    out.take(tmp);
}

/*
MatABR MatABR::operator+(const MatABR &A) {
    printf("Adim: %d selfdim: %d\n", A.dim(), _n);
    if (A.dim() != _n) {
        puts("**** NOT THE SAME *****!");
        MatABR blank;
        return blank;
    }

    else {
        MatABR *C = new MatABR(_n);
        MatABR tmp = *C;
        tmp._to_pass_up = C; 
        printf("TMPENEW %d\n", tmp.shallow());
        for (int i = 0; i < _n; ++i) {
            tmp[i] = _dat[i] + A[i];
        }
        return tmp;
    }
}

*/


void MatABR::expand(MatABR &result, FLOAT match, int expand_x_lt, int expand_x_rt, int expand_y_up, int expand_y_dn, int expand_diag_lt_up, int expand_diag_rt_up, int expand_diag_lt_dn, int expand_diag_rt_dn ) {
    int i;
    int m_len = this->dim1();
    int n_len = this->dim2();
    this->copy(result);
    for (int m = 0; m < m_len; ++m) {
        for (int n = 0; n < n_len; ++n) {
            if ((*this)(m,n) == match) {
                for (i = 1; i <= expand_x_lt; ++i) {
                    if (n-i >= 0) {
                        result(m,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_x_rt; ++i) {
                    if (n+i < n_len) {
                        result(m,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_y_up; ++i) {
                    if (m-i >= 0) {
                        result(m-i,n) = match;
                    }
                }
                for (i = 1; i <= expand_y_dn; ++i) {
                    if (m+i < m_len) {
                        result(m+i,n) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_up; ++i) {
                    if (n-i >= 0 && m-i >=0) {
                        result(m-i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_up; ++i) {
                    if (n+i < n_len && m-i >= 0) {
                        result(m-i,n+i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_lt_dn; ++i) {
                    if (n-i >= 0 && m+i < m_len) {
                        result(m+i,n-i) = match;
                    }
                }
                for (i = 1; i <= expand_diag_rt_dn; ++i) {
                    if (n+i < n_len && m+i < m_len) {
                        result(m+i,n+i) = match;
                    }
                }
            }
        }
    }
}


void MatABR::mask_as_vec(FLOAT return_val, MatI &mask, VecABR &out) {
    _dat.mask_as_vec(return_val, mask._dat, out);    
}


FLOAT MatABR::sum(int m) {
    FLOAT sum = 0;
    FLOAT *ptr = pointer(m);
    for (int i = 0; i < _n; ++i) {
        sum += ptr[i];
    }
    return sum;
}


void MatABR::print(bool without_axes) {
    MatABR tmp((*this),1); 
    if (!without_axes) {
        std::cout << _m << ' ' << _n << std::endl;
    }
    for (int m = 0; m < _m; ++m) {
        int n;
        for (n = 0; n < _n - 1; ++n) {
            std::cout << tmp(m,n) << " ";
        }
        std::cout << tmp(m,n);
        std::cout << std::endl;
    }
}

void MatABR::print(const char *filename, bool without_axes) {
    std::ofstream fh(filename);
    if (!fh) {
        std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_axes);
    fh.close();
}

void MatABR::print(std::ostream &fout, bool without_axes) {
    int m;
    if (!without_axes) {
        fout << _m << ' ' << _n << std::endl;
    }
    for (m = 0; m < _m; m++) {
        int n;
        for (n = 0; n < _n - 1; n++) {
            fout << _dat[(m*_n)+n] << " ";
        }
        fout << _dat[m*_n+n];
        fout << std::endl;
    }
}

void MatABR::write(const char *file) {
    if (file != NULL) {
        FILE *fh = fopen(file, "wb");
        fwrite(&_m, sizeof(int), 1, fh);
        fwrite(&_n, sizeof(int), 1, fh);
        fwrite((FLOAT*)(_dat), sizeof(FLOAT), _m*_n, fh);
        fclose(fh);
    }
    else {
        fwrite(&_m, sizeof(int), 1, stdout);
        fwrite(&_n, sizeof(int), 1, stdout);
        fwrite((FLOAT*)(_dat), sizeof(FLOAT), _m*_n, stdout);
    }
    //puts("**********************************************************");
    //puts("WRITING THE BINARY MAT:");
    //printf("1st val: %f\n", mptr[0]);
    //printf("last val: %f\n", mptr[rows_by_cols-1]);
    //puts("**********************************************************");
}

// END TEMPLATE

} // End namespace VEC


