
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include "vec.h"


#ifndef min
	#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#endif
#ifndef max
	#define max(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#endif

//#define JTP_BOUNDS_CHECK
//#define JTP_DEBUG

namespace VEC {

// BEGIN TEMPLATE


/****************************************************************
 * VecABR
 ***************************************************************/

// Constructors:
VecABR::VecABR() : _n(0), _shallow(true) {
#ifdef JTP_DEBUG
    puts("Creating DATA (NO ARGS)");
#endif
}

VecABR::VecABR(int n) : _n(n), _shallow(false) {
#ifdef JTP_BOUNDS_CHECK
    if (n < 0) { puts("n < 0, exiting"); exit(1); }
#endif
    _dat = new FLOAT[_n];
#ifdef JTP_DEBUG
    puts("Creating DATA(N)");
#endif
}

VecABR::VecABR(int n, const FLOAT &val) : _n(n), _shallow(false) {
    _dat = new FLOAT[_n];
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;            
    }
#ifdef JTP_DEBUG
    puts("Creating DATA(N,FLOAT)");
#endif
}

VecABR::VecABR(int n, FLOAT *arr, bool shallow) : _n(n), _dat(arr), _shallow(shallow) {
#ifdef JTP_DEBUG
    puts("SHALLOW, (N,*ARR)");
#endif
}

VecABR::VecABR(const VecABR &A, bool shallow) : _n(A._n), _shallow(shallow) { 
    if (!shallow) {
        _dat = new FLOAT[_n];
        for (int i = 0; i < _n; ++i) {
            _dat[i] = A._dat[i];
        }
    }
    else {
        _dat = A._dat;
    }
#ifdef JTP_DEBUG
    puts("created with VecABR(const VecABR &A)");
#endif
}

void VecABR::to_f(VecF &out) {
    VecF _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (float)(_dat[i]);
    }
    out.take(_tmp);
}

void VecABR::to_i(VecI &out) {
    VecI _tmp(_n);
    for (int i = 0; i < _n; ++i) {
        _tmp[i] = (int)(_dat[i]);
    }
    out.take(_tmp);
}


void VecABR::set(int n, FLOAT *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = true;
    _n = n;
}

void VecABR::take(int n, FLOAT *arr) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = arr;
    _shallow = false;
    _n = n;
}

void VecABR::set(VecABR &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    _dat = A._dat;
    _shallow = true;
    _n = A._n;
}


void VecABR::take(VecABR &A) {
    if (!_shallow) {
        delete[] _dat;
    }
    if (A._shallow) {
        puts("Can't take ownership of memory of a shallow Vec!");
        exit(1);
    }
    _dat = A._dat;
    A._shallow = true;
    _shallow = false;
    _n = A._n;
}


bool VecABR::operator==(const VecABR &A) {
    // We don't care if one is shallow and the A is not!
    if (A._n == _n) { // Same size
        if (A._dat == _dat) { return true; }  // Same data
        else {
            for (int i = 0; i < _n; ++i) {
                if (A._dat[i] != _dat[i]) { return false; } 
            }
            return true;
        }
    }
    else { 
        return false;
    }
}

void VecABR::copy(VecABR &receiver, bool shallow) const {
    if (!receiver._shallow) {
        delete[] receiver._dat;
    }
    if (shallow) {
        receiver._dat = _dat;
        receiver._n = _n;
        receiver._shallow = true;
    }
    else {
        receiver._n = _n;
        receiver._dat = new FLOAT[_n];
        _copy(receiver._dat, _dat, _n);
        receiver._shallow = false;
    }
}

void VecABR::_copy(FLOAT *p1, const FLOAT *p2, int len) const {
    // This is slightly faster on gcc Linux Mandrake
    for (int i = 0; i < len; ++i) {
        p1[i] = p2[i];
    }
    // Slightly slower
//    FLOAT *end = p1 + len;
//    while(p1 < end) {
//        *p1++ = *p2++;
//    }
}

VecABR & VecABR::operator=(const FLOAT &val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] = val;
    }
    return *this;
}

VecABR & VecABR::operator=(VecABR &A) {
#ifdef JTP_DEBUG
    puts("IN ASSIGNMENT OP tOP");
#endif
    if (this != &A) {
#ifdef JTP_DEBUG
        puts("IN ASSIGNMENT OP MID");
#endif
        if (!_shallow) {
            delete[] _dat;
        }
        _n = A._n;
        _dat = new FLOAT[_n];
        _copy(_dat, A._dat, _n);
        _shallow = false;
    }
    return *this;
}

VecABR::~VecABR( ) {
#ifdef JTP_DEBUG
    puts("DESTRUCTOR");
#endif
    if (!_shallow) {
#ifdef JTP_DEBUG
        puts("DELETING DATA");
#endif
        delete[] _dat;
    }
}

/*************************
 * MATH OPERATORS
 ************************/

void VecABR::operator+=(FLOAT val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] += val;
    }
}
void VecABR::operator-=(FLOAT val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= val;
    }
}
void VecABR::operator*=(FLOAT val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= val;
    }
}
void VecABR::operator/=(FLOAT val) {
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= val;
    }
}
void VecABR::operator+=(const VecABR &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] += A[i];
    }
}

void VecABR::operator-=(const VecABR &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] -= A[i];
    }
}

void VecABR::operator*=(const VecABR &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] *= A[i];
    }
}

void VecABR::operator/=(const VecABR &A) {
    int n = A.dim();
    if (n != _n) {
        return;
    }
    for (int i = 0; i < _n; ++i) {
        _dat[i] /= A[i];
    }
}

void VecABR::add(const VecABR &toadd, VecABR &out) {
    if (toadd._n == _n) {
        FLOAT *_tmparr = new FLOAT[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] + toadd[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecABR::sub(const VecABR &tosub, VecABR &out) {
    if (tosub._n == _n) {
        FLOAT *_tmparr = new FLOAT[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] - tosub[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecABR::mul(const VecABR &tomul, VecABR &out) {
    if (tomul._n == _n) {
        FLOAT *_tmparr = new FLOAT[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] * tomul[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}


void VecABR::div(const VecABR &todiv, VecABR &out) {
    if (todiv._n == _n) {
        FLOAT *_tmparr = new FLOAT[_n];
        for (int i = 0; i < _n; ++i) {
            _tmparr[i] = _dat[i] / todiv[i];
        }
        if (!out._shallow) {
            delete[] out._dat;
        }
        out._n = _n;
        out._shallow = false;
        out._dat = _tmparr;
    }
}

void VecABR::square_root() {
    FLOAT *me = (FLOAT*)(*this);
    for (int i = 0; i < _n; ++i) {
        me[i] = (FLOAT)sqrt((double)me[i]);      
    }
}

/*
VecABR VecABR::operator+(const VecABR &A) {
    printf("Adim: %d selfdim: %d\n", A.dim(), _n);
    if (A.dim() != _n) {
        puts("**** NOT THE SAME *****!");
        VecABR blank;
        return blank;
    }

    else {
        VecABR *C = new VecABR(_n);
        VecABR tmp = *C;
        tmp._to_pass_up = C; 
        printf("TMPENEW %d\n", tmp.shallow());
        for (int i = 0; i < _n; ++i) {
            tmp[i] = _dat[i] + A[i];
        }
        return tmp;
    }
}
*/

FLOAT VecABR::sum() {
    FLOAT *me = (FLOAT*)(*this);
    FLOAT sum = 0;
    for( int n = 0; n < _n; n++) {
        sum += me[n];
    }
    return sum;
}

char * VecABR::class_name() {
    char *name = new char[7]; 
    strcpy(name, "VecABR");
    return name;
}


void VecABR::abs_val() {
    for (int n = 0; n < _n; ++n) {
        if (_dat[n] < 0) { _dat[n] *= -1; }
    }
}


void VecABR::std_normal() {
    // @TODO: would like avg and stdev to stay double, even for floats!
    (*this) -= (FLOAT)this->avg();
    double mean, stdev;
    this->sample_stats(mean, stdev);
    (*this) /= (FLOAT)stdev;
}


void VecABR::remove(int index) {
    FLOAT *_tmp_arr = new FLOAT[_n - 1];
    int _cnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (i != index) {
            _tmp_arr[_cnt] = _dat[i];
            ++_cnt;
        }
    }
    if (!_shallow) {
        delete []_dat;
    }
    _n = _n - 1;
    _dat = _tmp_arr;
    _shallow = false;
}

// Could be faster for ints
int VecABR::FLOATCompare( const void *a, const void *b ) {
    FLOAT c = *(FLOAT *)a - *(FLOAT *)b;
    if ( c < 0 ) return -1;
    if ( c > 0 ) return 1;
    return 0;
}

void VecABR::sort() {
    qsort(_dat, _n, sizeof(FLOAT), FLOATCompare);
}

int VecABR::index(FLOAT val) { 
    for (int i = 0; i < _n; ++i) {
        if (val == _dat[i]) {
            return i;
        }
    }
    return -1;
}

double VecABR::avg() const {
    double total = 0;
    for( int n = 0; n < _n; ++n) {
        total += _dat[n];
    }
    return total/_n;
}

//double VecABR::prob_one_side_right(double x) {
//    double _mean;
//    double _sigma;
//    // zScore:
//    sample_stats( _mean, _sigma );
//    //printf("mean: %f sigma: %f\n", _mean, _sigma);
//
//    // @TODO: THIS IS BROKEN until I can figure out how to calc the normalCDF
//    // using public domain sources:
//    //return _normalCDF(_mean, _sigma, x); 
//}

void VecABR::sample_stats(double &mean, double &std_dev) {
    // Raw score method (below) is faster (~1.5X) than deviation score method
    // commonly used
    FLOAT *me = this->pointer();
    double _sum = 0.0;
    double _val;
    double _sumSq = 0.0;
    int _len = this->dim();
    for( int i=0; i<_len; ++i ) {
        _val = (double)me[i]; 
        _sum += _val;
        _sumSq += _val *_val;
    }
    double tmp = _sumSq - ((_sum * _sum)/_len);
    tmp /= _len>1 ? _len-1 : 1;
#ifdef WIN32
    std_dev = sqrt( tmp );
#else
    std_dev = std::sqrt( tmp );
#endif
    mean = _sum/_len;
}

double VecABR::_zScore(double mean, double sigma, double x) {
    return (x - mean)/(sigma == 0.0 ? 1E-20: sigma);
}

void VecABR::mask_as_vec(FLOAT return_val, VecI &mask, VecABR &out) {
    if (mask.size() != _n) { puts("mask.size() != this->length()"); exit(1); }
    FLOAT *me = (FLOAT*)(*this);
    int *maskptr = (int*)(mask);
    FLOAT *tmparr = new FLOAT[_n];
    int newcnt = 0;
    for (int i = 0; i < _n; ++i) {
        if (maskptr[i] == return_val) {
            tmparr[newcnt] = me[i];
            ++newcnt;
        }
    }
    out.take(newcnt, tmparr);
}

void VecABR::hist(int num_bins, VecD &bins, VecI &freqs) {
    int i;

    // Create the scaling factor
    FLOAT _min; FLOAT _max;
    min_max(_min, _max);
    double dmin = (double)_min;
    double conv = ((double)num_bins)/(double)(_max - _min);

    // initialize arrays
    VecD _bins(num_bins);
    VecI _freqs(num_bins, 0);
    int _len = this->dim();
    FLOAT *me = this->pointer();

    // Create the histogram:
    for (i = 0; i < _len; ++i) {
        int index = (int)((me[i]-_min)*conv);
        if (index == num_bins) {
            --index;
        }
        _freqs[index]++;
    }

    // Create the bins:
    double iconv = 1.0/conv;
    for (i = 0; i < num_bins; ++i) {
        _bins[i] = ((i+0.5) * iconv) + dmin;  // avg
        // _bins[i] = ((i+0.5) * iconv) + dmin; //min
    }

    bins.take(_bins);
    freqs.take(_freqs);
}

void VecABR::logarithm(double base) {
    FLOAT *me = (FLOAT*)(*this);
    for (int i = 0; i < _n; ++i) {
        //printf("ME: %f\n", me[i]);
        me[i] = (FLOAT)(log((double)(me[i]))/log(base));
        //printf("MELOGGED: %f\n", me[i]);
    }
}

void VecABR::min_max(FLOAT &mn, FLOAT &mx) {
    FLOAT *me = (FLOAT*)(*this);
    mn = me[0];
    mx = me[0];
    for (int n = 0; n < _n; ++n) {
        mn = min(mn, me[n]);
        mx = max(mx, me[n]);
    }
}



void VecABR::print(bool without_length ) {
    if (!without_length) {
        std::cout << _n << std::endl;
    }
    int i;
    for (i = 0; i < _n - 1; ++i) {
        std::cout << _dat[i] << " ";
    }
    std::cout << _dat[i]; // the last one
    std::cout << std::endl;
}

void VecABR::print(const char *filename, bool without_length) {
    std::ofstream fh(filename);
    if (!fh) {
        std::cout << "Error opening file " << filename << std::endl;
    }
    this->print(fh, without_length);
    fh.close();
}

void VecABR::print(std::ostream &fout, bool without_length) {
    int i;
    if (!without_length) {
        fout << _n << std::endl;
    }
    for (i = 0; i < _n - 1; ++i) {
        fout << _dat[i] << " ";
    }
    fout << _dat[i];
    fout << std::endl;
}


// Class functions:
// THIS MUST BE FOR FLOAT AND DOUBLE ONLY!!!
// This is a fairly precise Fortran->C translation of the SLATEC chim code
// Evaluate the deriv at each x point
// return 1 if less than 2 data points
// return 0 if no errors
// ASSUMES monotonicity of the X data points !!!!!
// ASSUMES that this->length() >= 2
// If length == 1 then derivs[0] is set to 0
// If length == 0 then prints message to STDERR and returns;
void VecABR::chim(VecABR &x, VecABR &y, VecABR &out_derivs) {
    FLOAT *tmp_derivs = new FLOAT[x.length()];
    // if they aren't the right type then warn
    
    //if (typeid(T) != typeid(float) || typeid(T) != typeid(double)) 
    //    printf("vec2 calling object must be of type float or double!\n");
    //    exit(2);
    // }

#ifdef JTP_BOUNDS_CHECK
    if (x.length() != y.length()) { puts("x.length() != y.length()"); exit(1); }
#endif
    int length = x.length();
    FLOAT del1;
    FLOAT del2;
    FLOAT h1;
    FLOAT h2;
    FLOAT hsum;
    FLOAT w1;
    FLOAT w2;
    FLOAT dmax;
    FLOAT dmin;
    FLOAT three = (FLOAT)3.0;
    FLOAT dsave;
    FLOAT drat1;
    FLOAT drat2;
    FLOAT hsumt3;

    int ierr = 0;
    int lengthLess1 = length - 1;
    
    if (length < 2) { 
        if (length == 1) {
            tmp_derivs[0] = 0;
            return;
        }
        else {
            std::cerr << "trying to chim with 0 data points!\n";
        }
    }

// THIS can be done BEFORE this routine if someone cares to...
//    // Check monotonicity
//    for (int i = 2; i < length; i++) {
//        if (x[i] <= x[i-1]) { 
//            return 2;
//        }
//    }

    h1 = x[1] - x[0];
    del1 = (y[1] - y[0]) / h1;
    dsave = del1;
    
    // special case length=2 --use linear interpolation
    if (lengthLess1 < 2) {
        tmp_derivs[0] = del1;
        tmp_derivs[1] = del1;
        out_derivs.take(3, tmp_derivs);
        return;
    }
    
    // Normal case (length >= 3)
// 10

    h2 = x[2] - x[1];
    del2 = (y[2] - y[1]) / h2;

// SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
//     SHAPE-PRESERVING.

    hsum = h1 + h2;
    w1 = (h1 + hsum)/hsum;
    w2 = -h1/hsum;
    tmp_derivs[0] = (w1*del1) + (w2*del2);
    if ( pchst(tmp_derivs[0],del1) <= 0 ) {
        tmp_derivs[0] = (FLOAT)0;
    }
    else if ( pchst(del1,del2) < 0 ) {
        // need to do this check only if monotonicity switches
        dmax = three * del1;
        if (fabs(tmp_derivs[0]) > fabs(dmax)) {
            tmp_derivs[0] = dmax;
        }
    }

    int pchstval;
    int ind;

    for (ind = 1; ind < lengthLess1; ind++) {
        if (ind != 1) {
            h1 = h2;
            h2 = x[ind+1] - x[ind];
            hsum = h1 + h2;
            del1 = del2;
            del2 = (y[ind+1] - y[ind])/h2;
        }
// 40
        tmp_derivs[ind] = (FLOAT)0;

        pchstval = pchst(del1,del2);

// 45
        if (pchstval > 0) {
            hsumt3 = hsum+hsum+hsum;
            w1 = (hsum + h1)/hsumt3;
            w2 = (hsum + h2)/hsumt3;
            dmax = (FLOAT)max( fabs(del1), fabs(del2) );
            dmin = (FLOAT)min( fabs(del1), fabs(del2) );
            drat1 = del1/dmax;
            drat2 = del2/dmax;
            tmp_derivs[ind] = dmin/(w1*drat1 + w2*drat2);
        }
// 42
        else if (pchstval < 0 ) {
            ierr = ierr + 1;
            dsave = del2;
            continue;
        }
// 41
        else {  // equal to zero
            if (del2 == (FLOAT)0) { continue; }
            if (VecABR::pchst(dsave,del2) < 0) { ierr = ierr + 1; }
            dsave = del2;
            continue;
        }
    }

// 50 
    w1 = -h2/hsum;
    w2 = (h2 + hsum)/hsum;
    tmp_derivs[ind] = w1*del1 + w2*del2;
    if ( VecABR::pchst(tmp_derivs[ind],del2) <= 0 ) {
        tmp_derivs[ind] = (FLOAT)0;
    }
    else if ( VecABR::pchst(del1, del2) < 0) {
        // NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
        dmax = three*del2;
        if (fabs(tmp_derivs[ind]) > fabs(dmax)) {
            tmp_derivs[ind] = dmax;      
        }
    }
    out_derivs.take(length, tmp_derivs);
    return;
}


void VecABR::xy_to_x(VecABR &x, VecABR &y) {
    FLOAT *_x = (FLOAT*)x;
    FLOAT *_y = (FLOAT*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] - _x[i]; 
    }
}

void VecABR::x_to_xy(VecABR &x, VecABR &y) {
    FLOAT *_x = (FLOAT*)x;
    FLOAT *_y = (FLOAT*)y;
    for (int i = 0; i < x.length(); i++) {
        _y[i] = _y[i] + _x[i]; 
    }
}


void VecABR::linear_derivs(VecABR &x, VecABR &y, VecABR &out_derivs) {
    VecABR tmp_d(x.size());
    for (int i = 0; i < x.size(); ++i) {
        tmp_d[i] = (y[i+1] - y[i]) / (x[i+1]-x[i]);
    }
    out_derivs.take(tmp_d);
}


void VecABR::linear_interp(VecABR &xin, VecABR &yin, VecABR &xe, VecABR &out_ye, int sorted) {

    if (out_ye.size() == 0) {
        FLOAT *to_take = new FLOAT[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecABR derivs;
    VecABR::linear_derivs(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;
    FLOAT dt;

    if (sorted) {
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input 
            out_ye[j] = yin[ifirst] + (dt*derivs[ifirst]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            dt = xe[j] - xin[ifirst];  // diff in x, eval to input 
            out_ye[j] = yin[ifirst] + (dt * ((yin[ir] - yin[ifirst]) / (xin[ir]-xin[ifirst])) );
        }
    }
}

FLOAT VecABR::sum_of_sq() {
    FLOAT *me = this->pointer();
    FLOAT total = 0;
    for( int n = 0; n < this->size(); n++) {
        total += me[n]*me[n];
    }
    return total;
}


double VecABR::pearsons_r(VecABR &x, VecABR &y) {
    
    // Preparation:
    double sum_xTy = VecABR::dot_product(x,y);
    double sum_x = x.sum();         
    double sum_y = y.sum();       
    // Could this step be sped up?
    double sum_x2 = x.sum_of_sq();
    double sum_y2 = y.sum_of_sq();
    int N = x.dim();
    
    // Here it is:
    // 'E' is Capital Sigma 
    // r = EXY - (EXEY/N) 
    //    -----------------
    //    sqrt( (EX^2 - (EX)^2/N) * (EY^2 - (EY)^2/N) )

    double top = sum_xTy - ((sum_x * sum_y)/N);
    double fbot = sum_x2 - ((sum_x*sum_x)/N);  //first part of bottom
    double sbot = sum_y2 - ((sum_y*sum_y)/N);  //second part of bottom
    return top / sqrt(fbot * sbot);

}

double VecABR::covariance(VecABR &x, VecABR &y) {
    int i;
    int len = x.size();
    double mean_x = 0;
    double mean_y = 0;
    // get the means and x * y
    for (i = 0; i < len; ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }
    mean_x /= len;
    mean_y /= len;
    double cov = 0;
    for (i = 0; i < len; ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov/len;
}

double VecABR::euclidean(VecABR &x, VecABR &y) {
    VecF diff(x.size());
    double sum_of_diffs = 0;
    for (int i = 0; i < x.size(); ++i) {
        sum_of_diffs += (x[i] - y[i]) * (x[i] - y[i]); 
    }
    return sqrt(sum_of_diffs);
}

FLOAT VecABR::dot_product(VecABR &x, VecABR &y) {
    //assert(x.dim() == y.dim());
    FLOAT sum = 0;
    for (int i = 0; i < x.dim(); i++) {
        sum += (x[i] * y[i]);
    }
    return sum;
}

void VecABR::chfe(VecABR &xin, VecABR &yin, VecABR &xe, VecABR &out_ye, int sorted) {
    //xin.print(); yin.print();

    if (out_ye.size() == 0) {
        FLOAT *to_take = new FLOAT[xe.size()];
        out_ye.take(xe.size(), to_take);
    }

    // Calc the derivs:
    VecABR derivs;
    VecABR::chim(xin,yin,derivs);
    int i,j,ir;  // i indexes xin, j indexes xnew
    int ifirst = 0;

    // find the bounding points in xin
    int istart;


    if (sorted) {
        VecABR c2(xin.size());
        VecABR c3(xin.size());
        calc_cubic_coeff(xin, yin, derivs, c2, c3);
        istart = 0;
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            istart = i;
            //printf("ifirst: %d ir %d i %d\n",ifirst, ir, i); printf("xin[ifirst]%f, xin[ir]%f, yin[ifirst]%f, yin[ir]%f, derivs[ifirst]%f, derivs[ir]%f, xe[j]%f\n", xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j]);
            chfev(xin[ifirst], yin[ifirst], derivs[ifirst], c2[ifirst], c3[ifirst], xe[j], out_ye[j]);
        }
    }
    else {

        // find the bounding points in xin
        for (j = 0; j < xe.size(); ++j) {
            ir = -1;
            istart = 0;
            // @TODO: This should be a binary search:
            for (i = istart; i < xin.size(); ++i) {
                // locate the interval
                if (xin[i] >= xe[j]) {
                    ir = i;
                    ifirst = i - 1;
                    break;
                }
            }
            if (ir == 0) { // left extrapolation
                ir = 1;
                ifirst = 0;
            }
            else if (ir == -1) { // right extrapolation
                ir = i - 1;
                ifirst = ir - 1;
            }
            //printf("ifirst: %d ir %d i %d\n",ifirst, ir, i); printf("xin[ifirst]%f, xin[ir]%f, yin[ifirst]%f, yin[ir]%f, derivs[ifirst]%f, derivs[ir]%f, xe[j]%f\n", xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j]);
            chfev_all(xin[ifirst], xin[ir], yin[ifirst], yin[ir], derivs[ifirst], derivs[ir], xe[j], out_ye[j]);
        }
    }
}




void VecABR::calc_cubic_coeff(VecABR &x, VecABR &y, VecABR &derivs, VecABR &c2, VecABR &c3) {

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    FLOAT DEL1, DEL2, DELTA, H;
    for (int i = 0; i < x.size() - 1; ++i) {
        H = x[i+1] - x[i];
        DELTA = (y[i+1] - y[i])/H;
        DEL1 = (derivs[i] - DELTA)/H;
        DEL2 = (derivs[i+1] - DELTA)/H;
        c2[i] = -(DEL1+DEL1 + DEL2);
        c3[i] = (DEL1 + DEL2)/H;
    }

}

void VecABR::chfev_all(FLOAT X1, FLOAT X2, FLOAT F1, FLOAT F2, FLOAT D1, FLOAT D2, FLOAT XE, FLOAT &FE) {
    FLOAT C2, C3, DEL1, DEL2, DELTA, H, X;

    H = X2 - X1;

    //  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
    DELTA = (F2 - F1)/H;
    DEL1 = (D1 - DELTA)/H;
    DEL2 = (D2 - DELTA)/H;
    C2 = -(DEL1+DEL1 + DEL2);
    C3 = (DEL1 + DEL2)/H;

    X = XE - X1;
    //printf("X: %f F1 %f D1 %f C2 %f C3 %f\n", X, F1, D1, C2, C3);
    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecABR::chfev(FLOAT X1, FLOAT F1, FLOAT D1, FLOAT C2, FLOAT C3, FLOAT XE, FLOAT &FE) {
    FLOAT X;
    X = XE - X1;
    //printf("X: %f F1 %f D1 %f C2 %f C3 %f\n", X, F1, D1, C2, C3);
    FE = F1 + X*(D1 + X*(C2 + X*C3));
}


void VecABR::chfe_xy(VecABR &x, VecABR &y, VecABR &new_x, VecABR &out_new_y, int sorted) {
    VecABR::xy_to_x(x,y);
    chfe(x,y,new_x, out_new_y, sorted);
    x_to_xy(new_x, out_new_y);
    VecABR::x_to_xy(x,y);
}

double VecABR::sum_sq_res_yeqx(VecABR &x, VecABR &y) {
    //// PDL way
    //return sum(0.5*(($y - $x)**2));
    double __sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        FLOAT diff = x[i] - y[i];
        __sum += 0.5*(diff*diff);
    }
    return __sum;
}

double VecABR::avg_sq_res_yeqx(VecABR &x, VecABR &y) {
    return (sum_sq_res_yeqx(x,y))/x.length();
}

double VecABR::avg_abs_diff(VecABR &x, VecABR &y) {
    double sum = 0.0;
    for (int i = 0; i < x.length(); ++i) {
        sum += fabs(x[i] - y[i]);
    }
    return sum/x.length();
}


void VecABR::rsq_slope_intercept(VecABR &x, VecABR &y, double &rsq, double &slope, double &y_intercept) {
    int i;
    double mean_x = x.avg();
    double mean_y = y.avg();
    double sum_sq_res_xx = 0.0;
    double sum_sq_res_yy = 0.0;
    double sum_sq_res_xy = 0.0;
//    double *sq_res_xx = new double[x.length()];
//    double *sq_res_yy = new double[y.length()];
//    double *sq_res_xy = new double[x.length()];
    //VecD sq_res_xx(x.length());
    //VecD sq_res_yy(x.length());
    //VecD sq_res_xy(x.length());
    for (i = 0; i < x.length(); ++i) {
        double x_minus_mean_i, y_minus_mean_i;
        x_minus_mean_i = ( (double)(x[i]) ) - mean_x;
        y_minus_mean_i = ( (double)(y[i]) ) - mean_y;
        sum_sq_res_xx += x_minus_mean_i*x_minus_mean_i;
        sum_sq_res_yy += y_minus_mean_i*y_minus_mean_i;
        sum_sq_res_xy += x_minus_mean_i*y_minus_mean_i;
    }
    slope = sum_sq_res_xy/sum_sq_res_xx;
    y_intercept = mean_y - (slope * mean_x); 
    rsq = (sum_sq_res_xy*sum_sq_res_xy)/(sum_sq_res_xx*sum_sq_res_yy);
//    delete[] sq_res_xx;
//    delete[] sq_res_yy;
//    delete[] sq_res_xy;
}


// END TEMPLATE

} // End namespace VEC


