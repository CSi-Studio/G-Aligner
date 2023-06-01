
#ifndef _PNGIO_H
#define _PNGIO_H

#include "mat.h"

using namespace VEC;

class PngIO {
    private:
        int _bw;

    public:
        PngIO(bool bw=0);
        void write(char *file, MatI &mat);
        //bool write(char *file, VecI vec);
        void write(char *file, MatF &mat);
};

#endif
