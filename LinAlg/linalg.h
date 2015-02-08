#ifndef LINALG_H
#define LINALG_H

#include "matrix.h"

namespace LinAlg {
    template<typename Type>
    Type Trace(const LinAlg::Matrix<Type>& mat);
}

#include "src/linalg.hpp"

#endif // LINALG_H




