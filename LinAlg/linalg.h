#ifndef LINALG_H
#define LINALG_H

#include "matrix.h"

namespace LinAlg {
    template<typename Type>
    Type Trace (const LinAlg::Matrix<Type>& mat);

    template<typename Type>
    void QR_Factorization (const LinAlg::Matrix<Type>& input_matrix,
                           LinAlg::Matrix<Type>& output_Q_matrix,
                           LinAlg::Matrix<Type>& output_R_matrix);

    //Simplified away to call QR_Factorization.
    template<typename Type>
    void QR (const LinAlg::Matrix<Type>& input_matrix,
             LinAlg::Matrix<Type>& output_Q_matrix,
             LinAlg::Matrix<Type>& output_R_matrix);

    template<typename Type>
    LinAlg::Matrix<Type> Hessemberg_Form (const LinAlg::Matrix<Type>& matrix_to_reduce);

    //Simplified away to call Hessemberg_Form.
    template<typename Type>
    LinAlg::Matrix<Type> Hess (const LinAlg::Matrix<Type>& matrix_to_reduce);

    template<typename Type>
    LinAlg::Matrix<Type> EigenValues(const LinAlg::Matrix<Type>& mat);

}

#include "src/linalg.hpp"

#endif // LINALG_H




