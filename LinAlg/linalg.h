#ifndef LINALG_H
#define LINALG_H

#define SWAP

#include "matrix.h"

#include <cfloat>

namespace LinAlg {
    //Not a linear algebra method. This is a computer to method to guarantee a matrix precision.
    template <typename  Type>
    void Balance (LinAlg::Matrix<Type>& matrix_to_balance);

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

    template<typename Type>
    void Hessemberg_Form (const LinAlg::Matrix<Type>& matrix_to_reduce,
                          LinAlg::Matrix<Type>& unitary_matrix,
                          LinAlg::Matrix<Type>& hessemberg_matrix);

    //Simplified away to call Hessemberg_Form.
    template<typename Type>
    LinAlg::Matrix<Type> Hess (const LinAlg::Matrix<Type>& matrix_to_reduce);

    template<typename Type>
    void Hess (const LinAlg::Matrix<Type>& matrix_to_reduce,
               LinAlg::Matrix<Type>& unitary_matrix,
               LinAlg::Matrix<Type>& hessemberg_matrix);

    template <typename Type>
    LinAlg::Matrix<Type> EigenValues(const LinAlg::Matrix<Type>& matrix_to_get_eigenvalues, unsigned iterations = 100); //Not well tested, is incomplete.

    template <typename Type>
    LinAlg::Matrix<Type> RotationMatrix2D(double angle);

    template <typename Type>
    LinAlg::Matrix<Type> RotationMatrix3D(double angle, char axis);

    template <typename Type>
    LinAlg::Matrix<Type> HomogeneousTransformation(LinAlg::Matrix<Type> BCoordinates, LinAlg::Matrix<Type> BCoordinatesInA, double angle, char axis);

    template <typename Type>
    LinAlg::Matrix<Type> FixedAngles(double angle);

}

#include "src/linalg.hpp"

#endif // LINALG_H




