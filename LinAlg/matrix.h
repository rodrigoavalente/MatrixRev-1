#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <iostream>

namespace LinAlg {
    template<typename Type>
    class Matrix
    {
        public:
            Matrix (std::string Mat);
            Matrix (unsigned row, unsigned column);
            Matrix (): rows(0), columns(0), mat(NULL){};
            Matrix (const LinAlg::Matrix<Type>& otherMatrix);
            virtual ~Matrix ();

            unsigned getNumberOfRows () const;
            unsigned getNumberOfColumns () const;

            Type& operator() (unsigned row, unsigned column);
            Type operator() (unsigned  row, unsigned column) const;
            void operator() (unsigned row, unsigned column, Type number);

            void operator= (std::string rhs);
            LinAlg::Matrix<Type>& operator= (const LinAlg::Matrix<Type>& otherMatrix);
            template<typename OtherMatrixType>
            LinAlg::Matrix<Type>& operator= (const LinAlg::Matrix<OtherMatrixType>& otherMatrix);


            LinAlg::Matrix<Type> operator- () const;

            LinAlg::Matrix<Type>& operator+= (const Type& rhs /*scalar*/);
            template<typename RightType>
            LinAlg::Matrix<Type>& operator+= (const LinAlg::Matrix<RightType>& rhs);

            LinAlg::Matrix<Type>& operator-= (const Type& rhs /*scalar*/);
            template<typename RightType>
            LinAlg::Matrix<Type>& operator-= (const LinAlg::Matrix<RightType>& rhs);


            LinAlg::Matrix<Type>& operator*= (const Type& rhs /*scalar*/);
            template<typename RightType>
            LinAlg::Matrix<Type>& operator*= (const LinAlg::Matrix<RightType>& rhs);


            LinAlg::Matrix<Type>& operator/= (const Type& rhs);
            friend LinAlg::Matrix<Type> operator/ (LinAlg::Matrix<Type> lhs, const Type& rhs) {return lhs /= rhs;};

            //LinAlg::Matrix<Type>& operator/= (const LinAlg::Matrix<Type>& rhs);
            //friend LinAlg::Matrix<Type> operator/ (const LinAlg::Matrix<Type> lhs, const LinAlg::Matrix<Type>& rhs) {return lhs /= rhs};

            template<typename OtherMatrixType>
            void swap (const LinAlg::Matrix<OtherMatrixType>& otherMatrix);
            template<typename OtherMatrixType>
            friend void swap (LinAlg::Matrix<Type>& lhs, LinAlg::Matrix<OtherMatrixType>& rhs) {lhs.swap(rhs);};

        private:
            void Init (std::string Mat);
            void Init (unsigned row, unsigned column);

            void ReInit (unsigned row, unsigned column);

            void Add (unsigned& row, unsigned& column, Type& number);

            template<typename OtherMatrixType>
            bool CheckDimensions(const LinAlg::Matrix<OtherMatrixType>& rhs, unsigned operation);

            unsigned rows, columns;
            Type** mat;
    };

    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator+ (LinAlg::Matrix<MatrixType> lhs, const ScalarType& rhs) {return lhs += rhs;}
    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator+ (const ScalarType& lhs, LinAlg::Matrix<MatrixType> rhs) {return rhs += lhs;}
    template<typename LeftType, typename RightType>
    LinAlg::Matrix<LeftType> operator+ (LinAlg::Matrix<LeftType> lhs, const LinAlg::Matrix<RightType>& rhs) {return lhs += rhs;}

    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator- (LinAlg::Matrix<MatrixType> lhs, const ScalarType& rhs) {return lhs -= rhs;}
    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator- (const ScalarType& lhs, LinAlg::Matrix<MatrixType> rhs) {return -rhs -= -lhs;}
    template<typename LeftType, typename RightType>
    LinAlg::Matrix<LeftType> operator- (LinAlg::Matrix<LeftType> lhs, const LinAlg::Matrix<RightType>& rhs) {return lhs -= rhs;}

    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator* (LinAlg::Matrix<MatrixType> lhs, const ScalarType& rhs) {return lhs *= rhs;}
    template<typename MatrixType, typename ScalarType>
    LinAlg::Matrix<MatrixType> operator* (const ScalarType& lhs, LinAlg::Matrix<MatrixType> rhs) {return rhs *= lhs;}
    template<typename LeftType, typename RightType>
    LinAlg::Matrix<LeftType> operator* (LinAlg::Matrix<LeftType> lhs, const LinAlg::Matrix<RightType>& rhs) {return lhs *= rhs;}


    template<typename Type>
    void Zeros (LinAlg::Matrix<Type>& Mat);

    template<typename Type>
    LinAlg::Matrix<Type> Zeros (unsigned rows, unsigned columns);

    template<typename Type>
    LinAlg::Matrix<Type> Eye (unsigned dimension);

    template<typename Type>
    Type Determinant (LinAlg::Matrix<Type>& mat);

    template<typename Type>
    LinAlg::Matrix<Type> Cofactor(LinAlg::Matrix<Type>& Mat);

    template<typename Type>
    void Print (const LinAlg::Matrix<Type>& Mat);
};

#include "MatrixHeaders/matrix.hpp"

#endif // MATRIX_H
