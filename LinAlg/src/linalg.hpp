#include "../linalg.h"

template<typename Type>
Type LinAlg::Trace (const LinAlg::Matrix<Type>& mat)
{
    Type ret = 0;

    if(mat.getNumberOfRows() != mat.getNumberOfColumns())
        std::cout << "O traco so e calculado para matrizes quadradas.";
    else
    {
        for(unsigned i = 1; i <= mat.getNumberOfRows(); i++)
            for(unsigned j = 1; j <= mat.getNumberOfColumns(); j++)
                if(i == j)
                    ret += mat(i, j);
    }

    return ret;
}

template<typename Type>
void LinAlg::QR_Factorization (const LinAlg::Matrix<Type>& input_matrix,
                               LinAlg::Matrix<Type>& output_Q_matrix,
                               LinAlg::Matrix<Type>& output_R_matrix)
{

    //Constants calculated and needed for the QR algorithm.
    unsigned R_rows = input_matrix.getNumberOfRows(), R_columns = input_matrix.getNumberOfColumns();
    Type tal, gama, sigma;

    output_Q_matrix = LinAlg::Eye<Type>(R_rows);
    output_R_matrix = input_matrix;

    for(unsigned j = 1; j <= R_columns; j++)
        for(unsigned i = R_rows; i >= j + 1; i--)
        {
            LinAlg::Matrix<Type> temp;

            temp = LinAlg::Eye<Type>(R_rows);

            if(output_R_matrix(i, j) != 0)
            {
                if(std::abs(output_R_matrix(i - 1, j)) > std::abs(output_R_matrix(i, j)))
                {
                    tal = output_R_matrix(i, j)/output_R_matrix(i - 1, j);
                    gama = 1/(std::sqrt(1 + std::pow(tal, 2)));
                    sigma = tal*gama;
                }
                else
                {
                    tal = output_R_matrix(i - 1, j)/output_R_matrix(i, j);
                    sigma = 1/(std::sqrt(1 + std::pow(tal, 2)));
                    gama = sigma*tal;
                }

                temp(i, i) = gama;
                temp(i, i - 1) = sigma;
                temp(i - 1, i) = -sigma;
                temp(i - 1, i - 1) = gama;

                output_R_matrix = (~temp)*output_R_matrix;
                output_Q_matrix *= temp;
            }


        }

}

//Simplified away to call QR_Factorization.
template<typename Type>
void LinAlg::QR (const LinAlg::Matrix<Type>& input_matrix,
                 LinAlg::Matrix<Type>& output_Q_matrix,
                 LinAlg::Matrix<Type>& output_R_matrix)
{
    LinAlg::QR_Factorization(input_matrix, output_Q_matrix, output_R_matrix);
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Hessemberg_Form (const LinAlg::Matrix<Type>& matrix_to_reduce)
{
    LinAlg::Matrix<Type> ret(matrix_to_reduce);

    if(ret.isSquare())
    {
        for(unsigned k = 1; k <= ret.getNumberOfRows() - 2; k++)
        {
            Type R, S = 0;
            LinAlg::Matrix<Type> X, W(ret.getNumberOfRows(), 1), V, c, Q;

            X = ret.GetColumn(k);

            for(unsigned i = k + 1; i <= X.getNumberOfRows(); i++)
                S += std::pow(X(i, 1), 2);

            if(X(k + 1,1) > 0)
                S = std::sqrt(S);
            else
                S = -1*std::sqrt(S);

            R = std::sqrt(2*S*(S + X(k + 1, 1)));

            for(unsigned i = k + 1; i <= W.getNumberOfRows(); i++)
            {
                if(i == k + 1)
                    W(i, 1) = X(i, 1) + S;
                else
                    W(i, 1) = X(i, 1);
            }

            W *= 1/R;
            V = ret * W;
            c = (~W) * V;
            Q = V - (c * W);
            ret -= ((2 * W * (~Q)) + (2 * Q * (~W)));
        }
    }
    else
        std::cout << "Funcao apenas para matrizes quadradas.";

    return ret;
}

//Simplified away to call Hessemberg_Form.
template<typename Type>
LinAlg::Matrix<Type> LinAlg::Hess (const LinAlg::Matrix<Type>& matrix_to_reduce)
{
    return LinAlg::Hessemberg_Form(matrix_to_reduce);
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::EigenValues(const LinAlg::Matrix<Type>& mat)
{
    LinAlg::Matrix<Type> ret(mat.getNumberOfRows(), 2);

    if(mat.isSquare())
    {
        LinAlg::Matrix<Type> H = LinAlg::Hess(mat), A;

        A = H;

        for(unsigned i = 0; i < 200; i++)
        {
            LinAlg::Matrix<Type> Q, R;

            LinAlg::QR(A, Q, R);
            A = R*Q;
        }

        ret = A;
    }
    else
        std::cout << "Funcao apenas para matrizes quadradas.";

    return ret;
}
