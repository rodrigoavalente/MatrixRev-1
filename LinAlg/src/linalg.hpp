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
void LinAlg::QR_Factoration (const LinAlg::Matrix<Type>& input_matrix,
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

