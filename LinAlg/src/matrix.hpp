#include "../matrix.h"

template<typename Type>
LinAlg::Matrix<Type>::Matrix (std::string Mat)
{
    this->Init(Mat);
}

template<typename Type>
LinAlg::Matrix<Type>::Matrix (unsigned row, unsigned column)
{
    this->Init(row, column);
}

template<typename Type>
LinAlg::Matrix<Type>::Matrix (const LinAlg::Matrix<Type>& otherMatrix)
{
    this->Init(otherMatrix.rows, otherMatrix.columns);

    for(unsigned i = 0; i < otherMatrix.rows; i++)
        for(unsigned j = 0; j < otherMatrix.columns; j++)
            this->mat[i][j] = otherMatrix.mat[i][j];
}

template<typename Type>
LinAlg::Matrix<Type>::~Matrix ()
{
    for(unsigned i = 0; i < this->rows; i++)
        delete this->mat[i];
    delete [] this->mat;

    this->rows = 0;
    this->columns = 0;

    this->mat = NULL;
}

template<typename Type>
void LinAlg::Matrix<Type>::Init (std::string Mat)
{
    unsigned commas = 1, semiColons = 1, row, column, lin = 0, col = 0;
    int posComma = 0, posSemiColon = 0;
    std::string temp;

    for(unsigned i = 0; i < Mat.length(); i++)
    {
        if(Mat[i] == ';')
            semiColons += 1;
        else if (Mat[i] == ',')
            commas += 1;
    }

    row = semiColons;
    column = (commas + semiColons - 1)/semiColons;

    this->Init(row, column);

    while(!(Mat.empty()))
    {
        posSemiColon = Mat.find(";");
        if(posSemiColon != -1)
            temp = Mat.substr(0, posSemiColon);
        else
        {
            temp = Mat;
            posSemiColon = Mat.length();
        }

        while(!(temp.empty()))
        {
            posComma = temp.find(",");
            if (posComma == -1)
                posComma = temp.length();

            std::string temp2 = temp.substr(0, posComma);
            Type number;

            if(temp2 == "")
                number = 0;
            else
                number = (Type)atof(temp2.c_str());

            this->mat[lin][col] =  number;
            temp.erase(0, posComma + 1);
            col++;
        }

        Mat.erase(0, posSemiColon + 1);
        col = 0;
        lin++;
    }
}

template<typename Type>
void LinAlg::Matrix<Type>::Init (unsigned row, unsigned column)
{
    if(row == 0)
        row = 1;
    if(column == 0)
        column = 1;

    this->rows = row;
    this->columns = column;

    this->mat = new Type*[row];
    for(unsigned i = 0; i < row; i++)
        this->mat[i] = new Type[column];

    LinAlg::Zeros(*this);
};

template<typename Type>
void LinAlg::Matrix<Type>::ReInit (unsigned row, unsigned column)
{
    LinAlg::Matrix<Type> temp(row, column);

    for(unsigned i = 0; i < this->rows; i++)
        for(unsigned j = 0; j < this->columns; j++)
            temp.mat[i][j] = this->mat[i][j];

    *this = temp;
}

template<typename Type>
void LinAlg::Matrix<Type>::Add (unsigned& row, unsigned& column, Type& number)
{
    unsigned greaterRow, greaterColumn;

    if(((row + 1) > this->rows) || ((column + 1) > this->columns))
    {
        if((row + 1) > this->rows)
            greaterRow = row + 1;
        else
            greaterRow = this->rows;

        if((column + 1) > this->columns)
            greaterColumn = column + 1;
        else
            greaterColumn = this->columns;

        this->ReInit(greaterRow, greaterColumn);
    }

    this->mat[row][column] = number;
}

template<typename Type> template<typename OtherMatrixType>
bool LinAlg::Matrix<Type>::CheckDimensions (const LinAlg::Matrix<OtherMatrixType>& rhs, unsigned operation)
{
    bool checked = false;

    switch(operation)
    {
    case 0:
        if((this->rows == rhs.getNumberOfRows()) && (this->columns == rhs.getNumberOfColumns()))
            checked = true;
        else
            std::cout << "As dimensoes nao batem. Impossivel realizar operacao." << std::endl;
        break;
    case 1:
        if(this->columns == rhs.getNumberOfRows())
            checked = true;
        else
            std::cout << "As dimensoes nao batem. Impossivel realizar operacao." << std::endl;
        break;
    }

    return checked;
}

template<typename Type> template<typename OtherMatrixType>
void LinAlg::Matrix<Type>::swap (const LinAlg::Matrix<OtherMatrixType>& otherMatrix)
{
    using std::swap;

    LinAlg::Matrix<Type> temp(otherMatrix.getNumberOfRows(), otherMatrix.getNumberOfColumns());

    for(unsigned i = 0; i < temp.rows; i++)
        for(unsigned j = 0; j < temp.columns; j++)
            temp.mat[i][j] = otherMatrix(i + 1, j + 1);

    swap (rows, temp.rows);
    swap (columns, temp.columns);

    swap (mat, temp.mat);
}

template<typename Type>
unsigned LinAlg::Matrix<Type>::getNumberOfRows () const
{
    return this->rows;
}

template<typename Type>
unsigned LinAlg::Matrix<Type>::getNumberOfColumns () const
{
    return this->columns;
}

template<typename Type>
bool LinAlg::Matrix<Type>::isNull ()
{
    bool ret = false;

    if(this->mat == NULL)
        ret = true;

    return ret;
}

template<typename Type>
bool LinAlg::Matrix<Type>::isSquare ()
{
    bool ret = false;

    if(this->rows == this->columns)
        ret = true;

    return ret;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Matrix<Type>::GetRow (unsigned number_of_the_row)
{
    LinAlg::Matrix<Type> ret(1, this->columns);

    for(unsigned j = 0; j < ret.columns; j++)
        ret.mat[0][j] = this->mat[number_of_the_row - 1][j];

    return ret;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Matrix<Type>::GetColumn (unsigned number_of_the_column)
{
    LinAlg::Matrix<Type> ret(this->rows, 1);

    for(unsigned i = 0; i < ret.rows; i++)
        ret.mat[i][0] = this->mat[i][number_of_the_column - 1];

    return ret;
}

template <typename Type>
void LinAlg::Matrix<Type>::SwapRows (unsigned row_to_be_swapped, unsigned row_to_take_place)
{
    LinAlg::Matrix<Type> aux1, aux2;

    aux1 = this->GetRow(row_to_be_swapped);
    aux2 = this->GetRow(row_to_take_place);

    for(unsigned j = 0; j < this->columns; j++)
    {
        this->mat[row_to_be_swapped - 1][j] = aux2.mat[0][j];
        this->mat[row_to_take_place - 1][j] = aux1.mat[0][j];
    }
}

template <typename Type>
void LinAlg::Matrix<Type>::SwapColumns (unsigned column_to_be_swapped, unsigned column_to_take_place)
{
    LinAlg::Matrix<Type>aux1, aux2;

    aux1 = this->GetColumn(column_to_be_swapped);
    aux2 = this->GetColumn(column_to_take_place);

    for(unsigned i = 0; i < this->rows; i++)
    {
        this->mat[i][column_to_be_swapped - 1] = aux2.mat[i][0];
        this->mat[i][column_to_take_place - 1] = aux1.mat[i][0];
    }
}

template <typename Type>
unsigned LinAlg::Matrix<Type>::Size ()
{
    unsigned ret;

    if(this->rows >= this->columns)
        ret = this->rows;
    else
        ret = this->columns;

    return ret;
}

template<typename Type>
Type& LinAlg::Matrix<Type>::operator() (unsigned row, unsigned column)
{
    return this->mat[row - 1][column - 1];
}

template<typename Type>
Type LinAlg::Matrix<Type>::operator() (unsigned row, unsigned column) const
{
    return this->mat[row - 1][column - 1];
}

template<typename Type>
void LinAlg::Matrix<Type>::operator() (unsigned row, unsigned column, Type number)
{
    this->Add(row, column, number);
}

template<typename Type>
void LinAlg::Matrix<Type>::operator= (std::string Mat)
{
    this->Init(Mat);
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator= (const LinAlg::Matrix<Type>& rhs)
{
    swap(rhs);

    return *this;
}

template<typename Type> template<typename OtherMatrixType>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator= (const LinAlg::Matrix<OtherMatrixType>& rhs)
{
    swap(rhs);

    return *this;
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator+= (const Type& rhs /*scalar*/)
{
    for(unsigned i = 0; i < this->rows; i++)
        for(unsigned j = 0; j < this-> columns; j++)
            this->mat[i][j] += rhs;

    return *this;
}

template<typename Type> template<typename RightType>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator+= (const LinAlg::Matrix<RightType>& rhs)
{
    if(CheckDimensions(rhs, 0))
    {
        for(unsigned i = 0; i < this->rows; i++)
            for(unsigned j = 0; j < this->columns; j++)
                this->mat[i][j] += rhs(i + 1, j + 1);
    }

    return *this;
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator-= (const Type& rhs /*scalar*/)
{
    return *this += -rhs;
}


template<typename Type> template<typename RightType>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator-= (const LinAlg::Matrix<RightType>& rhs)
{
    return *this += -rhs;
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator*= (const Type& rhs /*scalar*/)
{
    for(unsigned i = 0; i < this->rows; i++)
        for(unsigned j = 0; j < this->columns; j++)
            this->mat[i][j] *= rhs;

    return *this;
}

template<typename Type> template<typename RightType>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator*= (const LinAlg::Matrix<RightType>& rhs)
{

    if ((this->rows == 1) && (this->columns == 1))
    {
        *this = this->mat[0][0] * rhs;
    }
    else if(CheckDimensions(rhs, 1))
    {
        Type temp;
        LinAlg::Matrix<Type> tempMat(*this);
        this->Init(this->rows, rhs.getNumberOfColumns());

        for(unsigned i = 0; i < tempMat.rows; i++)
            for(unsigned col = 0; col < rhs.getNumberOfColumns(); col++)
            {
                temp = 0;
                for(unsigned j = 0; j < tempMat.columns; j++)
                    temp += tempMat.mat[i][j] * rhs(j + 1, col + 1);
                this->mat[i][col] = temp;
            }
    }

    return *this;
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator/= (const Type& rhs)
{
    return *this *= 1/rhs;
}

template<typename Type> template<typename RightType>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator/= (const LinAlg::Matrix<RightType>& rhs)
{
    return *this *= LinAlg::Inverse<RightType>(rhs);
}

template<typename Type>
LinAlg::Matrix<Type>& LinAlg::Matrix<Type>::operator^= (double exp)
{
    if(exp < 0)
    {
        *this = LinAlg::Inverse(*this);
        exp *= -1;
    }

    for(unsigned i = 0; i < this->rows; i++)
        for(unsigned j = 0; j < this->columns; j++)
            this->mat[i][j] = std::pow(this->mat[i][j], exp);

    return *this;
}

template<typename Type> template<typename RightType>
LinAlg::Matrix<Type> LinAlg::Matrix<Type>::operator| (LinAlg::Matrix<RightType>& rhs)
{
    LinAlg::Matrix<Type>ret;

    if(this->mat == NULL)
        ret = rhs;
    else
    {
        unsigned aux = this->columns;

        if(this->rows < rhs.getNumberOfRows())
            ret.Init(rhs.getNumberOfRows(), this->columns + rhs.getNumberOfRows());
        else
            ret.Init(this->rows, this->columns + rhs.getNumberOfRows());

        for(unsigned i = 0; i < this->rows; i++)
            for(unsigned j = 0; j < this->columns; j++)
                ret.mat[i][j] = this->mat[i][j];

        for(unsigned i = 1; i <= rhs.getNumberOfRows(); i++)
            for(unsigned j = 1; j <= rhs.getNumberOfColumns(); j++)
                ret(i, aux + j) = rhs(i, j);
    }

    return ret;
}

template<typename Type> template<typename RightType>
LinAlg::Matrix<Type> LinAlg::Matrix<Type>::operator|| (LinAlg::Matrix<RightType>& rhs)
{
    LinAlg::Matrix<Type>ret;

    if(this->mat == NULL)
        ret = rhs;
    else
    {
        unsigned aux = this->rows;

        if(this->columns < rhs.getNumberOfColumns())
            ret.Init(this->rows + rhs.getNumberOfRows(), rhs.getNumberOfColumns());
        else
            ret.Init(this->rows + rhs.getNumberOfRows(), this->columns);

        for(unsigned i = 0; i < this->rows; i++)
            for(unsigned j = 0; j < this->columns; j++)
                ret.mat[i][j] = this->mat[i][j];

        for(unsigned i = 1; i <= rhs.getNumberOfRows(); i++)
            for(unsigned j = 1; j <= rhs.getNumberOfColumns(); j++)
                ret(i + aux, j) = rhs(i, j);
    }

    return ret;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::operator- (const LinAlg::Matrix<Type>& mat)
{
    LinAlg::Matrix<Type> temp(mat);

    for(unsigned i = 1; i <= temp.getNumberOfRows(); i++)
        for(unsigned j = 1; j <= temp.getNumberOfColumns(); j++)
            temp(i, j) *= -1;

    return temp;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::operator~ (LinAlg::Matrix<Type>& mat)
{
    LinAlg::Matrix<Type> temp(mat.getNumberOfColumns(), mat.getNumberOfRows());

    for(unsigned i = 1; i <= mat.getNumberOfRows(); i++)
        for(unsigned j = 1; j <= mat.getNumberOfColumns(); j++)
            temp(j, i) = mat(i, j);

    return temp;
}

template<typename Type>
std::ostream& LinAlg::operator<< (std::ostream& output, const LinAlg::Matrix<Type>& mat)
{
    for(unsigned i = 1; i <= mat.getNumberOfRows(); i++)
    {
        for(unsigned j = 1; j <= mat.getNumberOfColumns(); j++)
            output << std::setw(10) << std::setprecision(5) << std::fixed << mat(i, j) << ' ';

        output << std::endl;
    }

    return output;
}

template<typename Type>
std::istream& LinAlg::operator>> (std::istream& input, LinAlg::Matrix<Type>& mat)
{
    std::string temp;

    input >> temp;
    mat = temp;

    return input;
}

template<typename Type>
bool LinAlg::operator== (const LinAlg::Matrix<Type>& lhs, const LinAlg::Matrix<Type>& rhs)
{
    bool ret = true;

    if((lhs.getNumberOfRows() == rhs.getNumberOfRows()) && (lhs.getNumberOfColumns() && rhs.getNumberOfColumns()))
    {
        for(unsigned i = 1; i <= lhs.getNumberOfRows(); i++)
            for(unsigned j = 1; j <= lhs.getNumberOfColumns(); j++)
                if(!(lhs(i, j) == rhs(i, j)))
                {
                    ret = false;
                    break;
                }
    }
    else
        ret = false;

    return ret;
}

template<typename Type>
void LinAlg::Zeros(Matrix<Type>& Mat)
{
    for(unsigned i = 1; i <= Mat.getNumberOfRows(); i++)
        for(unsigned j = 1; j <= Mat.getNumberOfColumns(); j++)
            Mat(i, j) = 0;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Zeros (unsigned rows, unsigned columns)
{

    LinAlg::Matrix<Type> Ret(rows, columns);

    return Ret;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Eye (unsigned dimension)
{
    LinAlg::Matrix<Type> Ret(dimension, dimension);

    for(unsigned i = 1; i <= dimension; i++)
        for(unsigned j = 1; j <= dimension; j++)
            {
                if( i == j)
                    Ret(i, j) = 1;
                else
                    Ret(i, j) = 0;
            }

    return Ret;
}

template<typename Type>
void LinAlg::Ones(LinAlg::Matrix<Type>& mat)
{
    for(unsigned i = 1; i <= mat.getNumberOfRows(); i++)
        for(unsigned j = 1; j <= mat.getNumberOfColumns(); j++)
            mat(i, j) = 1;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Ones(unsigned rows, unsigned columns)
{
    LinAlg::Matrix<Type> temp(rows, columns);

    LinAlg::Ones(temp);

    return temp;
}

template<typename Type>
Type LinAlg::Determinant(const LinAlg::Matrix<Type>& mat)
{
    Type determinant = 0;
    unsigned rows = mat.getNumberOfRows(), columns = mat.getNumberOfColumns(), aux1, aux2;
    LinAlg::Matrix<Type> temp(rows - 1, columns - 1);


    if(rows != columns)
    {
        determinant = 0;
        std::cout << "Operacao disponivel somente para matrizes quadradas.";
    }
    else if(rows == 1)
        determinant = mat(1, 1);
    else if(rows == 2)
        determinant = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1);
    else
    {
        for(unsigned k = 0; k < rows; k++)
        {
            aux1 = 0;
            aux2 = 0;
            for(unsigned i = 1; i < rows; i++)
            {
                for(unsigned j = 0; j < rows; j++)
                {
                    if(!(j == k))
                    {
                        temp(aux1 + 1, aux2 + 1) = mat(i + 1, j + 1);
                        aux2++;
                    }

                    if(aux2 == rows -1)
                    {
                        aux1++;
                        aux2 = 0;
                    }
                }
            }

            determinant += pow( -1, k)*mat(1, k + 1) * LinAlg::Determinant(temp);
        }
    }

    return determinant;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Cofactor(const LinAlg::Matrix<Type>& mat)
{
	unsigned rows = mat.getNumberOfRows(), columns = mat.getNumberOfColumns(), aux1, aux2;
	LinAlg::Matrix<Type> temp(rows - 1, columns - 1), ret(rows, columns);

	if(rows != columns)
    {
        LinAlg::Zeros(ret);
        std::cout << "Operacao disponivel somente para matrizes quadradas.";
    }
    else if(rows == 2)
    {
        ret(1, 1) = mat(2, 2);
        ret(2, 2) = mat(1, 1);
        ret(1, 2) = -mat(2, 1);
        ret(2, 1) = -mat(1, 2);
    }
    else
    {
        for(unsigned j = 1; j <= rows; j++)
            for(unsigned i = 1; i <= rows; i++)
            {
                aux1 = 1;

                for(unsigned m = 1; m <= rows; m++)
                {
                    if(!(m == i))
                    {
                        aux2 = 1;

                        for(unsigned n = 1; n <= rows; n++)
                        {
                            if(!(n == j))
                            {
                                temp(aux1, aux2) = mat(m, n);
                                aux2++;
                            }
                        }
                        aux1++;
                    }
                }

                ret(i, j) = pow(-1, i + j)*LinAlg::Determinant(temp);
            }
    }

    return ret;
}

template<typename Type>
LinAlg::Matrix<Type> LinAlg::Inverse(const LinAlg::Matrix<Type>& mat)
{
    Type determinant = LinAlg::Determinant(mat);
    unsigned rows = mat.getNumberOfRows(), columns = mat.getNumberOfColumns();
    LinAlg::Matrix<Type> ret(rows, columns);

    if(rows != columns)
        std::cout << "Operacao disponivel somente para matrizes quadradas.";
    else if( determinant == 0)
        std::cout << "Impossivel inverter, determinante igual a 0.";
    else
    {
        ret = LinAlg::Cofactor(mat);
        ret = (~ret)/LinAlg::Determinant(mat);
    }

    return ret;
}

template<typename Type>
void LinAlg::Print(const Matrix<Type>& mat)
{
  std::cout << std::endl;

  for(unsigned i = 1; i <= mat.getNumberOfRows(); i++)
  {
    for(unsigned j = 1; j <= mat.getNumberOfColumns(); j++)
        std::cout << std::setw(10) << mat(i, j) << ' ';

    std::cout << std::endl;
  }
}
