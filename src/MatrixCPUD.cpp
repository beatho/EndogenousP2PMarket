#include "../head/MatrixCPUD.h"
#include <cstring>

double MatrixCPUD::rand1()
{
    double a = (double)(rand()) / ((double)(RAND_MAX));
    return a;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MatrixCPUD::MatrixCPUD() {
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur appele" << std::endl;
#endif

}


MatrixCPUD::MatrixCPUD(int l, int c, double value)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur parametre appele" << std::endl;
    std::cout << _matrixCPU << std::endl;
#endif
    _row = l;
    _column = c;
    _matrixCPU = new double[l*c];
    for (int elem = 0; elem < l * c;elem++) {
        _matrixCPU[elem] = value; 
    }
#ifdef DEBUG_CONSTRUCTOR
    std::cout << _matrixCPU << std::endl;
#endif
}

MatrixCPUD::MatrixCPUD(const MatrixCPUD & m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur recopie appele" << std::endl;
#endif
    _row = m._row;
    _column = m._column;
    _matrixCPU = new double[_row * _column];
    memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
}


MatrixCPUD& MatrixCPUD::operator=(const MatrixCPUD& m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur operateur = appele" << std::endl;
#endif
    _row = m._row;
    _column = m._column;
    DELETEA(_matrixCPU);
    _matrixCPU = new double[_row * _column];
    memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
    
    
    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////
double MatrixCPUD::get(int i, int j) const
{
    if ((i >= _row) || ( j >= _column) || (i < 0) || ( j < 0)) {
        std::cout << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("get : index out of bounds");
    }
    return _matrixCPU[i*_column+j];
}



int MatrixCPUD::getNCol() const
{
    return _column;
}

int MatrixCPUD::getNLin() const
{
    return _row;
}

bool MatrixCPUD::dim(MatrixCPUD* m) const
{ 
    return ((_row == m->getNLin()) && (_column == m->getNCol()));
}



void MatrixCPUD::getLin(MatrixCPUD* vector, int i) const
{
    if ((_column != vector->getNCol()) || ( vector->getNLin() != 1)) 
    {
        throw std::invalid_argument("wrong dimension of the row vector");
    } 
    for (int j = 0; j < _column; j++) {
        vector->set(0, j, get(i, j));
    }
}

bool MatrixCPUD::isEqual(MatrixCPUD* m, double pre) const
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    else {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++) {
                if (fabs(get(i, j) - m->get(i, j)) > pre) {
                    return false;
                }
            }
        }
    }
    return true;
}

int MatrixCPUD::getNNull() const
{
    int n = 0;
    for (int i = 0;i < _row * _column;i++) {
        n = n + ((fabs(_matrixCPU[i]) > 0.000000001));
    }
    return n;
}

int MatrixCPUD::getNNullHalf() const
{
    int n = 0;
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j <= i; j++) {
            n = n + (fabs(get(i, j)) > 0.000000001);
        }
    }

    return n;
}
void MatrixCPUD::swap(MatrixCPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    
    double* temp = _matrixCPU;
    _matrixCPU = m->_matrixCPU;
    m->_matrixCPU = temp;
}

void MatrixCPUD::getBloc(MatrixCPUD* dest, int iBegin, int iEnd, int jBegin, int jEnd)
{

    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    } if (dest->getNLin() != (iEnd - iBegin) || dest->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension");
    }
    int row = 0;

    for (int i = iBegin; i < iEnd; i++) {
        int col = 0;
        for (int j = jBegin; j < jEnd;j++) {
            dest->set(row, col, get(i,j));
            col = col + 1;
        }
        row = row + 1;
    }
}


void MatrixCPUD::setSize(int l, int c)
{
    DELETEA(_matrixCPU);
    _matrixCPU = new double[l * c];
    _row = l;
    _column = c;
}

///////////////////////////////////////////////////////////////////////////////
// Setter
///////////////////////////////////////////////////////////////////////////////
void MatrixCPUD::set(int i, int j, double value)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << "err " << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds");
    }
    _matrixCPU[i * _column + j] = value;
}

void MatrixCPUD::set(MatrixCPUD* m)
{
    if (!dim(m)) {
        std::cout << _row << " * " << _column << " against " << m->_row << " * " << m->_column << std::endl;
        throw std::invalid_argument("not the same dimension");
    }
    memcpy(_matrixCPU, m->_matrixCPU, _row * _column * sizeof(double));
}

void MatrixCPUD::set(double value)
{
    for (int elem = 0; elem < _column * _row; elem++) {
        _matrixCPU[elem] = value;
    }
}

void MatrixCPUD::setTrans(MatrixCPUD* m)
{
    if (_column != m->getNLin() || _row != m->getNCol()) {
        throw std::invalid_argument("not the same transposed dimension");
    }

    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, m->get(j, i));
        }
    }
}



void MatrixCPUD::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        std::cout << " err out of bound " << iBegin << " " << jBegin << " " << iEnd << " " << jEnd << std::endl;
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin ) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension");
    }
    int row = 0;
    
    for (int i = iBegin; i < iEnd; i++) {
        int col = 0;
        for (int j = jBegin; j < jEnd;j++) {
            set(i, j, m->get(row, col));
            col = col + 1;
        }
        row = row + 1;
    }

}

void MatrixCPUD::setEyes(double v)
{
    int N = _row * (_row < _column) + _column * (_column <= _row);
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, v *(i==j));
        }
    }
 
}

void MatrixCPUD::setEyes(MatrixCPUD* vect)
{
    if (_column != _row) {
        throw std::invalid_argument("matrix must be square");
    }
    if (vect->_column == 1) {
        if (vect->_row != _row) {
            throw std::invalid_argument("matrix must have the same size as the vector");
        }
        else {
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _column; j++) {
                    set(i, j, vect->get(i, 0) * (i == j));
                }
            }
        }
    }
    else if (vect->_row == 1) {
        if (vect->_column != _row) {
            throw std::invalid_argument("matrix must have the same size as the vector");
        }
        else {
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _column; j++) {
                    set(i, j, vect->get(0, j) * (i == j));
                }
            }
        }
    }
    else {
        throw std::invalid_argument("argument must be a vector");
    }

}

void MatrixCPUD::setRand(int eps, int divide)
{
    //srand(time(nullptr));
    //exit(1);
    int N = _column * _row;
    for (int elem = 0; elem < N; elem++) {
        int r = rand() % eps;
        int signe = (rand() % 2) * 2 - 1;
        _matrixCPU[elem] = (double) r / divide * signe;
       
    }


}

void MatrixCPUD::setRand1(double eps)
{
    //srand(time(nullptr));

    int N = _column * _row;
    for (int elem = 0; elem < N; elem++) {
        _matrixCPU[elem] = 2 * (rand1() - 0.5) * eps;

    }
    exit(0);    

}

void MatrixCPUD::setFromFile(std::string filename, int entete) // aucune s�curit�, l'utilisateur doit v�rifier d'utiliser le bon fichier avec le bon nombre d'entr�e...
{
    std::ifstream myfile(filename, std::ios::in);
    //std::cout << filename << "or  " << "../../" + filename << std::endl;
    // si taille fichier < matrix : reste de la matrix sera remplis par le dernier terme
    // si taille fichier > matrix : seule les N premiers termes du fichier seront lus
    
    

    if (myfile)
    {
        for (int i = 0; i < _row; i++) {
            if (entete) {
                std::string s;
                myfile >> s;
            }
            for (int j = 0; j < _column;j++) {
                double v;
                myfile >> v;
                set(i, j, v);
                //std::cout << v << std::endl;
            }
        }
        myfile.close();
    }
    else {

        std::ifstream myfile2("../../" + filename, std::ios::in);
       
        //std::cout << filename << std::endl;
        // si taille fichier < matrix : reste de la matrix sera remplis par le dernier terme
        // si taille fichier > matrix : seule les N premiers termes du fichier seront lus



        if (myfile2)
        {
            for (int i = 0; i < _row; i++) {
                if (entete) {
                    std::string s;
                    myfile2 >> s;
                }
                for (int j = 0; j < _column; j++) {
                    double v;
                    myfile2 >> v;
                    set(i, j, v);
                    //std::cout << v << std::endl;
                }
            }
            myfile2.close();
        }
        else {
            std::cout << filename << " or  " << "../../" + filename << std::endl;
            throw std::invalid_argument("can't open this file");
        }
    }
       

}



///////////////////////////////////////////////////////////////////////////////
// Addition
///////////////////////////////////////////////////////////////////////////////
void MatrixCPUD::add(MatrixCPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) + m->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPUD::increment(int i, int j, double add)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds");
    }
    _matrixCPU[i * _column + j] += add;
}
void MatrixCPUD::addVector(MatrixCPUD* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector");
    }
    if (v->getNCol() == 1) {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) + v->get(i, 0);
                this->set(i, j, r);
            }
        }
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) + v->get(0, j);
                this->set(i, j, r);
            }
        }
    }
}
void MatrixCPUD::add(double c)
{
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) + c;
            this->set(i, j, r);
        }
    }
}
void MatrixCPUD::add(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = m1->get(i, j) + m2->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPUD::add(MatrixCPUD* m, double c)
{
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = m->get(i, j) + c;
            this->set(i, j, r);
        }
    }

}
void MatrixCPUD::addTrans(MatrixCPUD* m)
{
 
    if (_row != m->getNCol() && _column != m->getNLin()) 
    {
        throw std::invalid_argument("not the same dimension (transpose)");
    }
    for (int i = 0; i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) + m->get(j, i);
            set(i, j, r);
        }
    }
   
}
///////////////////////////////////////////////////////////////////////////////
// subtraction
///////////////////////////////////////////////////////////////////////////////
void MatrixCPUD::subtract(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = m1->get(i, j) - m2->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPUD::subtractRow(int row1, int row2, double factor)
{
    if ((row1 < _row) && (row2 < _row)) {
        for (int j = 0; j < _column; j++) {
            set(row1, j, get(row1, j) - factor * get(row2, j));
        }
    }
    else {
        throw std::invalid_argument("out of bound");
    }

}

void MatrixCPUD::subtractAbs(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if (!dim(m1) || !dim(m2)) {
        throw std::invalid_argument("not the same dimension");
    }
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            double r = fabs(m1->get(i, j)) - fabs(m2->get(i, j));
            set(i, j, r);
        }
    }
}

void MatrixCPUD::subtract(MatrixCPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) - m->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPUD::subtractVector(MatrixCPUD* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector");
    }
    if (v->getNCol() == 1) {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) - v->get(i, 0);
                this->set(i, j, r);
            }
        }
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) - v->get(0, j);
                this->set(i, j, r);
            }
        }
    }

}


void MatrixCPUD::subtractTrans(MatrixCPUD* m)
{
    MatrixCPUD temp(*this);
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("not the same dimension (transpose)");
    }
    for (int i = 0;i < _row;i++)
    {
        for (int j = 0;j < _column;j++)
        {
            double r = get(i, j) - m->get(j, i);
            temp.set(i, j, r);
        }
    }
    this->set(&temp);
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication
///////////////////////////////////////////////////////////////////////////////
void MatrixCPUD::multiply(MatrixCPUD* m1,MatrixCPUD* m2)
{
    if ((m1->getNCol() != m2->getNLin()) || (_row != m1->getNLin()) || (_column != m2->getNCol()) )
    {
        throw std::invalid_argument("not the good dimension");
    }
    double r = 0;
    int p = m1->getNCol();
    for (int i = 0; i < _row;++i)
    {
        for (int j = 0; j < _column;++j)
        {
            r = 0;
            for (int k = 0; k < p;++k)
            {
                r += m1->get(i, k) * m2->get(k, j);
            }
            this->set(i, j, r);
        }
    }
}

void MatrixCPUD::multiplyTrans(MatrixCPUD* m1, MatrixCPUD* m2, int numToTrans)
{
    if (numToTrans) {
        if ((m1->getNCol() != m2->getNCol()) || (_row != m1->getNLin()) || (_column != m2->getNLin()))
        {
            throw std::invalid_argument("not the good dimension");
        }
        double r = 0;
        int p = m1->getNCol();
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _column; ++j)
            {
                r = 0;
                for (int k = 0; k < p; ++k)
                {
                    r += m1->get(i, k) * m2->get(j, k);
                }
                this->set(i, j, r);
            }
        }
    }
    else {
        if ((m1->getNLin() != m2->getNLin()) || (_row != m1->getNCol()) || (_column != m2->getNCol()))
        {
            throw std::invalid_argument("not the good dimension");
        }
        double r = 0;
        int p = m1->getNLin();
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _column; ++j)
            {
                r = 0;
                for (int k = 0; k < p; ++k)
                {
                    r += m1->get(k, i) * m2->get(k, j);
                }
                this->set(i, j, r);
            }
        }
    }
    
}

void MatrixCPUD::multiplyDiag(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if((m1->getNCol() != m2->getNLin()) || (m2->getNCol() != m1->getNLin()) || (_row != m1->getNLin()) || (_column != 1))
    {
        throw std::invalid_argument("not the good dimension");
    }
    double r = 0;
    int p = m1->getNCol();
    for (int i = 0; i < _row; ++i)
    {
        r = 0;
        for (int k = 0; k < p; ++k)
        {
            r += m1->get(i, k) * m2->get(k, i);
        }
        this->set(i, 0, r);
        
    }
}

void MatrixCPUD::multiply(double c)
{
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) * c;
            this->set(i, j, r);
        }
    }
        
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication Terme � Terme
///////////////////////////////////////////////////////////////////////////////

void MatrixCPUD::multiplyT(MatrixCPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
       
    }

    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j) * m->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPUD::multiplyT(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
       
    }

    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = m1->get(i, j) * m2->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPUD::multiplyTVector(MatrixCPUD* m, MatrixCPUD* v, int sens)
{
    // vector v can be a row or a column vector it doesn't matter, it the sens which decide the calculation
    // ex : alpha.multiplyTVector(m,v); alpha : l*n , m : l*n, v : 1*n or n*1 -> sens = 0 
    double s = 0;
    if (!dim(m)) {
        throw std::invalid_argument("matrices must have the same size");
    } 
    if (sens) {
        if (v->getNCol() == 1 && v->getNLin() == _row)
        { // vecteur colonne
            for (int i = 0; i < _row; i++)
            {
                for (int j = 0; j < _column; j++)
                {
                    set(i, j, m->get(i, j) * v->get(i, 0));
                }

            }

        }
        else if (v->getNLin() == 1 && v->getNCol() == _row)
        { // vecteur ligne
            for (int i = 0; i < _row; i++)
            {
                for (int j = 0; j < _column; j++)
                {
                    set(i, j, m->get(i, j) * v->get(0, i));
                }

            }

        }
        else {
            throw std::invalid_argument("wrong size of the vector");
        }
    }
    else {
        if (v->getNCol() == 1 && v->getNLin() == _column)
        { // vecteur colonne
            for (int i = 0; i < _row; i++)
            {
                for (int j = 0; j < _column; j++)
                {
                    set(i, j, m->get(i, j) * v->get(j, 0));
                }

            }

        }
        else if (v->getNLin() == 1 && v->getNCol() == _column)
        { // vecteur ligne
            for (int i = 0; i < _row; i++)
            {
                for (int j = 0; j < _column; j++)
                {
                    set(i, j, m->get(i, j) * v->get(0, j));
                }

            }

        }
        else {
            throw std::invalid_argument("wrong size of the vector");
        }
    }
    
}


void MatrixCPUD::divide(double c)
{
    if (c == 0) {
        throw std::domain_error("divide by 0");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            double r = get(i, j)/c;
            this->set(i, j, r);
        }
    }
}

void MatrixCPUD::divideT(MatrixCPUD* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    MatrixCPUD temp(*this);
    double r = 0;
    double f = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            f = m->get(i, j);
            if (f == 0) {
                throw std::domain_error("divide by 0");
            }
            r = get(i, j) / f;
            temp.set(i, j, r);
        }
    }
    set(&temp);
}

void MatrixCPUD::invertGaussJordan(MatrixCPUD* mToInvert)
{
    MatrixCPUD m(*mToInvert);
    if (!dim(&m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_row != _column) {
        throw std::invalid_argument("must be a square matrix");
    }
    MatrixCPUD augmented(_row, _column);
    augmented.setEyes(1);
    MatrixCPUD indices(1, 2);
    int r = 0;
    for (int column = 0; column < _column; column++) {
        double pivotAbs = m.maxAbs(r, _row, column, column + 1, &indices);
        int k = indices.get(0, 0); // indice max de la colonne j
        double pivot = m.get(k, column);
        if (pivotAbs < 0.000001f) {
            throw std::invalid_argument("not invertible matrix");
        }
        else {

            for (int j = 0; j < _column; j++) {
                augmented.set(k, j, augmented.get(k, j) / pivot);
                m.set(k, j, m.get(k, j) / pivot);
            }
            if (k != r) {
                augmented.swapLine(k, r);
                m.swapLine(k, r);
            }
            for (int i = 0; i < _row; i++) {
                if (i != r) {
                    double local = m.get(i, column);
                    //std::cout << "substrat row i= " << i << " r =" << r << " column =" << column << " local =" << local << std::endl;
                    //m->display();
                    m.subtractRow(i, r, local);
                    
                    augmented.subtractRow(i, r, local);
                }
            }
            r++;
        }
    }
    
    set(&augmented);
}


void MatrixCPUD::LDLFactorization(MatrixCPUD* L, MatrixCPUD* D)
{
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("A must be symetrical and so must be square");
    }
    throw std::invalid_argument("Pas encore ecrite");
    
    
    int n = getNLin();
    *L = MatrixCPUD(n, n);
    *D = MatrixCPUD(n, n);

    L->set(0, 0, 1);
    D->set(0, 0, get(0, 0));

    for (int k = 1; k < n; k++) {
        // L1 : k - 1, 1 : k - 1 y = A1 : k - 1, k pour trouver y1 : k - 1
        // L k, 1 : k - 1 = (D - 11 :k - 1, 1 : k - 1 y1 : k - 1)T pour trouver L k, 1 : k - 1
        // lkk = 1
        //dkk = akk � Lk, 1 : k - 1 y1 : k - 1
        //dkk - 1 = 1 / dkk
    }




}

void MatrixCPUD::LUPFactorization(MatrixCPUD* A, MatrixCPUD* Po)
{
    double Tol = 0.0000001;
    int n = getNLin();
    A->set(this);
    
    // code from wikipedia adapted
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    if (Po->getNCol() != 1 || Po->getNLin() != (getNCol() + 1)) {
        throw std::invalid_argument("wrong size of P");
    }
    


    for (int i = 0; i < n; i++) {
        Po->set(i, 0, i); //Unit permutation matrix, P[N] initialized with N
    }
    
    double absA = 0;
    int j = 0;
    for (int col = 0; col < n; col++) {
        double maxA = 0.0;
        int imax = col;
        for (int k = col; k < n; k++){
            absA = fabs(A->get(k, col));
            if (absA  > maxA)
            {
                maxA = absA;
                imax = k;
            }
        }
        if (maxA < Tol) {
            throw std::invalid_argument("matrix is degenerate");
        }

        if (imax != col) { //le max pas sur la diagonal
            //pivoting P
            j = Po->get(col,0);
            Po->set(col, 0, Po->get(imax, 0));
            Po->set(imax, 0, j);

            //pivoting rows of A
            A->swapLine(col, imax);

            //counting pivots starting from N (for determinant)
            Po->increment(n, 0, 1);
        }

       

           
        for (int i = col + 1; i < n; i++) { 
            /*double lij = U.get(i, col) / U.get(col, col);
            L.set(i, col, 1);
            U.set(i, col, lij);*/ 
            A->set(i, col, A->get(i, col) / A->get(col, col)); //A[j][i] /= A[i][i];
            

            for (int k = col + 1; k < n; k++) {
                A->set(i, k, A->get(i, k) - A->get(i, col) * A->get(col, k)); //A[j][k] -= A[j][i] * A[i][k]; 
            }  
        }
    }
    // en vrai on peut tout stocker dans une matrice comme on sait que diag(L) = Id, Et donc on peut avoir A = (L-Id) + U
    
    


}

void MatrixCPUD::solveSysUpper(MatrixCPUD* U, MatrixCPUD* y)
{
    //n'utilise pas y !!!!! petit probleme...
    if (getNLin() != U->getNCol() || U->getNLin() != y->getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    int n = getNLin();
    for (int i = n - 1; i >= 0; i--)
    {
        
        for (int k = i + 1; k < n; k++) {
            increment(i, 0, -U->get(i, k) * get(k, 0)); // x[i] -= A[i][k] * x[k];
        }
            
        set(i, 0, get(i, 0) / U->get(i, i));
    }

}

void MatrixCPUD::solveSysLower(MatrixCPUD* L, MatrixCPUD* b, MatrixCPUD* P)
{
    if (getNLin() != L->getNCol() || L->getNLin() != b->getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    int n = getNLin();
    for (int i = 0; i < n; i++) {
        set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];

        for (int k = 0; k < i; k++) {
            increment(i, 0, - L->get(i, k) * get(k, 0));
        }
    }
}

void MatrixCPUD::solveSys(MatrixCPUD* A, MatrixCPUD* P, MatrixCPUD* b)
{
    if (A->getNCol() != A->getNLin() || A->getNLin() != b->getNLin()) {
        throw std::invalid_argument("wrong size of A");
    }
    if (b->getNCol() != 1) {
        throw std::invalid_argument("b must be a column vector");
    }
    if (P->getNLin() != (A->getNCol() + 1) || P->getNCol() != 1) {
        throw std::invalid_argument("wrong size of P");
    }
    /*std::cout << " A=" << std::endl;
    A->display();
    std::cout << " b=" << std::endl;
    b->display();*/
    int n = getNLin();
    for (int i = 0; i < n; i++) {
        set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];
    }
    /*std::cout << " permute b= or x" << std::endl;
    display();*/
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < i; k++) {
            increment(i, 0, -A->get(i, k) * get(k, 0));
        }
    }
    //std::cout << " y = " << std::endl;
    //display();
    for (int i = n - 1; i >= 0; i--)
    {

        for (int k = i + 1; k < n; k++) {
            increment(i, 0, -A->get(i, k) * get(k, 0)); // x[i] -= A[i][k] * x[k];
        }

        set(i, 0, get(i, 0) / A->get(i, i));
    }
    //std::cout << " sol = " << std::endl;
   // display();
    

    /*y.solveSysLower(A, b, P);
    solveSysUpper(A, &y);*/

}


void MatrixCPUD::solveSysGaussSeidel(MatrixCPUD* M, MatrixCPUD* b)
{
    // M doit �tre sym�trique definie positive... wait what ?

}

void MatrixCPUD::MultiplyMatVec(MatrixCPUD* m, MatrixCPUD* vect, int sens)
{
    if (getNLin() != m->getNLin() || getNCol() != 1) {
        throw std::invalid_argument("wrong dimension of this");
    }
    if (sens) {
        if (m->getNCol() != vect->getNCol() || vect->getNLin() != 1) {
            throw std::invalid_argument("wrong dimension of the row vector, (sens = 1) ");
        }
        for (int i = 0; i < getNLin(); i++) {
            double s = 0;
            for (int j = 0; j < m->getNCol(); j++) {
                s = s + m->get(i, j) * vect->get(0, j);
            }
            set(i, 0, s);
        }
    }
    else {
        if (m->getNCol() != vect->getNLin() || vect->getNCol() != 1) {
            throw std::invalid_argument("wrong dimension of the column vector, sens =0 ");
        }
        for (int i = 0; i < getNLin(); i++) {
            double s = 0;
            for (int j = 0; j < m->getNCol(); j++) {
                s = s + m->get(i, j) * vect->get(j, 0);
            }
            set(i, 0, s);
        }

    }

}

void MatrixCPUD::MultiplyMatMat(MatrixCPUD* m1, MatrixCPUD* m2)
{
    if (getNLin() != m1->getNLin() || getNCol() != m2->getNCol()) {
        throw std::invalid_argument("wrong dimension of this");
    }
      
    if (m1->getNCol() != m2->getNLin()) {
        throw std::invalid_argument("wrong dimension of the matrices ");
    }
    for (int i = 0; i < getNLin(); i++) {
        for (int j = 0; j < getNCol(); j++) {
            double s = 0;
            for (int p = 0; p < m1->getNCol(); p++) {
                s = s + m1->get(i, p) * m2->get(p, j);
            }
            set(i, j, s); 
        } 
    }

}


///////////////////////////////////////////////////////////////////////////////
// Operation complexe
///////////////////////////////////////////////////////////////////////////////




void MatrixCPUD::project(MatrixCPUD* Lb, MatrixCPUD* Ub)
{
    if (!dim(Lb) || !dim(Ub)) {
        throw std::invalid_argument("not the same dimension");
    }
    double ub = 0;
    double lb = 0;
    double r = 0;
    MatrixCPUD temp(*this);
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            r = get(i, j);
            ub = Ub->get(i, j);
            lb = Lb->get(i, j);
            if (ub < lb) {
                throw std::invalid_argument("impossible to have a value for the projection, ub>lb");
            }
            r = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r; // permet de ne pas faire de branchement if.
            temp.set(i, j, r);
        }
    }
    this->set(&temp);
}

void MatrixCPUD::projectNeg()
{
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            double r = get(i, j);
            r = (r < 0) * r;
            set(i, j, r);
        }
    }
}

void MatrixCPUD::projectPos()
{
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            double r = get(i, j);
            r = (r > 0) * r;
            set(i, j, r);
        }
    }
}




double MatrixCPUD::distance2(MatrixCPUD* m1) const
{
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    double d = 0;
    double r = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            r = get(i, j) - m1->get(i, j);
            d = d + r * r;
        }
    }
    return sqrtf(d);
}

double MatrixCPUD::distance2() const
{
    double d = 0;
    double r = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            r = get(i, j);
            d = d + r * r;
        }
    }
    return sqrtf(d);
    
}

double MatrixCPUD::max2() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    double M = fabs(get(0, 0));
    double m = 0;
    for (int i = 0; i < _row;++i)
    {
        for (int j = 0; j < _column;++j)
        {
            m = fabs(get(i, j));
            if (m > M) {
                M = m;
            }
        }
    }
    return M;
}

double MatrixCPUD::maxAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* indices)
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix");
    }
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    }
    double M = fabs(get(iBegin, jBegin));
    if (indices != nullptr) {
        indices->set(0, 0, iBegin);
        indices->set(0, 1, jBegin);
    }
    double m = 0;
    for (int i = iBegin;i < iEnd;++i)
    {
        for (int j = jBegin;j < jEnd;++j)
        {
            m = fabs(get(i, j));
            if (m > M) {
                M = m;
                if (indices != nullptr) {
                    indices->set(0, 0, i);
                    indices->set(0, 1, j);
                }
            }
        }
    }
    return M;
}

double MatrixCPUD::minAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* indices, bool Null)
{
    // Null == False : retourne le minimum non null ! si tout est nul retourne 0
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix");
    }
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    }
    double M = INFINITY;
    double m = 0;
    for (int i = iBegin;i < iEnd;++i)
    {
        for (int j = jBegin;j < jEnd;++j)
        {
            m = fabs(get(i, j));
            if (!Null && m==0) {
                
            }
            else {
                if (m < M) {
                    M = m;
                    if (indices != nullptr) {
                        indices->set(0, 0, i);
                        indices->set(0, 1, j);
                    }
                }
            }
        }
    }
    return (M!=INFINITY)*M; // retourne M si pas infini et 0 sinon
}

void MatrixCPUD::Moy(MatrixCPUD* m, MatrixCPUD* nb, int sens)
{   
    double s;
    int n;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol() ) || (_column != nb->getNCol()) || (nb->getNLin() !=1) ) 
        {
            throw std::invalid_argument("wrong dimension of the vector");
        }
        for (int j = 0; j < _column;j++)
        {
            n = nb->get(0, j);
            s = 0;
            if (n > 0) 
            {
                for (int i = 0; i < m->getNLin();i++)
                {
                    s = s + m->get(i, j);
                }
                s = s / n;
            }
            set(0, j, s);
        }

    }
    else { // on travaille sur les lignes 
        if ((_column != 1) || (_row != m->getNLin()) || (_row != nb->getNLin()) || (nb->getNCol() != 1)) {
            throw std::invalid_argument("wrong dimension of the vector");
        }
        for (int i = 0; i < _row;i++) 
        {
            n = nb->get(i, 0);
            s = 0;
            if (n > 0) {
                for (int j = 0; j < m->getNCol();j++) 
                {
                    s = s + m->get(i, j);
                }
                s = s / n;
            }            
            set(i, 0, s);
        }
    }


}

double MatrixCPUD::sum() const
{
    double d = 0;
    double r = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            r = get(i, j);
            d = d + r;
        }
    }
    return d;
}

double MatrixCPUD::sumD() const
{
    double d = 0;
    double r = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            r = get(i, j);
            d = d + r;
        }
    }
    return d;
}

void MatrixCPUD::sum(MatrixCPUD* m, int sens)
{
    double s = 0;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol()) )
        {
            throw std::invalid_argument("wrong dimension of the line vector, sens==1");
        }
        for (int j = 0; j < _column;j++)
        {
            s = 0;
            for (int i = 0; i < m->getNLin();i++)
            {
                s = s + m->get(i, j);
            }
            set(0, j, s);
        }

    }
    else { // on travaille sur les lignes 
        if ((_column != 1) || (_row != m->getNLin())) {
            throw std::invalid_argument("wrong dimension of the column vector (sens==0)");
        }
        for (int i = 0; i < _row;i++)
        {
            s = 0;
            for (int j = 0; j < m->getNCol();j++)
            {
                s = s + m->get(i, j);
            }
            set(i, 0, s);
        }
    }


}
void MatrixCPUD::sumT(MatrixCPUD* m, int sens)
{
    double s = 0;
    if (sens) { // on travaille sur les colonnes
        if ((_column != 1) || (_row != m->getNCol()))
        {
            throw std::invalid_argument("wrong dimension of the line vector, sens==1");
        }
        for (int i = 0; i < _row; i++)
        {
            s = 0;
            for (int j = 0; j < m->getNLin(); j++)
            {
                s = s + m->get(j, i);
            }
            set(i, 0, s);
        }

    }
    else { // on travaille sur les lignes 
        if ((_row != 1) || (_column != m->getNLin())){
            throw std::invalid_argument("wrong dimension of the column vector (sens==0)");
        }
        for (int i = 0; i < _column; i++)
        {
            s = 0;
            for (int j = 0; j < m->getNCol(); j++)
            {
                s = s + m->get(i, j);
            }
            set(0, i, s);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
// Tri de la matrix
///////////////////////////////////////////////////////////////////////////////



void MatrixCPUD::sort(int dim, int sens)
{
    if (sens) { // on travaille sur les colonnes
        if (dim >= _row) {
            throw std::invalid_argument("dim must be smaller then the matrix dimension");
        }

        int milieu = _column / 2;
        sortColumn(0, milieu, dim);
        sortColumn(milieu, _column, dim);
        fusionColumn(0, _column, dim);
    }
    else {
        if (dim >= _column) {
            throw std::invalid_argument("dim must be smaller then the matrix dimension");
        }
        int milieu = _row / 2;
        sortLine(0, milieu, dim);
        sortLine(milieu, _row, dim);
        fusionLine(0, _row, dim);
    }

}

void MatrixCPUD::sortLine(int line1, int line2, int dim)
{
    
    if ((line2 - line1) <= 1) {
    }
    else if ((line2 - line1) <= 2) {
        if (get(line1, dim) > get(line2 - 1, dim)) {
            swapLine(line1, line2 - 1);
        }
    }
    else {
        int milieu = (line1 + line2) / 2;
        sortLine(line1, milieu, dim);
        sortLine(milieu, line2, dim);
        fusionLine(line1, line2, dim);
    }
}

void MatrixCPUD::sortColumn(int col1, int col2, int dim)
{
    
    if ((col2 - col1) <= 1) {
        
    }
    else if ((col2 - col1) <= 2) {
        if (get(dim, col1) > get(dim, col2-1)) {
            swapColumn(col1, col2-1);
        }
    }else {
        int milieu = (col2 + col1) / 2;
        sortColumn(col1, milieu, dim);
        sortColumn(milieu, col2, dim);
        fusionColumn(col1, col2, dim);
    }
}

void MatrixCPUD::fusionLine(int line1, int line2, int dim)
{
    
    int milieu = (line2 - line1) / 2;
    int taille = line2 - line1;
    int* changement = new int[taille];
    int* reelPos = new int[taille];

    int k = 0;
    int l = milieu;
    int i = 0;
    while (k < milieu && l < taille) {
        if (get(k, dim) > get(l, dim)) {
            changement[i] = l;
            l++;
        }
        else {
            changement[i] = k;
            k++;
        }
        reelPos[i] = i;
        i++;
    }
    if (k == milieu) {
        for (int j = i; j < taille; j++) {
            changement[j] = l;
            l++;
            reelPos[j] = j;
        }
    } else {
        for (int j = i; j < taille; j++) {
            changement[j] = k;
            k++;
            reelPos[j] = j;
        }
    }
    int indice = 0;
    while (indice < taille) {
        if (changement[indice] != reelPos[indice]) {
            swapLine(indice + line1, changement[indice] + line1);
            reelPos[changement[indice]] = reelPos[indice];
        }
        indice++;
    }

    DELETEA(changement);
    DELETEA(reelPos);
}

void MatrixCPUD::fusionColumn(int col1, int col2, int dim)
{

    int milieu = (col2 - col1) / 2;
    int taille = col2 - col1;
    int* changement = new int[taille];
    int* reelPos = new int[taille];

    int k = 0;
    int l = milieu;
    int i = 0;
    while(k<milieu && l<taille){
        if (get(dim, k) > get(dim, l)) {
            changement[i] = l;
            
            l++;
        }
        else {
            changement[i] = k;
            k++;
        }
        reelPos[i] = i;
        i++;
    }
    if (k == milieu) {
        for (int j = i; j < taille; j++) {
            changement[j] = l;
            l++;
            reelPos[j] = j;
        }
    }
    else {
        for (int j = i; j < taille; j++) {
            changement[j] = k;
            k++;
            reelPos[j] = j;
        }
    }
    int indice = 0;
    while (indice < taille) {
        if (changement[indice] != reelPos[indice]) {
            swapColumn(indice+col1, changement[indice]+col1);
            reelPos[changement[indice]] = reelPos[indice];
        }
        indice++;
    }

    DELETEA(changement);
    DELETEA(reelPos);
}

void MatrixCPUD::swapLine(int line1, int line2)
{
    
    double temp = 0;
    for (int i = 0; i < _column; i++) {
        temp = get(line1, i);
        set(line1, i, get(line2, i));
        set(line2, i, temp);
    }
}

void MatrixCPUD::swapColumn(int col1, int col2)
{
    double temp = 0;
    for (int i = 0; i < _row; i++) {
        temp = get(i, col1);
        set(i, col1, get(i, col2));
        set(i, col2, temp);
    }
}

void MatrixCPUD::RelativeEror(MatrixCPUD* MatRef, MatrixCPUD* Mat)
{
    if (!dim(MatRef) || !dim(Mat)) {
        std::cout << _row << " " << _column << " " << MatRef->_row << " " << MatRef->_column 
            << " " << Mat->_row << " " << Mat->_column << " " << std::endl;
        throw std::invalid_argument("matrix must have the same size");
    }
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            double a = MatRef->get(i, j);
            double b = Mat->get(i, j);
            if (a != 0) {
                set(i, j, abs((a - b) / a));
            }
            else if (b != 0) {
                set(i, j, abs((b - a) / b));
            }
            else {
                set(i, j, abs(a - b));
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Display MatrixCPUD contents
///////////////////////////////////////////////////////////////////////////////
void MatrixCPUD::display() const
{   
    if (this) {
        if (_row == 0 || _column == 0)
        {
            std::cout << "matrix vide " << std::endl;
        }
        if (_column == 1) {
            std::cout << " transpose  : ";
            for (int i = 0;i < _row;++i)
            {
                for (int j = 0;j < _column;++j)
                {
                    double value = get(i, j);
                    std::cout << std::setprecision(6) << value;
                    //std::cout << std::fixed << std::setprecision(2) << value;
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }
        else {
            for (int i = 0;i < _row;++i)
            {
                for (int j = 0;j < _column;++j)
                {
                    double value = get(i, j);
                    std::cout << std::setprecision(6) << value;
                    //std::cout << std::fixed << std::setprecision(3) << value;
                    if (j != _column - 1) std::cout << " ";
                }

                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout << "matrix non definie " << std::endl;
    }

}

void MatrixCPUD::displayBloc(int iBegin, int iEnd, int jBegin, int jEnd) const
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    }
    if (this) {
        if (_row == 0 || _column == 0)
        {
            std::cout << "matrix vide " << std::endl;
        }
        if (jEnd - jBegin == 1) {
            std::cout << " transpose  : ";
            for (int i = iBegin; i < iEnd; ++i)
            {
                
                double value = get(i, jBegin);
                std::cout << std::setprecision(6) << value;
                //std::cout << std::fixed << std::setprecision(2) << value;
                std::cout << " ";
                
            }
            std::cout << std::endl;
        }
        else {
            for (int i = iBegin; i < iEnd; ++i)
            {
                for (int j = jBegin; j < jEnd; ++j)
                {
                    double value = get(i, j);
                    std::cout << std::setprecision(6) << value;
                    //std::cout << std::fixed << std::setprecision(3) << value;
                    if (j != jEnd - 1) std::cout << " ";
                }

                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout << "matrix non definie " << std::endl;
    }
}



void MatrixCPUD::saveCSV(const std::string& filename, std::ios_base::openmode mode, int trans) const
{
    std::ofstream myfile;
    myfile.open(filename, mode);
    myfile.precision(10);
    if (!trans) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column;j++) {
                myfile << get(i, j) << ";";
            }
            myfile << "\n";
        }
    }
    else {
        for (int j = 0; j < _column;j++) {
            for (int i = 0; i < _row; i++) {
                myfile << get(i, j) << ";";
            }
            myfile << "\n";
        }
    }
    
    myfile.close();
}

///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
MatrixCPUD::~MatrixCPUD()
{
    #ifdef DEBUG_DESTRUCTOR
        std::cout << "destruction matrix " << _matrixCPU << std::endl;
    #endif // DEBUG_DESTRUCTOR

    
    /*if (_matrixCPU != nullptr) {
        delete[] _matrixCPU;
        _matrixCPU = nullptr;
    }*/
   
   
    DELETEA(_matrixCPU);

}

///////////////////////////////////////////////////////////////////////////////
// Eigen
///////////////////////////////////////////////////////////////////////////////

#ifdef EIGEN
    void MatrixCPUD::set(Eigen::MatrixXd* eigenMatrix)
{
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
             set(i, j, (*eigenMatrix)(i, j));
        }
    }
}
   void MatrixCPUD::invertEigen(MatrixCPUD* mToInvert)
{
    if (!dim(mToInvert)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_row != _column) {
        throw std::invalid_argument("must be a square matrix");
    }

    Eigen::MatrixXd eigenMat(_row, _column);
    mToInvert->toEigenMatrix(&eigenMat);
    eigenMat = eigenMat.inverse();
    set(&eigenMat);
}
void MatrixCPUD::solveSysEigen(MatrixCPUD* M, MatrixCPUD* b)
{
    Eigen::MatrixXd eigenM(M->getNLin(), M->getNCol());
    M->toEigenMatrix(&eigenM);
    Eigen::MatrixXd eigenb(b->getNLin(), 1);
    b->toEigenMatrix(&eigenb);
   

    Eigen::MatrixXd x = (eigenM.partialPivLu()).solve(eigenb);
    
    set(&x);
}
void MatrixCPUD::toEigenMatrix(Eigen::MatrixXd* eigenMatrix)
{
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            (*eigenMatrix)(i, j) = get(i, j);
        }
    }
}


#endif



#ifdef OSQP
 c_float* MatrixCPUD::toCFloat()
{
    c_float* m = new c_float[_row * _column];

    for (int i = 0; i < _row * _column; i++) {
        m[i] = _matrixCPU[i];
    }

    return m;
}



c_int MatrixCPUD::toCSC(c_float* data, c_int* idx, c_int* ptr)
{
    int N = getNNull();
    int indiceData = 0;
 
    int indiceCol = 0;
    // cas particulier o� la/les premi�re colonnes sont nulles ????
    
    for (int j = 0; j < _column;j++) {
        ptr[j] = indiceCol;
        for (int i = 0; i < _row; i++) {
            if (fabs(get(i, j)) > 0.000000001) {
                data[indiceData] = get(i, j);
                idx[indiceData] =  i;
                indiceData++;
                indiceCol++;
            }
        }
    }
    ptr[_column] = N;
    return (c_int) N;
}

c_int MatrixCPUD::toCSCHalf(c_float* data, c_int* idx, c_int* ptr)
{
    if (_column != _row) {
        throw std::invalid_argument("the matrix must be square");
    }
    int N = getNNullHalf();

    int indiceData = 0;
    int indiceCol = 0;
    

    for (int j = 0; j < _column;j++) {
        ptr[j] = indiceCol;
        for (int i = 0; i <= j;i++) {
            if (fabs(get(i, j)) > 0.000000001) {
                data[indiceData] = get(i, j);
                idx[indiceData] =  i;
                indiceData++;
                indiceCol++;
            }
        }
    }
    ptr[_column] = N;
    return  N;
}


#endif
