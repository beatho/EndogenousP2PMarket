#include "../head/MatrixCPU.h"

float MatrixCPU::rand1()
{
    float a = (float)(rand()) / ((float)(RAND_MAX));
    return a;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MatrixCPU::MatrixCPU() {
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur appele" << std::endl;
#endif

}


MatrixCPU::MatrixCPU(int l, int c, float value)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur parametre appele" << std::endl;
    std::cout << _matrixCPU << std::endl;
#endif
    _row = l;
    _column = c;
    _matrixCPU = new float[l*c];
    for (int elem = 0; elem < l * c;elem++) {
        _matrixCPU[elem] = value; 
    }
#ifdef DEBUG_CONSTRUCTOR
    std::cout << _matrixCPU << std::endl;
#endif
}

MatrixCPU::MatrixCPU(const MatrixCPU & m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur recopie appele" << std::endl;
#endif
    _row = m._row;
    _column = m._column;
    _matrixCPU = new float[_row * _column];
    memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
}

MatrixCPU::MatrixCPU(const MatrixCPUD& m)
{
    _row = m.getNLin();
    _column = m.getNCol();
 
    if (_row * _column > 0) {
        _matrixCPU = new float[_row * _column];
    }
     

    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, (float)m.get(i, j));
        }
    }


}

MatrixCPU& MatrixCPU::operator=(const MatrixCPU& m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur operateur = appele" << std::endl;
#endif
    _row = m._row;
    _column = m._column;
    DELETEA(_matrixCPU);
    if (_row * _column > 0) {
        _matrixCPU = new float[_row * _column];
        memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
    }
    return *this;
}

MatrixCPU& MatrixCPU::operator=(const MatrixCPUD& m)
{
    _row = m.getNLin();
    _column = m.getNCol();
    DELETEA(_matrixCPU);
    if (_row * _column > 0) {
        _matrixCPU = new float[_row * _column];
    }
    

    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, (float)m.get(i, j));
        }
    }

    return *this;
}

void MatrixCPU::toMatCPUD(MatrixCPUD& m) const
{
    if (m.getNCol() != _column || m.getNLin() != _row) {
        m.setSize(_row, _column);
    }

    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++)
        {
            m.set(i, j, get(i, j));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////
inline float MatrixCPU::get(int i, int j) const
{
    if ((i >= _row) || ( j >= _column) || (i < 0) || ( j < 0)) {
        std::cout << "get " << _row << " " << _column << " " << i << " " << j << std::endl;
           throw std::out_of_range("index out of bounds (get)");
    }
    return _matrixCPU[i*_column+j];
}

double MatrixCPU::getD(int i, int j) const
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds (getD)");
    }
    return _matrixCPU[i * _column + j];
}

int MatrixCPU::getNCol() const
{
    return _column;
}

int MatrixCPU::getNLin() const
{
    return _row;
}

bool MatrixCPU::dim(MatrixCPU* m) const
{ 
    return ((_row == m->getNLin()) && (_column == m->getNCol()));
}



void MatrixCPU::getLin(MatrixCPU* vector, int i) const
{
    if ((_column != vector->getNCol()) || ( vector->getNLin() != 1)) 
    {
        std::cout << _column << " " << vector->getNCol() << " " << vector->getNLin() << std::endl;
        throw std::invalid_argument("wrong dimension of the row vector, (getLin)");
    } 
    for (int j = 0; j < _column; j++) {
        vector->set(0, j, get(i, j));
    }
}

bool MatrixCPU::isEqual(MatrixCPU* m, float pre) const
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension, (isEqual)");
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

int MatrixCPU::getNNull() const
{
    int n = 0;
    for (int i = 0;i < _row * _column;i++) {
        n = n + ((fabs(_matrixCPU[i]) > 0.000000001));
    }
    return n;
}

int MatrixCPU::getNNullHalf() const
{
    int n = 0;
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j <= i; j++) {
            n = n + (fabs(get(i, j)) > 0.000000001);
        }
    }

    return n;
}
void MatrixCPU::swap(MatrixCPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension (swap)");
    }
    
    float* temp = _matrixCPU;
    _matrixCPU = m->_matrixCPU;
    m->_matrixCPU = temp;
}

void MatrixCPU::getBloc(MatrixCPU* dest, int iBegin, int iEnd, int jBegin, int jEnd)
{

    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds (getBloc)");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd (getBloc)");
    } if (dest->getNLin() != (iEnd - iBegin) || dest->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension (getBloc)");
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


void MatrixCPU::setSize(int l, int c)
{
    DELETEA(_matrixCPU);
    _matrixCPU = new float[l * c];
    _row = l;
    _column = c;
}

///////////////////////////////////////////////////////////////////////////////
// Setter
///////////////////////////////////////////////////////////////////////////////
inline void MatrixCPU::set(int i, int j, float value)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << "err " << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds, (set)");
    }
    _matrixCPU[i * _column + j] = value;
}

void MatrixCPU::set(float value)
{
    for (int elem = 0; elem < _column * _row; elem++) {
        _matrixCPU[elem] = value;
    }
}

void MatrixCPU::set(MatrixCPU* m)
{
    if (!dim(m)) {
        std::cout << _row << " * " << _column << " against " << m->_row << " * " << m->_column << std::endl;
        throw std::invalid_argument("not the same dimension, (set)");
    }
    memcpy(_matrixCPU, m->_matrixCPU, _row * _column * sizeof(float));
}

void MatrixCPU::setTrans(MatrixCPU* m)
{
    if (_column != m->getNLin() || _row != m->getNCol()) {
        throw std::invalid_argument("not the same transposed dimension, (setTrans)");
    }

    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, m->get(j, i));
        }
    }
}

void MatrixCPU::set(Eigen::MatrixXd* eigenMatrix)
{
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
             set(i, j, (*eigenMatrix)(i, j));
        }
    }
}

void MatrixCPU::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        std::cout << " err out of bound " << iBegin << " " << jBegin << " " << iEnd << " " << jEnd << std::endl;
        throw std::out_of_range("index out of bounds, (setBloc)");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd (setBloc)");
    } if (m->getNLin() != (iEnd - iBegin ) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension, (setBloc)");
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

void MatrixCPU::setEyes(float v)
{
    int N = _row * (_row < _column) + _column * (_column <= _row);
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            set(i, j, v *(i==j));
        }
    }
 
}

void MatrixCPU::setEyes(MatrixCPU* vect)
{
    if (_column != _row) {
        throw std::invalid_argument("matrix must be square (setEyes)");
    }
    if (vect->_column == 1) {
        if (vect->_row != _row) {
            throw std::invalid_argument("matrix must have the same size as the vector (setEyes)");
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
            throw std::invalid_argument("matrix must have the same size as the vector (setEyes)");
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
        throw std::invalid_argument("argument must be a vector (setEyes)");
    }

}

void MatrixCPU::setRand(int eps, int divide)
{
    //srand(time(nullptr));
    //exit(1);
    int N = _column * _row;
    for (int elem = 0; elem < N; elem++) {
        int r = rand() % eps;
        int signe = (rand() % 2) * 2 - 1;
        _matrixCPU[elem] = (float) r / divide * signe;
       
    }


}

void MatrixCPU::setRand1(float eps)
{
    //srand(time(nullptr));

    int N = _column * _row;
    for (int elem = 0; elem < N; elem++) {
        _matrixCPU[elem] = 2 * (rand1() - 0.5) * eps;

    }
    exit(0);    

}

void MatrixCPU::setFromFile(std::string filename, int entete) // aucune sécurité, l'utilisateur doit vérifier d'utiliser le bon fichier avec le bon nombre d'entrée...
{
    std::ifstream myfile(filename, std::ios::in);
    //std::cout << filename << std::endl;
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
                float v;
                myfile >> v;
                set(i, j, v);
                //std::cout << v << std::endl;
            }
        }
        myfile.close();
    }
    else {
        std::ifstream myfile2("../../" + filename, std::ios::in);

        if (myfile2)
        {
            for (int i = 0; i < _row; i++) {
                if (entete) {
                    std::string s;
                    myfile2 >> s;
                }
                for (int j = 0; j < _column; j++) {
                    float v;
                    myfile2 >> v;
                    set(i, j, v);
                    //std::cout << v << std::endl;
                }
            }
            myfile2.close();
        }
        else {
            std::cout << filename << " or  " << "../../" + filename << std::endl;

            throw std::invalid_argument("can't open this file (setFromFile)");
        }
    }
       

}



///////////////////////////////////////////////////////////////////////////////
// Addition
///////////////////////////////////////////////////////////////////////////////
void MatrixCPU::add(MatrixCPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension (add)");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) + m->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPU::increment(int i, int j, float add)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds (increment)");
    }
    _matrixCPU[i * _column + j] += add;
}
void MatrixCPU::addVector(MatrixCPU* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector (addVector)");
    }
    if (v->getNCol() == 1) {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) + v->get(i, 0);
                this->set(i, j, r);
            }
        }
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) + v->get(0, j);
                this->set(i, j, r);
            }
        }
    }
}
void MatrixCPU::add(float c)
{
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) + c;
            this->set(i, j, r);
        }
    }
}
void MatrixCPU::add(MatrixCPU* m1, MatrixCPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension, m1 with m2 (add)" );
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension m with m1 (add)");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = m1->get(i, j) + m2->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPU::add(MatrixCPU* m, float c)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension this with m (add)");
    }
    for (int i = 0; i < _row;++i)
    {
        for (int j = 0; j < _column;++j)
        {
            float r = m->get(i, j) + c;
            this->set(i, j, r);
        }
    }

}
void MatrixCPU::addTrans(MatrixCPU* m)
{
 
    if (_row != m->getNCol() && _column != m->getNLin()) 
    {
        throw std::invalid_argument("not the same dimension (addTrans)");
    }
    for (int i = 0; i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) + m->get(j, i);
            set(i, j, r);
        }
    }
   
}
///////////////////////////////////////////////////////////////////////////////
// subtraction
///////////////////////////////////////////////////////////////////////////////
void MatrixCPU::subtract(MatrixCPU* m1, MatrixCPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension, (substract)");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension, (substract)");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = m1->get(i, j) - m2->get(i, j);
            this->set(i, j, r);
        }
    }
}
void MatrixCPU::subtractRow(int row1, int row2, float factor)
{
    if ((row1 < _row) && (row2 < _row)) {
        for (int j = 0; j < _column; j++) {
            set(row1, j, get(row1, j) - factor * get(row2, j));
        }
    }
    else {
        throw std::invalid_argument("out of bound , (substractRow)");
    }

}

void MatrixCPU::subtractAbs(MatrixCPU* m1, MatrixCPU* m2)
{
    if (!dim(m1) || !dim(m2)) {
        throw std::invalid_argument("not the same dimension , (substractAbs)");
    }
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            float r = fabs(m1->get(i, j)) - fabs(m2->get(i, j));
            set(i, j, r);
        }
    }
}

void MatrixCPU::subtract(MatrixCPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension , (substract)");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) - m->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPU::subtractVector(MatrixCPU* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector , (substractVector)");
    }
    if (v->getNCol() == 1) {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) - v->get(i, 0);
                this->set(i, j, r);
            }
        }
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) - v->get(0, j);
                this->set(i, j, r);
            }
        }
    }

}


void MatrixCPU::subtractTrans(MatrixCPU* m)
{
    MatrixCPU temp(*this);
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("not the same dimension (substractTrans)");
    }
    for (int i = 0;i < _row;i++)
    {
        for (int j = 0;j < _column;j++)
        {
            float r = get(i, j) - m->get(j, i);
            temp.set(i, j, r);
        }
    }
    this->set(&temp);
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication
///////////////////////////////////////////////////////////////////////////////
void MatrixCPU::multiply(MatrixCPU* m1, MatrixCPU* m2)
{
    if ((m1->getNCol() != m2->getNLin()) || (_row != m1->getNLin()) || (_column != m2->getNCol()) )
    {
        std::cout << m1->getNCol() << " " << m2->getNLin() << " " << _row << " " <<  m1->getNLin() << " " << _column << " " <<  m2->getNCol() << std::endl;
        throw std::invalid_argument("not the good dimension (multiply)");
    }
    float r = 0;
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

void MatrixCPU::multiplyTrans(MatrixCPU* m1, MatrixCPU* m2, int secondToTrans)
{
    if (secondToTrans) {
        if ((m1->getNCol() != m2->getNCol()) || (_row != m1->getNLin()) || (_column != m2->getNLin()))
        {
            throw std::invalid_argument("not the good dimension (multiplyTrans)");
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
            throw std::invalid_argument("not the good dimension (multiplyTrans)");
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

void MatrixCPU::multiplyDiag(MatrixCPU* m1, MatrixCPU* m2)
{
    if((m1->getNCol() != m2->getNLin()) || (m2->getNCol() != m1->getNLin()) || (_row != m1->getNLin()) || (_column != 1))
    {
        throw std::invalid_argument("not the good dimension (multiplyDiag)");
    }
    float r = 0;
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

void MatrixCPU::multiply(float c)
{
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) * c;
            this->set(i, j, r);
        }
    }
        
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication Terme à Terme
///////////////////////////////////////////////////////////////////////////////

void MatrixCPU::multiplyT(MatrixCPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension (multiplyT)");
       
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j) * m->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPU::multiplyT(MatrixCPU* m1, MatrixCPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension (multiplyT)");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension (multiplyT)");
       
    }

    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = m1->get(i, j) * m2->get(i, j);
            this->set(i, j, r);
        }
    }
}

void MatrixCPU::multiplyTVector(MatrixCPU* m, MatrixCPU* v, int sens)
{
    // vector v can be a row or a column vector it doesn't matter, it the sens which decide the calculation
    // ex : alpha.multiplyTVector(m,v); alpha : l*n , m : l*n, v : 1*n or n*1 -> sens = 0 
    float s = 0;
    if (!dim(m)) {
        throw std::invalid_argument("matrices must have the same size (multiplyTVector)");
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
            throw std::invalid_argument("wrong size of the vector (multiplyTVector)");
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
            throw std::invalid_argument("wrong size of the vector (multiplyTVector)");
        }
    }
    
}


void MatrixCPU::divide(float c)
{
    if (c == 0) {
        throw std::domain_error("divide by 0 (divide)");
    }
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            float r = get(i, j)/c;
            this->set(i, j, r);
        }
    }
}

void MatrixCPU::divideT(MatrixCPU* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension, (divideT)");
    }
    MatrixCPU temp(*this);
    float r = 0;
    float f = 0;
    for (int i = 0;i < _row;++i)
    {
        for (int j = 0;j < _column;++j)
        {
            f = m->get(i, j);
            if (f == 0) {
                throw std::domain_error("divide by 0, (divideT)");
            }
            r = get(i, j) / f;
            temp.set(i, j, r);
        }
    }
    set(&temp);
}

void MatrixCPU::invertGaussJordan(MatrixCPU* mToInvert)
{
    MatrixCPU m(*mToInvert);
    if (!dim(&m)) {
        throw std::invalid_argument("not the same dimension (invertGaussJordan)");
    }
    if (_row != _column) {
        throw std::invalid_argument("must be a square matrix (invertGaussJordan)");
    }
    MatrixCPU augmented(_row, _column);
    augmented.setEyes(1);
    MatrixCPU indices(1, 2);
    int r = 0;
    for (int column = 0; column < _column; column++) {
        float pivotAbs = m.maxAbs(r, _row, column, column + 1, &indices);
        int k = indices.get(0, 0); // indice max de la colonne j
        float pivot = m.get(k, column);
        if (pivotAbs < 0.000001f) {
            throw std::invalid_argument("not invertible matrix (invertGaussJordan)");
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
                    float local = m.get(i, column);
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
void MatrixCPU::invertEigen(MatrixCPU* mToInvert)
{
    if (!dim(mToInvert)) {
        throw std::invalid_argument("not the same dimension (invertEigen)");
    }
    if (_row != _column) {
        throw std::invalid_argument("must be a square matrix (invertEigen)");
    }

    Eigen::MatrixXd eigenMat(_row, _column);
    mToInvert->toEigenMatrix(&eigenMat);
    eigenMat = eigenMat.inverse();
    set(&eigenMat);
}

void MatrixCPU::LDLFactorization(MatrixCPU* L, MatrixCPU* D)
{
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("A must be symetrical and so must be square (LDLFactorization)");
    }
    throw std::invalid_argument("Pas encore ecrite (LDLFactorization)");
    
    
    int n = getNLin();
    *L = MatrixCPU(n, n);
    *D = MatrixCPU(n, n);

    L->set(0, 0, 1);
    D->set(0, 0, get(0, 0));

    for (int k = 1; k < n; k++) {
        // L1 : k - 1, 1 : k - 1 y = A1 : k - 1, k pour trouver y1 : k - 1
        // L k, 1 : k - 1 = (D - 11 :k - 1, 1 : k - 1 y1 : k - 1)T pour trouver L k, 1 : k - 1
        // lkk = 1
        //dkk = akk – Lk, 1 : k - 1 y1 : k - 1
        //dkk - 1 = 1 / dkk
    }




}

void MatrixCPU::LUPFactorization(MatrixCPU* A, MatrixCPU* Po)
{
    float Tol = 0.000001;
    int n = getNLin();
    A->set(this);
    
    // code from wikipedia adapted
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("A must be square (LUPFactorization)");
    }
    if (Po->getNCol() != 1 || Po->getNLin() != (getNCol() + 1)) {
        throw std::invalid_argument("wrong size of P (LUPFactorization)");
    }
    


    for (int i = 0; i < n; i++) {
        Po->set(i, 0, i); //Unit permutation matrix, P[N] initialized with N
    }
    
    float absA = 0;
    int j = 0;
    float maxA = 0.0;
    int imax = -1;
    for (int col = 0; col < n; col++) {
        //A->display();
        maxA = 0.0;
        imax = col;
        for (int row = col; row < n; row++){
            absA = fabs(A->get(row, col));
            
            if (absA  > maxA)
            {
                maxA = absA;
                imax = row;
            }
        }
        if (maxA < Tol) {
            throw std::invalid_argument("matrix is degenerate (LUPFactorization)");
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
            /*float lij = U.get(i, col) / U.get(col, col);
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

void MatrixCPU::solveSysUpper(MatrixCPU* U, MatrixCPU* y)
{
    //n'utilise pas y !!!!! petit probleme...
    if (getNLin() != U->getNCol() || U->getNLin() != y->getNLin()) {
        throw std::invalid_argument("A must be square (solveSysUpper)");
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

void MatrixCPU::solveSysLower(MatrixCPU* L, MatrixCPU* b, MatrixCPU* P)
{
    if (getNLin() != L->getNCol() || L->getNLin() != b->getNLin()) {
        throw std::invalid_argument("A must be square (solveSysLower)");
    }
    int n = getNLin();
    for (int i = 0; i < n; i++) {
        set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];

        for (int k = 0; k < i; k++) {
            increment(i, 0, - L->get(i, k) * get(k, 0));
        }
    }
}

void MatrixCPU::solveSys(MatrixCPU* A, MatrixCPU* P, MatrixCPU* b)
{
    if (A->getNCol() != A->getNLin() || A->getNLin() != b->getNLin()) {
        throw std::invalid_argument("wrong size of A (solveSys)");
    }
    if (b->getNCol() != 1) {
        throw std::invalid_argument("b must be a column vector (solveSys)");
    }
    if (P->getNLin() != (A->getNCol() + 1) || P->getNCol() != 1) {
        throw std::invalid_argument("wrong size of P (solveSys)");
    }

    int n = getNLin();
    for (int i = 0; i < n; i++) {
        set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];

        for (int k = 0; k < i; k++) {
            increment(i, 0, -A->get(i, k) * get(k, 0));
        }
    }
    for (int i = n - 1; i >= 0; i--)
    {

        for (int k = i + 1; k < n; k++) {
            increment(i, 0, -A->get(i, k) * get(k, 0)); // x[i] -= A[i][k] * x[k];
        }

        set(i, 0, get(i, 0) / A->get(i, i));
    }
    MatrixCPU y(*this);

    /*y.solveSysLower(A, b, P);
    solveSysUpper(A, &y);*/

}

void MatrixCPU::solveSysEigen(MatrixCPU* M, MatrixCPU* b)
{
    Eigen::MatrixXd eigenM(M->getNLin(), M->getNCol());
    M->toEigenMatrix(&eigenM);
    Eigen::MatrixXd eigenb(b->getNLin(), 1);
    b->toEigenMatrix(&eigenb);
   

    Eigen::MatrixXd x = (eigenM.partialPivLu()).solve(eigenb);
    
    set(&x);
}




void MatrixCPU::MultiplyMatVec(MatrixCPU* m, MatrixCPU* vect, int sens)
{
    if (getNLin() != m->getNLin() || getNCol() != 1) {
        throw std::invalid_argument("wrong dimension of this (MultiplyMatVec)");
    }
    if (sens) {
        if (m->getNCol() != vect->getNCol() || vect->getNLin() != 1) {
            throw std::invalid_argument("wrong dimension of the row vector, (sens = 1) (MultiplyMatVec)");
        }
        for (int i = 0; i < getNLin(); i++) {
            float s = 0;
            for (int j = 0; j < m->getNCol(); j++) {
                s = s + m->get(i, j) * vect->get(0, j);
            }
            set(i, 0, s);
        }
    }
    else {
        if (m->getNCol() != vect->getNLin() || vect->getNCol() != 1) {
            throw std::invalid_argument("wrong dimension of the column vector, sens =0 (MultiplyMatVec)");
        }
        for (int i = 0; i < getNLin(); i++) {
            float s = 0;
            for (int j = 0; j < m->getNCol(); j++) {
                s = s + m->get(i, j) * vect->get(j, 0);
            }
            set(i, 0, s);
        }
    }
}

void MatrixCPU::MultiplyMatTransVec(MatrixCPU* mToTrans, MatrixCPU* vect, int sensRow)
{
    if (getNLin() != mToTrans->getNCol() || getNCol() != 1) {
        throw std::invalid_argument("wrong dimension of this, (MultiplyMatTransVec)");
    }
    if (sensRow) {
        if (mToTrans->getNLin() != vect->getNCol() || vect->getNLin() != 1) {
            throw std::invalid_argument("wrong dimension of the row vector, (sens = 1) (MultiplyMatTransVec) ");
        }
        for (int i = 0; i < getNLin(); i++) {
            float s = 0;
            for (int j = 0; j < mToTrans->getNLin(); j++) {
                s = s + mToTrans->get(j, i) * vect->get(0, j);
            }
            set(i, 0, s);
        }
    }
    else {
        if (mToTrans->getNLin() != vect->getNLin() || vect->getNCol() != 1) {
            throw std::invalid_argument("wrong dimension of the column vector, sens =0 (MultiplyMatTransVec)");
        }
        for (int i = 0; i < getNLin(); i++) {
            float s = 0;
            for (int j = 0; j < mToTrans->getNLin(); j++) {
                s = s + mToTrans->get(j, i) * vect->get(j, 0);
            }
            set(i, 0, s);
        }
    }
}

void MatrixCPU::MultiplyMatMat(MatrixCPU* m1, MatrixCPU* m2)
{
    if (getNLin() != m1->getNLin() || getNCol() != m2->getNCol()) {
        std::cout << getNLin() << " " << m1->getNLin() << " " << getNCol() << "  " << m2->getNCol() << std::endl;
        throw std::invalid_argument("wrong dimension of this (MultiplyMatMat)");
    }
      
    if (m1->getNCol() != m2->getNLin()) {
        throw std::invalid_argument("wrong dimension of the matrices (MultiplyMatMat) ");
    }
    for (int i = 0; i < getNLin(); i++) {
        for (int j = 0; j < getNCol(); j++) {
            float s = 0;
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




void MatrixCPU::project(MatrixCPU* Lb, MatrixCPU* Ub)
{
    if (!dim(Lb) || !dim(Ub)) {
        throw std::invalid_argument("not the same dimension (project)");
    }
    float ub = 0;
    float lb = 0;
    float r = 0;
    MatrixCPU temp(*this);
    for (int i = 0; i < _row; i++)
    {
        for (int j = 0; j < _column; j++)
        {
            r = get(i, j);
            ub = Ub->get(i, j);
            lb = Lb->get(i, j);
            if (ub < lb) {
                std::cout << "'lb  ub : " << lb << " " << ub << std::endl;
                throw std::invalid_argument("impossible to have a value for the projection, ub>lb (project)");
            }
            r = ub * (ub < r) + lb * (lb > r) + r * (r >= lb) * (r <= ub); // permet de ne pas faire de branchement if.
            temp.set(i, j, r);
        }
    }
    this->set(&temp);
}

void MatrixCPU::projectNeg()
{
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            float r = get(i, j);
            r = (r < 0) * r;
            set(i, j, r);
        }
    }
}

void MatrixCPU::projectPos()
{
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            float r = get(i, j);
            r = (r > 0) * r;
            set(i, j, r);
        }
    }
}




float MatrixCPU::distance2(MatrixCPU* m1) const
{
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension (distance2)");
    }
    float d = 0;
    float r = 0;
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

float MatrixCPU::distance2() const
{
    float d = 0;
    float r = 0;
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

float MatrixCPU::min2() const
{
    if (_row == 0 || _column == 0) {
        return 0;
    }
    float M = fabs(get(0, 0));
    float m = 0;
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            m = fabs(get(i, j));
            if (m < M ) {
                M = m;
            }
        }
    }
    return M;
}

float MatrixCPU::min2Nnull(float eps) const
{
    if (_row == 0 || _column == 0) {
        throw std::invalid_argument("Empty Matrix");
    }
    float M = fabs(get(0, 0));
    float m = 0;
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            m = fabs(get(i, j));
            if ((m < M && m > eps) || (m>M && M <= eps)) {
                M = m;
            }
        }
    }
    return M;
}

float MatrixCPU::max2() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    float M = fabs(get(0, 0));
    float m = 0;
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

float MatrixCPU::max2(MatrixCPU* m1)
{
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension (max2)");
    }

    if (_row == 0 || _column == 0) {
        return 0;
    }
    float M = fabs(get(0, 0) - m1->get(0,0));
    float m = 0;
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < _column; ++j)
        {
            m = fabs(get(i, j) - m1->get(i, j));
            if (m > M) {
                M = m;
            }
        }
    }
    return M;
}

float MatrixCPU::maxAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* indices)
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix (maxAbs)");
    }
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds (maxAbs)");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd (maxAbs)");
    }
    float M = fabs(get(iBegin, jBegin));
    if (indices != nullptr) {
        indices->set(0, 0, iBegin);
        indices->set(0, 1, jBegin);
    }
    float m = 0;
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

float MatrixCPU::minAbs(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* indices, bool Null)
{
    // Null == False : retourne le minimum non null ! si tout est nul retourne 0
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix (minAbs)");
    }
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds (minAbs)");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd (minAbs)");
    }
    float M = INFINITY;
    float m = 0;
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

void MatrixCPU::Moy(MatrixCPU* m, MatrixCPU* nb, int sens)
{   
    float s;
    int n;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol() ) || (_column != nb->getNCol()) || (nb->getNLin() !=1) ) 
        {
            throw std::invalid_argument("wrong dimension of the vector (moy)");
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
            throw std::invalid_argument("wrong dimension of the vector (moy)");
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

float MatrixCPU::sum() const
{
    float d = 0;
    float r = 0;
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

float MatrixCPU::sumD() const
{
    double d = 0;
    float r = 0;
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

void MatrixCPU::sum(MatrixCPU* m, int sens)
{
    float s = 0;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol()) )
        {
            throw std::invalid_argument("wrong dimension of the row vector, sens==1 (sum)");
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
            throw std::invalid_argument("wrong dimension of the column vector (sens==0) (sum)");
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
void MatrixCPU::sumT(MatrixCPU* m, int sens)
{
    float s = 0;
    if (sens) { // on travaille sur les colonnes
        if ((_column != 1) || (_row != m->getNCol()))
        {
            throw std::invalid_argument("wrong dimension of the line vector, sens==1 (sumT)");
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
            throw std::invalid_argument("wrong dimension of the column vector (sens==0) (sumT)");
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



void MatrixCPU::sort(int dim, int sens)
{
    if (sens) { // on travaille sur les colonnes
        if (dim >= _row) {
            throw std::invalid_argument("dim must be smaller then the matrix dimension (sort)");
        }

        int milieu = _column / 2;
        sortColumn(0, milieu, dim);
        sortColumn(milieu, _column, dim);
        fusionColumn(0, _column, dim);
    }
    else {
        if (dim >= _column) {
            throw std::invalid_argument("dim must be smaller then the matrix dimension (sort)");
        }
        int milieu = _row / 2;
        sortLine(0, milieu, dim);
        sortLine(milieu, _row, dim);
        fusionLine(0, _row, dim);
    }

}

void MatrixCPU::sortLine(int line1, int line2, int dim)
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

void MatrixCPU::sortColumn(int col1, int col2, int dim)
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

void MatrixCPU::fusionLine(int line1, int line2, int dim)
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

void MatrixCPU::fusionColumn(int col1, int col2, int dim)
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

void MatrixCPU::swapLine(int line1, int line2)
{
    
    float temp = 0;
    for (int i = 0; i < _column; i++) {
        temp = get(line1, i);
        set(line1, i, get(line2, i));
        set(line2, i, temp);
    }
}

void MatrixCPU::swapColumn(int col1, int col2)
{
    float temp = 0;
    for (int i = 0; i < _row; i++) {
        temp = get(i, col1);
        set(i, col1, get(i, col2));
        set(i, col2, temp);
    }
}

void MatrixCPU::RelativeEror(MatrixCPU* MatRef, MatrixCPU* Mat)
{
    if (!dim(MatRef) || !dim(Mat)) {
        std::cout << _row << " " << _column << " " << MatRef->_row << " " << MatRef->_column 
            << " " << Mat->_row << " " << Mat->_column << " " << std::endl;
        throw std::invalid_argument("matrix must have the same size (RelativeEror) ");
    }
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            float a = MatRef->get(i, j);
            float b = Mat->get(i, j);
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
// Conversion
///////////////////////////////////////////////////////////////////////////////


c_float* MatrixCPU::toCFloat()
{
    c_float* m = new c_float[_row * _column];

    for (int i = 0; i < _row * _column; i++) {
        m[i] = _matrixCPU[i];
    }

    return m;
}



c_int MatrixCPU::toCSC(c_float* data, c_int* idx, c_int* ptr)
{
    int N = getNNull();
    int indiceData = 0;
 
    int indiceCol = 0;
    // cas particulier où la/les première colonnes sont nulles ????
    
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

c_int MatrixCPU::toCSCHalf(c_float* data, c_int* idx, c_int* ptr)
{
    if (_column != _row) {
        throw std::invalid_argument("the matrix must be square (toCSCHalf)");
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

void MatrixCPU::toEigenMatrix(Eigen::MatrixXd* eigenMatrix)
{
    for (int i = 0; i < _row; i++) {
        for (int j = 0; j < _column; j++) {
            (*eigenMatrix)(i, j) = get(i, j);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Display MatrixCPU contents
///////////////////////////////////////////////////////////////////////////////
void MatrixCPU::display() const
{   
    if (this) {
        if (_row == 0 || _column == 0)
        {
            std::cout << "matrix vide " << std::endl;
        }
        else if (_column == 1) {
            std::cout << " transpose  : ";
            for (int i = 0;i < _row;++i)
            {
                for (int j = 0;j < _column;++j)
                {
                    float value = get(i, j);
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
                    float value = get(i, j);
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

void MatrixCPU::displayBloc(int iBegin, int iEnd, int jBegin, int jEnd) const
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds (displayBloc)");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd (displayBloc)");
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
                
                float value = get(i, jBegin);
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
                    float value = get(i, j);
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



void MatrixCPU::saveCSV(const std::string& filename, std::ios_base::openmode mode, int trans, std::string separator) const
{
    //std::cout << "save " << filename << std::endl;
    std::ofstream myfile;
    myfile.open(filename, mode);
    myfile.precision(10);
    if (!trans) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column;j++) {
                myfile << get(i, j) << separator;
            }
            myfile << "\n";
        }
    }
    else {
        for (int j = 0; j < _column;j++) {
            for (int i = 0; i < _row; i++) {
                myfile << get(i, j) << separator;
            }
            myfile << "\n";
        }
    }
    
    myfile.close();
    //std::cout << "fin save " << std::endl;
}

void MatrixCPU::saveTabMatCSV(MatrixCPU* tab, int sizeTab, const std::string& filename, std::ios_base::openmode mode, int trans, std::string separator)
{
    //std::cout << "save " << filename << std::endl;
    std::ofstream myfile;
    myfile.open(filename, mode);
    myfile.precision(10);

    int rowLoc = tab[0].getNLin();
    int columnLoc = tab[0].getNCol();


    for (int n = 0; n < sizeTab; n++)
    {
        if ((tab[n].getNLin() != rowLoc) && (trans == 0)) {
            throw std::invalid_argument("saveTabMatCSV : tab must be a MatrixCPU aay of the same number of row ");
        }
        if ((tab[n].getNCol() != columnLoc) && (trans == 1)) {
            throw std::invalid_argument("saveTabMatCSV : tab must be a MatrixCPU aay of the same number of column, trans  ");
        }
    }

    if (!trans) {
        
        for (int i = 0; i < rowLoc; i++) {
            for (int n = 0; n < sizeTab; n++)
            {
                int columnn = tab[n].getNCol();
                for (int j = 0; j < columnn; j++) {
                    myfile << tab[n].get(i, j) << separator;
                }
            myfile << "\n";
            }
        }
    }
        
    else {
        for (int j = 0; j < columnLoc; j++) {
            for (int n = 0; n < sizeTab; n++)
            {
                int row = tab[n].getNLin();
                for (int i = 0; i < row; i++) {
                    myfile << tab[n].get(i, j) << separator;
                }
               
            }
            myfile << "\n";
        }
    }

    myfile.close();
    //std::cout << "fin save " << std::endl;


}



///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
MatrixCPU::~MatrixCPU()
{
    #ifdef DEBUG_DESTRUCTOR
        std::cout << "destruction matrix " << _matrixCPU << std::endl;
    #endif // DEBUG_DESTRUCTOR

    
    /*if (_matrixCPU != nullptr) {
        delete[] _matrixCPU;
        _matrixCPU = nullptr;
    }*/
   
    //std::cout << _matrixCPU << std::endl;
    DELETEA(_matrixCPU);
    
}


