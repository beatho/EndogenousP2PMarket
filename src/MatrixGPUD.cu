#include "../head/MatrixGPUD.cuh" 

static const int warpSizeD = 32;



double MatrixGPUD::rand1()
{
    double a = (double)(rand()) / ((double)(RAND_MAX));
    return a;
}


///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MatrixGPUD::MatrixGPUD() {
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur appele" << std::endl;
#endif
    _row = 0;
    _column = 0;
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
}

MatrixGPUD::MatrixGPUD(int l, int c, double value, bool pos)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur parametre appele" << std::endl;
    std::cout << _matrixCPU << std::endl;
#endif
    _row = l;
    _column = c;
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
    if (pos) {
        cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, value, _N);
        _GPU = true;
    }
    else {
        _matrixCPU = new double[l * c];
        for (int elem = 0; elem < l * c; elem++) {
            _matrixCPU[elem] = value;
        }
    }
#ifdef DEBUG_CONSTRUCTOR
    std::cout << _matrixGPU << std::endl;
#endif
}

MatrixGPUD::MatrixGPUD(const MatrixCPUD& m, bool pos)
{
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);

    if (pos) {
        _GPU = true;
        cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
    }
    else {
        _matrixCPU = new double[_row * _column];
        memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
    }
    
}

MatrixGPUD::MatrixGPUD(const MatrixCPU& m, bool pos)
{
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
    _matrixCPU = new double[_row * _column];
    for (int elem = 0; elem < _row * _column; elem++) {
        _matrixCPU[elem] = (double)m._matrixCPU[elem];
    }

    if (pos) {
        _GPU = true;
        cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        cudaMemcpy(_matrixGPU, _matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
    }
    
}

MatrixGPUD::MatrixGPUD(const MatrixGPUD & m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur recopie appele" << std::endl;
#endif
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);

    if (m.getPos()) {
        cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        setGPU <<<_numBlocks, _blockSize >>> (_matrixGPU, m._matrixGPU, _N);
        _GPU = true;
    }
    else {
        _matrixCPU = new double[_row * _column];
        memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
    }
}

MatrixGPUD& MatrixGPUD::operator=(const MatrixGPUD& m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur operateur = appele" << std::endl;
#endif
    if (_row == m.getNLin() && _column == m.getNCol()) {
        //matrix already has the good size no free needed
        if (getPos()) {
            if (m.getPos()) {
                setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
            }
            else {
                cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
            }
        }
        else {
            if (m.getPos()) {
                cudaMemcpy(_matrixCPU, m._matrixGPU, sizeof(double) * _row * _column, cudaMemcpyDeviceToHost);
            }
            else {
                memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
            }
        }
    }
    else {
        _row = m.getNLin();
        _column = m.getNCol();
        _N = _row * _column;
        _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
        if (getPos()) {
            _GPU = false;
            cudaFree(_matrixGPU);
        }
        else
        {
            DELETEA(_matrixCPU);
        }
        if (m.getPos()) {
            cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
            setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
            _GPU = true;
        }
        else {
            _matrixCPU = new double[_row * _column];
            memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
        }
    }
   
    return *this;
}

MatrixGPUD& MatrixGPUD::operator=(const MatrixCPUD& m)
{
    if (_row == m.getNLin() && _column == m.getNCol()) {
        //matrix already has the good size no free needed
        if (getPos()) {
            cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
        }
        else {
            memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
        }
    }
    else {
        _row = m.getNLin();
        _column = m.getNCol();
        _N = _row * _column;
        _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
        if (getPos()) {
            cudaFree(_matrixGPU);
            cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
            cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
            _GPU = true;
        }
        else
        {
            DELETEA(_matrixCPU);
            _matrixCPU = new double[_row * _column];
            memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(double));
        } 
    }
    return *this;
}

void MatrixGPUD::preallocateReduction()
{
    if (preallocation) {
        cudaFreeHost(_preallocationFloat);
        cudaFree(_preallocation);
        preallocation = false;
    }
    cudaError_t c;
    do
    {
        c = cudaHostAlloc(&_preallocationFloat, sizeof(double), cudaHostAllocDefault);
     
        /*if (_preallocationFloat == nullptr) {
            
            std::cout << "prealocation echouer ? " << c << std::endl;
            std::cout <<  cudaGetErrorName(c) << std::endl;
        }*/
    } while (c == 700);

    
    c = cudaMalloc((void**)&_preallocation, sizeof(double) * _numBlocks);
    /*if (_preallocation == nullptr) {
        std::cout << _row << " " << _column << " " << _blockSize << std::endl;
        std::cout << "prealocation echouer ? " << c << " " << _numBlocks  <<std::endl;
        std::cout << cudaGetErrorName(c) << std::endl;
    }*/
    if (c == 700) {
        exit(-1);
    }
    
    preallocation = true;
    setGPU <<<_numBlocks, _blockSize >>> (_preallocation, 0.0f, _numBlocks);
}

void MatrixGPUD::transferGPU()
{
    if (!_GPU) {
        if (!_matrixGPU) {
            cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        }
        cudaMemcpy(_matrixGPU, _matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
        //DELETEA(_matrixCPU);
        _GPU = true;
    }
    else {
        throw std::domain_error("already in the GPU");
    }
    
}

void MatrixGPUD::transferCPU()
{
    
    if (_GPU) {
        
        if (!_matrixCPU) {
            
            _matrixCPU = new double[_row * _column];
        }
        cudaMemcpy(_matrixCPU, _matrixGPU, sizeof(double) * _row * _column, cudaMemcpyDeviceToHost);
        //cudaFree(_matrixGPU);
        //_matrixGPU = nullptr;
        _GPU = false;
    }
    else {
        std::cout << "already in the CPU" << _GPU <<std::endl;
        throw std::domain_error("already in the CPU");
    }

}

///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////
 double MatrixGPUD::get(int i, int j, bool verbose) const
{
    //std::cout << "hey de taille " << _row << " " << _column << "pos "<< i <<" "<< j << std::endl;
    if ((i >= _row) || ( j >= _column) || (i < 0) || ( j < 0)) {
        throw std::out_of_range("index out of bounds");
    }
    if (_GPU) {
        double value;
        cudaMemcpy(&value, _matrixGPU + i*_column+j, sizeof(double), cudaMemcpyDeviceToHost);
        if (verbose) {
            std::cout << " Warning matrix on GPU" << std::endl;
        }
        return value;
        //throw std::invalid_argument("Matrix on GPU");
    }
    else {
        return _matrixCPU[i * _column + j];
    }
}

int MatrixGPUD::getNCol() const
{
    return _column;
}

int MatrixGPUD::getNLin() const
{
    return _row;
}

void MatrixGPUD::getCol(MatrixGPUD* col, int numCol, int offset)
{
    if (numCol < 0 || numCol >= _column) {
        throw std::out_of_range("index out of bounds");
    }
    if (offset < 0 || offset >= _row) {
        throw std::out_of_range("index out of bounds");
    }
    if (col->getNLin() != _row) {
        throw std::invalid_argument("not the same dimension");
    }
    if (col->getNCol() != 1) {
        throw std::invalid_argument("must be a column vector");
    }

    if (!_GPU && !col->getPos()) {
        for (int i = 0; i < offset; i++) {
            col->set(i, 0, 0);
        }
        for (int i = offset; i < _row; i++) {
            col->set(i, 0, get(i, numCol));
        }
    }
    else if (_GPU && col->getPos()) {
        setColGPU <<< _numBlocks, _blockSize >>> (col->_matrixGPU, _matrixGPU, numCol, _column, _row, offset);
    }

}

bool  MatrixGPUD::getPos() const
{
    return _GPU;
}
bool MatrixGPUD::dim(MatrixGPUD* m) const
{ 
    return ((_row == m->getNLin()) && (_column == m->getNCol()));
}


bool MatrixGPUD::isEqual(MatrixGPUD* m, double pre) const
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    else {
        if (_GPU || m->getPos()) {
            throw std::invalid_argument("Matrix on GPU");
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
    }
    return true;
}

void MatrixGPUD::toMatCPU(MatrixCPUD& m) const // passer m en paramètre
{
    if (m.getNCol() != _column || m.getNLin() != _row) {
        m.setSize(_row, _column);
    }
    if (_GPU) {
        cudaMemcpy(m._matrixCPU, _matrixGPU, sizeof(double) * _row * _column, cudaMemcpyDeviceToHost);
    }
    else {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++) 
            {
                m.set(i, j, get(i, j));
            }
        }
    }
}



void MatrixGPUD::setSize(int row, int column)
{
    _row = row;
    _column = column;
    _N = _row * _column;
    _numBlocks = ceil((_N + _blockSize - 1) / _blockSize);
    if (getPos()) {
        cudaFree(_matrixGPU);
        cudaMalloc((void**)&_matrixGPU, sizeof(double) * _row * _column);
        _GPU = true;
    }
    else
    {
        DELETEA(_matrixCPU);
        _matrixCPU = new double[_row * _column];
        memset(_matrixCPU, 0, _N * sizeof(double));
    }
}

///////////////////////////////////////////////////////////////////////////////
// Setter
///////////////////////////////////////////////////////////////////////////////
 void MatrixGPUD::set(int i, int j, double value, bool force)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        throw std::out_of_range("index out of bounds");
    }
     if (_GPU && !force) {
        throw std::invalid_argument("Matrix on GPU");
     }
     else if (_GPU && force) {
         setGPUunique <<< 1, 1 >>> (_matrixGPU, value, i * _column + j);
     }
     else {
         //std::cout << "changement de valeur " << value << " en " << i << " " << j << std::endl;
         _matrixCPU[i * _column + j] = value;
     }
}

 void MatrixGPUD::setEyes(double value)
 {
     if (!_GPU) {
         int N = _row * (_row < _column) + _column * (_column <= _row);

         for (int i = 0; i < _row; i++) {
             for (int j = 0; j < _column; j++)
             {
                 if (i == j) {
                     set(i, j, value);
                 }
                 else {
                     set(i, j, 0);
                 }  
             }
         }
     }
     else {
         setEyesGPU<<< _numBlocks, _blockSize >>>(_matrixGPU, value, _column, _row);
     }
     
 }

 void MatrixGPUD::setEyes(MatrixGPUD* m)
 {
     if (m->getNLin() != _row || _row != _column || m->getNCol() != 1) 
     {
         throw std::invalid_argument("not the good dimension");
     }
     if (_GPU && m->getPos()) {
         setEyesGPU <<< _numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, _column, _row);

     }
     else if (!_GPU && !(m->getPos()))
     {
         for (int i = 0; i < _row; ++i)
         {
            this->set(i, i, m->get(i, 0));
         }
     }
     else {
         throw std::invalid_argument("Matrix not at the same place");
     }



 }

void MatrixGPUD::set(MatrixGPUD* m, bool synchrone, cudaStream_t stream)
{
   
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos()) {
        if (synchrone) {
            setGPU <<<_numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, _N);
        }
        else {
            setGPU <<< _numBlocks, _blockSize, 0, stream>>> (_matrixGPU, m->_matrixGPU, _N);
        }
        
    }
    else if (!_GPU && !(m->getPos())) 
    {
        for (int i = 0; i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            { 
                this->set(i, j, m->get(i, j));
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }

}

void MatrixGPUD::set(MatrixCPUD* m)
{
    if (m->getNCol() != _column || m->getNLin() != _row) {
        throw std::invalid_argument("not the same dimension");
    }
    
    if (getPos()) {
        cudaMemcpy(_matrixGPU, m->_matrixCPU, sizeof(double) * _row * _column, cudaMemcpyHostToDevice);
    }
    else {
        memcpy(_matrixCPU, m->_matrixCPU, _row * _column * sizeof(double));
    }
   
    
}

void MatrixGPUD::setTrans(MatrixGPUD* m)
{
    if (_column != m->getNLin() || _row != m->getNCol()) {
        throw std::invalid_argument("not the same transposed dimension");
    }
    if (getPos() && m->getPos()) {
        setTransGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, _column, _row);
    }
    else if (!getPos() && !(m->getPos()))
    {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++) {
                set(i, j, m->get(j, i));
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }

}


void MatrixGPUD::setRand(double eps)
{
    //exit(1);
    if (_GPU) {

        /*curandGenerator_t gen;
        curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
        curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);
        curandGenerateUniform(gen, _matrixGPU, _N);
        curandDestroyGenerator(gen);*/
        
        curandState* state = nullptr;
        cudaMalloc((void**)&state, _N * sizeof(curandState));
        setup_kernelD <<<_numBlocks, _blockSize >>> (state);
        generate_kernel << <_numBlocks, _blockSize >> > (state, _matrixGPU, eps, _N);
        //throw std::invalid_argument("Matrix on GPU");
    }
    else {
        int N = _column * _row;
        for (int elem = 0; elem < N; elem++) {
            _matrixCPU[elem] = 2 * (rand1() - 0.5) * eps;

        }
    }
    
}

void MatrixGPUD::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPUD* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!_GPU && !(m->getPos())) {
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
    else if (getPos() && (m->getPos())) {
        /*const int nThread = 16;
         const int bx = (jEnd - jBegin + nThread - 1) / nThread;
         const int by = (iEnd - iBegin + nThread - 1) / nThread;
         dim3 gridBlock(bx, by);
         dim3 dimBlock(nThread, nThread);*/
        SetBlocGPU << <m->_numBlocks, m->_blockSize >> > (_matrixGPU, m->_matrixGPU, iBegin, iEnd, jBegin, jEnd, _column);
    }
    else {
       
        throw std::invalid_argument("Matrix not at the same place");
    }
}



void MatrixGPUD::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPUD* m, double factor)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!_GPU && !(m->getPos())) {
        int row = 0;

        for (int i = iBegin; i < iEnd; i++) {
            int col = 0;
            for (int j = jBegin; j < jEnd; j++) {
                set(i, j, factor * m->get(row, col));
                col = col + 1;
            }
            row = row + 1;
        }
    }
    else if (getPos() && (m->getPos())) {
        /*const int nThread = 16;
        const int bx = (jEnd - jBegin + nThread - 1) / nThread;
        const int by = (iEnd - iBegin + nThread - 1) / nThread;
        dim3 gridBlock(bx, by);
        dim3 dimBlock(nThread, nThread);*/
        SetBlocGPU <<<m->_numBlocks, m->_blockSize >> > (_matrixGPU, m->_matrixGPU, iBegin, iEnd, jBegin, jEnd, _column, factor);
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}
void MatrixGPUD::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPUD* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!_GPU) {
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
    else {
        throw std::domain_error("Matrix on GPU");
    }
}


void MatrixGPUD::swap(MatrixGPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos()) {
        double* temp = _matrixGPU;
        _matrixGPU = m->_matrixGPU;
        m->_matrixGPU = temp;

    }
    else if (!_GPU && !(m->getPos())) {
        double* temp = _matrixCPU;
        _matrixCPU = m->_matrixCPU;
        m->_matrixCPU = temp;
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }

    
    
}

void MatrixGPUD::replace(double previous, double newValue)
{
    if (_GPU) {
        replaceGPU <<<_numBlocks, _blockSize >> > (_matrixGPU, previous, newValue, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                if (get(i, j) == previous) {
                    this->set(i, j, newValue);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// Addition
///////////////////////////////////////////////////////////////////////////////
void MatrixGPUD::add(MatrixGPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos()) 
    {
        addGPU<<<_numBlocks, _blockSize>>>(_matrixGPU,m->_matrixGPU,_N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) + m->get(i, j);
                this->set(i, j, r);
            }
        }
    } else {
        throw std::invalid_argument("Matrix not at the same place");
    } 
}

void MatrixGPUD::addVector(MatrixGPUD* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector");
    }
    if (v->getNCol() == 1) {
        if (_GPU && v->getPos()) 
        {
            addVectorGPU1<<<_numBlocks, _blockSize >>>(_matrixGPU, v->_matrixGPU, _column, _N);
        }
        else if ((!_GPU) && !(v->getPos())) {
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
            throw std::invalid_argument("Matrix not at the same place");
        } 
    }
    else {
        if (_GPU && v->getPos())
        {
            addVectorGPU2<<<_numBlocks, _blockSize >>>(_matrixGPU, v->_matrixGPU, _column, _N);
        }
        else if ((!_GPU) && !(v->getPos())) {
            for (int i = 0;i < _row;++i)
            {
                for (int j = 0;j < _column;++j)
                {
                    double r = get(i, j) + v->get(0, j);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("Matrix not at the same place");
        }
    }
}
void MatrixGPUD::add(double c)
{
    if (_GPU) {
        addGPU<<<_numBlocks,_blockSize >>>(_matrixGPU,c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) + c;
                this->set(i, j, r);
            }
        }
    }
}

void MatrixGPUD::add(MatrixGPUD* m1, MatrixGPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m1->getPos() && m2->getPos()) 
    {
        addGPU<<<_numBlocks, _blockSize >> > (_matrixGPU, m1->_matrixGPU,m2->_matrixGPU, _N);
    }
    else if (!_GPU && !(m1->getPos()) && !(m2->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = m1->get(i, j) + m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    
}
void MatrixGPUD::add(MatrixGPUD* m, double c)
{
    if (_GPU && m->getPos()) 
    {
        addGPU<<<_numBlocks, _blockSize >>> (_matrixGPU,m->_matrixGPU, c, _N);
    }
    else if ((!_GPU) && !(m->getPos())) 
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
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    

}
void MatrixGPUD::addTrans(MatrixGPUD* m)
{
    MatrixGPUD temp(*this);
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("not the same dimension (transpose)");
    }
    if (_GPU && m->getPos())
    {
        addTransGPU<<<_numBlocks, _blockSize >>>(temp._matrixGPU, _matrixGPU, m->_matrixGPU,_column,_row,_N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) + m->get(j, i);
                temp.set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    this->set(&temp);
    
}
///////////////////////////////////////////////////////////////////////////////
// subtraction
///////////////////////////////////////////////////////////////////////////////
void MatrixGPUD::subtract(MatrixGPUD* m1, MatrixGPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m1->getPos() && m2->getPos())
    {
        substractGPU<<<_numBlocks, _blockSize >> > (_matrixGPU, m1->_matrixGPU, m2->_matrixGPU, _N);
    }
    else if (!_GPU && !(m1->getPos()) && !(m2->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = m1->get(i, j) - m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    
}
void MatrixGPUD::subtract(MatrixGPUD* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos())
    {
        substractGPU <<<_numBlocks, _blockSize >> > ( _matrixGPU, m->_matrixGPU, _N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) - m->get(i, j);
                set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
   
    
}
void MatrixGPUD::subtractVector(MatrixGPUD* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("wrong dimension of the vector");
    }
    if (v->getNCol() == 1) {
        if (_GPU && v->getPos())
        {
            substractVectorGPU1 <<<_numBlocks, _blockSize >>>(_matrixGPU, v->_matrixGPU, _column, _N);
        }
        else if ((!_GPU) && !(v->getPos())) {
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
            throw std::invalid_argument("Matrix not at the same place");
        }
    }
    else {
        if (_GPU && v->getPos())
        {
            substractVectorGPU2 <<<_numBlocks, _blockSize >>>(_matrixGPU, v->_matrixGPU, _column, _N);
        }
        else if ((!_GPU) && !(v->getPos())) {
            for (int i = 0;i < _row;++i)
            {
                for (int j = 0;j < _column;++j)
                {
                    double r = get(i, j) - v->get(0, j);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("Matrix not at the same place");
        }
    }

}
void MatrixGPUD::subtractTrans(MatrixGPUD* m)
{
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("not the same dimension (transpose)");
    }
    MatrixGPUD temp(*this);
    if (_GPU && m->getPos())
    {
        substractTransGPU <<<_numBlocks, _blockSize >>>(temp._matrixGPU, _matrixGPU, m->_matrixGPU, _column, _row, _N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) - m->get(j, i);
                temp.set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    this->set(&temp);
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication
///////////////////////////////////////////////////////////////////////////////


void MatrixGPUD::multiply(double c)
{
    if (_GPU) {
        multiplyGPU<<<_numBlocks, _blockSize >>> (_matrixGPU, c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) * c;
                this->set(i, j, r);
            }
        }
    }
        
}

void MatrixGPUD::multiplyMat(MatrixGPUD* A, MatrixGPUD* B)
{
    if (A->getNLin() != getNLin()) {
        throw std::invalid_argument("result must be compatible with A (row)");
    }
    else if (A->getNCol() != B->getNLin()) {
        throw std::invalid_argument("A must be compatible with B (column with row)");
    }
    else if (getNCol() != B->getNCol()) {
        throw std::invalid_argument("result must be compatible with B (column)");

    }
    if (_GPU && A->getPos() && B->getPos()) { // solution temporaire
        transferCPU();
        A->transferCPU();
        B->transferCPU();
        double r = 0;
        int p = A->getNCol();
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _column; ++j)
            {
                r = 0;
                for (int k = 0; k < p; ++k)
                {
                    r += A->get(i, k) * B->get(k, j);
                }
                this->set(i, j, r);
            }
        }
        transferGPU();
        A->transferGPU();
        B->transferGPU();
    }
    else if (!_GPU && !(A->getPos()) && !B->getPos()) {
        double r = 0;
        int p = A->getNCol();
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _column; ++j)
            {
                r = 0;
                for (int k = 0; k < p; ++k)
                {
                    r += A->get(i, k) * B->get(k, j);
                }
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

void MatrixGPUD::multiply(MatrixGPUD* Mat, MatrixGPUD* vect, bool trans)
{
    // result = Mat*vect  nLine*nCol = nLine*Taille *Taille * Ncol 
    if (trans) {
        if (getNLin() != 1) {
            throw std::invalid_argument("result must be a row vector ");
        }
        else if (getNCol() != Mat->getNLin()) {
            throw std::invalid_argument("result must be compatible with Mat");
        }
        else if (vect->getNLin() != 1) {
            throw std::invalid_argument("vect must be a row vector ");
        }
        else if (vect->getNCol() != Mat->getNCol()) {
            throw std::invalid_argument("vect must be compatible with Mat");
        }
    }
    else {
        if (getNCol() != 1) {
            throw std::invalid_argument("result must be a column vector ");
        }
        else if (getNLin() != Mat->getNLin()) {
            throw std::invalid_argument("result must have the same row number as the Mat");
        }
        else if (vect->getNCol() != 1) {
            throw std::invalid_argument("vect must be a column vector ");
        }
        else if (vect->getNLin() != Mat->getNCol()) {
            throw std::invalid_argument("vect must be compatible with Mat");
        }
    }
   
    if (_GPU && Mat->getPos() && vect->getPos())
    {
        int numBlock = Mat->getNLin();
        switch (_blockSize) {
        case 512:
            multiplyGPU<512> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case 256:
            multiplyGPU<256> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case 128:
            multiplyGPU<128> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case 64:
            multiplyGPU< 64> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case 32:
            multiplyGPU< 32> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case 16:
            multiplyGPU< 16> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case  8:
            multiplyGPU<  8> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case  4:
            multiplyGPU<  4> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case  2:
            multiplyGPU<  2> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        case  1:
            multiplyGPU<  1> << <numBlock, _blockSize >> > (_matrixGPU, Mat->_matrixGPU, vect->_matrixGPU, Mat->getNCol());
            break;
        }
        
    }
    else if (!_GPU && !(Mat->getPos()) && !vect->getPos())
    {
        if (trans) {
            for (int i = 0; i < Mat->getNLin(); ++i)
            {
                double sum = 0;
                for (int j = 0; j < Mat->getNCol(); ++j)
                {
                    sum += Mat->get(i, j) * vect->get(0, j);
                }
                set(0, i, sum);
            }
        }
        else {
            for (int i = 0; i < _row; ++i)
            {
                double sum = 0;
                for (int j = 0; j < Mat->getNCol(); ++j)
                {
                    sum += Mat->get(i, j) * vect->get(j, 0);
                }
                set(i, 0, sum);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

void MatrixGPUD::linearOperation(MatrixGPUD* A, MatrixGPUD* x, MatrixGPUD* b, bool trans)
{
    if (trans) {
        if (getNLin() != 1) {
            throw std::invalid_argument("result must be a row vector ");
        }
        else if (getNCol() != A->getNLin()) {
            throw std::invalid_argument("result must be compatible with A");
        }
        else if (x->getNLin() != 1 || b->getNLin() != 1) {
            throw std::invalid_argument("x and b must be a row vector ");
        }
        else if (x->getNCol() != A->getNCol()) {
            throw std::invalid_argument("x must be compatible with A");
        }
    }
    else {
        if (getNCol() != 1) {
            throw std::invalid_argument("result must be a column vector ");
        }
        else if (getNLin() != A->getNLin()) {
            throw std::invalid_argument("result must have the same row number as A");
        }
        else if (x->getNCol() != 1 || b->getNCol() != 1) {
            throw std::invalid_argument("x and b must be a column vector ");
        }
        else if (x->getNLin() != A->getNCol()) {
            throw std::invalid_argument("x must be compatible with Mat");
        }
    }

    if (_GPU && A->getPos() && b->getPos() && x->getPos())
    {
        int numBlock = A->getNLin();
        switch (_blockSize) {
        case 512:
            linearOpGPU<512> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case 256:
            linearOpGPU<256> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case 128:
            linearOpGPU<128> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case 64:
            linearOpGPU< 64> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case 32:
            linearOpGPU< 32> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case 16:
            linearOpGPU< 16> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case  8:
            linearOpGPU<  8> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case  4:
            linearOpGPU<  4> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case  2:
            linearOpGPU<  2> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        case  1:
            linearOpGPU<  1> << <numBlock, _blockSize >> > (_matrixGPU, A->_matrixGPU, x->_matrixGPU, b->_matrixGPU, A->getNCol());
            break;
        }

    }
    else if (!_GPU && !(A->getPos()) && !x->getPos() && !b->getPos())
    {
        if (trans) {
            for (int i = 0; i < A->getNLin(); ++i)
            {
                double sum = 0;
                for (int j = 0; j < A->getNCol(); ++j)
                {
                    sum += A->get(i, j) * x->get(0, j);
                }
                set(0, i, sum + b->get(0,i));
            }
        }
        else {
            for (int i = 0; i < _row; ++i)
            {
                double sum = 0;
                for (int j = 0; j < A->getNCol(); ++j)
                {
                    sum += A->get(i, j) * x->get(j, 0);
                }
                set(i, 0, sum + b->get(i,0));
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication Terme � Terme
///////////////////////////////////////////////////////////////////////////////

void MatrixGPUD::multiplyT(MatrixGPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos())
    {
        multiplyTGPU <<<_numBlocks, _blockSize >>> (_matrixGPU, m->_matrixGPU, _N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) * m->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

void MatrixGPUD::multiplyT(MatrixGPUD* m1, MatrixGPUD* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (!dim(m1)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m1->getPos() && m2->getPos())
    {
        multiplyTGPU<<<_numBlocks, _blockSize >>>(_matrixGPU, m1->_matrixGPU, m2->_matrixGPU, _N);
    }
    else if (!_GPU && !(m1->getPos()) && !(m2->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = m1->get(i, j) * m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}


void MatrixGPUD::divide(double c)
{
    if (c == 0) {
        throw std::domain_error("divide by 0");
    }
    if (_GPU) {
        divideGPU <<<_numBlocks, _blockSize >>> (_matrixGPU, c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                double r = get(i, j) / c;
                this->set(i, j, r);
            }
        }
    }
    
}

void MatrixGPUD::divideT(MatrixGPUD* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && m->getPos())
    {
        divideGPU<<<_numBlocks, _blockSize >>>(_matrixGPU, m->_matrixGPU, _N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        MatrixGPUD temp(*this);
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
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    
}

void MatrixGPUD::invertGaussJordan(MatrixGPUD* mToInvert)
{
    
    if (!dim(mToInvert)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_row != _column) {
        throw std::invalid_argument("must be a square matrix");
    }
   
    if (!_GPU && !mToInvert->getPos()) {
        MatrixCPUD m;
        mToInvert->toMatCPU(m);
       
        MatrixCPUD augmented(_row, _column);
        augmented.setEyes(1);
        MatrixCPUD indices(1, 2);
        int r = 0;
        for (int column = 0; column < _column; column++) {
            
            double pivotabs = m.maxAbs(r, _row, column, column + 1, &indices);
            int k = indices.get(0, 0); // indice max de la colonne j
            double pivot = m.get(k, column);
            if (pivotabs < 0.000001f) {
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
                        m.subtractRow(i, r, local);
                        augmented.subtractRow(i, r, local);
                    }
                }
                r++;
            }
        }
        
        set(&augmented);
       
    }
   
    else if (_GPU && mToInvert->getPos()) {
        MatrixGPUD m(*mToInvert);
        setEyes(1);
        dim3 threadsPerBlock(32, 32);
        int r = 0;
        MatrixGPUD matCol(_row, 1, 0, 1);
        int k = 0;
        
        for (int column = 0; column < _column; column++) {
           
            //transferCPU();
            ///display();
            //transferGPU();

            m.getCol(&matCol, column, r);
            double pivotAbs = matCol.max2(&k); // comme matCol est un vecteur colonne, la position du maximum correspond directement à la ligne !
            double pivot = matCol.get(k, 0, false);
            
            
            if (pivotAbs < 0.000001f) {
                std::cout << "not invertible " << column <<" " << pivotAbs << std::endl;
                if (_N < 100) {
                    mToInvert->display(true);
                }
                throw std::invalid_argument("not invertible matrix");
            }
            else {
                normalisationGJ << <_numBlocks, _blockSize >> > (_matrixGPU, k, _column, pivot); // normalisation pour toute la ligne k
                normalisationGJ <<<_numBlocks, _blockSize >>> (m._matrixGPU, k, _column, pivot);
                
              
                if (k != r) {
                    swapLineGJ << <_numBlocks, _blockSize >> > (_matrixGPU, k, r, _column);// swap des lignes k et r 
                    swapLineGJ << <_numBlocks, _blockSize >> > (m._matrixGPU, k, r, _column);// swap des lignes k et r 
                }
                
                // soustration de la ligne sauf pour la r
                eliminationGJ <<<_row, _blockSize >> > (m._matrixGPU, _matrixGPU, r, _row, _column);
                r++;
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }


}

void MatrixGPUD::LUPFactorization(MatrixGPUD* A, MatrixGPUD* P)
{
    double Tol = 0.0000001;
    int n = getNLin();
    A->set(this);

    // code from wikipedia adapted
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    if (P->getNCol() != 1 || P->getNLin() != (getNCol() + 1)) {
        throw std::invalid_argument("wrong size of P");
    }

    if (!_GPU && !A->getPos() && !P->getPos()) {
        
        for (int i = 0; i < n; i++) {
            P->set(i, 0, i); //Unit permutation matrix, P[N] initialized with N
        }

        double absA = 0;
        int j = 0;
        for (int col = 0; col < n; col++) {
            double maxA = 0.0;
            int imax = col;
            for (int k = col; k < n; k++) {
                absA = fabs(get(k, col));
                if (absA > maxA)
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
                j = P->get(col, 0);
                P->set(col, 0, P->get(imax, 0));
                P->set(imax, 0, j);

                //pivoting rows of A
                A->swapLine(col, imax);

                //counting pivots starting from N (for determinant)
                P->set(n, 0, P->get(n, 0) + 1);
            }


            for (int i = col + 1; i < n; i++) {

                A->set(i, col, A->get(i, col) / A->get(col, col)); //A[j][i] /= A[i][i];


                for (int k = col + 1; k < n; k++) {
                    A->set(i, k, A->get(i, k) - A->get(i, col) * A->get(col, k)); //A[j][k] -= A[j][i] * A[i][k];

                }


            }
        }
    }
    else if (_GPU && A->getPos() && P->getPos()) {

        MatrixGPUD matCol(_row, 1, 0, 1);
        initPermMatr << <_numBlocks, _blockSize >> > (P->_matrixGPU, n);
        int k = 0;
        for (int col = 0; col < n; col++) {

            A->getCol(&matCol, col, col);
            double pivotAbs = matCol.max2(&k); // comme matCol est un vecteur colonne, la position du maximum correspond directement à la ligne !
            //double pivot = matCol.get(k, 0, false);


            if (pivotAbs < Tol) {
                //std::cout << "failure, matrix is degenerate" << std::endl;
                throw std::invalid_argument("matrix is degenerate");
            }
            else {
                if (k != col) { //le max pas sur la diagonal
                    //pivoting P and counting pivots starting from N (for determinant)
                    updatePermMatr << <1, 1 >> > (P->_matrixGPU, k, col, n);

                    //pivoting rows of A
                    A->swapLine(col, k);
                }

                updateLUPFactorization <<<n, _blockSize >>> (A->_matrixGPU, col, n);
            }
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }

    // en vrai on peut tout stocker dans une matrice comme on sait que diag(L) = Id, Et donc on peut avoir A = (L-Id) + U


}

void MatrixGPUD::solveSysUpper(MatrixGPUD* U)
{
    if (getNLin() != U->getNCol() || U->getNLin() != getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    int n = getNLin();
    if (!_GPU && !U->getPos()) {

        for (int i = n - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < n; k++) {
                set(i, 0, get(i, 0) - U->get(i, k) * get(k, 0));// x[i] -= A[i][k] * x[k];
            }

            set(i, 0, get(i, 0) / U->get(i, i));
        }
    }
    else if (_GPU && U->getPos()) {
        for (int i = n - 1; i >= 0; i--)
        {
            solveUpSys << < 1, _blockSize >> > (U->_matrixGPU, _matrixGPU, i, n);
        }

    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

void MatrixGPUD::solveSysLower(MatrixGPUD* L, MatrixGPUD* b, MatrixGPUD* P) // element diag equal to 1
{
    if (getNLin() != L->getNCol() || L->getNLin() != b->getNLin()) {
        throw std::invalid_argument("A must be square");
    }
    int n = getNLin();
    if (!_GPU && !L->getPos() && !P->getPos() && !b->getPos()) {
        for (int i = 0; i < n; i++) {
            set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];

            for (int k = 0; k < i; k++) {
                set(i, 0, get(i, 0) - L->get(i, k) * get(k, 0));
            }
        }
    }
    else if (_GPU && L->getPos() && P->getPos() && b->getPos()) {
        setPermute << <_numBlocks, _blockSize >> > (_matrixGPU, b->_matrixGPU, P->_matrixGPU, n);
        for (int i = 0; i < n; i++) {
            solveLowSys << < 1, _blockSize >> > (L->_matrixGPU, _matrixGPU, i, n);
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }

}

void MatrixGPUD::solveSys(MatrixGPUD* A, MatrixGPUD* P, MatrixGPUD* b)
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

    if (!_GPU && !A->getPos() && !P->getPos() && !b->getPos()) {
       
        int n = getNLin();
       
        for (int i = 0; i < n; i++) {
                set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];
        }
       
        for (int iter = 0; iter < n; iter++) {
            for (int k = iter + 1; k < n; k++) { // en parallele
                set(k, 0, get(k, 0) - A->get(k, iter) * get(iter, 0)); // x[k] = x[k] - A[k][iter] *x[iter] ; avec k>iter
            }
        }
       
        for (int iter = n - 1; iter >= 0; iter--) {
            set(iter, 0, get(iter, 0) / A->get(iter, iter));
            for (int k = 0; k < iter; k++) { // en parallele
                set(k, 0, get(k, 0) - A->get(k, iter) * get(iter, 0)); // x[k] = x[k] - A[k][iter] *x[iter] ; avec k<iter
            }

        }
  
    }
    else if (_GPU && A->getPos() && P->getPos() && b->getPos()) {
        int n = getNLin();
       
        setPermute << <_numBlocks, _blockSize >> > (_matrixGPU, b->_matrixGPU, P->_matrixGPU, n);
       
        
        solveSysGPU << <1, _blockSize, n * sizeof(double) >> > (A->_matrixGPU, _matrixGPU, n);

    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }


}


///////////////////////////////////////////////////////////////////////////////
// Fonction autres
///////////////////////////////////////////////////////////////////////////////

double MatrixGPUD::max2() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    if (!_GPU) {
        double M = fabs(get(0, 0));
        double m = 0;
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                m = fabs(get(i, j));
                if (m > M) {
                    M = m;
                }
            }
        }
        return M;
    } 
    else {    
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double odata = 0;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            std::cout << "allocation !!!" << std::endl;
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }
        //std::cout << _numBlocks << std::endl;

        switch (_blockSize) {
        case 512:
            maxMultiBlock<512> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<512> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 256:
            maxMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 128:
            maxMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 64:
            maxMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            maxMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            maxMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            maxMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            maxMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            maxMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            maxMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            maxMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        }
        //cudaDeviceSynchronize();
         if (preallocation) {
            cudaMemcpy(_preallocationFloat, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
            return sqrt(*_preallocationFloat);
        }
        else
        {
            cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
            std::cout << "free !!!" << std::endl;
            cudaFree(d_odata);
            return sqrt(odata);
        }
    }
}

double MatrixGPUD::max2(int* indice)
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix");
    }
    if (!_GPU) {
        double M = fabs(get(0, 0));
        double m = 0;
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _column; ++j)
            {
                m = fabs(get(i, j));
                if (m > M) {
                    *indice = i * _column + j;
                    M = m;
                }
            }
        }
        return M;
    }
    else {

        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double odata = 0;
        double* d_odata;
        int* d_pos;
        int pos = 0;
        cudaMalloc((void**)&d_pos, sizeof(int) * numBlocks);
        if (preallocation) {

            d_odata = _preallocation;
        }
        else {

            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }
        //std::cout << _numBlocks << std::endl;

        switch (_blockSize) {
        case 512:
            maxMultiBlock<512> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<512> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case 256:
            maxMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case 128:
            maxMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case 64:
            maxMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case 32:
            maxMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case 16:
            maxMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case  8:
            maxMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case  4:
            maxMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case  2:
            maxMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        case  1:
            maxMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n, d_pos);
            maxMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks, d_pos);
            break;
        }
        //cudaDeviceSynchronize();
        cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&pos, d_pos, sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(d_pos);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        *indice = pos;
        return sqrt(odata);
    }
}
double MatrixGPUD::max2(MatrixGPUD* m) const
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("Empty Matrix");
    }
    if (!_GPU && !(m->getPos())) 
    {
        double M = fabs(get(0, 0));
        double f = 0;
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                f = fabs(get(i, j)- m->get(i,j));
                if (f > M) {
                    M = f;
                }
            }
        }
        return M;
    }
    else if (_GPU && m->getPos()) {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double odata;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            std::cout << "allocation !!!" << std::endl;
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }
        //std::cout << _numBlocks << std::endl;

        switch (_blockSize) {
        case 512:
            maxMultiBlock<512> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<512> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 256:
            maxMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 128:
            maxMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 64:
            maxMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            maxMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            maxMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            maxMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            maxMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            maxMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            maxMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            maxMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        }
        //cudaDeviceSynchronize();
        if (preallocation) {
            cudaMemcpy(_preallocationFloat, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
            return sqrt(*_preallocationFloat);
        }
        else
        {
            cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
            std::cout << "free !!!" << std::endl;
            cudaFree(d_odata);
            return sqrt(odata);
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

double MatrixGPUD::distance2(MatrixGPUD* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("not the same size");
    }
    if (_GPU && m->getPos())
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double odata = 0;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }
        //std::cout << _numBlocks << std::endl;

        switch (_blockSize) {
        case 512:
            distanceMultiBlock<512> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<512> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 256:
            distanceMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 128:
            distanceMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 64:
            distanceMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            distanceMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            distanceMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            distanceMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            distanceMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            distanceMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            distanceMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, d_odata, n);
            sumMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        }
        //cudaDeviceSynchronize();
        cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }


        return sqrtf(odata);
    }
    else if (!_GPU && !(m->getPos()))
    {
        double d = 0;
        double r = 0;
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                r = get(i, j) - m->get(i, j);
                d = d + r * r;
            }
        }
        return sqrtf(d);
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

void MatrixGPUD::Moy(MatrixGPUD* m, MatrixGPUD* nb, int sens)
{
    double s;
    int n;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol()) || (_column != nb->getNCol()) || (nb->getNLin() != 1))
        {
            throw std::invalid_argument("wrong dimension of the vector");
        }
        if (_GPU && nb->getPos() && m->getPos())
        {
            
            moyGPU1 <<<_numBlocks, _blockSize >>> (_matrixGPU, m->_matrixGPU, nb->_matrixGPU, m->getNLin(), _column);
        }
        else if ((!_GPU) && !(nb->getPos()) && !(m->getPos())) {
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
        else {
            throw std::invalid_argument("Matrix not at the same place");
        }
        

    }
    else { // on travaille sur les lignes 
        if ((_column != 1) || (_row != m->getNLin()) || (_row != nb->getNLin()) || (nb->getNCol() != 1)) {
            throw std::invalid_argument("wrong dimension of the vector");
        }
        if (_GPU && nb->getPos())
        {
            
            moyGPU2<<<_numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, nb->_matrixGPU, _row , m->getNCol());
        }
        else if ((!_GPU) && !(nb->getPos())) {
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
        else {
            throw std::invalid_argument("Matrix not at the same place");
        }

    }
}

void MatrixGPUD::project(MatrixGPUD* Lb, MatrixGPUD* Ub)
{
    if (!dim(Lb) || !dim(Ub)) {
        throw std::invalid_argument("not the same dimension");
    }
    if (_GPU && Lb->getPos() && Ub->getPos())
    {
        projectGPU<<<_numBlocks, _blockSize >>>(_matrixGPU, Lb->_matrixGPU, Ub->_matrixGPU, _N);
    }
    else if (!_GPU && !(Lb->getPos()) && !(Ub->getPos()))
    {
        double ub = 0;
        double lb = 0;
        double r = 0;
        MatrixGPUD temp(*this);
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
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
    
}
void MatrixGPUD::projectNeg()
{
    if (_GPU)
    {
        projectGPUNeg << <_numBlocks, _blockSize >> > (_matrixGPU, _N);
    }
    else if (!_GPU)
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

}
void MatrixGPUD::projectPos()
{
   
    if (_GPU)
    {
        projectGPUPos <<<_numBlocks, _blockSize >> > (_matrixGPU, _N);
    }
    else if (!_GPU)
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
}



double MatrixGPUD::sum() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    if (_GPU) 
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double odata = 0;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }


        switch (_blockSize) {
        case 512:
            SumMultiBlock<512> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<512> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 256:
            SumMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 128:
            SumMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 64:
            SumMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            SumMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            SumMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            SumMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            SumMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            SumMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            SumMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        }
        //cudaDeviceSynchronize();
        cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        //std::cout << "sum " << odata << " " <<_blockSize << " " << numBlocks << std::endl;
        return odata;
    }
    else 
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
}

double MatrixGPUD::sum(int begin, int end)
{
    if (begin < 0 || end < 0) {
        throw std::invalid_argument("indice must be positve");
    }
    if (begin > end) {
        throw std::invalid_argument("begin must be smaller than end");
    }
    if (begin > _N || end > _N) {
        throw std::out_of_range("indice must smaller than N");
    }
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }

    if (_GPU)
    {
        int numBlocks = _numBlocks;
        double odata = 0;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }


        switch (_blockSize) {
        case 512:
            SumMultiBlock<512> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<512> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 256:
            SumMultiBlock<256> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<256> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 128:
            SumMultiBlock<128> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<128> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 64:
            SumMultiBlock< 64> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock< 64> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            SumMultiBlock< 32> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock< 32> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            SumMultiBlock< 16> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock< 16> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            SumMultiBlock<  8> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock< 8> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            SumMultiBlock<  4> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<  4> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            SumMultiBlock<  2> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<  2> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            SumMultiBlock<  1> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, begin, end);
            sumMonoBlock<  1> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        }
        //cudaDeviceSynchronize();
        cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        //std::cout << "sum " << odata << " " <<_blockSize << " " << numBlocks << std::endl;
        return odata;
    }
    else if (!_GPU)
    {
        double d = 0;
        double r = 0;
        for (int elem = begin; elem < end; ++elem)
        {

            r = _matrixCPU[elem];
            d = d + r;
        }
        return d;
    }
}

void MatrixGPUD::sum(MatrixGPUD* m)
{
    double s = 0;
     // on travaille sur les lignes 
    if ((_column != 1) || (_row != m->getNLin())) {
        throw std::invalid_argument("wrong dimension of the column vector ");
    }
    int col = m->getNCol();
    if (_GPU && m->getPos())
    {
        int numBlocks = _row;
        switch (_blockSize) {
        case 512:
            SumEachRow<512> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case 256:
            SumEachRow<256> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case 128:
            SumEachRow<128> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case 64:
            SumEachRow< 64> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case 32:
            SumEachRow< 32> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case 16:
            SumEachRow< 16> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case  8:
            SumEachRow<  8> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case  4:
            SumEachRow<  4> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case  2:
            SumEachRow<  2> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        case  1:
            SumEachRow<  1> << <numBlocks, _blockSize >> > (m->_matrixGPU, _matrixGPU, col);
            break;
        }
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0; i < _row;i++)
        {
            s = 0;
            for (int j = 0; j < col;j++)
            {
                s = s + m->get(i, j);
            }
            set(i, 0, s);
        }
    }
    else {
        throw std::invalid_argument("Matrix not at the same place");
    }
}

double MatrixGPUD::distance2() {

    if (_GPU ) //&& m->getPos())
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        double* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(double) * numBlocks);
        }
        double odata = 0;
        
        
        //std::cout << _numBlocks << std::endl;
 
        switch (_blockSize) {
        case 512:
            distanceMultiBlock<512> << <numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<512> << < 1, _blockSize >> > (d_odata, d_odata, numBlocks);
        break;
        case 256:
            distanceMultiBlock<256> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<256> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
        break;
        case 128:
            distanceMultiBlock<128> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<128> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
        break;
        case 64:
            distanceMultiBlock< 64> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 64> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 32:
            distanceMultiBlock< 32> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 32> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case 16:
            distanceMultiBlock< 16> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 16> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  8:
            distanceMultiBlock<  8> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock< 8> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  4:
            distanceMultiBlock<  4> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<  4> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  2:
            distanceMultiBlock<  2> <<<numBlocks, _blockSize >> > (_matrixGPU, d_odata, n);
            sumMonoBlock<  2> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break;
        case  1:
            distanceMultiBlock<  1> <<<numBlocks, _blockSize >>> (_matrixGPU, d_odata, n);
            sumMonoBlock<  1> <<< 1, _blockSize >> > (d_odata, d_odata, numBlocks);
            break; 
        }
        //cudaDeviceSynchronize();
        cudaMemcpy(&odata, d_odata, sizeof(double), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        return sqrtf(odata);
    }
    else if (!_GPU)// && !(m->getPos()))
    {
        double d = 0;
        double r = 0;
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                r = get(i, j);// -m->get(i, j);
                d = d + r * r;
            }
        }
        return sqrtf(d); 
    }
}



///////////////////////////////////////////////////////////////////////////////
// Display MatrixGPUD contents
///////////////////////////////////////////////////////////////////////////////
void MatrixGPUD::display(bool force) 
{   
    bool transfert = false;
    if (this) {
        if (_GPU && !force ) {
           std::cout << " Matrix stockee sur GPU, faire le transfertCPU avant d'afficher " << std::endl;
        }
        if (_row == 0 || _column == 0)
        {
            std::cout << "matrix vide " << std::endl;
            return;
        }
        else {
            
            if (_GPU) {
                transferCPU();
                transfert = true;
            }
            if (_column == 1) {
                std::cout << " transpose  : ";
                for (int i = 0;i < _row;++i)
                {
                    for (int j = 0;j < _column;++j)
                    {
                        double value = get(i, j);
                        std::cout << std::setprecision(7) << value;
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
                        std::cout << std::setprecision(7) << value;
                        //std::cout << std::fixed << std::setprecision(3) << value;
                        if (j != _column - 1) std::cout << " ";
                    }

                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            if (transfert) {
                transferGPU();
            }
        }
    }
    else 
    {
        std::cout << "matrix non definie " << std::endl;
    }
}

void MatrixGPUD::displayBloc(int iBegin, int iEnd, int jBegin, int jEnd, bool force)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("xBegin must be smaller than xEnd");
    }
    bool transfert = false;
    if (this) {
        if (_GPU && !force) {
            std::cout << " Matrix stockee sur GPU, faire le transfertCPU avant d'afficher " << std::endl;
        }
        else {
            if (_row == 0 || _column == 0)
            {
                std::cout << "matrix vide " << std::endl;
            }
            if (_GPU) {
                transferCPU();
                transfert = true;
            }
            if (jEnd - jBegin == 1 ) {
                std::cout << " transpose  : ";
                for (int i = iBegin; i < iEnd; ++i)
                {
                    double value = get(i, jBegin);
                    std::cout << std::setprecision(7) << value;
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
                        std::cout << std::setprecision(7) << value;
                        if (j != jEnd - 1) {
                            std::cout << " ";
                        }
                    }

                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            if (transfert) {
                transferGPU();
            }
        }
    }
    else
    {
        std::cout << "matrix non definie " << std::endl;
    }
}

void MatrixGPUD::swapLine(int line1, int line2)
{
    if (_GPU) {
        swapLineGJ << <_numBlocks, _blockSize >> > (_matrixGPU, line1, line2, _column);// swap des lignes
    }
    else {
        double temp = 0;
        for (int i = 0; i < _column; i++) {
            temp = get(line1, i);
            set(line1, i, get(line2, i));
            set(line2, i, temp);
        }
    }
    
}




///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
MatrixGPUD::~MatrixGPUD()
{
    #ifdef DEBUG_DESTRUCTOR
        std::cout << "destruction matrix " << _matrixGPU << std::endl;
    #endif // DEBUG_DESTRUCTOR
    if (_preallocationFloat != nullptr) {
        cudaFreeHost(_preallocationFloat);
        _preallocationFloat = nullptr;
    }
    if (_preallocation != nullptr) {
        cudaFree(_preallocation);
        _preallocation = nullptr;
    }
    if (_matrixGPU) {
        cudaFree(_matrixGPU);
        _matrixGPU = nullptr;
    }
     DELETEA(_matrixCPU);
    
}



void MatrixGPUD::saveCSV(const std::string& filename, std::ios_base::openmode mode, int trans) const
{
    if (_GPU) {
        throw std::domain_error("Matrix on GPU");
    }
    std::ofstream myfile;
    myfile.open(filename, mode);
    myfile.precision(50);
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

void MatrixGPUD::saveCSVForce(const std::string& filename, std::ios_base::openmode mode, int trans)
{
    int transfert = 0;
    if (_GPU) {
        transfert = 1;
        transferCPU();
    }
    std::ofstream myfile;
    myfile.open(filename, mode);
    myfile.precision(10);
    if (!trans) {
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++) {
                myfile << get(i, j) << ";";
            }
            myfile << "\n";
        }
    }
    else {
        for (int j = 0; j < _column; j++) {
            for (int i = 0; i < _row; i++) {
                myfile << get(i, j) << ";";
            }
            myfile << "\n";
        }
    }

    myfile.close();
    if (transfert) {
        transferGPU();
    }
}


///////////////////////////////////////////////////////////////////////////////
// Fonction globale
///////////////////////////////////////////////////////////////////////////////


__global__ void setup_kernelD(curandState* state) {

    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(1234, idx, 0, &state[idx]);
}


__global__ void generate_kernel(curandState* my_curandstate, double* result, double eps, const unsigned int N) {

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        result[i] = (2*curand_uniform(my_curandstate + i)-1) * eps;
    }
}





__global__ void setGPU(double* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i];
    }
}

__global__ void setGPU(double* mat1, const double value, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = value;
    }
}
__global__ void setGPUunique(double* mat1, const double value, int pos) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    if (index == 0) {
        mat1[pos] = value;
    }

}

__global__ void setTransGPU(double* mat1, double* matToTrans, const int column, const int row) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    const int N = column * row;

    for (int e = index; e < N; e += step)
    {
        int i = e / column;
        int j = e % column;
        mat1[e] = matToTrans[j * row + i];

    }
}

__global__ void setColGPU(double* mat1, double* mat2, const int numCol, const int column, const int row, const int offset) {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < row; i += step)
    {
        mat1[i] = i < offset ? 0 : mat2[i*column+numCol];
    }

}

__global__ void setEyesGPU(double* mat2, const double value, const int col, const int row) 
{
    int N = row * col;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int l = index; l < N; l+=step) {
        int i = l / col;
        int j = l % col;
        mat2[l] = (i==j) ? value : 0; // pas coalescent, mais bon...
    }
}
__global__ void setEyesGPU(double* mat2, double* mat1, const int col, const int row)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    int N = row * col;
    for (int l = index; l < N; l += step) {
        int i = l / col;
        int j = l % col;
        mat2[l] = (i == j) ? mat1[i] : 0; // pas coalescent, mais bon...
    }
}


__global__ void SetBlocGPU(double* out, double* in, int ibegin, int iend, int jbegin, int jend, int col)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int step = blockDim.x * gridDim.x;
    int offset = jbegin + ibegin * col;
    int N = (jend - jbegin) * (iend - ibegin);

    for (int j = index; j < N; j += step)
    {
        int rowLoc = j / (jend - jbegin);
        int colLoc = j % (jend - jbegin);
        int GlobalInd = offset + rowLoc *col + colLoc;
        out[GlobalInd] = in[j];
    }
}

__global__ void SetBlocGPU(double* out, double* in, int ibegin, int iend, int jbegin, int jend, int col, double factor)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int step = blockDim.x * gridDim.x;
    int offset = jbegin + ibegin * col;
    int N = (jend - jbegin) * (iend - ibegin);

    for (int j = index; j < N; j += step)
    {
        int rowLoc = j / (jend - jbegin);
        int colLoc = j % (jend - jbegin);
        int GlobalInd = offset + rowLoc * col + colLoc;
        out[GlobalInd] = factor * in[j];
    }
}

/*__global__ void SetBlocGPU(double* out, double* in, int ibegin, int iend, int jbegin, int jend, int col, double factor) // fait que la première ligne
{
    int indexX = threadIdx.x + blockIdx.x * blockDim.x;
    int stepX = blockDim.x * gridDim.x;
    int indexY = threadIdx.y + blockIdx.y * blockDim.y;
    int stepY = blockDim.y * gridDim.y;

    for (int j = indexX + jbegin; j < jend; j += stepX)
    {
        for (int i = indexY + ibegin; i < iend; i += stepY)
        {
            out[j + i * col] = factor * in[indexX + indexY * col];
        }
    }
}*/




__global__ void replaceGPU(double* mat,const double previous, const double newValue,const int N) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = (mat[i] == previous) * (newValue-mat[i]) + mat[i];
    }
}




__global__ void addGPU(double* mat, double c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] + c;
    }
}
__global__ void addGPU(double* mat1, double* mat2, double c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] + c;
    }
}
__global__ void addGPU(double* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] + mat2[i];
    }
}
__global__ void addGPU(double* mat1, double* mat2, double* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] + mat3[i];
    }
}

__global__ void addVectorGPU1(double* mat1, double* vect, const int n, int N) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i / n; // division entière
        mat1[i] = mat1[i] + vect[k];
    }

}
__global__ void addVectorGPU2(double* mat1, double* vect, const int n, int N) // vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i % n; // modulo
        mat1[i] = mat1[i] + vect[k];
    }


}

__global__ void addTransGPU(double* out, double* mat1, double* mat2, const int col, const int line, int N) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int l = index; l < N; l += step)
    {
        int i = l / col;
        int j = l % col;
        int k = i + j * line;
        out[l] = mat1[l] + mat2[k];
    }
}

__global__ void substractGPU(double* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] - mat2[i];
    }
}
__global__ void substractGPU(double* mat1, double* mat2, double* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] - mat3[i];
    }
}

__global__ void substractVectorGPU1(double* mat1, double* vect, const int n, int N) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i / n; // division entière
        mat1[i] = mat1[i] - vect[k];
    }

}
__global__ void substractVectorGPU2(double* mat1, double* vect, const int n, int N) // vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i % n; // modulo
        mat1[i] = mat1[i] - vect[k];
    }

}

__global__ void substractTransGPU(double* out, double* mat1, double* mat2, const int col, const int line, int N)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int l = index; l < N; l += step)
    {
        int i = l / col;
        int j = l % col;
        int k = i + j * line;
        out[l] = mat1[l] - mat2[k];
    }
}

__global__ void multiplyGPU(double* mat, const double c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] * c;
    }
}

__global__ void multiplyTGPU(double* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] * mat2[i];
    }
}
__global__ void multiplyTGPU(double* mat1, double* mat2, double* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] * mat3[i];
    }
}

__global__ void divideGPU(double* mat, const double c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] / c;
    }
}
__global__ void divideGPU(double* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] / mat2[i];
    }
}

__global__ void moyGPU1(double* res, double* mat1, double* nb, const int line, const int column) //vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
     
    for (int i = index; i < column; i += step)
    {
        double s = 0.0;
        for (int j = 0; j < line; j++)
        {
            s += mat1[i + column *j];
        }
        res[i] = s / nb[i];
    }

}
__global__ void moyGPU2(double* res, double* mat1, double* nb, const int line, const int column) // vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    
    for (int i = index; i < line; i += step)
    {
        double s = 0.0;
        for (int j = 0; j < column; j++)
        {
            s +=  mat1[i*column + j];
        }
        res[i] = s /nb[i];
    }
}

__global__ void projectGPU(double* mat, double* Lb, double* Ub, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        double r = mat[i];
        double ub = Ub[i];
        double lb = Lb[i];
        r = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
        mat[i] = r;//(Ub[i] - mat[i])* (mat[i] > Ub[i]) + (Lb[i] - mat[i]) * (mat[i] < Lb[i]) + mat[i];
    }
}

__global__ void projectGPUPos(double* mat, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        double r = mat[i];
        mat[i] = (r > 0) * r;
    }
}

__global__ void projectGPUNeg(double* mat, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        double r = mat[i];
        mat[i] = (r < 0) * r;
    }
}


__global__ void sumGPU(double* res, double* mat1, const int line, const int column) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < line; i += step)
    {
        double s = 0.0;
        for (int j = 0; j < column; j++)
        {
            s += mat1[i*column + j];
        }
        res[i] = s;
    }
}

__global__ void sumGPU2(double* res, double* mat1, const int line) //vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index==0)
    {
        double s = 0.0;
        for (int j = 0; j < line; j++)
        {
            s += mat1[j];
        }
        
        *res = s ;
    }
}



__device__ int sumCommSingleWarp(volatile double* shArr) {
    int idx = threadIdx.x % warpSizeD; //the lane index in the warp
    if (idx < 16) {
        shArr[idx] += shArr[idx + 16];
        shArr[idx] += shArr[idx + 8];
        shArr[idx] += shArr[idx + 4];
        shArr[idx] += shArr[idx + 2];
        shArr[idx] += shArr[idx + 1];
    }
    return shArr[0];
}

template <unsigned int blockSize>
__global__ void sumMonoBlock(double* g_idata, double* g_odata, unsigned int n) {
    
    int idx = threadIdx.x;
    double sum = 0;
    for (int i = idx; i < n; i += blockSize)
        sum += g_idata[i];
    __shared__ double r[blockSize];
    r[idx] = sum;
    sumCommSingleWarp(&r[idx & ~(warpSizeD - 1)]);
    __syncthreads();
    if (idx < warpSizeD) { //first warp only
        r[idx] = idx * warpSizeD < blockSize ? r[idx * warpSizeD] : 0;
        sumCommSingleWarp(r);
        if (idx == 0)
            *g_odata = r[0];
    }
}


template <unsigned int blockSize>
__global__ void SumMultiBlock(double* g_idata, double* g_odata, unsigned int n) {

    __shared__ double shArr[blockSize];
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double sum = 0;
    for (int i = gthIdx; i < n; i += gridSize)
       sum += g_idata[i];
    
    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] += shArr[thIdx + size];
        __syncthreads();
    }
    if (thIdx == 0) {
        g_odata[blockIdx.x] = shArr[0];
    }
}

template <unsigned int blockSize>
__global__ void SumMultiBlock(double* g_idata, double* g_odata, unsigned int begin, unsigned int end) {
    __shared__ float shArr[blockSize];
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double sum = 0;
    for (int i = gthIdx + begin; i < end; i += gridSize)
        sum += g_idata[i];

    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] += shArr[thIdx + size];
        __syncthreads();
    }
    if (thIdx == 0) {
        g_odata[blockIdx.x] = shArr[0];
    }
}


template <unsigned int blockSize>
__global__ void SumEachRow(double* g_idata, double* g_odata, const int nCol) {
    __shared__ double shArr[blockSize];
    int thIdx = threadIdx.x;
    int row = blockIdx.x;
    int idBegin = thIdx + row * nCol;
    int idEnd = (row + 1) * nCol;
    int step = blockDim.x;

    double sum = 0;
    for (int i = idBegin; i < idEnd; i += step) {
        sum += g_idata[i]; 
    }
        

    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] += shArr[thIdx + size];
        __syncthreads();
    }
    if (thIdx == 0) {
        g_odata[blockIdx.x] = shArr[0];
    }
}

template <unsigned int blockSize>
__global__ void distanceMultiBlock(double* g_idata, double* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double sum = 0;
    for (int i = gthIdx; i < n; i += gridSize)
        sum += (g_idata[i] * g_idata[i]);
    __shared__ double shArr[blockSize];
    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] += shArr[thIdx + size];
        __syncthreads();
    }
    if (thIdx == 0)
        g_odata[blockIdx.x] = shArr[0];
   
}



template <unsigned int blockSize>
__global__ void distanceMultiBlock(double* g_idata, double* g_idata2, double* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double sum = 0;
    for (int i = gthIdx; i < n; i += gridSize)
        sum += ((g_idata[i]- g_idata2[i]) * (g_idata[i] - g_idata2[i]));
    __shared__ double shArr[blockSize];
    shArr[thIdx] = sum;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] += shArr[thIdx + size];
        __syncthreads();
    }
    if (thIdx == 0)
        g_odata[blockIdx.x] = shArr[0];

}

__device__ double warpReduceMax(volatile double* r) {
    int idx = threadIdx.x % warpSizeD; //the lane index in the warp
    if (idx < 16) {
        r[idx] = r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] * (r[idx + 16] > r[idx]) + r[idx] * (r[idx] <= r[idx + 16]);
        r[idx] = r[idx + 8] > r[idx] ? r[idx + 8] : r[idx];//r[idx +  8] > r[idx] ? r[idx +  8] : r[idx];//r[idx +  8] * (r[idx +  8] > r[idx]) + r[idx] * (r[idx] <= r[idx +  8]);
        r[idx] = r[idx + 4] > r[idx] ? r[idx + 4] : r[idx];//r[idx +  4] > r[idx] ? r[idx +  4] : r[idx];//r[idx +  4] * (r[idx +  4] > r[idx]) + r[idx] * (r[idx] <= r[idx +  4]);
        r[idx] = r[idx + 2] > r[idx] ? r[idx + 2] : r[idx];//r[idx +  2] > r[idx] ? r[idx +  2] : r[idx];//r[idx +  2] * (r[idx +  2] > r[idx]) + r[idx] * (r[idx] <= r[idx +  2]);
        r[idx] = r[idx + 1] > r[idx] ? r[idx + 1] : r[idx];//r[idx +  1] > r[idx] ? r[idx +  1] : r[idx];//r[idx +  1] * (r[idx +  1] > r[idx]) + r[idx] * (r[idx] <= r[idx +  1]);
    }
    return r[0];
}

__device__ void warpReduceMax(volatile double* r, volatile int* pos){
    int idx = threadIdx.x % warpSizeD; //the lane index in the warp
    if (idx < 16) {
        pos[idx] = r[idx + 16] > r[idx] ? pos[idx + 16] : pos[idx];
        r[idx] = r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] * (r[idx + 16] > r[idx]) + r[idx] * (r[idx] <= r[idx + 16]);
        
        pos[idx] = r[idx + 8] > r[idx] ? pos[idx + 8] : pos[idx];
        r[idx] = r[idx + 8] > r[idx] ? r[idx + 8] : r[idx];//r[idx +  8] > r[idx] ? r[idx +  8] : r[idx];//r[idx +  8] * (r[idx +  8] > r[idx]) + r[idx] * (r[idx] <= r[idx +  8]);
        
        pos[idx] = r[idx + 4] > r[idx] ? pos[idx + 4] : pos[idx];
        r[idx] = r[idx + 4] > r[idx] ? r[idx + 4] : r[idx];//r[idx +  4] > r[idx] ? r[idx +  4] : r[idx];//r[idx +  4] * (r[idx +  4] > r[idx]) + r[idx] * (r[idx] <= r[idx +  4]);
        
        pos[idx] = r[idx + 2] > r[idx] ? pos[idx + 2] : pos[idx];
        r[idx] = r[idx + 2] > r[idx] ? r[idx + 2] : r[idx];//r[idx +  2] > r[idx] ? r[idx +  2] : r[idx];//r[idx +  2] * (r[idx +  2] > r[idx]) + r[idx] * (r[idx] <= r[idx +  2]);
        
        pos[idx] = r[idx + 1] > r[idx] ? pos[idx + 1] : pos[idx];
        r[idx] = r[idx + 1] > r[idx] ? r[idx + 1] : r[idx];//r[idx +  1] > r[idx] ? r[idx +  1] : r[idx];//r[idx +  1] * (r[idx +  1] > r[idx]) + r[idx] * (r[idx] <= r[idx +  1]);
    }
}

template <unsigned int blockSize>
__global__ void maxMonoBlock(double* g_idata, double* g_odata, unsigned int n) {
    int idx = threadIdx.x;
    double max = 0;
    for (int i = idx; i < n; i += blockSize) {
        double s = g_idata[i];
        max = s > max ? s : max;// s>max ? s:max;//s * (s > max) + max * (max <= s);
    }
    __shared__ double r[blockSize];
    r[idx] = max;
    warpReduceMax(&r[idx & ~(warpSizeD - 1)]);
    __syncthreads();
    if (idx < warpSizeD) { //first warp only
        r[idx] = idx * warpSizeD < blockSize ? r[idx * warpSizeD] : 0;
        warpReduceMax(r);
        if (idx == 0)
            *g_odata = r[0];
    }
}
template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double max = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        double s = (g_idata[i] * g_idata[i]);
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ double shArr[blockSize];
    shArr[thIdx] = max;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size)
            shArr[thIdx] = shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 
        __syncthreads();
    }
    if (thIdx == 0)
        g_odata[blockIdx.x] = shArr[0];
}

template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_odata, unsigned int n, int* pos) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double max = 0;
    int indice = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        double s = (g_idata[i] * g_idata[i]);
        indice = s > max ? i : indice;
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ double shArr[blockSize];
    __shared__ double shPos[blockSize];
    shArr[thIdx] = max;
    shPos[thIdx] = indice;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //can unroll the for loop
        if (thIdx < size) {
            shPos[thIdx] = shArr[thIdx + size] > shArr[thIdx] ? shPos[thIdx + size] : shPos[thIdx];
            shArr[thIdx] = shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx]; 
            //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 
        }
           
        __syncthreads();
    }
    if (thIdx == 0) {
        g_odata[blockIdx.x] = shArr[0];
        pos[blockIdx.x] = shPos[0];

    }
}
       

template <unsigned int blockSize>
__global__ void maxMonoBlock(double* g_idata, double* g_odata, unsigned int n, int* pos) {
    int idx = threadIdx.x;
    double max = 0;
    int indice = 0;
    for (int i = idx; i < n; i += blockSize) {
        double s = g_idata[i];
        indice = s > max ? pos[i] : indice;
        max = s > max ? s : max;// s>max ? s:max;//s * (s > max) + max * (max <= s);
    }
    __shared__ double r[blockSize];
    __shared__ int shPos[blockSize];
    r[idx] = max;
    shPos[idx] = indice;
    warpReduceMax(&r[idx & ~(warpSizeD - 1)], &shPos[idx & ~(warpSizeD - 1)]);
    __syncthreads();
    if (idx < warpSizeD) { //first warp only
        r[idx] = idx * warpSizeD < blockSize ? r[idx * warpSizeD] : 0;
        warpReduceMax(r, shPos);
        if (idx == 0) {
            *g_odata = r[0];
            *pos = shPos[0];
        }
            
        
    }
}


template <unsigned int blockSize>
__global__ void maxMultiBlock(double* g_idata, double* g_idata2, double* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    double max = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        double s = (g_idata[i] - g_idata2[i]);
        s = s*s;
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ double shArr[blockSize];
    shArr[thIdx] = max;
    __syncthreads();
    for (int size = blockSize / 2; size > 0; size /= 2) { //uniform
        if (thIdx < size) {
            shArr[thIdx] = shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 
        }
            
        __syncthreads();
    }
    if (thIdx == 0) {
        g_odata[blockIdx.x] = shArr[0];
    }
}


__global__ void normalisationGJ(double* mat, const int row, const int nCol, const double factor) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < nCol; i += step)
    {
        mat[i + row * nCol] = mat[i + row * nCol] / factor;
    }


}

__global__ void swapLineGJ(double* mat, const int row1, const int row2, const int nCol) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < nCol; i += step)
    {
        double temp = mat[i + row1 * nCol];
        double temp2 = mat[i + row2 * nCol]; 
        mat[i + row1 * nCol] = temp2; // or mat[i + row1 * nCol] = mat[i + row2 * nCol];
        mat[i + row2 * nCol] = temp;
    }
}

__global__ void eliminationGJ(double* mat, double* matAug, const int r, const int nRow, const int nCol) {

    // un bloc = une ligne, 
    int index = threadIdx.x;
    int row = blockIdx.x;
    int step = blockDim.x;
    __shared__ double shFactor;
    if (row != r) { // le bloc r ne fait rien... bah...
        if (index == 0) {
            shFactor = mat[row * nCol + r];
        }
        __syncthreads();
        for (int j = index; j < nCol; j+=step) {
            double value1 = mat[r * nCol + j];
            double oldvalue1 = mat[row * nCol + j];
            double oldvalue2 = matAug[row * nCol + j];
            double value2 = matAug[r * nCol + j];

            mat[row * nCol + j] = oldvalue1 - shFactor * value1;
            //matAug[row * nCol + j] -=  shFactor * value2;
            matAug[row * nCol + j] = oldvalue2 - shFactor * value2;
        }
    }
    

    /*int indexX = blockIdx.x * blockDim.x + threadIdx.x;
    int indexY = blockIdx.y * blockDim.y + threadIdx.y;
    int stepX = blockDim.x * gridDim.x;
    int stepY = blockDim.y * gridDim.y;


    for (int i = indexY; i < nRow; i += stepY)
    {
        if (i != r) {
            double factor = mat[i * nCol + r]; // ne doit pas changer tant que la ligne n'est pas fini
            for (int j = indexX; j < nCol; j += stepX)
            {
                if (j != r) {
                    double value1 = mat[r * nCol + j];
                    double value2 = matAug[r * nCol + j];

                    mat[i * nCol + j] = mat[i * nCol + j] - factor * value1;
                    matAug[i * nCol + j] = matAug[i * nCol + j] - factor * value2;
                }
            }
        }
    }
    __syncthreads();
    for (int i = indexY; i < nRow; i += stepY)
    {
        double factor = mat[i * nCol + r]; // ne doit pas changer tant que la ligne n'est pas fini
        if (i != r) {
            if (indexX == 1) {
                mat[i * nCol + r] = 0;
                matAug[i * nCol + r] = matAug[i * nCol + r] - factor * matAug[r * nCol + r];
            }
        }
    }*/
}



__global__ void initPermMatr(double* P, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < (N + 1); i += step)
    {
        P[i] = i * (i < N);
    }
}


__global__ void updatePermMatr(double* P, const int line1, const int line2, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    if (index == 0) {
        int inter = P[line1];
        P[line1] = P[line2];
        P[line2] = inter;
        P[N] = P[N] + 1;
    }
}


__global__ void updateLUPFactorization(double* A, const int col, const int N) {
    // un bloc par ligne i ?
    int index = threadIdx.x;
    int i = blockIdx.x;
    int step = blockDim.x;

    __shared__ double Aicol;

    if (i > col) { // les blocs trop petits ne font rien, en s'en fout ?
        if (index == 0) {
            Aicol = A[i * N + col] / A[col * N + col];
            A[i * N + col] = Aicol;
        }
        __syncthreads();
        for (int k = index + col + 1; k < N; k += step) {
            A[i * N + k] = A[i * N + k] - Aicol * A[col * N + k];
        }
    }
}





__global__ void setPermute(double* y, double* b, double* P, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < N; i += step)
    {
        int indice = P[i];
        y[i] = b[indice]; // c'est absolument moche...
    }


}
__global__ void solveLowSys(double* A, double* y, const int iter, const int N) {
    int index = threadIdx.x;
    int step = blockDim.x;
    __shared__ double yiter;

    if (index == 0) {
        yiter = y[iter];

    }
    __syncthreads();
    for (int i = index + iter + 1; i < N; i += step)
    {
        y[i] = y[i] - y[iter] * A[i * N + iter]; // moche ne faudrait-il pas stocker A^T ?

    }


}

__global__ void solveUpSys(double* A, double* y, const int iter, const int N) {
    int index = threadIdx.x;
    int step = blockDim.x;
    __shared__ double yiter;

    if (index == 0) {
        yiter = y[iter] / A[iter * N + iter];
        y[iter] = yiter;
    }
    __syncthreads();

    for (int i = index; i < iter; i += step)
    {
        y[i] = y[i] - yiter * A[i * N + iter]; // moche ne faudrait-il pas stocker A^T ?

    }


}

__global__ void solveSysGPU(double* A, double* y, const int N) {


    int index = threadIdx.x;
    int step = blockDim.x;
    extern __shared__ double ytemp[];


    for (int n = index; n < N; n += step)
    {
        ytemp[n] = y[n];
    }
    __syncthreads();
    for (int iter = 0; iter < N; iter++) {
        for (int i = index + iter + 1; i < N; i += step)
        {
            ytemp[i] = ytemp[i] - ytemp[iter] * A[i * N + iter]; // moche ne faudrait-il pas stocker A^T ?

        }
        __syncthreads();
    }
    for (int iter = N - 1; iter >= 0; iter--) {
        if (index == 0) {
            ytemp[iter] = ytemp[iter] / A[iter * N + iter];
            y[iter] = ytemp[iter];
        }
        __syncthreads();

        for (int i = index; i < iter; i += step)
        {
            ytemp[i] = ytemp[i] - ytemp[iter] * A[i * N + iter]; // moche ne faudrait-il pas stocker A^T ?

        }
        __syncthreads();
    }
}