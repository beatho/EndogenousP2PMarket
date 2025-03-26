#include "../head/MatrixGPU.cuh" 

#define MAX(X, Y) X * (X >= Y) + Y * (Y > X)
const int warpSize = 32;

#define CHECK_CUDA_ERROR_MAT(val) checkMat((val), #val, __FILE__, __LINE__);
#define CHECK_LAST_CUDA_ERROR_MAT() checkMatLast(__FILE__, __LINE__);

template <typename T>
void checkMat(T err, char const* const func, char const* const file,
    int const line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
            << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        std::exit(EXIT_FAILURE);
    }
}
void checkMatLast(char const* const file, int const line)
{
    cudaError_t err{ cudaGetLastError() };
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
            << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

float MatrixGPU::rand1()
{
    float a = (float)(rand()) / ((float)(RAND_MAX));
    return a;
}


///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MatrixGPU::MatrixGPU() {
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur appele" << std::endl;
#endif
    _row = 0;
    _column = 0;
    _N = _row * _column;
     _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize),1);
}

MatrixGPU::MatrixGPU(int l, int c, float value, bool pos)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur parametre appele" << std::endl;
    std::cout << _matrixCPU << std::endl;
#endif
    _row = l;
    _column = c;
    _N = _row * _column;
    _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);
    if (pos) {
        if (_N > 0) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _N);
            setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, value, _N);
        }
        _GPU = true;
    }
    else {
        if (_N > 0) {
            _matrixCPU = new float[_N];
        }
        for (int elem = 0; elem < _N; elem++) {
            _matrixCPU[elem] = value;
        }
    }
#ifdef DEBUG_CONSTRUCTOR
    std::cout << _matrixGPU << std::endl;
#endif
}

MatrixGPU::MatrixGPU(const MatrixCPU& m, bool pos)
{
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);

    if (pos) {
        _GPU = true;
        if (_N > 0) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
            cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
        }
    }
    else {
        if (_N > 0) {
            _matrixCPU = new float[_row * _column];
            memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
        }
    }
    
}

MatrixGPU::MatrixGPU(const MatrixGPU & m)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "contructeur recopie appele" << std::endl;
#endif
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);

    if (m.getPos()) {
        if (_N > 0) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
            setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
        }
        _GPU = true;
    }
    else {
        if (_N > 0) {
            _matrixCPU = new float[_row * _column];
            memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
        }
    }
}

MatrixGPU::MatrixGPU(const MatrixGPUD& m)
{
    _row = m.getNLin();
    _column = m.getNCol();
    _N = _row * _column;
    _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);

    if (m.getPos()) {
        if (_N > 0) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
            setGPUFD << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
        }
        _GPU = true;
}
    else {
        if (_N > 0) {
            _matrixCPU = new float[_row * _column];
        }
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++)
            {
                set(i, j, m.get(i, j));
            }
        }
    }
}

MatrixGPU& MatrixGPU::operator=(const MatrixGPU& m)
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
                if (_N > 0) {
                    cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
                }
            }
        }
        else {
            if (m.getPos()) {
                if (_N > 0) {
                    cudaMemcpy(_matrixCPU, m._matrixGPU, sizeof(float) * _row * _column, cudaMemcpyDeviceToHost);
                }
            }
            else {
                if (_N > 0) {
                    memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
                }
            }
        }
    }
    else {
        _row = m.getNLin();
        _column = m.getNCol();
        _N = _row * _column;
        _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);
        _GPU = false;
        if (_matrixGPU) {
            cudaFree(_matrixGPU);
            _matrixGPU = nullptr;
        }
        DELETEA(_matrixCPU);
        if (m.getPos()) {
            if (_N > 0) {
                cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
                setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
            }
            _GPU = true;
        }
        else {
            if (_N > 0) {
                _matrixCPU = new float[_row * _column];
                memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
            }
        }
    }
   
    return *this;
}

MatrixGPU& MatrixGPU::operator=(const MatrixGPUD& m)
{
    if (_row == m.getNLin() && _column == m.getNCol()) {
        //matrix already has the good size no free needed
        if (getPos()) {
            if (m.getPos()) {
                setGPUFD << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
            }
            else {
                if (_matrixCPU == nullptr) {
                    _matrixCPU = new float[_row * _column];
                }
                for (int i = 0; i < _row; i++) {
                    for (int j = 0; j < _column; j++)
                    {
                        set(i, j, m.get(i, j));
                    }
                }
                transferGPU();
            }
        }
        else {
            if (m.getPos()) {
                if (_matrixGPU == nullptr) {
                    cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
                }
                setGPUFD << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
                _GPU = true;
                transferCPU();
            }
            else {
                for (int i = 0; i < _row; i++) {
                    for (int j = 0; j < _column; j++)
                    {
                        set(i, j, m.get(i, j));
                    }
                }
            }
        }
    }
    else {
        _row = m.getNLin();
        _column = m.getNCol();
        _N = _row * _column;
        _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize), 1);
        if (_matrixGPU) {
            cudaFree(_matrixGPU);
            _matrixGPU = nullptr;
        }
        DELETEA(_matrixCPU);
        if (m.getPos()) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
            setGPUFD << <_numBlocks, _blockSize >> > (_matrixGPU, m._matrixGPU, _N);
            _GPU = true;
        }
        else {
            _matrixCPU = new float[_row * _column];
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _column; j++)
                {
                    set(i, j, m.get(i, j));
                }
            }
        }
    }
    return *this;
}

MatrixGPU& MatrixGPU::operator=(const MatrixCPU& m)
{
    if (_row == m.getNLin() && _column == m.getNCol()) {
        //matrix already has the good size no free needed
        if (getPos()) {
            if (_N > 0) {
                cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
            }
        }
        else {
            if (_N > 0) {
                memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
            }
        }
    }
    else {
        _row = m.getNLin();
        _column = m.getNCol();
        _N = _row * _column;
        _numBlocks = MAX(ceil((_N + _blockSize - 1) / _blockSize),1);
        if (_matrixGPU) {
            cudaFree(_matrixGPU);
            _matrixGPU = nullptr;
        }
        DELETEA(_matrixCPU);

        if (getPos()) {
            if (_N > 0) {
                cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
                cudaMemcpy(_matrixGPU, m._matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
            }
            _GPU = true;
        }
        else
        {
            if (_N > 0) {
                _matrixCPU = new float[_row * _column];
                memcpy(_matrixCPU, m._matrixCPU, _row * _column * sizeof(float));
            }
        } 
    }
    return *this;
}

void MatrixGPU::preallocateReduction()
{
    CHECK_LAST_CUDA_ERROR();
    if (preallocation) {
        cudaFreeHost(_preallocationFloat);
        cudaFree(_preallocation);
        preallocation = false;
    }
    cudaError_t c;
    int counter = 0;
     /* do
    {
        c = cudaHostAlloc(&_preallocationFloat, sizeof(float), cudaHostAllocDefault);
        counter++;
    } while (_preallocationFloat == nullptr && counter < 10);
   
     
    if (_preallocationFloat == nullptr) {
            
        std::cout << "prealocation echouer ? " << c << std::endl;
        std::cout <<  cudaGetErrorName(c) << std::endl;
    }/
     if (c == 700) {
        std::cout << "c=700" << std::endl;
        exit(-1);
     }*/
     do
     {
         c = cudaMalloc((void**)&_preallocation, sizeof(float) * _numBlocks);
         counter++;
         cudaDeviceSynchronize();
     } while (_preallocation == nullptr && counter < 20);

    if (_preallocation == nullptr) {
        std::cout << _row << " " << _column << " " << _blockSize << std::endl;
        std::cout << "prealocation echouer ? " << c << " " << _numBlocks  <<std::endl;
        std::cout << cudaGetErrorName(c) << std::endl;
    }
    if (c == 700) {
        exit(-1);
    }/**/
    
    preallocation = true;
    setGPU <<<_numBlocks, _blockSize >>> (_preallocation, 0.0f, _numBlocks);
}

void MatrixGPU::transferGPU()
{
    if (!_GPU) {
        if (!_matrixGPU) {
            cudaMalloc((void**)&_matrixGPU, sizeof(float) * _row * _column);
        }
        cudaMemcpy(_matrixGPU, _matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
        //DELETEA(_matrixCPU);
        _GPU = true;
    }
    else {
        throw std::domain_error("already in the GPU");
    }
    
}

void MatrixGPU::transferCPU()
{
    
    if (_GPU) {
        
        if (!_matrixCPU) {
            
            _matrixCPU = new float[_row * _column];
        }
        cudaMemcpy(_matrixCPU, _matrixGPU, sizeof(float) * _row * _column, cudaMemcpyDeviceToHost);
        //cudaFree(_matrixGPU);
        //_matrixGPU = nullptr;
        _GPU = false;
    }
    else {
        std::cout << "transferCPU : already in the CPU " << _GPU <<std::endl;
        throw std::domain_error("already in the CPU");
    }

}

///////////////////////////////////////////////////////////////////////////////
// Getter
///////////////////////////////////////////////////////////////////////////////
 float MatrixGPU::get(int i, int j, bool verbose) const
{
    //std::cout << "hey de taille " << _row << " " << _column << "pos "<< i <<" "<< j << std::endl;
    if ((i >= _row) || ( j >= _column) || (i < 0) || ( j < 0)) {
        std::cout << "get" << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("index out of bounds");
    }
    if (_GPU) {
        float value;
        cudaMemcpy(&value, _matrixGPU + i*_column+j, sizeof(float), cudaMemcpyDeviceToHost);
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

int MatrixGPU::getNCol() const
{
    return _column;
}

int MatrixGPU::getNLin() const
{
    return _row;
}

void MatrixGPU::getCol(MatrixGPU* col, int numCol, int offset)
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
    } else {
        throw std::invalid_argument("getCol Matrix not at the same place");
    }

}


bool  MatrixGPU::getPos() const
{
    return _GPU;
}
bool MatrixGPU::dim(MatrixGPU* m) const
{ 
    return ((_row == m->getNLin()) && (_column == m->getNCol()));
}


bool MatrixGPU::isEqual(MatrixGPU* m, float pre) const
{
    if (!dim(m)) {
        throw std::invalid_argument("is Equal : not the same dimension");
    }
    else {
        if (_GPU || m->getPos()) {
            throw std::invalid_argument("is Equal : Matrix on GPU");
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

void MatrixGPU::toMatCPU(MatrixCPU& m) const // passer m en paramètre
{
    if (m.getNCol() != _column || m.getNLin() != _row) {
        m.setSize(_row, _column);
    }
    if (_GPU) {
        cudaMemcpy(m._matrixCPU, _matrixGPU, sizeof(float) * _row * _column, cudaMemcpyDeviceToHost);
    }
    else {
        memcpy(m._matrixCPU, _matrixCPU, sizeof(float) * _row * _column );
        /*for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++) 
            {
                m.set(i, j, get(i, j));
            }
        }*/
    }
}

void MatrixGPU::toMatGPUD(MatrixGPUD& m) const
{
    if (m.getNCol() != _column || m.getNLin() != _row) {
        std::cout << "pas de bonne taille" << std::endl;
        m.setSize(_row, _column);
    }
    if (_GPU) {
        if (!m.getPos()) {
            std::cout << "transfert GPU" << std::endl;
            m.transferGPU();
        }
       
        setGPUDF << <_numBlocks, _blockSize >> > (m._matrixCPU, _matrixGPU, _N);
        m.transferCPU();
        m.display();
        m.transferGPU();

    }
    else {
        if (m.getPos()) {
            m.transferCPU();
        }
        for (int i = 0; i < _row; i++) {
            for (int j = 0; j < _column; j++)
            {
                m.set(i, j, get(i, j));
            }
        }
    }

}



///////////////////////////////////////////////////////////////////////////////
// Setter
///////////////////////////////////////////////////////////////////////////////
 void MatrixGPU::set(int i, int j, float value, bool force)
{
    if ((i >= _row) || (j >= _column) || (i < 0) || (j < 0)) {
        std::cout << _row << " " << _column << " " << i << " " << j << std::endl;
        throw std::out_of_range("set : index out of bounds");
    }
    if (_GPU && !force) {
        throw std::invalid_argument("set : Matrix on GPU");
    }
    else if (_GPU && force) {
        setGPUunique <<< 1, 1 >>> (_matrixGPU, value, i * _column + j);
    }
    else {
        //std::cout << "changement de valeur " << value << " en " << i << " " << j << std::endl;
        _matrixCPU[i * _column + j] = value;
    }
}

 void MatrixGPU::setEyes(float value)
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

 void MatrixGPU::setEyes(MatrixGPU* m)
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
         throw std::invalid_argument("setEyes Matrix not at the same place");
     }



 }

void MatrixGPU::set(MatrixGPU* m, bool synchrone, cudaStream_t stream)
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
        throw std::invalid_argument("set Matrix not at the same place");
    }

}

void MatrixGPU::set(MatrixCPU* m)
{
    if (m->getNCol() != _column || m->getNLin() != _row) {
        throw std::invalid_argument("not the same dimension");
    }
    
    if (getPos()) {
        cudaMemcpy(_matrixGPU, m->_matrixCPU, sizeof(float) * _row * _column, cudaMemcpyHostToDevice);
    }
    else {
        memcpy(_matrixCPU, m->_matrixCPU, _row * _column * sizeof(float));
    }
}

void MatrixGPU::set(double value)
{
   
    if (_GPU) {
        setGPU << <_numBlocks, _blockSize >> > (_matrixGPU, value, _N);
    
    }
    else {
        for (int elem = 0; elem < _N; elem++) {
            _matrixCPU[elem] = value;
        }
    }
}

void MatrixGPU::setTrans(MatrixGPU* m)
{
    if (_column != m->getNLin() || _row != m->getNCol()) {
        std::cout << _row << " " << _column << " " << m->getNLin() <<" " << m->getNCol() <<std::endl;
        throw std::invalid_argument(" setTans : not the same transposed dimension");
    }
    if (getPos() && m->getPos()) {
        if (_N > 0) {
            setTransGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, _column, _row);
        }
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
        throw std::invalid_argument("setTrans Matrix not at the same place");
    }

}


void MatrixGPU::setRand(float eps)
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
        setup_kernel <<<_numBlocks, _blockSize >>> (state);
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

void MatrixGPU::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPU* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        std::cout << iBegin << " " << iEnd  << " " << jBegin << " " << jEnd << " " << _row << " " << _column << std::endl;
        throw std::out_of_range(" setBloc : index out of bounds"); 
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("setBloc : xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("setBloc : not the same dimension");
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
       
        throw std::invalid_argument("setBloc : Matrix not at the same place");
    }
}
void MatrixGPU::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixGPU* m, float factor)
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
        throw std::invalid_argument("setBloc Matrix not at the same place");
    }
}
void MatrixGPU::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, MatrixCPU* m)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        std::cout << _row << " " << _column << " " << iEnd << " " << jEnd << std::endl;
        throw std::out_of_range("setBloc: index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("setBloc : xBegin must be smaller than xEnd");
    } if (m->getNLin() != (iEnd - iBegin) || m->getNCol() != (jEnd - jBegin)) {
        throw std::invalid_argument("setBloc : not the same dimension");
    }
    
    if (!_GPU) {
        int row = 0;
        
        for (int i = iBegin; i < iEnd; i++) {
            int col = 0;
            for (int j = jBegin; j < jEnd;j++) {
               // if (iBegin == 44597 && row > 911 && col > 700) {  
                set(i, j, m->get(row, col));
                col = col + 1;
            }
          
            row = row + 1;
        }
    }
    else {
        throw std::domain_error("setBloc : Matrix on GPU");
    }
}


void MatrixGPU::setBloc(int iBegin, int iEnd, int jBegin, int jEnd, float value)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        std::cout << iBegin << " " << iEnd << " " << jBegin << " " << jEnd << " " << _row << " " << _column << std::endl;
        throw std::out_of_range(" setBloc : index out of bounds");
    } if ((iBegin >= iEnd) || (jBegin >= jEnd)) {
        throw std::invalid_argument("setBloc : xBegin must be smaller than xEnd");
    } 
    if (!_GPU) {
        int row = 0;

        for (int i = iBegin; i < iEnd; i++) {
            for (int j = jBegin; j < jEnd; j++) {
                set(i, j, value);
            }
        }
    }
    else if (getPos()) {
       
        SetBlocGPU << <_numBlocks, _blockSize >> > (_matrixGPU, value, iBegin, iEnd, jBegin, jEnd, _column);
    }
    else {

        throw std::invalid_argument("setBloc :Matrix not at the same place");
    }
}



void MatrixGPU::swap(MatrixGPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("swap : not the same dimension");
    }
    if (_GPU && m->getPos()) {
        float* temp = _matrixGPU;
        _matrixGPU = m->_matrixGPU;
        m->_matrixGPU = temp;

    }
    else if (!_GPU && !(m->getPos())) {
        float* temp = _matrixCPU;
        _matrixCPU = m->_matrixCPU;
        m->_matrixCPU = temp;
    }
    else {
        throw std::invalid_argument("swap : Matrix not at the same place");
    } 
}

void MatrixGPU::replace(float previous, float newValue)
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
void MatrixGPU::add(MatrixGPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("add : not the same dimension");
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
                float r = get(i, j) + m->get(i, j);
                this->set(i, j, r);
            }
        }
    } else {
        throw std::invalid_argument("add : Matrix not at the same place");
    } 
}

void MatrixGPU::addVector(MatrixGPU* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument("addVector : wrong dimension of the vector");
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
                    float r = get(i, j) + v->get(i, 0);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("addVector : Matrix not at the same place");
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
                    float r = get(i, j) + v->get(0, j);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("addVector : Matrix not at the same place");
        }
    }
}
void MatrixGPU::add(float c)
{
    if (_GPU) {
        addGPU<<<_numBlocks,_blockSize >>>(_matrixGPU,c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) + c;
                this->set(i, j, r);
            }
        }
    }
}

void MatrixGPU::add(MatrixGPU* m1, MatrixGPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("not the same dimension, fct add, m1 with m2");
    }
    if (!dim(m1)) {
        std::cout << _row << " " << m1->_row << " " << _column << " " << m1->_column << std::endl;
        throw std::invalid_argument("not the same dimension, fct add, this with m1");
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
                float r = m1->get(i, j) + m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("add Matrix not at the same place");
    }
    
}
void MatrixGPU::add(MatrixGPU* m, float c)
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
                float r = m->get(i, j) + c;
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("add Matrix not at the same place");
    }
    

}
void MatrixGPU::addTrans(MatrixGPU* m)
{
    MatrixGPU temp(*this);
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("addTrans not the same dimension (transpose)");
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
                float r = get(i, j) + m->get(j, i);
                temp.set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("addTrans Matrix not at the same place");
    }
    this->set(&temp);
    
}
///////////////////////////////////////////////////////////////////////////////
// subtraction
///////////////////////////////////////////////////////////////////////////////
void MatrixGPU::subtract(MatrixGPU* m1, MatrixGPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("subtract not the same dimension m1 with m2");
        
    }
    if (!dim(m1)) {
        throw std::invalid_argument("subtract not the same dimension m1 with this");
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
                float r = m1->get(i, j) - m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument(" subtract Matrix not at the same place");
    }
    
}
void MatrixGPU::subtract(MatrixGPU* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("subtract not the same dimension");
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
                float r = get(i, j) - m->get(i, j);
                set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("subtract Matrix not at the same place");
    }
   
    
}
void MatrixGPU::subtractVector(MatrixGPU* v)
{
    if (((v->getNCol() != 1) || (v->getNLin() != _row)) && ((v->getNLin() != 1) || (v->getNCol() != _column))) {
        throw std::invalid_argument( " subtractVector wrong dimension of the vector");
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
                    float r = get(i, j) - v->get(i, 0);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("subtractVector Matrix not at the same place");
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
                    float r = get(i, j) - v->get(0, j);
                    this->set(i, j, r);
                }
            }
        }
        else {
            throw std::invalid_argument("subtractVector Matrix not at the same place");
        }
    }

}
void MatrixGPU::subtractTrans(MatrixGPU* m)
{
    if (_row != m->getNCol() && _column != m->getNLin())
    {
        throw std::invalid_argument("subtractTrans not the same dimension (transpose)");
    }
    MatrixGPU temp(*this);
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
                float r = get(i, j) - m->get(j, i);
                temp.set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("subtractTrans Matrix not at the same place");
    }
    this->set(&temp);
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication
///////////////////////////////////////////////////////////////////////////////


void MatrixGPU::multiply(float c)
{
    if (_GPU) {
        multiplyGPU<<<_numBlocks, _blockSize >>> (_matrixGPU, c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) * c;
                this->set(i, j, r);
            }
        }
    }
        
}

void MatrixGPU::multiplyMat(MatrixGPU* A, MatrixGPU* B)
{
    if (A->getNLin() != getNLin()) {
        throw std::invalid_argument("multiplyMat result must be compatible with A (row)");
    }
    else if (A->getNCol() != B->getNLin()) {
        throw std::invalid_argument("multiplyMat A must be compatible with B (column with row)");
    }
    else if (getNCol() != B->getNCol()) {
        throw std::invalid_argument("multiplyMat result must be compatible with B (column)");

    }
    if (_GPU && A->getPos() && B->getPos()) { // solution temporaire
        transferCPU();
        A->transferCPU();
        B->transferCPU();
        float r = 0;
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
        float r = 0;
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
        throw std::invalid_argument("multiplyMat Matrix not at the same place");
    }
}

void MatrixGPU::multiply(MatrixGPU* Mat, MatrixGPU* vect, bool trans)
{
    // result = Mat*vect  nLine*nCol = nLine*Taille *Taille * Ncol 
    if (trans) {
        if (getNLin() != 1) {
            throw std::invalid_argument("multiply result must be a row vector ");
        }
        else if (getNCol() != Mat->getNLin()) {
            throw std::invalid_argument("multiply result must be compatible with Mat");
        }
        else if (vect->getNLin() != 1) {
            throw std::invalid_argument("multiply vect must be a row vector ");
        }
        else if (vect->getNCol() != Mat->getNCol()) {
            throw std::invalid_argument("multiply vect must be compatible with Mat");
        }
    }
    else {
        if (getNCol() != 1) {
            throw std::invalid_argument("multiply result must be a column vector ");
        }
        else if (getNLin() != Mat->getNLin()) {
            throw std::invalid_argument("multiply result must have the same row number as the Mat");
        }
        else if (vect->getNCol() != 1) {
            throw std::invalid_argument("multiply vect must be a column vector ");
        }
        else if (vect->getNLin() != Mat->getNCol()) {
            throw std::invalid_argument("multiply vect must be compatible with Mat");
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
                float sum = 0;
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
                float sum = 0;
                for (int j = 0; j < Mat->getNCol(); ++j)
                {
                    sum += Mat->get(i, j) * vect->get(j, 0);
                }
                set(i, 0, sum);
            }
        }
    }
    else {
        throw std::invalid_argument("multiply Matrix not at the same place");
    }
}

void MatrixGPU::MultiplyMatTransVec(MatrixGPU* MatToTrans, MatrixGPU* vect, bool rowVector)
{
    // methode tres peu efficace, acces memoire pas coalescent, mais on n'y peut rien...
    if (rowVector) {
        if (getNLin() != 1) {
            throw std::invalid_argument("MultiplyMatTransVec result must be a row vector ");
        }
        else if (getNCol() != MatToTrans->getNCol()) {
            throw std::invalid_argument("MultiplyMatTransVec result must be compatible with Mat");
        }
        else if (vect->getNLin() != 1) {
            throw std::invalid_argument("MultiplyMatTransVec vect must be a row vector ");
        }
        else if (vect->getNCol() != MatToTrans->getNLin()) {
            throw std::invalid_argument("MultiplyMatTransVec vect must be compatible with Mat");
        }
    }
    else {
        if (getNCol() != 1) {
            throw std::invalid_argument("MultiplyMatTransVec result must be a column vector ");
        }
        else if (getNLin() != MatToTrans->getNCol()) {
            throw std::invalid_argument("MultiplyMatTransVec result must have the same row number as the Mat");
        }
        else if (vect->getNCol() != 1) {
            throw std::invalid_argument("MultiplyMatTransVec vect must be a column vector ");
        }
        else if (vect->getNLin() != MatToTrans->getNLin()) {
            throw std::invalid_argument("MultiplyMatTransVec vect must be compatible with Mat");
        }
    }

    if (_GPU && MatToTrans->getPos() && vect->getPos())
    {
        int numBlock = MatToTrans->getNCol();
        switch (_blockSize) {
        case 512:
            multiplyGPUMatVectTrans<512> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case 256:
            multiplyGPUMatVectTrans<256> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case 128:
            multiplyGPUMatVectTrans<128> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case 64:
            multiplyGPUMatVectTrans< 64> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case 32:
            multiplyGPUMatVectTrans< 32> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case 16:
            multiplyGPUMatVectTrans< 16> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case  8:
            multiplyGPUMatVectTrans<  8> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case  4:
            multiplyGPUMatVectTrans<  4> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case  2:
            multiplyGPUMatVectTrans<  2> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        case  1:
            multiplyGPUMatVectTrans<  1> << <numBlock, _blockSize >> > (_matrixGPU, MatToTrans->_matrixGPU, vect->_matrixGPU, MatToTrans->getNCol());
            break;
        }

    }
    else if (!_GPU && !(MatToTrans->getPos()) && !vect->getPos())
    {
        if (rowVector) {
            for (int i = 0; i < MatToTrans->getNLin(); ++i)
            {
                float sum = 0;
                for (int j = 0; j < MatToTrans->getNCol(); ++j)
                {
                    sum += MatToTrans->get(j, i) * vect->get(0, j);
                }
                set(0, i, sum);
            }
        }
        else {
            for (int i = 0; i < _row; ++i)
            {
                float sum = 0;
                for (int j = 0; j < MatToTrans->getNCol(); ++j)
                {
                    sum += MatToTrans->get(j, i) * vect->get(j, 0);
                }
                set(i, 0, sum);
            }
        }
    }
    else {
        throw std::invalid_argument("MultiplyMatTransVec Matrix not at the same place");
    }

}

void MatrixGPU::linearOperation(MatrixGPU* A, MatrixGPU* x, MatrixGPU* b, bool trans)
{
    if (trans) {
        if (getNLin() != 1) {
            throw std::invalid_argument("linearOperation result must be a row vector ");
        }
        else if (getNCol() != A->getNLin()) {
            throw std::invalid_argument("linearOperation result must be compatible with A");
        }
        else if (x->getNLin() != 1 || b->getNLin() != 1) {
            throw std::invalid_argument("linearOperation x and b must be a row vector ");
        }
        else if (x->getNCol() != A->getNCol()) {
            throw std::invalid_argument("linearOperation x must be compatible with A");
        }
    }
    else {
        if (getNCol() != 1) {
            throw std::invalid_argument("linearOperation result must be a column vector ");
        }
        else if (getNLin() != A->getNLin()) {
            throw std::invalid_argument("linearOperation result must have the same row number as A");
        }
        else if (x->getNCol() != 1 || b->getNCol() != 1) {
            throw std::invalid_argument("linearOperation x and b must be a column vector ");
        }
        else if (x->getNLin() != A->getNCol()) {
            throw std::invalid_argument("linearOperation x must be compatible with Mat");
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
                float sum = 0;
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
                float sum = 0;
                for (int j = 0; j < A->getNCol(); ++j)
                {
                    sum += A->get(i, j) * x->get(j, 0);
                }
                set(i, 0, sum + b->get(i,0));
            }
        }
    }
    else {
        throw std::invalid_argument("linearOperation Matrix not at the same place");
    }
}

///////////////////////////////////////////////////////////////////////////////
// Multiplication Terme � Terme
///////////////////////////////////////////////////////////////////////////////

void MatrixGPU::multiplyT(MatrixGPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("multiplyT not the same dimension");
    }
   
    if (_GPU && m->getPos())
    {
        if (_N > 0) {
            multiplyTGPU << <_numBlocks, _blockSize >> > (_matrixGPU, m->_matrixGPU, _N);
        }
    }
    else if (!_GPU && !(m->getPos()))
    {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) * m->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("multiplyT Matrix not at the same place");
    }
}

void MatrixGPU::multiplyT(MatrixGPU* m1, MatrixGPU* m2)
{
    if (!m1->dim(m2)) {
        throw std::invalid_argument("multiplyT not the same dimension");
    }
    if (!dim(m1)) {
        throw std::invalid_argument("multiplyT not the same dimension");
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
                float r = m1->get(i, j) * m2->get(i, j);
                this->set(i, j, r);
            }
        }
    }
    else {
        throw std::invalid_argument("multiplyT Matrix not at the same place");
    }
}


void MatrixGPU::divide(float c)
{
    if (c == 0) {
        throw std::domain_error("divide : divide by 0");
    }
    if (_GPU) {
        divideGPU <<<_numBlocks, _blockSize >>> (_matrixGPU, c, _N);
    }
    else {
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                float r = get(i, j) / c;
                this->set(i, j, r);
            }
        }
    }
    
}

void MatrixGPU::divideT(MatrixGPU* m)
{
    
    if (!dim(m)) {
        throw std::invalid_argument("divideT not the same dimension");
    }
    if (_GPU && m->getPos())
    {
        divideGPU<<<_numBlocks, _blockSize >>>(_matrixGPU, m->_matrixGPU, _N);
    }
    else if (!_GPU && !(m->getPos()))
    {
        MatrixGPU temp(*this);
        float r = 0;
        float f = 0;
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                f = m->get(i, j);
                if (f == 0) {
                    throw std::domain_error("divideT divide by 0");
                }
                r = get(i, j) / f;
                temp.set(i, j, r);
            }
        }
        set(&temp);
    }
    else {
        throw std::invalid_argument("divideT Matrix not at the same place");
    }
    
}

void MatrixGPU::invertGaussJordan(MatrixGPU* mToInvert)
{
    
    if (!dim(mToInvert)) {
        throw std::invalid_argument("invertGaussJordan not the same dimension");
    }
    if (_row != _column) {
        throw std::invalid_argument("invertGaussJordan must be a square matrix");
    }
   
    if (!_GPU && !mToInvert->getPos()) {
        MatrixCPU m;
        mToInvert->toMatCPU(m);
       
        MatrixCPU augmented(_row, _column);
        augmented.setEyes(1);
        MatrixCPU indices(1, 2);
        int r = 0;
        for (int column = 0; column < _column; column++) {
            
            float pivotabs = m.maxAbs(r, _row, column, column + 1, &indices);
            int k = indices.get(0, 0); // indice max de la colonne j
            float pivot = m.get(k, column);
            if (pivotabs < 0.000001f) {
                throw std::invalid_argument("invertGaussJordan not invertible matrix");
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
        MatrixGPU m(*mToInvert);
        setEyes(1);
        dim3 threadsPerBlock(32, 32);
        int r = 0;
        MatrixGPU matCol(_row, 1, 0, 1);
        int k = 0;
        
        for (int column = 0; column < _column; column++) {
           
            //transferCPU();
            ///display();
            //transferGPU();

            m.getCol(&matCol, column, r);
            float pivotAbs = matCol.max2(&k); // comme matCol est un vecteur colonne, la position du maximum correspond directement à la ligne !
            float pivot = matCol.get(k, 0, false);
            
            
            if (pivotAbs < 0.000001f) {
                std::cout << "not invertible " << column <<" " << pivotAbs << std::endl;
                if (_N < 100) {
                    mToInvert->display(true);
                }
                throw std::invalid_argument("invertGaussJordan not invertible matrix");
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
        throw std::invalid_argument("invertGaussJordan Matrix not at the same place");
    }


}

void MatrixGPU::LUPFactorization(MatrixGPU* A, MatrixGPU* P)
{
    float Tol = 0.0000001;
    int n = getNLin();
    A->set(this);

    // code from wikipedia adapted
    if (getNCol() != getNLin()) {
        throw std::invalid_argument("LUPFactorization A must be square");
    }
    if (P->getNCol() != 1 || P->getNLin() != (getNCol() + 1)) {
        throw std::invalid_argument("LUPFactorization wrong size of P");
    }

    if (!_GPU && !A->getPos() && !P->getPos()) {
        for (int i = 0; i < n; i++) {
            P->set(i, 0, i); //Unit permutation matrix, P[N] initialized with N
        }

        float absA = 0;
        int j = 0;
        for (int col = 0; col < n; col++) {
            float maxA = 0.0;
            int imax = col;
            for (int k = col; k < n; k++) {
                absA = fabs(A->get(k, col));
                if (absA > maxA)
                {
                    maxA = absA;
                    imax = k;
                }
            }
            //std::cout << "max de " << maxA << "en position " << imax << std::endl;
            if (maxA < Tol) {
                //std::cout << "failure, matrix is degenerate" << std::endl;
                throw std::invalid_argument("LUPFactorization matrix is degenerate");
                return; //failure, matrix is degenerate
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

        MatrixGPU matCol(_row, 1, 0, 1);
        initPermMatr <<<_numBlocks, _blockSize >> > (P->_matrixGPU, n);
       
        int k = 0;
        for (int col = 0; col < n; col++) {

            A->getCol(&matCol, col, col);
            float pivotAbs = matCol.max2(&k); // comme matCol est un vecteur colonne, la position du maximum correspond directement à la ligne !
            float pivot = matCol.get(k, 0, false);

            //std::cout << "max de " << pivot << "en position " << k << std::endl;

            if (pivotAbs < Tol) {
                //std::cout << "failure, matrix is degenerate" << std::endl;
             
                throw std::invalid_argument("LUPFactorization matrix is degenerate");
            }
            else {
                if (k != col) { //le max pas sur la diagonal
                    //pivoting P and counting pivots starting from N (for determinant)
                    updatePermMatr <<<1, 1 >>> (P->_matrixGPU, k, col, n);
                 
                    //pivoting rows of A
                    A->swapLine(col, k);
                }
                
                updateLUPFactorization << <n, _blockSize >> > (A->_matrixGPU, col, n);
            }
        }
    }
    else {
        throw std::invalid_argument("LUPFactorization Matrix not at the same place");
    }
    
    // en vrai on peut tout stocker dans une matrice comme on sait que diag(L) = Id, Et donc on peut avoir A = (L-Id) + U


}

void MatrixGPU::solveSysUpper(MatrixGPU* U)
{
    if (getNLin() != U->getNCol() || U->getNLin() != getNLin()) {
        throw std::invalid_argument("solveSysUpper A must be square");
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
        throw std::invalid_argument("solveSysUpper Matrix not at the same place");
    }
}

void MatrixGPU::solveSysLower(MatrixGPU* L, MatrixGPU* b, MatrixGPU* P) // element diag equal to 1
{
    if (getNLin() != L->getNCol() || L->getNLin() != b->getNLin()) {
        throw std::invalid_argument("solveSysLower A must be square");
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
        throw std::invalid_argument("solveSysLower Matrix not at the same place");
    }
    
}

void MatrixGPU::solveSys(MatrixGPU* A, MatrixGPU* P, MatrixGPU* b)
{
    if (A->getNCol() != A->getNLin() || A->getNLin() != b->getNLin()) {
        throw std::invalid_argument("solveSys wrong size of A");
    }
    if (b->getNCol() != 1) {
        throw std::invalid_argument("solveSys b must be a column vector");
    }
    if (P->getNLin() != (A->getNCol() + 1) || P->getNCol() != 1) {
        throw std::invalid_argument("solveSys wrong size of P");
    }

    if (!_GPU && !A->getPos() && !P->getPos() && !b->getPos()) {
        int n = getNLin();
        for (int i = 0; i < n; i++) {
            set(i, 0, b->get(P->get(i, 0), 0)); // x[i] = b[P[i]];

            for (int k = 0; k < i; k++) {
                set(i, 0, get(i, 0) - A->get(i, k) * get(k, 0));
            }
        }
        for (int i = n - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < n; k++) {
                set(i, 0, get(i, 0) - A->get(i, k) * get(k, 0));
               
            }

            set(i, 0, get(i, 0) / A->get(i, i));
        }
    }
    else if (_GPU && A->getPos() && P->getPos() && b->getPos()) {
        int n = getNLin();
        setPermute << <_numBlocks, _blockSize >> > (_matrixGPU, b->_matrixGPU, P->_matrixGPU, n);
        
        //CHECK_CUDA_ERROR_MAT(cudaFuncSetAttribute(solveSysGPU, cudaFuncAttributeMaxDynamicSharedMemorySize, n * sizeof(float)));

        solveSysGPU<<<1, _blockSize, n*sizeof(float)>>> (A->_matrixGPU, _matrixGPU, n);
        
        
        ////CHECK_LAST_CUDA_ERROR_MAT();

        //CHECK_CUDA_ERROR_MAT(cudaDeviceSynchronize());/**/
        /*for (int i = 0; i < n; i++) {
            solveLowSys <<< 1, _blockSize >>> (A->_matrixGPU, _matrixGPU, i, n);
        }
       
        for (int i = n - 1; i >= 0; i--)
        {
            solveUpSys << < 1, _blockSize >> > (A->_matrixGPU, _matrixGPU, i, n);
        }*/
    }
    else {
        throw std::invalid_argument("solveSys Matrix not at the same place");
    }
    
   
}

///////////////////////////////////////////////////////////////////////////////
// Fonction autres
///////////////////////////////////////////////////////////////////////////////

float MatrixGPU::max2() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    if (!_GPU) {
        float M = fabs(get(0, 0));
        float m = 0;
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
        float odata = 0;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            std::cout << "allocation !!!" << std::endl;
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
             cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
             //cudaMemcpy(_preallocationFloat, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
             return sqrt(odata);//sqrt(*_preallocationFloat);
        }
        else
        {
            cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
            std::cout << "free !!!" << std::endl;
            cudaFree(d_odata);
            return sqrt(odata);
        }
    }
}

float MatrixGPU::max2(int* indice)
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("max2 Empty Matrix");
    }
    if (!_GPU) {
        float M = fabs(get(0, 0));
        float m = 0;
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
        float odata = 0;
        float* d_odata;
        int* d_pos;
        int pos = 0;
        cudaMalloc((void**)&d_pos, sizeof(int) * numBlocks);
        if (preallocation) {

            d_odata = _preallocation;
        }
        else {

            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
        cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&pos, d_pos, sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(d_pos);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        *indice = pos;
        return sqrt(odata);
    }
}
float MatrixGPU::max2(MatrixGPU* m) const
{
    if (_row == 0 || _column == 0) {
        throw std::out_of_range("max2 Empty Matrix");
    }
    if (!_GPU && !(m->getPos())) 
    {
        float M = 0;
        float f = 0;
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
        float odata;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            std::cout << "allocation !!!" << std::endl;
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
            cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
            //cudaMemcpy(_preallocationFloat, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
            return  sqrt(odata);//sqrt(*_preallocationFloat);
        }
        else
        {
            cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
            std::cout << "free !!!" << std::endl;
            cudaFree(d_odata);
            return sqrt(odata);
        }
    }
    else {
        throw std::invalid_argument("max2 Matrix not at the same place");
    }
}

float MatrixGPU::distance2(MatrixGPU* m)
{
    if (!dim(m)) {
        throw std::invalid_argument("distance2 not the same size");
    }
    if (_GPU && m->getPos())
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        float odata = 0;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
        cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }


        return sqrtf(odata);
    }
    else if (!_GPU && !(m->getPos()))
    {
        double d = 0;
        double r = 0;
        for (int i = 0; i < _row;++i)
        {
            for (int j = 0; j < _column;++j)
            {
                r = get(i, j) - m->get(i, j);
                d = d + r * r;
            }
        }
        return sqrtf(d);
    }
    else {
        throw std::invalid_argument("distance2 Matrix not at the same place");
    }
}

void MatrixGPU::Moy(MatrixGPU* m, MatrixGPU* nb, int sens)
{
    float s;
    int n;
    if (sens) { // on travaille sur les colonnes
        if ((_row != 1) || (_column != m->getNCol()) || (_column != nb->getNCol()) || (nb->getNLin() != 1))
        {
            throw std::invalid_argument("Moy wrong dimension of the vector");
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
            throw std::invalid_argument("Moy Matrix not at the same place");
        }
        

    }
    else { // on travaille sur les lignes 
        if ((_column != 1) || (_row != m->getNLin()) || (_row != nb->getNLin()) || (nb->getNCol() != 1)) {
            throw std::invalid_argument("Moy wrong dimension of the vector");
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
            throw std::invalid_argument("Moy Matrix not at the same place");
        }

    }
}

void MatrixGPU::project(MatrixGPU* Lb, MatrixGPU* Ub)
{
    if (!dim(Lb) || !dim(Ub)) {
        throw std::invalid_argument("project not the same dimension");
    }
    if (_GPU && Lb->getPos() && Ub->getPos())
    {
        projectGPU<<<_numBlocks, _blockSize >>>(_matrixGPU, Lb->_matrixGPU, Ub->_matrixGPU, _N);
    }
    else if (!_GPU && !(Lb->getPos()) && !(Ub->getPos()))
    {
        float ub = 0;
        float lb = 0;
        float r = 0;
        MatrixGPU temp(*this);
        for (int i = 0;i < _row;++i)
        {
            for (int j = 0;j < _column;++j)
            {
                r = get(i, j);
                ub = Ub->get(i, j);
                lb = Lb->get(i, j);
                if (ub < lb) {
                    throw std::invalid_argument("project impossible to have a value for the projection, ub>lb");
                }
                r = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r; // permet de ne pas faire de branchement if.
                temp.set(i, j, r);
            }
        }
        this->set(&temp);
    }
    else {
        throw std::invalid_argument("project Matrix not at the same place");
    }
    
}
void MatrixGPU::projectNeg()
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
                float r = get(i, j);
                r = (r > 0) * r;
                set(i, j, r);
            }
        }
    }

}
void MatrixGPU::projectPos()
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
                float r = get(i, j);
                r = (r > 0) * r;
                set(i, j, r);
            }
        }
    }
}



float MatrixGPU::sum() const
{
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }
    if (_GPU) 
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        float odata = 0;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
        cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        //std::cout << "sum " << odata << " " <<_blockSize << " " << numBlocks << std::endl;
        return odata;
    }
    else if (!_GPU)
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
}

float MatrixGPU::sum(int begin, int end)
{
    if (begin < 0 || end < 0) {
        throw std::invalid_argument("sum indice must be positve");
    }
    if (begin > end ) {
        throw std::invalid_argument("sum begin must be smaller than end");
    }
    if (begin > _N || end > _N) {
        throw std::out_of_range("sum indice must smaller than N");
    }
    if (_row == 0 || _column == 0) {
        return 0;
        //throw std::out_of_range("Empty Matrix");
    }

    if (_GPU)
    {
        int numBlocks = _numBlocks;
        float odata = 0;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
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
        cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        //std::cout << "sum " << odata << " " <<_blockSize << " " << numBlocks << std::endl;
        return odata;
    }
    else if (!_GPU)
    {
        float d = 0;
        float r = 0;
        for (int elem = begin; elem < end; ++elem)
        {
                r = _matrixCPU[elem];
                d = d + r;
        }
        return d;
    }
}

void MatrixGPU::sum(MatrixGPU* m)
{
    float s = 0;
     // on travaille sur les lignes 
    if ((_column != 1) || (_row != m->getNLin())) {
        throw std::invalid_argument("sum wrong dimension of the column vector ");
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
        throw std::invalid_argument("sum Matrix not at the same place");
    }
}

float MatrixGPU::distance2() {

    if (_GPU ) //&& m->getPos())
    {
        int numBlocks = _numBlocks;
        unsigned int n = _N;
        float* d_odata;
        if (preallocation) {
            d_odata = _preallocation;
        }
        else {
            cudaMalloc((void**)&d_odata, sizeof(float) * numBlocks);
        }
        float odata = 0;
        
        
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
        cudaMemcpy(&odata, d_odata, sizeof(float), cudaMemcpyDeviceToHost);
        if (!preallocation) {
            cudaFree(d_odata);
        }
        return sqrtf(odata);
    }
    else if (!_GPU)// && !(m->getPos()))
    {
        float d = 0;
        float r = 0;
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
// Display MatrixGPU contents
///////////////////////////////////////////////////////////////////////////////
void MatrixGPU::display(bool force) 
{   
    bool transfert = false;
    if (this) {
        if (_GPU && !force ) {
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
            if (_column == 1) {
                std::cout << " transpose  : ";
                for (int i = 0;i < _row;++i)
                {
                    for (int j = 0;j < _column;++j)
                    {
                        float value = get(i, j);
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
                        float value = get(i, j);
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

void MatrixGPU::displayBloc(int iBegin, int iEnd, int jBegin, int jEnd, bool force)
{
    if ((iBegin < 0) || (jBegin < 0) || iEnd > _row || jEnd > _column) {
        throw std::out_of_range("displayBloc index out of bounds");
    } if ((iBegin > iEnd) || (jBegin > jEnd)) {
        throw std::invalid_argument("displayBloc xBegin must be smaller than xEnd");
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
                    float value = get(i, jBegin);
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
                        float value = get(i, j);
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

void MatrixGPU::swapLine(int line1, int line2)
{
    if (_GPU) {
        swapLineGJ << <_numBlocks, _blockSize >> > (_matrixGPU, line1, line2, _column);// swap des lignes
    }
    else {
        float temp = 0;
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
MatrixGPU::~MatrixGPU()
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
    if (_matrixGPU != nullptr) {
        cudaFree(_matrixGPU);
        _matrixGPU = nullptr;
    }
     DELETEA(_matrixCPU);
    
}



void MatrixGPU::saveCSV(const std::string& filename, std::ios_base::openmode mode, int trans) const
{
    if (_GPU) {
        throw std::domain_error("saveCSV : Matrix on GPU");
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

void MatrixGPU::saveCSVForce(const std::string& filename, std::ios_base::openmode mode, int trans)
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


__global__ void setup_kernel(curandState* state) {

    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(1234, idx, 0, &state[idx]);
}


__global__ void generate_kernel(curandState* my_curandstate, float* result, float eps, const unsigned int N) {

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        result[i] = (2*curand_uniform(my_curandstate + i)-1) * eps;
    }
}





__global__ void setGPU(float* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i];
    }
}
__global__ void setGPUFD(float* mat1, double* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i];
    }
}
__global__ void setGPUDF(double* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = (double) mat2[i];
    }
}


__global__ void setGPU(float* mat1, const float value, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = value;
    }
}
__global__ void setGPUunique(float* mat1, const float value, int pos) {
    int index = threadIdx.x;
    if (index == 0) {
        mat1[pos] = value;
    }

}

__global__ void setTransGPU(float* mat1, float* matToTrans, const int column, const int row) {
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

__global__ void setColGPU(float* mat1, float* mat2, const int numCol, const int column, const int row, const int offset) {
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < row; i += step)
    {
        mat1[i] = i < offset ? 0 : mat2[i*column + numCol];
    }

}

__global__ void setEyesGPU(float* mat2, const float value, const int col, const int row) 
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
__global__ void setEyesGPU(float* mat2, float* mat1, const int col, const int row)
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


__global__ void SetBlocGPU(float* out, float* in, int ibegin, int iend, int jbegin, int jend, int col)
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

__global__ void SetBlocGPU(float* out, float value, int ibegin, int iend, int jbegin, int jend, int col)
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
        out[GlobalInd] = value;
    }
}

__global__ void SetBlocGPU(float* out, float* in, int ibegin, int iend, int jbegin, int jend, int col, float factor)
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

/*__global__ void SetBlocGPU(float* out, float* in, int ibegin, int iend, int jbegin, int jend, int col, float factor) // fait que la première ligne
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




__global__ void replaceGPU(float* mat,const float previous, const float newValue,const int N) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = (mat[i] == previous) * (newValue-mat[i]) + mat[i];
    }
}




__global__ void addGPU(float* mat, float c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] + c;
    }
}
__global__ void addGPU(float* mat1, float* mat2, float c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] + c;
    }
}
__global__ void addGPU(float* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] + mat2[i];
    }
}
__global__ void addGPU(float* mat1, float* mat2, float* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] + mat3[i];
    }
}

__global__ void addVectorGPU1(float* mat1, float* vect, const int n, int N) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i / n; // division entière
        mat1[i] = mat1[i] + vect[k];
    }

}
__global__ void addVectorGPU2(float* mat1, float* vect, const int n, int N) // vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < N; i += step)
    {
        int k = i % n; // modulo
        mat1[i] = mat1[i] + vect[k];
    }


}

__global__ void addTransGPU(float* out, float* mat1, float* mat2, const int col, const int line, int N) 
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

__global__ void substractGPU(float* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] - mat2[i];
    }
}
__global__ void substractGPU(float* mat1, float* mat2, float* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] - mat3[i];
    }
}

__global__ void substractVectorGPU1(float* mat1, float* vect, const int n, int N) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i / n; // division entière
        mat1[i] = mat1[i] - vect[k];
    }

}
__global__ void substractVectorGPU2(float* mat1, float* vect, const int n, int N) // vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        int k = i % n; // modulo
        mat1[i] = mat1[i] - vect[k];
    }

}

__global__ void substractTransGPU(float* out, float* mat1, float* mat2, const int col, const int line, int N)
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

__global__ void multiplyGPU(float* mat, const float c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] * c;
    }
}

__global__ void multiplyTGPU(float* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] * mat2[i];
    }
}
__global__ void multiplyTGPU(float* mat1, float* mat2, float* mat3, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat2[i] * mat3[i];
    }
}

__global__ void divideGPU(float* mat, const float c, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat[i] = mat[i] / c;
    }
}
__global__ void divideGPU(float* mat1, float* mat2, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        mat1[i] = mat1[i] / mat2[i];
    }
}

__global__ void moyGPU1(float* res, float* mat1, float* nb, const int line, const int column) //vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
     
    for (int i = index; i < column; i += step)
    {
        float s = 0.0;
        for (int j = 0; j < line; j++)
        {
            s += mat1[i + column *j];
        } 
        res[i] = s / nb[i];
    }

}
__global__ void moyGPU2(float* res, float* mat1, float* nb, const int line, const int column) // vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    
    for (int i = index; i < line; i += step)
    {
        float s = 0.0;
        for (int j = 0; j < column; j++)
        {
            s +=  mat1[i*column + j];
        }
        res[i] = s /nb[i];
    }
}

__global__ void projectGPU(float* mat, float* Lb, float* Ub, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        float r = mat[i];
        float ub = Ub[i];
        float lb = Lb[i];
        r = (ub - r) * (r > ub) + (lb - r) * (r < lb) + r;
        mat[i] = r;//(Ub[i] - mat[i])* (mat[i] > Ub[i]) + (Lb[i] - mat[i]) * (mat[i] < Lb[i]) + mat[i];
    }
}

__global__ void projectGPUPos(float* mat, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        float r = mat[i];
        mat[i] = (r > 0) * r;
    }
}

__global__ void projectGPUNeg(float* mat, int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += step)
    {
        float r = mat[i];
        mat[i] = (r < 0) * r;
    }
}


__global__ void sumGPU(float* res, float* mat1, const int line, const int column) //vecteur colonne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < line; i += step)
    {
        float s = 0.0;
        for (int j = 0; j < column; j++)
        {
            s += mat1[i*column + j];
        }
        res[i] = s;
    }
}

__global__ void sumGPU2(float* res, float* mat1, const int line) //vecteur ligne
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index==0)
    {
        float s = 0.0;
        for (int j = 0; j < line; j++)
        {
            s += mat1[j];
        }
        
        *res = s ;
    }
}



__device__ int sumCommSingleWarp(volatile float* shArr) {
    int idx = threadIdx.x % warpSize; //the lane index in the warp
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
__global__ void sumMonoBlock(float* g_idata, float* g_odata, unsigned int n) {
    
    int thIdx = threadIdx.x;
    float sum = 0;
    for (int i = thIdx; i < n; i += blockSize)
        sum += g_idata[i];
    __shared__ float r[blockSize];
    r[thIdx] = sum;
    __syncthreads();
    if (blockSize >= 512) {
        if (thIdx < 256) {
            r[thIdx] += r[thIdx + 256];
        }
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (thIdx < 128) {
            r[thIdx] += r[thIdx + 128];
        }
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (thIdx < 64) {
            r[thIdx] += r[thIdx + 64];
        }
        __syncthreads();
    }
    //if (blockSize >= 64) {
    //    if (thIdx < 32) {
            warpReduce<blockSize>(r, thIdx);
    //    }
    // }
    //else if (blockSize >= 32) { // cas blockSize = 32
       // warpReduce<blockSize>(r, thIdx);
   // }
    __syncthreads;
    if (thIdx == 0) {
         *g_odata = r[0];
    }
       
    
}


template <unsigned int blockSize>
__global__ void SumMultiBlock(float* g_idata, float* g_odata, unsigned int n) {

    __shared__ float shArr[blockSize];
    int thIdx = threadIdx.x;

    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float sum = 0;
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
__global__ void SumMultiBlock(float* g_idata, float* g_odata, unsigned int begin, unsigned int end) {
    __shared__ float shArr[blockSize];
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float sum = 0;

    for (int i = gthIdx + begin; i < end; i += gridSize) {
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
__global__ void SumEachRow(float* g_idata, float* g_odata, const int nCol) {
    __shared__ float shArr[blockSize];
    int thIdx = threadIdx.x;
    int row = blockIdx.x;
    int idBegin = thIdx + row * nCol;
    int idEnd = (row + 1) * nCol;
    int step = blockDim.x;

    float sum = 0;
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
__global__ void distanceMultiBlock(float* g_idata, float* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float sum = 0;
    for (int i = gthIdx; i < n; i += gridSize)
        sum += (g_idata[i] * g_idata[i]);
    __shared__ float shArr[blockSize];
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
__global__ void distanceMultiBlock(float* g_idata, float* g_idata2, float* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float sum = 0;
    for (int i = gthIdx; i < n; i += gridSize)
        sum += ((g_idata[i]- g_idata2[i]) * (g_idata[i] - g_idata2[i]));
    __shared__ float shArr[blockSize];
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
__device__ void warpReduceMaxPos(volatile float* r, volatile int* pos){
    int idx = threadIdx.x % warpSize; //the lane index in the warp

    if (idx < 32 && blockSize >= 64) {
        pos[idx] = r[idx + 32] > r[idx] ? pos[idx + 32] : pos[idx];
        r[idx] = r[idx + 32] > r[idx] ? r[idx + 32] : r[idx];//r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] * (r[idx + 16] > r[idx]) + r[idx] * (r[idx] <= r[idx + 16]);
    }
    __syncwarp();
    if (idx < 16 && blockSize >= 32) {
        pos[idx] = r[idx + 16] > r[idx] ? pos[idx + 16] : pos[idx];
        r[idx] = r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] > r[idx] ? r[idx + 16] : r[idx];//r[idx + 16] * (r[idx + 16] > r[idx]) + r[idx] * (r[idx] <= r[idx + 16]);
    }
    __syncwarp();
    if (idx <  8 && blockSize >= 16) {
        pos[idx] = r[idx + 8] > r[idx] ? pos[idx + 8] : pos[idx];
        r[idx] = r[idx + 8] > r[idx] ? r[idx + 8] : r[idx];//r[idx +  8] > r[idx] ? r[idx +  8] : r[idx];//r[idx +  8] * (r[idx +  8] > r[idx]) + r[idx] * (r[idx] <= r[idx +  8]);
    }
    __syncwarp();
    if (idx < 4 && blockSize >=  8) {
        pos[idx] = r[idx + 4] > r[idx] ? pos[idx + 4] : pos[idx];
        r[idx] = r[idx + 4] > r[idx] ? r[idx + 4] : r[idx];//r[idx +  4] > r[idx] ? r[idx +  4] : r[idx];//r[idx +  4] * (r[idx +  4] > r[idx]) + r[idx] * (r[idx] <= r[idx +  4]);
    }
    __syncwarp();
    if (idx < 2 && blockSize >=  4) {
        pos[idx] = r[idx + 2] > r[idx] ? pos[idx + 2] : pos[idx];
        r[idx] = r[idx + 2] > r[idx] ? r[idx + 2] : r[idx];//r[idx +  2] > r[idx] ? r[idx +  2] : r[idx];//r[idx +  2] * (r[idx +  2] > r[idx]) + r[idx] * (r[idx] <= r[idx +  2]);
    }
    __syncwarp();
    if (idx < 1 && blockSize >=  2) {
        pos[idx] = r[idx + 1] > r[idx] ? pos[idx + 1] : pos[idx];
        r[idx] = r[idx + 1] > r[idx] ? r[idx + 1] : r[idx];//r[idx +  1] > r[idx] ? r[idx +  1] : r[idx];//r[idx +  1] * (r[idx +  1] > r[idx]) + r[idx] * (r[idx] <= r[idx +  1]);
    }
}


template <unsigned int blockSize>
__global__ void maxMonoBlock(float* g_idata, float* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    float max = 0;

    for (int i = thIdx; i < n; i += blockSize) {
        float s = g_idata[i];
        max = s > max ? s : max;// s>max ? s:max;//s * (s > max) + max * (max <= s);
    }
    __shared__ float shArr[blockSize];
    shArr[thIdx] = max;
    __syncthreads();
    if (blockSize >= 512) {
        if (thIdx < 256) {
            shArr[thIdx] = shArr[thIdx + 256] > shArr[thIdx] ? shArr[thIdx + 256] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 

        }
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (thIdx < 128) {
            shArr[thIdx] = shArr[thIdx + 128] > shArr[thIdx] ? shArr[thIdx + 128] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 

        }
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (thIdx < 64) {
            shArr[thIdx] = shArr[thIdx + 64] > shArr[thIdx] ? shArr[thIdx + 64] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 

        }
        __syncthreads();
    }
    if (blockSize >= 64) {
        if (thIdx < 32) {
            warpReduceMax<blockSize>(shArr, thIdx);
        }
    }
    if (blockSize >= 32 && blockSize < 64) {
       warpReduceMax<blockSize>(shArr, thIdx);
    }
    __syncthreads;

    if (thIdx == 0)
        *g_odata = shArr[0];
    
}

template <unsigned int blockSize>
__global__ void maxMultiBlock(float* g_idata, float* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float max = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        float s = (g_idata[i] * g_idata[i]);
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ float shArr[blockSize];
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
__global__ void maxMultiBlock(float* g_idata, float* g_odata, unsigned int n, int* pos) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float max = 0;
    int indice = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        float s = (g_idata[i] * g_idata[i]);
        indice = s > max ? i : indice;
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ float shArr[blockSize];
    __shared__ float shPos[blockSize];
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
__global__ void maxMonoBlock(float* g_idata, float* g_odata, unsigned int n, int* pos) {
    int thIdx = threadIdx.x;
    float max = 0;
    int indice = 0;
    for (int i = thIdx; i < n; i += blockSize) {
        float s = g_idata[i];
        indice = s > max ? pos[i] : indice;
        max = s > max ? s : max;// s>max ? s:max;//s * (s > max) + max * (max <= s);
    }
    __shared__ float shArr[blockSize];
    __shared__ int shPos[blockSize];
    shArr[thIdx] = max;
    shPos[thIdx] = indice;

    __syncthreads();
    if (blockSize >= 512) {
        if (thIdx < 256) {
            shPos[thIdx] = shArr[thIdx + 256] > shArr[thIdx] ? shPos[thIdx + 256] : shPos[thIdx];
            shArr[thIdx] = shArr[thIdx + 256] > shArr[thIdx] ? shArr[thIdx + 256] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 
        }
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (thIdx < 128) {
            shPos[thIdx] = shArr[thIdx + 128] > shArr[thIdx] ? shPos[thIdx + 128] : shPos[thIdx];
            shArr[thIdx] = shArr[thIdx + 128] > shArr[thIdx] ? shArr[thIdx + 128] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 

        }
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (thIdx < 64) {
            shPos[thIdx] = shArr[thIdx + 64] > shArr[thIdx] ? shPos[thIdx + 64] : shPos[thIdx];
            shArr[thIdx] = shArr[thIdx + 64] > shArr[thIdx] ? shArr[thIdx + 64] : shArr[thIdx]; //shArr[thIdx + size] > shArr[thIdx] ? shArr[thIdx + size] : shArr[thIdx];//shArr[thIdx + size] * (shArr[thIdx + size] > shArr[thIdx]) + shArr[thIdx] * (shArr[thIdx] <= shArr[thIdx + size]); 
        }
        __syncthreads();
    }
    if (blockSize >= 64) {
        if (thIdx < 32) {
            warpReduceMaxPos<blockSize>(shArr, shPos);
        }
    }
    if (blockSize >= 32 && blockSize < 64) {
        warpReduceMaxPos<blockSize>(shArr, shPos);
    }
    __syncthreads();
   
    if (thIdx == 0) {
        *g_odata = shArr[0];
        *pos = shPos[0];
    }
       
}


template <unsigned int blockSize>
__global__ void maxMultiBlock(float* g_idata, float* g_idata2, float* g_odata, unsigned int n) {
    int thIdx = threadIdx.x;
    int gthIdx = thIdx + blockIdx.x * blockSize;
    const int gridSize = blockSize * gridDim.x;
    float max = 0;
    for (int i = gthIdx; i < n; i += gridSize) {
        float s = (g_idata[i] - g_idata2[i]);
        s = s * s;
        max = s > max ? s : max;//s > max ? s : max; //s * (s > max) + max * (max <= s);
    }
    __shared__ float shArr[blockSize];
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


__global__ void normalisationGJ(float* mat, const int row, const int nCol, const float factor) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < nCol; i += step)
    {
        mat[i + row * nCol] = mat[i + row * nCol] / factor;
    }


}

__global__ void swapLineGJ(float* mat, const int row1, const int row2, const int nCol) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < nCol; i += step)
    {
        float temp = mat[i + row1 * nCol];
        float temp2 = mat[i + row2 * nCol]; 
        mat[i + row1 * nCol] = temp2; // or mat[i + row1 * nCol] = mat[i + row2 * nCol];
        mat[i + row2 * nCol] = temp;
    }
}

__global__ void eliminationGJ(float* mat, float* matAug, const int r, const int nRow, const int nCol) {

    // un bloc = une ligne, 
    int index = threadIdx.x;
    int row = blockIdx.x;
    int step = blockDim.x;
    __shared__ float shFactor;
    if (row != r) { // le bloc r ne fait rien... bah...
        if (index == 0) {
            shFactor = mat[row * nCol + r];
        }
        __syncthreads();
        for (int j = index; j < nCol; j+=step) {
            float value1 = mat[r * nCol + j];
            float oldvalue1 = mat[row * nCol + j];
            float oldvalue2 = matAug[row * nCol + j];
            float value2 = matAug[r * nCol + j];

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
            float factor = mat[i * nCol + r]; // ne doit pas changer tant que la ligne n'est pas fini
            for (int j = indexX; j < nCol; j += stepX)
            {
                if (j != r) {
                    float value1 = mat[r * nCol + j];
                    float value2 = matAug[r * nCol + j];

                    mat[i * nCol + j] = mat[i * nCol + j] - factor * value1;
                    matAug[i * nCol + j] = matAug[i * nCol + j] - factor * value2;
                }
            }
        }
    }
    __syncthreads();
    for (int i = indexY; i < nRow; i += stepY)
    {
        float factor = mat[i * nCol + r]; // ne doit pas changer tant que la ligne n'est pas fini
        if (i != r) {
            if (indexX == 1) {
                mat[i * nCol + r] = 0;
                matAug[i * nCol + r] = matAug[i * nCol + r] - factor * matAug[r * nCol + r];
            }
        }
    }*/
}



__global__ void initPermMatr(float* P, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    for (int i = index; i < (N+1); i += step)
    {
        P[i] = i*(i<N);
    }
}


__global__ void updatePermMatr(float* P, const int line1, const int line2, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;
    if (index == 0) {
        int inter = P[line1];
        P[line1] = P[line2];
        P[line2] = inter;
        P[N] = P[N] + 1;
    }
}


__global__ void updateLUPFactorization(float* A, const int col, const int N) {
    // un bloc par ligne i ?
    int index = threadIdx.x;
    int i = blockIdx.x;
    int step = blockDim.x;

    __shared__ float Aicol;

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





__global__ void setPermute(float* y, float* b, float* P, const int N) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int step = blockDim.x * gridDim.x;

    for (int i = index; i < N; i += step)
    {
        int indice = P[i];
        y[i] = b[indice]; // c'est absolument moche...
    }


}
__global__ void solveLowSys(float* A, float* y, const int iter, const int N) {
    int index = threadIdx.x;
    int step = blockDim.x;
    __shared__ float yiter;

    if (index == 0) {
        yiter = y[iter];
        
    }
    __syncthreads();
    for (int i = index + iter+1; i < N; i += step)
    {
        y[i] = y[i] - yiter * A[i * N + iter]; // moche ne faudrait-il pas stocker A^T ?

    }


}

__global__ void solveUpSys(float* A, float* y, const int iter, const int N) {
    int index = threadIdx.x;
    int step = blockDim.x;
    __shared__ float yiter;

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


/* 
 int n = getNLin();
        setPermute << <_numBlocks, _blockSize >> > (_matrixGPU, b->_matrixGPU, P->_matrixGPU, n);
        for (int i = 0; i < n; i++) {
            solveLowSys <<< 1, _blockSize >>> (A->_matrixGPU, _matrixGPU, i, n);
        }


        for (int i = n - 1; i >= 0; i--)
        {
            solveUpSys << < 1, _blockSize >> > (A->_matrixGPU, _matrixGPU, i, n);
        }
*/

__global__ void solveSysGPU(float* A, float* y, const int N) {


    int index = threadIdx.x;
    int step = blockDim.x;
    extern __shared__ float ytemp[];


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
