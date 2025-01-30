from pathlib import Path
from setuptools import setup
from setuptools_cuda_cpp import CUDAExtension, BuildExtension, fix_dll
import os
cuda_ext_path = Path('./src')
cuda_ext = CUDAExtension(
    name='EndoCuda',
    include_dirs=[cuda_ext_path / 'include'],
    sources=list(map(lambda x: cuda_ext_path / x, os.listdir("src"))),
    libraries=fix_dll(['cudart']),  # Use fix_dll() only for Windows compatibility (check documentation for more info).
    extra_compile_args={
        'cxx': ['-g'],  # cpp compiler flags
        'nvcc': ['-O2'],  # nvcc flags
    },
)

setup(
    name='my-cuda-package',
    version='0.0.1',
    install_requires=['numpy', ],
    extras_require={'cython': ['cython'], },
    ext_modules=[cuda_ext],
    cmdclass={'build_ext': BuildExtension},
)