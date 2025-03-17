
from pathlib import Path

from setuptools import setup, Extension
from cuda_extension import BuildCudaFiles, CustomCUDABuild, CudaBuildExt
from sysconfig import get_paths

import os
cuda_ext_path = Path('./src')
print(get_paths()["include"] )

CudaFlags = ['--gpu-architecture=sm_52', '--device-c', '-dlink' ]

import sys

extra_linker_options = []
extra_compiler_options = []
if sys.platform == "win32":
	extra_linker_options += []
	extra_compiler_options += ["/MT"]

cuda_ext = Extension(
	name="EndoCuda",
	include_dirs=[cuda_ext_path / 'include', os.environ["CUDA_PATH"] + "/include", get_paths()["include"] ],
	library_dirs=[os.environ["CUDA_PATH"] + "/lib/x64", "./build/cuda/"],
	sources=list(map(lambda x: str(cuda_ext_path / x), os.listdir("src"))),
	libraries=["cudart", "cudadevrt"],  # Use fix_dll() only for Windows compatibility (check documentation for more info).
	extra_compile_args=[] + extra_compiler_options, #
	 extra_link_args=["build/cuda/dlink.lib", ] + extra_linker_options,
)




setup(
	name='EndoCuda',
	version='0.0.1',
	install_requires=['numpy', ],
	extras_require={'cython': ['cython'], },
	ext_modules=[cuda_ext],
	cmdclass={
		"build" : CustomCUDABuild,
		"build_cu" : BuildCudaFiles,
		"build_ext" : CudaBuildExt
	},
	exclude_package_data={
		"EndoCuda": ["*.cu"],  # Files to ignore
	},
)