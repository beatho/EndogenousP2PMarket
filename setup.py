
from pathlib import Path

from setuptools import setup, Extension, find_packages
from cuda_extension import BuildCudaFiles, CustomCUDABuild, CudaBuildExt
from sysconfig import get_paths
import os
cuda_ext_path = Path('./src')
CudaFlags = ['--gpu-architecture=sm_52', '--device-c', '-dlink' ]

#print(list(map(lambda x: str(cuda_ext_path / x), os.listdir("src"))))
print(get_paths()["include"] )
import sys

out_ext_lib = ".a"
extra_linker_options = []
extra_compiler_options = []
optional_libraries = []
optional_library_dirs = []
data_files = []
if sys.platform == "win32":
	out_ext_lib = ".lib"
	extra_linker_options += ["build\\cuda\\dlink" + out_ext_lib, "/EXPORT:PyInit_EndoCuda"]
	extra_compiler_options += ["/MT", ]
	optional_libraries += ["user32", "kernel32"]
	optional_library_dirs += ["C:/Windows/System32"]
else:
	try:
		os.environ["CUDA_PATH"] + "42"
	except KeyError:
		os.environ["CUDA_PATH"] = "/usr/lib/cuda"
	#extra_linker_options = list(map(lambda x: "build/cuda/obj/" + x.replace(".cu", ".o"), filter(lambda x : x.endswith(".cu"),os.listdir("src"))))  + ["build/cuda/dlink.o"]
	extra_linker_options += ["-Wl,-rpath=$ORIGIN", "-Lbuild/cuda","-ldlink", ]
	extra_compiler_options = ["-fpermissive", "-fpic"]
	data_files = ["build/cuda/libdlink.so"]
cuda_ext = Extension(
	name='EndoCuda',
	include_dirs=[cuda_ext_path / 'include', os.environ["CUDA_PATH"] + "/include", get_paths()["include"]],
	library_dirs=[os.environ["CUDA_PATH"] + "/lib/x64", "./build/cuda/"] + optional_library_dirs,
	sources=list(map(lambda x: str(cuda_ext_path / x), os.listdir("src"))),
	libraries=["cudart_static"] + optional_libraries,  # Use fix_dll() only for Windows compatibility (check documentation for more info).
	extra_compile_args=["-O2"] + extra_compiler_options, #
	 extra_link_args=[] + extra_linker_options,
)




setup(
	name='EndoCuda',
	version='0.0.1',
	extras_require={'cython': ['cython'], },
	ext_modules=[cuda_ext],
	cmdclass={
		"build" : CustomCUDABuild,
		"build_cu" : BuildCudaFiles,
		"build_ext" : CudaBuildExt
	},
	data_files=data_files

)