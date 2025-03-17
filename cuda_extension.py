from setuptools import Command
import os
from setuptools.command.build import build
from setuptools.command.build_ext import build_ext
from pathlib import Path
from sysconfig import get_paths
import sys
import subprocess
class BuildCudaFiles(Command):
	"""a build process for the custom files"""
	user_options = []

	def initialize_options(self):
		pass
	def finalize_options(self):
		pass
	def run(self):
		from pathlib import Path
		Path("./build/cuda/obj").mkdir(exist_ok=True, parents=True)
		extra_compiler_options = []
		extra_linker_options = []
		out_ext_o = ".o"
		out_ext_lib = ".a"
		print(sys.platform)
		if sys.platform == "win32":
			out_ext_o = ".obj"
			out_ext_lib = ".lib"
		else:
			extra_compiler_options += ["-std=c++14", "--compiler-options" , "-fPIC", "-rdc=true"]
			extra_linker_options = ["-rdc=true"]
		print("building CUDA files")
		built = []
		for root, _, files in os.walk("src"):
			for file in files:
				if file.endswith(".cu"):
					in_path = os.path.join(root, file)
					out_path = os.path.join("./build/cuda/obj", file.split(".cu")[0] + out_ext_o)
					subprocess.run(["nvcc", "-O2", ] + extra_compiler_options + [ in_path, "-o", out_path, "--device-c", "-I", get_paths()["include"] ])
					built.append(out_path)
		if sys.platform == "win32":
			print("performing dlink step")
			subprocess.run(["nvcc","-dlink"] + extra_linker_options + built + [ "-o", "./build/cuda/dlink" + out_ext_o ])
			print("making a library")
			subprocess.run(["nvcc","--lib" ] + extra_linker_options + built + [ "./build/cuda/dlink" + out_ext_o , "-o", "./build/cuda/dlink" + out_ext_lib ])
		else:
			print("making a shared library")
			subprocess.run(["nvcc", "--shared", "-o", "./build/cuda/libdlink.so"] + built)

		

class CustomCUDABuild(build):
	"""Extend the build command to include custom file processing"""

	def run(self):
		"""Run the original build command and then process custom files"""
		self.run_command("build_cu")  # Run our custom command
		#print(zip(self))
		build.run(self)

class CudaBuildExt(build_ext):
    """Custom build_ext to ignore unwanted files"""

    def run(self):
        """Run the normal build_ext but filter out unwanted files"""
        for ext in self.extensions:
            ext.sources = [s for s in ext.sources if not s.endswith(".cu")]
		
        super().run() 