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
		out_ext = ".o"
		print(sys.platform)
		if sys.platform == "win32":
			out_ext = ".obj"
		print("building CUDA files")
		built = []
		for root, _, files in os.walk("src"):
			for file in files:
				if file.endswith(".cu"):
					in_path = os.path.join(root, file)
					out_path = os.path.join("./build/cuda/obj", file.split(".cu")[0] + out_ext)
					subprocess.run(["nvcc", in_path, "-o", out_path, "--device-c", "-I", get_paths()["include"]])
					built.append(out_path)
		print("performing dlink step")
		subprocess.run(["nvcc","-dlink"] + built + [ "-o", "./build/cuda/dlink.obj" ])
		subprocess.run(["nvcc","--lib" ] + built + [ "./build/cuda/dlink.obj" , "-o", "./build/cuda/dlink.lib" ])		
		

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