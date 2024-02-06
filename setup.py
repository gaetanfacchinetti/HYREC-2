import subprocess
import os

# go to where the project needs to be compiled
os.chdir("./src/pyhyrec/src/")

try:
    subprocess.run(['make'], check=True)
    print("make executed successfully.")

    # Run the 'make' command with the desired target
    subprocess.run(['make', 'clean'], check=True)
    print("make clean executed successfully.")

except subprocess.CalledProcessError as e:
    print("Error executing Makefile:", e)

# go back to the main project folder
os.chdir("../../../")

from setuptools import setup, Extension
from Cython.Build import cythonize

cythonize("./src/pyhyrec/wrapperhyrec.pyx", force = True)

extensions = [Extension("pyhyrec.wrapperhyrec", sources = ["./src/pyhyrec/wrapperhyrec.c"], libraries = ["hyrec"], library_dirs=["./src/pyhyrec/src/"], build_temp="src/pyhyrec/")]

setup(ext_modules=extensions, package_dir={'': 'src'},  packages=['pyhyrec'], package_data= {'pyhyrec' :['data/*']}, include_package_data=True)