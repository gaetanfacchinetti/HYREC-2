import subprocess

try:
    subprocess.run(['make'], check=True, cwd="./src/")
    print("make executed successfully.")

    # Run the 'make' command with the desired target
    subprocess.run(['make', 'clean'], check=True, cwd="./src/")
    print("make clean executed successfully.")

except subprocess.CalledProcessError as e:
    print("Error executing Makefile:", e)


from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension("pyhyrec", sources = ["pyhyrec.pyx"], libraries = ["hyrec"], library_dirs=["./src/"])

setup(ext_modules=cythonize(ext))