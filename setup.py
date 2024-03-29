import subprocess
import os

from setuptools import setup, Extension
from Cython.Build import cythonize

def main():

    # go to where the project needs to be compiled
    
    os.chdir("./src/pyhyrec/src/")

    try:
        
        if os.path.exists("libhyrec.so"):
            os.remove("libhyrec.so")

        subprocess.run(['make'], check=True)
        print("make executed successfully.")

        # Run the 'make' command with the desired target
        subprocess.run(['make', 'clean'], check=True)
        print("make clean executed successfully.")

    except subprocess.CalledProcessError as e:
        print("Error executing Makefile:", e)

    # go back to the main project folder
    os.chdir("../../../")
   
    #extensions = [Extension("pyhyrec.src.hyrec", sources = ["./src/pyhyrec/src/hyrec.c"])]
    extensions = [Extension("pyhyrec.wrapperhyrec", sources = ["./src/pyhyrec/wrapperhyrec.c"], libraries = ["hyrec"], library_dirs=["./src/pyhyrec/src/"], build_temp="src/pyhyrec/")]
    cythonize("./src/pyhyrec/wrapperhyrec.pyx", force = True)
    
    setup(ext_modules=extensions, package_dir={'' : 'src'},  packages=['pyhyrec'], package_data= {'pyhyrec' :['data/*']}, include_package_data=True)


if __name__ == "__main__":
    main()