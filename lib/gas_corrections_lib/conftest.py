import os
import sys
directory = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, directory + os.sep + "gas_corrections_cpp" + os.sep + "src")

# This will run the build.py script in the src/ directory
# The build script will compile the gas transmittance library and use cppimport to make the python module available for import
# Doing this in conftest.py will ensure that the C++ library is built and up-to-date prior to running the tests
import build