import shutil
import glob
import os
from pathlib import Path
import cppimport

print('Building gas transmittance lib...')
previous_cwd = os.getcwd()
this_script_location = Path(__file__).resolve().parent
os.chdir(this_script_location)

funcs = cppimport.imp("gas_transmittance")
for file in glob.glob(r'*.so'):
    dest = '../bin/' + file
    print('Copying ', file, ' to ', dest)
    shutil.copy(file, dest)  # Using shutil.move() caused a "free(): invalid pointer" error at the end of the program. Using copy() then remove() solves this issue.
    if os.path.exists(file):
        os.remove(file)
    
os.chdir(previous_cwd)