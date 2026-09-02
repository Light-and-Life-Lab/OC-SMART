# OC-SMART
Ocean Color - Simultaneous Marine and Aerosol Retrieval Tool (OC-SMART) is a multi-sensor data analysis platform designed to retrieve aerosol and ocean color products from satellite remote sensing images.

# Installation
## Linux and MacOS

### 1. Download and Extract OC-SMART

The [OC-SMART releases page](https://github.com/Light-and-Life-Lab/OC-SMART/releases) lists the versions of OC-SMART available for installation. 
Select the version you would like to download and then copy the link to the tar.gz file to your clipboard. 
Navigate to the location where you would like to download OC-SMART on your machine. Then you can download and extract the archive using the command:

`curl -L https://github.com/Light-and-Life-Lab/OC-SMART/archive/refs/tags/v2.6.5.tar.gz | tar -xzf -`

### 2. (Optional) Install the Conda Package Manager

Some steps in the install process are simpler if you use the Conda package manager to handle the dependencies. However, Conda is not strictly required so instructions are included for both Conda and the default Python package manager, pip.

If you are able to use Conda, the recommended distribution is Miniconda. Instructions for installing Miniconda can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### 3. Set Up and Activate Python Virtual Environment

It is recommended that you create a fresh virtual environment in which to install OC-SMART to avoid potential version conflicts with existing Python package installations on your machine.

#### Conda

To set up a Conda virtual environment, run the command:

`conda create -n <envname>`

Then activate it using:

`conda activate <envname>`

Even if you are using Conda, there are certain steps in the process that require installation of Python wheels via URL. This cannot be done directly using the `conda install` command, so it is required that Conda users install pip within their Conda virtual environment. Once you have activated your virtual environment using `conda activate <envname>`, run the command:

`conda install pip`

This will install an instance of the pip package manager that lives inside the Conda virtual environment. 

**Note:** if you already have a different installation of Python (common for most Linux distributions, which ship with Python) then simply using `pip install <package_name>` will typically default to the version of pip that was originally installed with the OS. But for our purposes we need to use the newly installed version of pip that lives inside the Conda virtual environment. **To avoid ambiguity, all `pip install <package_name>` commands used throughout this installation process should instead be `python -m pip install <package_name>`.** The `python -m` part of the command will ensure that the version of pip used for installing packages is the one installed to the currently active environment. So as long as this command is run while the Conda virtual environment is active, then the correct version of pip will be used.

#### Pip

To set up a regular Python virtual environment, navigate to the root of the project directory (i.e. OC-SMART) and run the command:

`python -m venv <envname>`

Then activate it using:

`source <envname>/bin/activate`

### 4. Install Regular Python Dependencies

OC-SMART has a number of dependencies, most of which can be pip installed via the requirements.txt file. There are two additional dependencies that are not available on the pip package index and so must be installed manually via Python wheels.

#### Conda and Pip

For both Conda and pip this step is the same. A requirements.txt file is provided in the root directory of the OC-SMART repository. From this directory, install these dependencies using the following command:

`python -m pip install -r requirements.txt`

### 5. Install GDAL

#### Conda

GDAL is registered to the conda-forge channel, so installing is straightforward. Simply run the command:

`conda install conda-forge::gdal`

#### Pip

a. The GDAL library can be installed via your OS's package manager. Please follow the instructions for your OS/Distribution:

Debian/Ubuntu (including WSL):

`sudo apt update`

`sudo apt install -y gdal-bin libgdal-dev`

RHEL/Fedora/CentOS:

`sudo dnf install -y gdal gdal-devel`

macOS:

`brew install gdal`

b. Install the matching Python GDAL bindings:

`pip install "GDAL==$(gdal-config --version)"`

c. Troubleshooting

If installed successfully, the command

`gdal-config --version`

will print the GDAL version number. Some versions of GDAL are not available on the Python Package Index (PyPI) and so the above steps will fail. If this is the case for your version, check the available versions at https://pypi.org/project/GDAL/#history and install the closest match to the version produced by `gdal-config --version`. For example:

`pip install "GDAL==3.8.4"`

### 6. Install l8angles

The l8angles library is a Python wrapper around a C-based USGS Landsat 8 tool for computing per-pixel solar and sensor azimuth and zenith angles from Angle Coefficient Files. A set of Python wheels has been generated for the platforms supported by OC-SMART, which may be found here: **[l8-angles](https://github.com/Light-and-Life-Lab/l8-angles)**. Follow the installation instructions on the linked page in order to install the appropriate Python wheel for your platform.

### 6. Install Atmospheric Gas Correction Library

A C++-based library for computing gas transmittance values is also used by OC-SMART. This library provides a Python interface and may be imported and used as a Python module for seamless integration with Python workflows. A set of Python wheels has been generated for the platforms supported by OC-SMART, which may be found here: **[Atmospheric Gas Correction Library](https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library)**. Follow the installation instructions on the linked page in order to install the appropriate Python wheel for your platform.

## Windows
Native Windows wheels are not built or supported. Windows users should install and run this library from within **[WSL (Windows Subsystem for Linux)](https://learn.microsoft.com/en-us/windows/wsl/install)**, then follow the Linux installation instructions above from inside your WSL environment.


### 7. Consult User Guide

The attached User Guide (UserGuide_Python_Linux.pdf) contains these same install instructions in Section 1. In Section 2, there is additional information on one-time setup of data download credentials. The remaining sections of the guide provide additional information on configuration settings, input/output files, etc.