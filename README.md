# OC-SMART
Ocean Color - Simultaneous Marine and Aerosol Retrieval Tool (OC-SMART) is a multi-sensor data analysis platform designed to retrieve aerosol and ocean color products from satellite remote sensing images.

# Installation
## Linux and MacOS

OC-SMART has a number of dependencies, most of which can be pip installed via the requirements.txt file. There are two additional dependencies that are not available on the pip package index and so must be installed manually via Python wheels.

### 1. Install Regular Python Dependencies

A requirements.txt file is provided in the root directory of the OC-SMART repository. Install these dependencies using the following command:

pip install -r requirements.txt

### 2. Install l8angles

The l8angles library is a Python wrapper around a C-based USGS Landsat 8 tool for computing per-pixel solar and sensor azimuth and zenith angles from Angle Coefficient Files. A set of Python wheels has been generated for the platforms supported by OC-SMART, which may be found here: **[l8-angles](https://github.com/Light-and-Life-Lab/l8-angles)**. Follow the installation instructions on the linked page in order to install the appropriate Python wheel for your platform.

### 3. Install Atmospheric Gas Correction Library

A C++-based library for computing gas transmittance values is also used by OC-SMART. This library provides a Python interface and may be imported and used as a Python module for seamless integration with Python workflows. A set of Python wheels has been generated for the platforms supported by OC-SMART, which may be found here: **[Atmospheric Gas Correction Library](https://github.com/Light-and-Life-Lab/Atmospheric_Gas_Correction_Library)**. Follow the installation instructions on the linked page in order to install the appropriate Python wheel for your platform.

## Windows
Native Windows wheels are not built or supported. Windows users should install and run this library from within **[WSL (Windows Subsystem for Linux)](https://learn.microsoft.com/en-us/windows/wsl/install)**, then follow the Linux installation instructions above from inside your WSL environment.
