import h5py
import numpy as np
from scipy import interpolate

from ..gas_corrections_cpp.bin import gas_transmittance


def L1_Data():
    return gas_transmittance.L1_Data()
    

def Ancillary_Data():
    return gas_transmittance.Ancillary_Data()


def Oxygen_A_Band_Option():
    return gas_transmittance.Oxygen_A_Band_Option


def validate_keyword_args(args):
    for key, value in args.items():
        if value is None:
            raise ValueError(f"Keyword argument {key} cannot have value {value}. Please populate the {key} variable and pass it as a keyword argument.")
        

def check_for_invalid_keyword_args(args, valid_kwargs):
    for key in args.keys():
        if key not in valid_kwargs:
            raise ValueError(f"{key} is not a valid keyword argument in this context.")


class Gas_Correction_Manager:
    def __init__(self):
        self.l1_data = None
        self.l1_filename = None
        self.gas_transmittance_table = None
        self.gas_transmittance_table_filename = None
            

    def read_PACE_data(self, l1_filename, **kwargs):
        """
        Reads PACE data from a file in the NetCDF file format and stores it into the l1_data member variable.
        The file contents are cached, so if this function is called for a second time on the same file name then the file contents
        that were previously read will simply be returned again to prevent redundant (and potentially time consuming) file reads.
        Optional keyword arguments can be provided to specify subimage (via start/end lines and start/end pixels).

        Args: 
            l1_filename (str): Name of the PACE NetCDF file to be read (with path included if necessary).

        Keyword Args:
            The following keyword args may be used to define a subimage of the PACE data. The keyword args are used as indices to
                slice a numpy array, i.e. arr[start_line:end_line, start_pixel:end_pixel]. If no keyword args are specified,
                then the whole image is used, i.e. arr[0:, 0:]
            start_line (int): Used to define the first line (i.e. row) of a subimage. Default value is 0.
            end_line (int): Used to define the last line (i.e. row) of a subimage. Default value is None.
            start_pixel (int): Used to define the first line (i.e. column) of a subimage. Default value is 0.
            end_line (int): Used to define the last line (i.e. column) of a subimage. Default value is None.

        Returns (gas_transmittance.L1_Data):
            An instance of the L1_Data class (available in the gas corrections library), which contains the PACE data that was read
            by this function. It is the same instance of the class that is stored in the self.l1_data member variable.

        Example Usage:
            gas_correction_manager = gas_corrections.Gas_Correction_Manager()
            gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
        """
        if self.l1_data is not None and l1_filename == self.l1_filename:
            # Avoid redundant reading from a file we already read data from
            return self.l1_data
        else:
            l1_data = gas_transmittance.L1_Data()
            start_line = kwargs.get("start_line", 0)
            end_line = kwargs.get("end_line", None)
            start_pixel = kwargs.get("start_pixel", 0)
            end_pixel = kwargs.get("end_pixel", None)

            with h5py.File(l1_filename, 'r') as f:
                solar_zenith = 0.01*np.array(f['/geolocation_data/solar_zenith'][start_line:end_line, start_pixel:end_pixel])
                sensor_zenith = 0.01*np.array(f['/geolocation_data/sensor_zenith'][start_line:end_line, start_pixel:end_pixel])
                l1_data.cos_solar_zenith = np.cos(np.deg2rad(solar_zenith))
                l1_data.cos_sensor_zenith = np.cos(np.deg2rad(sensor_zenith))

                l1_data.latitude = np.flip(np.array(f['/geolocation_data/latitude']), 0)
                l1_data.longitude = np.flip(np.array(f['/geolocation_data/longitude']), 0)

                blue_wavelengths = np.array(f['/sensor_band_parameters/blue_wavelength'][1:])
                red_wavelengths = np.array(f['/sensor_band_parameters/red_wavelength'][3:])
                sensor_wavelengths = np.zeros(len(blue_wavelengths) + len(red_wavelengths))
                sensor_wavelengths[0:len(blue_wavelengths)] = blue_wavelengths
                sensor_wavelengths[len(blue_wavelengths):] = red_wavelengths

                blue_rhot = np.array(f['/observation_data/rhot_blue'][1:, start_line:end_line, start_pixel:end_pixel])
                red_rhot = np.array(f['/observation_data/rhot_red'][3:, start_line:end_line, start_pixel:end_pixel])
                assert(blue_rhot.shape[1] == red_rhot.shape[1])
                assert(blue_rhot.shape[2] == red_rhot.shape[2])

                rhot = np.zeros((blue_rhot.shape[0] + red_rhot.shape[0], blue_rhot.shape[1], blue_rhot.shape[2]))
                rhot[0:len(blue_wavelengths), :, :] = blue_rhot
                rhot[len(blue_wavelengths):, :, :] = red_rhot
                rhot = np.rollaxis(rhot, 0, 3)
                l1_data.reflectance = rhot/np.pi*l1_data.cos_sensor_zenith[:, :, None]

                assert(len(sensor_wavelengths) == l1_data.reflectance.shape[2])

                l1_data.wavelengths = sensor_wavelengths
                l1_data.num_pixels = len(l1_data.cos_solar_zenith.flatten())
                l1_data.num_wavelengths = len(l1_data.wavelengths)

            self.l1_filename = l1_filename
            self.l1_data = l1_data
            return self.l1_data


    def read_HICO_data(self, l1_filename, **kwargs):
        """
        Reads HICO data from a file in the NetCDF file format and stores it into the l1_data member variable.
        The file contents are cached, so if this function is called for a second time on the same file name then the file contents
        that were previously read will simply be returned again to prevent redundant (and potentially time consuming) file reads.
        Optional keyword arguments can be provided to specify subimage (via start/end lines and start/end pixels).

        Args: 
            l1_filename (str): Name of the HICO NetCDF file to be read (with path included if necessary).

        Keyword Args:
            The following keyword args may be used to define a subimage of the HICO data. The keyword args are used as indices to
                slice a numpy array, i.e. arr[start_line:end_line, start_pixel:end_pixel]. If no keyword args are specified,
                then the whole image is used, i.e. arr[0:, 0:]
            start_line (int): Used to define the first line (i.e. row) of a subimage. Default value is 0.
            end_line (int): Used to define the last line (i.e. row) of a subimage. Default value is None.
            start_pixel (int): Used to define the first line (i.e. column) of a subimage. Default value is 0.
            end_line (int): Used to define the last line (i.e. column) of a subimage. Default value is None.

        Returns (gas_transmittance.L1_Data):
            An instance of the L1_Data class (available in the gas corrections library), which contains the PACE data that was read
            by this function. It is the same instance of the class that is stored in the self.l1_data member variable.

        Example Usage:
            gas_correction_manager = gas_corrections.Gas_Correction_Manager()
            gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
        """
        if self.l1_data is not None and l1_filename == self.l1_filename:
            # Avoid redundant reading from a file we already read data from
            return self.l1_data
        else:
            l1_data = gas_transmittance.L1_Data()
            start_line = kwargs.get("start_line", 0)
            end_line = kwargs.get("end_line", None)
            start_pixel = kwargs.get("start_pixel", 0)
            end_pixel = kwargs.get("end_pixel", None)

            with h5py.File(l1_filename, 'r') as f:
                solar_zenith = 0.01*np.array(f['/navigation/solar_zenith'][start_line:end_line, start_pixel:end_pixel])
                sensor_zenith = 0.01*np.array(f['/navigation/sensor_zenith'][start_line:end_line, start_pixel:end_pixel])
                l1_data.cos_solar_zenith = np.cos(np.deg2rad(solar_zenith))
                l1_data.cos_sensor_zenith = np.cos(np.deg2rad(sensor_zenith))

                l1_data.latitude = np.flip(np.array(f['/navigation/latitudes']), 0)
                l1_data.longitude = np.flip(np.array(f['/navigation/longitudes']), 0)

                # blue_wavelengths = np.array(f['/sensor_band_parameters/blue_wavelength'][1:])
                # red_wavelengths = np.array(f['/sensor_band_parameters/red_wavelength'][3:])
                # sensor_wavelengths = np.zeros(len(blue_wavelengths) + len(red_wavelengths))
                # sensor_wavelengths[0:len(blue_wavelengths)] = blue_wavelengths
                # sensor_wavelengths[len(blue_wavelengths):] = red_wavelengths

                # blue_rhot = np.array(f['/observation_data/rhot_blue'][1:, start_line:end_line, start_pixel:end_pixel])
                # red_rhot = np.array(f['/observation_data/rhot_red'][3:, start_line:end_line, start_pixel:end_pixel])
                # assert(blue_rhot.shape[1] == red_rhot.shape[1])
                # assert(blue_rhot.shape[2] == red_rhot.shape[2])

                # rhot = np.zeros((blue_rhot.shape[0] + red_rhot.shape[0], blue_rhot.shape[1], blue_rhot.shape[2]))
                # rhot[0:len(blue_wavelengths), :, :] = blue_rhot
                # rhot[len(blue_wavelengths):, :, :] = red_rhot
                # rhot = np.rollaxis(rhot, 0, 3)
                # l1_data.reflectance = rhot/np.pi*l1_data.cos_sensor_zenith[:, :, None]

                # assert(len(sensor_wavelengths) == l1_data.reflectance.shape[2])

                # l1_data.wavelengths = sensor_wavelengths
                l1_data.num_pixels = len(l1_data.cos_solar_zenith.flatten())
                # l1_data.num_wavelengths = len(l1_data.wavelengths)

            self.l1_filename = l1_filename
            self.l1_data = l1_data
            return self.l1_data


    def read_gas_transmittance_table(self, gas_transmittance_table_filename):
        """
        Reads gas transmittance lookup tables from a NetCDF file and stores them in the gas_transmittance_table member variable.
        The file contents are cached, so if this function is called for a second time on the same file name then the file contents
        that were previously read will simply be returned again to prevent redundant (and potentially time consuming) file reads.

        Args:
            gas_transmittance_table_filename (str): Name of the lookup table NetCDF file to be read (with path included if necessary).

        Returns (gas_transmittance.Gas_Transmittance_Lookup_Table):
            An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains the 
            lookup tables read in by this function. It is the same instance of the class that is stored in the self.gas_transmittance_table member variable.

        Example Usage:
            gas_correction_manager = gas_corrections.Gas_Correction_Manager()
            gas_correction_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
        """
        if self.gas_transmittance_table is not None and gas_transmittance_table_filename == self.gas_transmittance_table_filename:
            # Avoid redundant reading from a file we already read data from
            return self.gas_transmittance_table
        else:
            gas_transmittance_table = gas_transmittance.Gas_Transmittance_Lookup_Table()

            netcdf_variables = {
                    'air_mass_factor_mixed' : "air_mass_factor_mixed_gases",
                    'air_mass_factor_wv' : "air_mass_factor_water_vapor",
                    'carbon_dioxide_transmittance' : "co2_transmittance",
                    'carbon_monoxide_transmittance' : "co_transmittance",
                    'methane_transmittance' : "ch4_transmittance",
                    'nitrous_oxide_transmittance' : "n2o_transmittance",
                    'oxygen_transmittance' : "o2_transmittance",
                    'water_vapor' : "water_vapor_concentration",
                    'water_vapor_transmittance' : "h2o_transmittance",
                    'wavelength' : "wavelengths",
                }
            
            netcdf_dimensions = {
                'n_air_mass_factor' : "num_amf_grid_points",
                'n_water_vapor' : "num_water_vapor_concentrations",
                'nmodels' : 'num_models',
                'nwavelengths' : 'num_wavelengths',
            }

            model_map = {
                "Tropical": 0,
                "MidLatSummer" : 1,
                "MidLatWinter" : 2,
                "SubarcticSummer" : 3,
                "SubarcticWinter" : 4,
                "USstandard62" : 5
            }
                
            with h5py.File(gas_transmittance_table_filename, 'r') as f:
                for netcdf_node, var_name in netcdf_variables.items():
                    if netcdf_node in f:
                        var_value = np.array(f[netcdf_node])
                        setattr(gas_transmittance_table, var_name, var_value)

                for netcdf_node, var_name in netcdf_dimensions.items():
                    if netcdf_node in f:
                        var_value = np.array(f[netcdf_node])
                        setattr(gas_transmittance_table, var_name, len(var_value))

                setattr(gas_transmittance_table, 'model', model_map["USstandard62"])

            self.gas_transmittance_table_filename = gas_transmittance_table_filename
            self.gas_transmittance_table = gas_transmittance_table
            return self.gas_transmittance_table


def ozone_transmittance(**kwargs):
    """
    Computes the ozone transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which
            contains the image for which transmittance data is to be applied. Default value is None.
        **ancillary_data (gas_transmittance.Ancillary_Data)**: An instance of the Ancillary_Data class (available in the gas corrections library),
            which contains ozone cross section and ozone concentration data that has been interpolated to the L1 Data grid. Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in ozone data from appropriate ancillary file, e.g. ozone_climatology_v2014.hdf. Assume ozone lat and lon grids are
        stored in oz_lat and oz_lon, and that ozone concentration values are stored in an array ozmap with the same dimensions.
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = cos_solar_zenith
        l1_data.cos_sensor_zenith = cos_sensor_zenith
        l1_data.num_pixels = len(cos_solar_zenith)
        l1_data.num_wavelengths = len(ozone_absorption_cross_section)

        # Interpolate ozone map to the L1B grid
        func = interpolate.RegularGridInterpolator((np.flip(oz_lat), oz_lon), np.flip(ozmap, 0))
        ozone_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())
        ozone_concentration = ozone_concentration

        ancillary_data = gas_corrections.Ancillary_Data()
        ancillary_data.ozone_absorption_cross_section = ozone_absorption_cross_section
        ancillary_data.ozone_concentration = ozone_concentration

        gas_transmittances = gas_corrections.ozone_transmittance(l1_data=l1_data, \
                                                                 ancillary_data=ancillary_data)
    """
    valid_kwargs = ['l1_data', 'ancillary_data']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    ancillary_data = args['ancillary_data'] = kwargs.get("ancillary_data", None)

    validate_keyword_args(args)

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(ancillary_data.ozone_concentration.size != 0)
    assert(ancillary_data.ozone_absorption_cross_section.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.cos_solar_zenith.shape == ancillary_data.ozone_concentration.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_pixels == ancillary_data.ozone_concentration.size)
    assert(l1_data.num_wavelengths == ancillary_data.ozone_absorption_cross_section.size)

    return gas_transmittance.ozone_transmittance(l1_data, ancillary_data)


def co2_transmittance(**kwargs):
    """
    Computes the Carbon Dioxide transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **gas_transmittance_table (gas_transmittance.Gas_Transmittance_Lookup_Table)**: An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains transmittance lookup tables read in from a NetCDF file. 
            Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
        gas_transmittance_manager.read_gas_transmittance_table("oci_gas_transmittance_cia_amf_v3.2.nc")

        gas_transmittances = gas_corrections.co2_transmittance(l1_data=l1_data, \
                                                               gas_transmittance_table=gas_transmittance_table)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    l1_data.num_wavelengths = len(l1_data.wavelengths)

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.co2_transmittance, axis = 0)
    co2_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.co2_transmittance = co2_transmittance_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.co2_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.co2_transmittance.shape[0])

    return gas_transmittance.co2_transmittance(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension)


def co_transmittance(**kwargs):
    """
    Computes the Carbon Monoxide transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **gas_transmittance_table (gas_transmittance.Gas_Transmittance_Lookup_Table)**: An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains transmittance lookup tables read in from a NetCDF file. 
            Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
        gas_transmittance_manager.read_gas_transmittance_table("oci_gas_transmittance_cia_amf_v3.2.nc")

        gas_transmittances = gas_corrections.co_transmittance(l1_data=l1_data, \
                                                              gas_transmittance_table=gas_transmittance_table)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    l1_data.num_wavelengths = len(l1_data.wavelengths)

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.co_transmittance, axis = 0)
    co_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.co_transmittance = co_transmittance_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.co_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.co_transmittance.shape[0])

    return gas_transmittance.co_transmittance(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension)


def ch4_transmittance(**kwargs):
    """
    Computes the Methane transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **gas_transmittance_table (gas_transmittance.Gas_Transmittance_Lookup_Table)**: An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains transmittance lookup tables read in from a NetCDF file. 
            Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
        gas_transmittance_manager.read_gas_transmittance_table("oci_gas_transmittance_cia_amf_v3.2.nc")

        gas_transmittances = gas_corrections.ch4_transmittance(l1_data=l1_data, \
                                                               gas_transmittance_table=gas_transmittance_table)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    l1_data.num_wavelengths = len(l1_data.wavelengths)

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.ch4_transmittance, axis = 0)
    ch4_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.ch4_transmittance = ch4_transmittance_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.ch4_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.ch4_transmittance.shape[0])

    return gas_transmittance.ch4_transmittance(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension)


def n2o_transmittance(**kwargs):
    """
    Computes the Nitrous Oxide transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **ancillary_data (gas_transmittance.Ancillary_Data)**: An instance of the Ancillary_Data class (available in the gas corrections library),
            which contains ozone cross section and ozone concentration data that has been interpolated to the L1 Data grid. Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(l1_data.wavelengths)

        gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
        gas_transmittance_manager.read_gas_transmittance_table("oci_gas_transmittance_cia_amf_v3.2.nc")

        gas_transmittances = gas_corrections.n2o_transmittance(l1_data=l1_data, \
                                                               gas_transmittance_table=gas_transmittance_table)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    l1_data.num_wavelengths = len(l1_data.wavelengths)

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.n2o_transmittance, axis = 0)
    n2o_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.n2o_transmittance = n2o_transmittance_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.n2o_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.n2o_transmittance.shape[0])

    return gas_transmittance.n2o_transmittance(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension)


def no2_transmittance(**kwargs):
    """
    Computes the Nitrogen Dioxide transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **ancillary_data (gas_transmittance.Ancillary_Data)**: An instance of the Ancillary_Data class (available in the gas corrections library),
            which contains no2 absorption cross section, stratospheric and tropospheric no2 concentrations, and fraction of no2 that lies above 200m.
            All of these quantities must be interpolated to the L1 Data grid. Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in no2 data from appropriate ancillary files, e.g. 'no2_climatology_v2013.hdf' and 'trop_f_no2_200m.hdf'. Assume no2 lat and lon grids are
        stored in no2_lat and no2_lon, and that no2 concentration values in  the troposphere and stratosphere are stored in an arrays no2_tropo and no2_strat with the same dimensions.
        Similarly, the fraction of no2 above 200 m is stored in an array no2_frac with associated lat/lon grids stored in no2_frac_lat and no2_frac_lon.
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_solar_zenith = cos_solar_zenith
        l1_data.cos_sensor_zenith = cos_sensor_zenith
        l1_data.num_pixels = len(cos_solar_zenith)
        l1_data.num_wavelengths = len(no2_absorption_cross_section)
        
        # Interpolate no2 map to the L1B grid
        func = interpolate.RegularGridInterpolator((np.flip(no2_frac_lat), no2_frac_lon), np.flip(no2_frac, 0))
        fraction_tropospheric_no2_above_200m = func(np.array([l1b_lat, l1b_lon]).transpose())

        no2_strat = no2_strat[int(month)-1, :, :]
        func = interpolate.RegularGridInterpolator((np.flip(no2_lat), no2_lon), np.flip(no2_strat, 0))
        stratospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

        no2_tropo = no2_tropo[int(month)-1, :, :]
        func = interpolate.RegularGridInterpolator((np.flip(no2_lat), no2_lon), np.flip(no2_tropo, 0))
        tropospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

        ancillary_data = gas_corrections.Ancillary_Data()
        ancillary_data.no2_absorption_cross_section = no2_absorption_cross_section
        ancillary_data.fraction_tropospheric_no2_above_200m = fraction_tropospheric_no2_above_200m
        ancillary_data.tropospheric_no2_concentration = tropospheric_no2_concentration
        ancillary_data.stratospheric_no2_concentration = stratospheric_no2_concentration

        gas_transmittances = gas_corrections.no2_transmittance(l1_data=l1_data, \
                                                               ancillary_data=ancillary_data)
    """
    valid_kwargs = ['l1_data', 'ancillary_data']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    ancillary_data = args['ancillary_data'] = kwargs.get("ancillary_data", None)

    validate_keyword_args(args)

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(ancillary_data.no2_absorption_cross_section.size != 0)
    assert(ancillary_data.fraction_tropospheric_no2_above_200m.size != 0)
    assert(ancillary_data.tropospheric_no2_concentration.size != 0)
    assert(ancillary_data.stratospheric_no2_concentration.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.cos_solar_zenith.shape == ancillary_data.fraction_tropospheric_no2_above_200m.shape)
    assert(l1_data.cos_solar_zenith.shape == ancillary_data.tropospheric_no2_concentration.shape)
    assert(l1_data.cos_solar_zenith.shape == ancillary_data.stratospheric_no2_concentration.shape)

    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_pixels == ancillary_data.fraction_tropospheric_no2_above_200m.size)
    assert(l1_data.num_pixels == ancillary_data.tropospheric_no2_concentration.size)
    assert(l1_data.num_pixels == ancillary_data.stratospheric_no2_concentration.size)
    assert(l1_data.num_wavelengths == ancillary_data.no2_absorption_cross_section.size)

    return gas_transmittance.no2_transmittance(l1_data, ancillary_data)


def o2_transmittance(**kwargs):
    """
    Computes the Oxygen transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **gas_transmittance_table (gas_transmittance.Gas_Transmittance_Lookup_Table)**: An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains transmittance lookup tables read in from a NetCDF file. 
            Default value is None.
        **oxygen_A_band_option (gas_corrections.Oxygen_A_Band_Option)**: An set of enumerated values describing the different options for computing the oxygen transmittance. The options are:
            NO_CORRECTION: Do not do any oxygen corrections
            DING_GORDON: Apply Ding and Gordon (1995) correction
            TRANSMITTANCE_TABLE: Apply oxygen transmittance from gas transmittance table
            SURROUNDING_WINDOW_BANDS: Compute oxygen transmittance from A-band and surrounding window bands (requires AMF gas trasmittance table)
            Default value is gas_corrections.Oxygen_A_Band_Option.TRANSMITTANCE_TABLE.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.reflectance = l1b_reflectance
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(sensor_wavelengths)

        gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table(gas_transmittance_filepath)

        gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, \
                                                              gas_transmittance_table=gas_transmittance_table, \
                                                              oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option.TRANSMITTANCE_TABLE)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table', 'oxygen_A_band_option']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)

    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    oxygen_A_band_option = kwargs.get("oxygen_A_band_option", gas_transmittance.Oxygen_A_Band_Option.TRANSMITTANCE_TABLE)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.o2_transmittance, axis = 0)
    o2_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.o2_transmittance = o2_transmittance_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.o2_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.o2_transmittance.shape[0])

    return gas_transmittance.o2_transmittance(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension, oxygen_A_band_option)


def h2o_transmittance(**kwargs):
    """
    Computes the Oxygen transmittance pixel-by-pixel for the image stored in the l1_data member variable.

    Keyword Args:
        **l1_data (gas_transmittance.L1_Data)**: An instance of the L1_Data class (available in the gas corrections library), which contains the image for which transmittance data is to be applied.
            Default value is None.
        **gas_transmittance_table (gas_transmittance.Gas_Transmittance_Lookup_Table)**: An instance of the Gas_Transmittance_Lookup_Table class (available in the gas corrections library), which contains transmittance lookup tables read in from a NetCDF file. 
            Default value is None.

    Returns:
        Returns a dataclass with three members: (i) solar_zenith, (ii) sensor_zenith, and (iii) total. These contain arrays of 
        the transmittance values computed along the slant paths at the solar and sensor zenith angles. The total
        transmittance is the product of the solar zenith and sensor zenith transmittances. Each element of these matrices is 
        a transmittance value corresponding to the pixel at the same index in the input L1 Data. For example, sensor_zenith[10, 50]
        is the sensor zenith transmittance corresponding to the pixel located at l1_data[10, 50].

    Example Usage:
        Read in solar and sensor zenith arrays from L1 input file and store them in l1b_solz and l1b_senz.

        l1_data = gas_corrections.L1_Data()
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(l1b_senz))
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(l1b_solz))
        l1_data.num_pixels = len(l1_data.cos_solar_zenith)
        l1_data.reflectance = l1b_reflectance
        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_wavelengths = len(sensor_wavelengths)

        gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table(gas_transmittance_filepath)

        gas_transmittances = gas_corrections.h2o_transmittance(l1_data=l1_data, \
                                                              gas_transmittance_table=gas_transmittance_table)
    """
    valid_kwargs = ['l1_data', 'gas_transmittance_table']
    check_for_invalid_keyword_args(kwargs, valid_kwargs)
    
    args = dict()
    l1_data = args['l1_data'] = kwargs.get("l1_data", None)
    gas_transmittance_table = args['gas_transmittance_table'] = kwargs.get("gas_transmittance_table", None)
    lookup_table_has_amf_dimension = bool(gas_transmittance_table.num_amf_grid_points)

    validate_keyword_args(args)

    ancillary_data = gas_transmittance.Ancillary_Data()

    ancillary_data.precipitable_water = np.zeros(l1_data.cos_solar_zenith.size)
    ancillary_data.water_vapor_bands = np.array([782, 817, 857], dtype=np.float64)
    ancillary_data.num_water_vapor_bands = ancillary_data.water_vapor_bands.size

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.h2o_transmittance, axis = 1)
    h2o_transmittance_at_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.num_wavelengths = len(l1_data.wavelengths)
    gas_transmittance_table.h2o_transmittance = h2o_transmittance_at_sensor_wavelengths

    assert(l1_data.cos_solar_zenith.size != 0)
    assert(l1_data.cos_sensor_zenith.size != 0)
    assert(l1_data.wavelengths.size != 0)
    assert(gas_transmittance_table.h2o_transmittance.size != 0)

    assert(l1_data.cos_solar_zenith.shape == l1_data.cos_sensor_zenith.shape)
    assert(l1_data.num_pixels == l1_data.cos_solar_zenith.size)
    assert(l1_data.num_pixels == l1_data.cos_sensor_zenith.size)
    assert(l1_data.num_wavelengths == l1_data.wavelengths.size)
    assert(l1_data.num_wavelengths == gas_transmittance_table.h2o_transmittance.shape[1])    

    return gas_transmittance.h2o_transmittance(l1_data, ancillary_data, gas_transmittance_table, lookup_table_has_amf_dimension)