#ifndef GAS_TRANSMITTANCE_H
#define GAS_TRANSMITTANCE_H

enum Atmosphere_Model
{
    TROPICAL = 0,
    MID_LATITUDE_SUMMER = 1,
    MID_LATITUDE_WINTER = 2,
    SUBARCTIC_SUMMER = 3,
    SUBARCTIC_WINTER = 4,
    US_STANDARD_62 = 5,
};

struct Ancillary_Data
{
    // See https://oceancolor.gsfc.nasa.gov/resources/docs/ancillary/ for sources of Ancillary Data
    // The algorithms for retrieval of ocean color radiometry and sea surface temperature from satellite radiometry make use of a number of ancillary data sources, 
    // in addition to the sensor observed radiometry. Examples include meteorological data (windspeed, surface pressure, relative humidity) 
    // and concentrations of atmospheric gases (water vapor, ozone, nitrogen dioxide). 
    // In general, the optimal ancillary data is not available when the satellite radiometry is first acquired, so the NASA OBPG processing uses best available ancillary data sources 
    // for processing satellite data in near realtime, and then reprocesses the data about one week later to refine the products using the optimal ancillary data.

    // From e.g. ozone_climatology_v2014.hdf
    double* ozone_absorption_cross_section{}; // k_oz
    double* ozone_concentration{};

    // From e.g. no2_climatology_v2013.hdf and trop_f_no2_200m.hdf
    double* no2_absorption_cross_section{}; // k_no2
    double* fraction_tropospheric_no2_above_200m{};
    double* tropospheric_no2_concentration{};
    double* stratospheric_no2_concentration{};

    // From e.g. GMAO_MERRA2.20240411T180000.MET.nc
    double* precipitable_water{};

    double* water_vapor_bands{};
    int num_water_vapor_bands{};
};


struct Gas_Transmittance_Lookup_Table
{
    // From e.g. oci_gas_transmittance_cia_amf_v3.2.nc
    double* co2_transmittance{};
    double* co_transmittance{};
    double* ch4_transmittance{};
    double* o2_transmittance{};
    double* n2o_transmittance{};
    double* h2o_transmittance{};

    int model{};
    double* wavelengths{};
    double* air_mass_factor_mixed_gases{};
    double* air_mass_factor_water_vapor{};
    double* water_vapor_concentration{};
    int num_models{};
    int num_wavelengths{};
    int num_amf_grid_points{}; // Length of both mixed gases and water vapor air mass factor tables
    int num_water_vapor_concentrations{};
};

struct L1_Data
{
    double* reflectance{};

    double* cos_solar_zenith{}; 
    double* cos_sensor_zenith{};

    double* latitude{};
    double* longitude{};
    
    double* wavelengths{};

    int num_pixels{};
    int num_wavelengths{};

};

struct Gas_Transmittances
{
    double* solar_zenith{};
    double* sensor_zenith{};
    double* total{};
};

enum Oxygen_A_Band_Option
{
    // Enum values are set to match those found in OCSSW
    TRANSMITTANCE_TABLE = 2, // Apply oxygen transmittance from gas transmittance table
    SURROUNDING_WINDOW_BANDS = 3, // Compute oxygen transmittance from A-band and surrounding window bands. Requires AMF gas trasmittance table.
};

void ozone_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittances* gas_transmittances);
void co2_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension);
void co_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension);
void ch4_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension);
void o2_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension, Oxygen_A_Band_Option oxygen_A_band_option);
void n2o_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension);
void no2_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittances* gas_transmittances);
void h2o_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension);

#endif // GAS_TRANSMITTANCE_H