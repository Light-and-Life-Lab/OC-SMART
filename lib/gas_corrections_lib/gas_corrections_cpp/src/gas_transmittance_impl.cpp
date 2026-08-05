#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <omp.h>
#include <iostream>

#include "gas_transmittance.h"

int32_t get_index_lowerbound(double* table_val, int num_val, float val) 
{
    int32_t index;
    int32_t i;

    for (i = 0; i < num_val; i++)
        if (val < table_val[i])
            break;
    index = std::max(i-1, 0);
    index = std::min(index, num_val-2);

    return index;
}


int32_t get_index_upperbound(float *table_val, int32_t num_val, float val) {
    int32_t index;

    for (index = 0; index < num_val; index++)
        if (val >= table_val[index])
            break;

    return index;
}


std::pair<int, double> get_amf_index_and_ratio(Gas_Transmittance_Lookup_Table* gas_transmittance_table, double amf_value)
{
    int index_amf = get_index_lowerbound(gas_transmittance_table->air_mass_factor_mixed_gases, gas_transmittance_table->num_amf_grid_points, amf_value);
    double ratio = (amf_value - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf]) /
                    (gas_transmittance_table->air_mass_factor_mixed_gases[index_amf + 1] - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf]);

    return std::make_pair(index_amf, ratio);
}


double interpolate_transmittance_to_amf(double* transmittance_table, int32_t index, double ratio)
{
    // In the case where amf correction is performed, the gas transmittance input file will 
    // have transmittance values as a 2D matrix. The transmittance values in this matrix are
    // a function of both wavelength and air mass factor (amf)

    // This desired amf value may fall between two of the points in the amf grid, so in general
    // it is necessary to interpolate the values in the transmittance table to obtain the 
    // transmittance at the desired amf value

    double transmittance_interpolated_to_amf = transmittance_table[index]*(1-ratio)
                                             + transmittance_table[index+1]*ratio;
    return transmittance_interpolated_to_amf;
}


int windex(float wave, double twave[], int ntwave) {
    int iw, index;
    double wdiff;
    double wdiffmin = 99999.;

    for (iw = 0; iw < ntwave; iw++) {

        /* break on exact match */
        if (twave[iw] == wave) {
            index = iw;
            break;
        }

        /* look for closest */
        wdiff = fabs(twave[iw] - wave);
        if (wdiff < wdiffmin) {
            wdiffmin = wdiff;
            index = iw;
        }
    }

    return (index);
}


float get_airmass_oxygen(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, int32_t ip, double window1, double absorp_band, double window2)
{
    int32_t i;

    int32_t num_wavelengths = l1_data->num_wavelengths;
    int32_t row_offset = ip*num_wavelengths;
    double* wavelength_array = l1_data->wavelengths;
    double reflectances[3];

    int absorption_window_lower_wavelength_index = windex(window1, wavelength_array, num_wavelengths);
    reflectances[0] = l1_data->reflectance[absorption_window_lower_wavelength_index];

    int absorption_window_upper_wavelength_index = windex(window2, wavelength_array, num_wavelengths);
    reflectances[1] = l1_data->reflectance[absorption_window_upper_wavelength_index];

    int band_absorp = windex(absorp_band, wavelength_array, num_wavelengths);
    reflectances[2] = l1_data->reflectance[band_absorp];

    double reflectances_interpolated = reflectances[0]+(absorp_band-window1)*(reflectances[1]-reflectances[0])/(window2-window1);

    double trans_o2_true = reflectances[2]/reflectances_interpolated;

    int num_airmass = gas_transmittance_table->num_amf_grid_points;
    int gas_transmittance_table_row_offset = band_absorp*num_airmass;
    
    for (i = 0; i < num_airmass; i++) 
    {
        if (trans_o2_true >= gas_transmittance_table->o2_transmittance[gas_transmittance_table_row_offset + i])
            break;
    }
    if (i == 0)
        i = 1;
    if (i == num_airmass)
        i = num_airmass - 1;

    double amf_interp = gas_transmittance_table->air_mass_factor_mixed_gases[i] + (trans_o2_true - gas_transmittance_table->o2_transmittance[gas_transmittance_table_row_offset+i]) 
                 * (gas_transmittance_table->air_mass_factor_mixed_gases[i] - gas_transmittance_table->air_mass_factor_mixed_gases[i - 1]) 
                 / (gas_transmittance_table->o2_transmittance[gas_transmittance_table_row_offset+i] - gas_transmittance_table->o2_transmittance[gas_transmittance_table_row_offset+i - 1]);

    return (amf_interp);
}


float get_wv_band_ratio(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, int32_t ip, float window1, float absorp_band, float window2, double amf_total)
{
    std::vector<float> interpolated_transmittances(gas_transmittance_table->num_water_vapor_concentrations);

    // As far as I can tell, what this function does is the following:
    // Take as input 3 wavelength values. 
    // These values are a water vapor absorption band (e.g. 940 nm) and a pair of wavelengths on either side of the absorption window (e.g. 880 nm and 1038 nm)
    // The ToA radiance (Lt) values at these wavelengths (from the L1 file) are then converted to reflectances (rhot) at the same wavelengths
    // Then a linear interpolation is performed beteween the reflectances at the two endpoints
    // The absorption band reflectance is then divided by this interpolated reflectance
    // The idea here appears to be that the endpoint reflectances are outside of the absorption band, and so will have transmittance close to 1.0
    // The interpolated reflectance will therefore also have transmittance close to 1.0
    // The absorption band reflectance will be significantly smaller because it is affected by the absorption window
    // Therefore, if you think of the absorption band reflectance as lying in a "valley", then the interpolated reflectance would be "on a bridge over the valley"
    // The absorption band reflectance divided by the interpolated reflectance therefore gives the "depth of the valley", i.e. the transmittance at the absorption band
    // This is a way of computing the "true" transmittance value at the absorption band using the actual data, and does not require any a priori theoretical knowledge of water vapor transmittance

    int32_t num_L1_wavelengths = l1_data->num_wavelengths;
    int32_t row_offset = ip*num_L1_wavelengths;
    double* wavelength_array = l1_data->wavelengths;
    double reflectances[3];

    // derive a transmittance using a line height (or in this case, depth) approach
    int absorption_window_lower_wavelength_index = windex(window1, wavelength_array, num_L1_wavelengths);
    reflectances[0] = l1_data->reflectance[row_offset + absorption_window_lower_wavelength_index];

    int absorption_window_upper_wavelength_index = windex(window2, wavelength_array, num_L1_wavelengths);
    reflectances[1] = l1_data->reflectance[row_offset + absorption_window_upper_wavelength_index];

    int absorption_band_index = windex(absorp_band, wavelength_array, num_L1_wavelengths);
    reflectances[2] =l1_data->reflectance[row_offset + absorption_band_index];

    double interpolated_reflectance = reflectances[0] + ((absorp_band - window1) / (window2 - window1)) * (reflectances[1] - reflectances[0]);
    double true_water_vapor_transmittance = reflectances[2] / interpolated_reflectance;

    // Once the "true" water vapor transmittance value at the absorption band is calculated, a lookup table of water vapor values is used
    // This lookup table is from the amf NetCDF file corresponding to the sensor of interest (e.g. for OCI the file oci_gas_transmittance_cia_amf_v3.2.nc is used)
    // The lookup table has water vapor transmittance values at various wavelength, water vapor concentration, and amf values
    // 

    // For the given absorption band and pixel air mass factor, interpolate the water vapor transmittance table for the
    // each tabular water vapor

    double* amf_wv = gas_transmittance_table->air_mass_factor_water_vapor;
    int model = gas_transmittance_table->model;
    int num_gas_transmittance_wavelengths = gas_transmittance_table->num_wavelengths;
    int num_airmass = gas_transmittance_table->num_amf_grid_points;
    int num_water_vapors = gas_transmittance_table->num_water_vapor_concentrations;

    
    int amf_index = get_index_lowerbound(amf_wv, num_airmass, amf_total);
    double amf_ratio = (amf_total - amf_wv[amf_index]) / (amf_wv[amf_index + 1] - amf_wv[amf_index]);
    int water_vapor_transmittance_table_index = (model * num_gas_transmittance_wavelengths * num_airmass * num_water_vapors) +
             (absorption_band_index * num_airmass * num_water_vapors) +
             (amf_index * num_water_vapors);

    for (int i = 0; i < num_water_vapors; i++) {
        interpolated_transmittances[i] = gas_transmittance_table->h2o_transmittance[water_vapor_transmittance_table_index + i] * (1 - amf_ratio) +
                            gas_transmittance_table->h2o_transmittance[water_vapor_transmittance_table_index + num_water_vapors + i] *  amf_ratio;
    }

    // Find the bounding transmittance index matching the "true" (computed) transmittance
    int columnar_water_vapor_table_index = get_index_upperbound(interpolated_transmittances.data(), num_water_vapors, true_water_vapor_transmittance);

    // retrieve water vapor by interpolating the tabular column water vapor assocaited with the "true" transmittance
    double columnar_water_vapor_interpolated_to_true_transmittance = 
    gas_transmittance_table->water_vapor_concentration[columnar_water_vapor_table_index] 
    + (true_water_vapor_transmittance - interpolated_transmittances[columnar_water_vapor_table_index]) * (gas_transmittance_table->water_vapor_concentration[columnar_water_vapor_table_index] - gas_transmittance_table->water_vapor_concentration[columnar_water_vapor_table_index - 1]) 
    / (interpolated_transmittances[columnar_water_vapor_table_index] - interpolated_transmittances[columnar_water_vapor_table_index - 1]);

    return (columnar_water_vapor_interpolated_to_true_transmittance);
}


void ozone_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittances* gas_transmittances)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            double ozone_optical_depth = ancillary_data->ozone_concentration[ip] * ancillary_data->ozone_absorption_cross_section[iw];
            gas_transmittances->solar_zenith[row_offset + iw] = exp(-(ozone_optical_depth / l1_data->cos_solar_zenith[ip]));

            // if (lookup_table_has_amf_dimension) 
            // {
            //     gas_transmittances->total[row_offset + iw] = exp(-ozone_optical_depth * (1.0/l1_data->cos_solar_zenith[ip] + 1.0/l1_data->cos_sensor_zenith[ip]));
            //     gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            // } 
            // else 
            // {
                gas_transmittances->sensor_zenith[row_offset + iw] = exp(-(ozone_optical_depth / l1_data->cos_sensor_zenith[ip]));
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            // }
        }
    }
}


void co2_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(gas_transmittance_table, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(gas_transmittance_table, amf_total);

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            if (lookup_table_has_amf_dimension)
            {
                int32_t row_index = iw*gas_transmittance_table->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                gas_transmittances->solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->co2_transmittance, table_index_solz, ratio_solz);
                gas_transmittances->total[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->co2_transmittance, table_index_total, ratio_total);
                gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            }
            else
            {
                gas_transmittances->solar_zenith[row_offset + iw] = pow(gas_transmittance_table->co2_transmittance[iw], amf_solar_zenith);
                gas_transmittances->sensor_zenith[row_offset + iw] = pow(gas_transmittance_table->co2_transmittance[iw], amf_sensor_zenith);
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            }
        }
    }
}


void co_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension) 
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(gas_transmittance_table, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(gas_transmittance_table, amf_total);

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            if (lookup_table_has_amf_dimension)
            {
                int32_t row_index = iw*gas_transmittance_table->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                gas_transmittances->solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->co_transmittance, table_index_solz, ratio_solz);
                gas_transmittances->total[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->co_transmittance, table_index_total, ratio_total);
                gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            }
            else
            {
                gas_transmittances->solar_zenith[row_offset + iw] = pow(gas_transmittance_table->co_transmittance[iw], amf_solar_zenith);
                gas_transmittances->sensor_zenith[row_offset + iw] = pow(gas_transmittance_table->co_transmittance[iw], amf_sensor_zenith);
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            }
        }
    }
}

void ch4_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(gas_transmittance_table, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(gas_transmittance_table, amf_total);

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            if (lookup_table_has_amf_dimension)
            {
                int32_t row_index = iw*gas_transmittance_table->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                gas_transmittances->solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->ch4_transmittance, table_index_solz, ratio_solz);
                gas_transmittances->total[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->ch4_transmittance, table_index_total, ratio_total);
                gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            }
            else
            {
                gas_transmittances->solar_zenith[row_offset + iw] = pow(gas_transmittance_table->ch4_transmittance[iw], amf_solar_zenith);
                gas_transmittances->sensor_zenith[row_offset + iw] = pow(gas_transmittance_table->ch4_transmittance[iw], amf_sensor_zenith);
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            }
        }
    }
}


void o2_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension, Oxygen_A_Band_Option oxygen_A_band_option)  
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int32_t nwave = l1_data->num_wavelengths;
        int32_t ipb = ip*nwave;

        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        int32_t index_amf_solz;
        int32_t index_amf_total;

        float ratio_solz;
        float ratio_total;

        if (lookup_table_has_amf_dimension) 
        {
            index_amf_solz = get_index_lowerbound(gas_transmittance_table->air_mass_factor_mixed_gases, gas_transmittance_table->num_amf_grid_points, amf_solar_zenith);
            index_amf_total = get_index_lowerbound(gas_transmittance_table->air_mass_factor_mixed_gases, gas_transmittance_table->num_amf_grid_points, amf_total);

            ratio_solz = (amf_solar_zenith - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz]) /
                            (gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz + 1] - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz]);
            ratio_total = (amf_total - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total]) /
                            (gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total + 1] - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total]);
        }

        int index_amf_solz_o2;
        float ratio_solz_o2;
        float ratio_total_o2;
        int index_amf_total_o2;

        if (lookup_table_has_amf_dimension && oxygen_A_band_option == Oxygen_A_Band_Option::SURROUNDING_WINDOW_BANDS) 
        {
            float amf_total_o2 = get_airmass_oxygen(l1_data, gas_transmittance_table, ip, 753.0221, 761.7891, 776.81335);
            float scaling_factor = amf_total_o2 / amf_total;

            index_amf_solz_o2 = get_index_lowerbound(gas_transmittance_table->air_mass_factor_mixed_gases, gas_transmittance_table->num_amf_grid_points, amf_solar_zenith * scaling_factor);
            index_amf_total_o2 = get_index_lowerbound(gas_transmittance_table->air_mass_factor_mixed_gases, gas_transmittance_table->num_amf_grid_points, amf_total * scaling_factor);

            ratio_solz_o2 = (amf_solar_zenith * scaling_factor - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz_o2]) /
                            (gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz_o2 + 1] - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_solz_o2]);
            ratio_total_o2 = (amf_total * scaling_factor - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total_o2]) /
                            (gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total_o2 + 1] - gas_transmittance_table->air_mass_factor_mixed_gases[index_amf_total_o2]);
        }

        double* t_o2 = gas_transmittance_table->o2_transmittance;

        for (int iw = 0; iw < nwave; iw++) {
            if (lookup_table_has_amf_dimension) 
            {
                int32_t index=iw*gas_transmittance_table->num_amf_grid_points;
                float t_o2_interp;

                if (oxygen_A_band_option == Oxygen_A_Band_Option::SURROUNDING_WINDOW_BANDS) 
                {
                    t_o2_interp = t_o2[index + index_amf_solz_o2] * (1 - ratio_solz_o2) +
                                t_o2[index + index_amf_solz_o2 + 1] * ratio_solz_o2;
                    gas_transmittances->solar_zenith[ipb + iw] = t_o2_interp;

                    t_o2_interp = t_o2[index + index_amf_total_o2] * (1 - ratio_total_o2) +
                                t_o2[index + index_amf_total_o2 + 1] * ratio_total_o2;
                    gas_transmittances->total[ipb + iw] = t_o2_interp;
                }
                else
                {
                    t_o2_interp = t_o2[index + index_amf_solz] * (1 - ratio_solz) +
                                t_o2[index + index_amf_solz + 1] * ratio_solz;
                    gas_transmittances->solar_zenith[ipb + iw] = t_o2_interp;

                    t_o2_interp = t_o2[index + index_amf_total] * (1 - ratio_total) +
                                t_o2[index + index_amf_total + 1] * ratio_total;
                    gas_transmittances->total[ipb + iw] = t_o2_interp;
                }

                gas_transmittances->sensor_zenith[ipb + iw] = gas_transmittances->total[ipb + iw] / gas_transmittances->solar_zenith[ipb + iw];
            }
            else
            {
                int t_o2_at_amf_1_index = iw*gas_transmittance_table->num_amf_grid_points;
                
                gas_transmittances->solar_zenith[ipb + iw] = pow(t_o2[t_o2_at_amf_1_index], amf_solar_zenith);
                gas_transmittances->sensor_zenith[ipb + iw] = pow(t_o2[t_o2_at_amf_1_index], amf_sensor_zenith);
                gas_transmittances->total[ipb + iw] = gas_transmittances->sensor_zenith[ipb + iw] * gas_transmittances->solar_zenith[ipb + iw];
            }
        }
    }
}


void n2o_transmittance(L1_Data* l1_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(gas_transmittance_table, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(gas_transmittance_table, amf_total);

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            if (lookup_table_has_amf_dimension)
            {
                int32_t row_index = iw*gas_transmittance_table->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                gas_transmittances->solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->n2o_transmittance, table_index_solz, ratio_solz);
                gas_transmittances->total[row_offset + iw] = interpolate_transmittance_to_amf(gas_transmittance_table->n2o_transmittance, table_index_total, ratio_total);
                gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            }
            else
            {
                gas_transmittances->solar_zenith[row_offset + iw] = pow(gas_transmittance_table->n2o_transmittance[iw], amf_solar_zenith);
                gas_transmittances->sensor_zenith[row_offset + iw] = pow(gas_transmittance_table->n2o_transmittance[iw], amf_sensor_zenith);
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            }
        }
    }
}

void no2_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittances* gas_transmittances)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        double sec0 = 1.0 / l1_data->cos_solar_zenith[ip];
        double sec = 1.0 / l1_data->cos_sensor_zenith[ip];
        double tropospheric_no2_concentration_above_200m{0.0};

        if (ancillary_data->tropospheric_no2_concentration[ip] > 0.0)
        {
            /* compute tropo no2 above 200m (Z.Ahmad)
            tropospheric_no2_concentration_above_200m = exp(12.6615 + 0.61676*log(no2_tropo));
            new, location-dependent method */
            tropospheric_no2_concentration_above_200m = ancillary_data->fraction_tropospheric_no2_above_200m[ip] * ancillary_data->tropospheric_no2_concentration[ip];
        }

        int row_offset = ip*l1_data->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            if (ancillary_data->no2_absorption_cross_section[iw] > 0.0) 
            {
                double a_285 = ancillary_data->no2_absorption_cross_section[iw] * (1.0 - 0.003 * (285.0 - 294.0));
                double a_225 = ancillary_data->no2_absorption_cross_section[iw] * (1.0 - 0.003 * (225.0 - 294.0));

                double no2_optical_depth_to_200m = a_285 * tropospheric_no2_concentration_above_200m 
                                                 + a_225 * ancillary_data->stratospheric_no2_concentration[ip];

                gas_transmittances->solar_zenith[row_offset + iw] = exp(-(no2_optical_depth_to_200m * sec0));

                // if (lookup_table_has_amf_dimension) 
                // {
                //     gas_transmittances->total[row_offset + iw] = exp(-(no2_optical_depth_to_200m * (sec + sec0)));
                //     gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
                // }
                // else
                // {
                    gas_transmittances->sensor_zenith[row_offset + iw] = exp(-(no2_optical_depth_to_200m * sec));
                    gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
                // }
            }
        }
    }
}




void h2o_transmittance(L1_Data* l1_data, Ancillary_Data* ancillary_data, Gas_Transmittance_Lookup_Table* gas_transmittance_table, Gas_Transmittances* gas_transmittances, bool lookup_table_has_amf_dimension) 
{
    // #pragma omp parallel for
    for (int ip = 0; ip < l1_data->num_pixels; ip++)
    {
        double amf_solar_zenith = 1.0/l1_data->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_data->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        // wv variable is a table of values quantifying the water vapor present in each pixel of the image
        // amf table has table of columnar (i.e. integrated) water vapor in cm, whereas MET files have total precipitable water vapor in kg/m^2
        // Note: from this source (https://remss.com/measurements/atmospheric-water-vapor/) we have that kg/m2 = 1 mm due to density of water being 1000 kg/m^3


        // watervapor_bands is an array of wavelengths specified in the mls12 input config file for OCSSW
        // For example, OCI has several input config files:
        //      msl12_defaults.par: watervapor_bands=[679,719,749,782,817,859]
        //      msl12_defaults_SFREFL.par: watervapor_bands=[880,940,1038,880,940,1038]
        //      msl12_defaults_LANDVI.par: watervapor_bands=[880,940,1038,880,940,1038]
        //      msl12_defaults_CLD.par: watervapor_bands=0

        if (lookup_table_has_amf_dimension && ancillary_data->water_vapor_bands) 
        {
            double total_columnar_water_vapor = 0;
            for (int iw = 0; iw < ancillary_data->num_water_vapor_bands;) 
            {
                total_columnar_water_vapor += get_wv_band_ratio(l1_data, gas_transmittance_table, ip, ancillary_data->water_vapor_bands[iw], ancillary_data->water_vapor_bands[iw + 1],
                                        ancillary_data->water_vapor_bands[iw + 2], amf_total);
                iw += 3;
            }
            total_columnar_water_vapor /= (ancillary_data->num_water_vapor_bands / 3);
            ancillary_data->precipitable_water[ip] = total_columnar_water_vapor;
        }


        float tempratio;

        int index_amf_wv_solz{};
        int index_amf_wv_total{};
        double ratio_amf_solz{};
        double ratio_amf_total{};

        if (lookup_table_has_amf_dimension) 
        {
            index_amf_wv_solz = get_index_lowerbound(gas_transmittance_table->air_mass_factor_water_vapor, gas_transmittance_table->num_amf_grid_points, amf_solar_zenith);
            index_amf_wv_total = get_index_lowerbound(gas_transmittance_table->air_mass_factor_water_vapor, gas_transmittance_table->num_amf_grid_points, amf_total);

            ratio_amf_solz = (amf_solar_zenith - gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_solz]) /
                            (gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_solz + 1] - gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_solz]);
            ratio_amf_total = (amf_total - gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_total]) /
                            (gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_total + 1] - gas_transmittance_table->air_mass_factor_water_vapor[index_amf_wv_total]);
        }

        const double wv = ancillary_data->precipitable_water[ip];

        int ja = get_index_lowerbound(gas_transmittance_table->water_vapor_concentration, gas_transmittance_table->num_water_vapor_concentrations, wv );
        int ja_sen = get_index_lowerbound(gas_transmittance_table->water_vapor_concentration, gas_transmittance_table->num_water_vapor_concentrations, wv*amf_sensor_zenith );
        int ja_sol = get_index_lowerbound(gas_transmittance_table->water_vapor_concentration, gas_transmittance_table->num_water_vapor_concentrations, wv*amf_solar_zenith );

        double water_vapor_concentration_interpolated = (wv - gas_transmittance_table->water_vapor_concentration[ja])/(gas_transmittance_table->water_vapor_concentration[ja+1] - gas_transmittance_table->water_vapor_concentration[ja]);
        

        for (int iw = 0; iw < l1_data->num_wavelengths; iw++) 
        {
            int32_t row_offset = ip*l1_data->num_wavelengths;

            if (lookup_table_has_amf_dimension) 
            {
                // (model, num_wavelengths, num_airmass, numwatervapors) is the set of dimensions of the
                // gas transmittance table NetCDF file (aka the number of wavelength rows in the water vapor transmittance table)
                // (nmodels, nwavelengths, n_air_mass_factor, and n_water_vapor) respectively
                int index = gas_transmittance_table->model
                            *gas_transmittance_table->num_wavelengths
                            *gas_transmittance_table->num_amf_grid_points
                            *gas_transmittance_table->num_water_vapor_concentrations
                            + iw*gas_transmittance_table->num_amf_grid_points*gas_transmittance_table->num_water_vapor_concentrations;

                double f00 = gas_transmittance_table->h2o_transmittance[index + index_amf_wv_solz*gas_transmittance_table->num_water_vapor_concentrations + ja];
                double f10 = gas_transmittance_table->h2o_transmittance[index + (index_amf_wv_solz+1)*gas_transmittance_table->num_water_vapor_concentrations + ja];
                double f01 = gas_transmittance_table->h2o_transmittance[index + index_amf_wv_solz*gas_transmittance_table->num_water_vapor_concentrations + ja + 1];
                double f11 = gas_transmittance_table->h2o_transmittance[index + (index_amf_wv_solz+1)*gas_transmittance_table->num_water_vapor_concentrations + ja + 1];

                double water_vapor_transmittance_solar_zenith = (1. - ratio_amf_solz)*(1. - water_vapor_concentration_interpolated) * f00 
                                                                + ratio_amf_solz * water_vapor_concentration_interpolated * f11 
                                                                + ratio_amf_solz * (1. - water_vapor_concentration_interpolated) * f10 
                                                                + water_vapor_concentration_interpolated * (1. - ratio_amf_solz) * f01;

                gas_transmittances->solar_zenith[row_offset + iw] = water_vapor_transmittance_solar_zenith;

                f00 = gas_transmittance_table->h2o_transmittance[index + index_amf_wv_total*gas_transmittance_table->num_water_vapor_concentrations + ja];
                f10 = gas_transmittance_table->h2o_transmittance[index + (index_amf_wv_total+1)*gas_transmittance_table->num_water_vapor_concentrations + ja];
                f01 = gas_transmittance_table->h2o_transmittance[index + index_amf_wv_total*gas_transmittance_table->num_water_vapor_concentrations + ja + 1];
                f11 = gas_transmittance_table->h2o_transmittance[index + (index_amf_wv_total+1)*gas_transmittance_table->num_water_vapor_concentrations + ja + 1];

                double water_vapor_transmittance_total = (1. - ratio_amf_total)*(1. - water_vapor_concentration_interpolated) * f00 
                                                        + ratio_amf_total * water_vapor_concentration_interpolated * f11 
                                                        + ratio_amf_total * (1. - water_vapor_concentration_interpolated) * f10 
                                                        + water_vapor_concentration_interpolated * (1. - ratio_amf_total) * f01;
                gas_transmittances->total[row_offset + iw] = water_vapor_transmittance_total;
                gas_transmittances->sensor_zenith[row_offset + iw] = gas_transmittances->total[row_offset + iw] / gas_transmittances->solar_zenith[row_offset + iw];
            }
            else
            {
                int index = gas_transmittance_table->model
                            *gas_transmittance_table->num_wavelengths
                            *gas_transmittance_table->num_water_vapor_concentrations
                            + iw*gas_transmittance_table->num_water_vapor_concentrations;

                std::cout << "gas_transmittance_table->model: " << gas_transmittance_table->model << std::endl;
                std::cout << "gas_transmittance_table->num_wavelengths: " << gas_transmittance_table->num_wavelengths << std::endl;
                std::cout << "gas_transmittance_table->num_water_vapor_concentrations: " << gas_transmittance_table->num_water_vapor_concentrations << std::endl;
                std::cout << "gas_transmittance_table->num_water_vapor_concentrations: " << gas_transmittance_table->num_water_vapor_concentrations << std::endl;

                tempratio = (wv*amf_solar_zenith - gas_transmittance_table->water_vapor_concentration[ja_sol])/(gas_transmittance_table->water_vapor_concentration[ja_sol+1] - gas_transmittance_table->water_vapor_concentration[ja_sol]);
                double water_vapor_transmittance_solar_zenith = gas_transmittance_table->h2o_transmittance[index+ja_sol]*(1-tempratio) + gas_transmittance_table->h2o_transmittance[index+ja_sol+1]*tempratio;

                std::cout << "wv: " << wv << std::endl;
                std::cout << "amf_solar_zenith: " << amf_solar_zenith << std::endl;

                std::cout << "tempratio numerator: " << (wv*amf_solar_zenith - gas_transmittance_table->water_vapor_concentration[ja_sol]) << std::endl;
                std::cout << "tempration denom: " << (gas_transmittance_table->water_vapor_concentration[ja_sol+1] - gas_transmittance_table->water_vapor_concentration[ja_sol]) << std::endl;
                


                std::cout << "tempratio: " << tempratio << std::endl;

                std::cout << "gas_transmittance_table->water_vapor_concentration[ja_sol]: " << gas_transmittance_table->water_vapor_concentration[ja_sol] << std::endl;
                std::cout << "gas_transmittance_table->water_vapor_concentration[ja_sol+1]: " << gas_transmittance_table->water_vapor_concentration[ja_sol+1] << std::endl;
                
                std::cout << "index: " << index << std::endl;

                std::cout << "gas_transmittance_table->h2o_transmittance[index+ja_sol]: " << gas_transmittance_table->h2o_transmittance[index+ja_sol] << std::endl;
                std::cout << "gas_transmittance_table->h2o_transmittance[index+ja_sol+1]: " << gas_transmittance_table->h2o_transmittance[index+ja_sol+1] << std::endl;
                
                std::cout << "water_vapor_transmittance_solar_zenith: " << water_vapor_transmittance_solar_zenith << std::endl;
                gas_transmittances->solar_zenith[row_offset + iw] = water_vapor_transmittance_solar_zenith;

                tempratio = (wv*amf_sensor_zenith -gas_transmittance_table->water_vapor_concentration[ja_sen])/(gas_transmittance_table->water_vapor_concentration[ja_sen+1]-gas_transmittance_table->water_vapor_concentration[ja_sen]);
                double water_vapor_transmittance_sensor_zenith = gas_transmittance_table->h2o_transmittance[index+ja_sen]*(1-tempratio) + gas_transmittance_table->h2o_transmittance[index+ja_sen+1]*tempratio;
                gas_transmittances->sensor_zenith[row_offset + iw] = water_vapor_transmittance_sensor_zenith;
                gas_transmittances->total[row_offset + iw] = gas_transmittances->sensor_zenith[row_offset + iw] * gas_transmittances->solar_zenith[row_offset + iw];
            }
        }
    }
}