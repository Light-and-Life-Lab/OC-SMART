<%
cfg['compiler_args'] = ['-std=c++2a', '-O3', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
cfg['sources'] = ['gas_transmittance_impl.cpp']
setup_pybind11(cfg)
%>

#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "allocate_output_array.hpp"
#include "pybind_interface_types.hpp"
#include "gas_transmittance.h"

namespace py = pybind11;

Gas_Transmittances_PY ozone_transmittance(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data) 
{
    Ancillary_Data ancillary_data_c{};
    
    ancillary_data_c.ozone_absorption_cross_section = static_cast<double*>(ancillary_data.ozone_absorption_cross_section.request().ptr);
    ancillary_data_c.ozone_concentration = static_cast<double*>(ancillary_data.ozone_concentration.request().ptr);
    
    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    ozone_transmittance(&l1_data_c, &ancillary_data_c, &gas_transmittances_c);

    return gas_transmittances;
}


Gas_Transmittances_PY co2_transmittance(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.co2_transmittance = static_cast<double*>(gas_transmittance_table.co2_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    co2_transmittance(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    return gas_transmittances;
}


Gas_Transmittances_PY co_transmittance(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.co_transmittance = static_cast<double*>(gas_transmittance_table.co_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    
    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    co_transmittance(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    return gas_transmittances;
}


Gas_Transmittances_PY ch4_transmittance(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.ch4_transmittance = static_cast<double*>(gas_transmittance_table.ch4_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    ch4_transmittance(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    return gas_transmittances;
}


Gas_Transmittances_PY o2_transmittance(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension, Oxygen_A_Band_Option oxygen_A_band_option) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.o2_transmittance = static_cast<double*>(gas_transmittance_table.o2_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    
    L1_Data l1_data_c{};

    l1_data_c.reflectance = static_cast<double*>(l1_data.reflectance.request().ptr);
    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;
    l1_data_c.wavelengths = static_cast<double*>(l1_data.wavelengths.request().ptr);

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    o2_transmittance(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension, oxygen_A_band_option);

    return gas_transmittances;
}


Gas_Transmittances_PY n2o_transmittance(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.n2o_transmittance = static_cast<double*>(gas_transmittance_table.n2o_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    
    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    n2o_transmittance(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    return gas_transmittances;
}


Gas_Transmittances_PY no2_transmittance(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data) 
{
    Ancillary_Data ancillary_data_c{};

    ancillary_data_c.no2_absorption_cross_section = static_cast<double*>(ancillary_data.no2_absorption_cross_section.request().ptr);
    ancillary_data_c.fraction_tropospheric_no2_above_200m = static_cast<double*>(ancillary_data.fraction_tropospheric_no2_above_200m.request().ptr);
    ancillary_data_c.tropospheric_no2_concentration = static_cast<double*>(ancillary_data.tropospheric_no2_concentration.request().ptr);
    ancillary_data_c.stratospheric_no2_concentration = static_cast<double*>(ancillary_data.stratospheric_no2_concentration.request().ptr);
    
    L1_Data l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    no2_transmittance(&l1_data_c, &ancillary_data_c, &gas_transmittances_c);

    return gas_transmittances;
}


Gas_Transmittances_PY h2o_transmittance(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Ancillary_Data ancillary_data_c{};

    ancillary_data_c.no2_absorption_cross_section = static_cast<double*>(ancillary_data.no2_absorption_cross_section.request().ptr);
    ancillary_data_c.fraction_tropospheric_no2_above_200m = static_cast<double*>(ancillary_data.fraction_tropospheric_no2_above_200m.request().ptr);
    ancillary_data_c.tropospheric_no2_concentration = static_cast<double*>(ancillary_data.tropospheric_no2_concentration.request().ptr);
    ancillary_data_c.stratospheric_no2_concentration = static_cast<double*>(ancillary_data.stratospheric_no2_concentration.request().ptr);
    ancillary_data_c.precipitable_water = static_cast<double*>(ancillary_data.precipitable_water.request().ptr);
    ancillary_data_c.water_vapor_bands = static_cast<double*>(ancillary_data.water_vapor_bands.request().ptr);
    ancillary_data_c.num_water_vapor_bands = ancillary_data.num_water_vapor_bands;

    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};

    gas_transmittance_table_c.h2o_transmittance = static_cast<double*>(gas_transmittance_table.h2o_transmittance.request().ptr);
    gas_transmittance_table_c.model = gas_transmittance_table.model;
    gas_transmittance_table_c.air_mass_factor_water_vapor = static_cast<double*>(gas_transmittance_table.air_mass_factor_water_vapor.request().ptr);
    gas_transmittance_table_c.wavelengths = static_cast<double*>(gas_transmittance_table.wavelengths.request().ptr);
    gas_transmittance_table_c.water_vapor_concentration = static_cast<double*>(gas_transmittance_table.water_vapor_concentration.request().ptr);
    gas_transmittance_table_c.num_models = gas_transmittance_table.num_models;
    gas_transmittance_table_c.num_wavelengths = gas_transmittance_table.num_wavelengths;
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    gas_transmittance_table_c.num_water_vapor_concentrations = gas_transmittance_table.num_water_vapor_concentrations;
    
    L1_Data l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<double*>(l1_data.cos_solar_zenith.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<double*>(l1_data.cos_sensor_zenith.request().ptr);
    l1_data_c.reflectance = static_cast<double*>(l1_data.reflectance.request().ptr);
    l1_data_c.wavelengths = static_cast<double*>(l1_data.wavelengths.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    Gas_Transmittances_PY gas_transmittances{};

    gas_transmittances.solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    gas_transmittances.total = allocate_output_array<double>(n_rows, n_cols);

    Gas_Transmittances gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<double*>(gas_transmittances.solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<double*>(gas_transmittances.sensor_zenith.request().ptr);
    gas_transmittances_c.total = static_cast<double*>(gas_transmittances.total.request().ptr);

    h2o_transmittance(&l1_data_c, &ancillary_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    return gas_transmittances;
}


PYBIND11_MODULE(gas_transmittance, m) 
{
    py::class_<Ancillary_Data_PY>(m, "Ancillary_Data", py::module_local())
        .def(py::init<>())
        .def_readwrite("ozone_absorption_cross_section", &Ancillary_Data_PY::ozone_absorption_cross_section)
        .def_readwrite("ozone_concentration", &Ancillary_Data_PY::ozone_concentration)
        .def_readwrite("no2_absorption_cross_section", &Ancillary_Data_PY::no2_absorption_cross_section)
        .def_readwrite("fraction_tropospheric_no2_above_200m", &Ancillary_Data_PY::fraction_tropospheric_no2_above_200m)
        .def_readwrite("tropospheric_no2_concentration", &Ancillary_Data_PY::tropospheric_no2_concentration)
        .def_readwrite("stratospheric_no2_concentration", &Ancillary_Data_PY::stratospheric_no2_concentration)
        .def_readwrite("precipitable_water", &Ancillary_Data_PY::precipitable_water)
        .def_readwrite("water_vapor_bands", &Ancillary_Data_PY::water_vapor_bands)
        .def_readwrite("num_water_vapor_bands", &Ancillary_Data_PY::num_water_vapor_bands);

    py::class_<Gas_Transmittance_Lookup_Table_PY>(m, "Gas_Transmittance_Lookup_Table", py::module_local())
        .def(py::init<>())
        .def_readwrite("co2_transmittance", &Gas_Transmittance_Lookup_Table_PY::co2_transmittance)
        .def_readwrite("co_transmittance", &Gas_Transmittance_Lookup_Table_PY::co_transmittance)
        .def_readwrite("ch4_transmittance", &Gas_Transmittance_Lookup_Table_PY::ch4_transmittance)
        .def_readwrite("o2_transmittance", &Gas_Transmittance_Lookup_Table_PY::o2_transmittance)
        .def_readwrite("n2o_transmittance", &Gas_Transmittance_Lookup_Table_PY::n2o_transmittance)
        .def_readwrite("h2o_transmittance", &Gas_Transmittance_Lookup_Table_PY::h2o_transmittance)
        .def_readwrite("model", &Gas_Transmittance_Lookup_Table_PY::model)
        .def_readwrite("wavelengths", &Gas_Transmittance_Lookup_Table_PY::wavelengths)
        .def_readwrite("air_mass_factor_mixed_gases", &Gas_Transmittance_Lookup_Table_PY::air_mass_factor_mixed_gases)
        .def_readwrite("air_mass_factor_water_vapor", &Gas_Transmittance_Lookup_Table_PY::air_mass_factor_water_vapor)
        .def_readwrite("water_vapor_concentration", &Gas_Transmittance_Lookup_Table_PY::water_vapor_concentration)
        .def_readwrite("num_models", &Gas_Transmittance_Lookup_Table_PY::num_models)
        .def_readwrite("num_wavelengths", &Gas_Transmittance_Lookup_Table_PY::num_wavelengths)
        .def_readwrite("num_amf_grid_points", &Gas_Transmittance_Lookup_Table_PY::num_amf_grid_points)
        .def_readwrite("num_water_vapor_concentrations", &Gas_Transmittance_Lookup_Table_PY::num_water_vapor_concentrations);

    py::class_<L1_Data_PY>(m, "L1_Data", py::module_local())
        .def(py::init<>())
        .def_readwrite("reflectance", &L1_Data_PY::reflectance)
        .def_readwrite("cos_solar_zenith", &L1_Data_PY::cos_solar_zenith)
        .def_readwrite("cos_sensor_zenith", &L1_Data_PY::cos_sensor_zenith)
        .def_readwrite("latitude", &L1_Data_PY::latitude)
        .def_readwrite("longitude", &L1_Data_PY::longitude)
        .def_readwrite("num_pixels", &L1_Data_PY::num_pixels)
        .def_readwrite("num_wavelengths", &L1_Data_PY::num_wavelengths)
        .def_readwrite("wavelengths", &L1_Data_PY::wavelengths);

    py::class_<Gas_Transmittances_PY>(m, "Gas_Transmittances", py::module_local())
        .def(py::init<>())
        .def_readwrite("solar_zenith", &Gas_Transmittances_PY::solar_zenith)
        .def_readwrite("sensor_zenith", &Gas_Transmittances_PY::sensor_zenith)
        .def_readwrite("total", &Gas_Transmittances_PY::total);

    py::enum_<Oxygen_A_Band_Option>(m, "Oxygen_A_Band_Option", py::module_local())
        .value("TRANSMITTANCE_TABLE", Oxygen_A_Band_Option::TRANSMITTANCE_TABLE)
        .value("SURROUNDING_WINDOW_BANDS", Oxygen_A_Band_Option::SURROUNDING_WINDOW_BANDS);

    m.def("ozone_transmittance", py::overload_cast<const L1_Data_PY&, const Ancillary_Data_PY&>(&ozone_transmittance));
    m.def("co2_transmittance", py::overload_cast<const L1_Data_PY&, const Gas_Transmittance_Lookup_Table_PY&, bool>(&co2_transmittance));
    m.def("co_transmittance", py::overload_cast<const L1_Data_PY&, const Gas_Transmittance_Lookup_Table_PY&, bool>(&co_transmittance));
    m.def("ch4_transmittance", py::overload_cast<const L1_Data_PY&, const Gas_Transmittance_Lookup_Table_PY&, bool>(&ch4_transmittance));
    m.def("o2_transmittance", py::overload_cast<const L1_Data_PY&, const Gas_Transmittance_Lookup_Table_PY& ,bool, Oxygen_A_Band_Option>(&o2_transmittance));
    m.def("n2o_transmittance", py::overload_cast<const L1_Data_PY&, const Gas_Transmittance_Lookup_Table_PY&, bool>(&n2o_transmittance));
    m.def("no2_transmittance", py::overload_cast<const L1_Data_PY&, const Ancillary_Data_PY&>(&no2_transmittance));
    m.def("h2o_transmittance", py::overload_cast<const L1_Data_PY&, const Ancillary_Data_PY&, const Gas_Transmittance_Lookup_Table_PY&, bool>(&h2o_transmittance));
}