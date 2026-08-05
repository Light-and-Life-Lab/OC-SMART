import numpy as np
import matplotlib.pyplot as plt
import h5py
from pyhdf.V import *
from pyhdf.HDF import *
from pyhdf.SD import *
from scipy import interpolate
import time
import pytest

from gas_corrections_lib.src import gas_corrections

save_transmittances = False

@pytest.fixture(scope="session")
def read_HICO_geometry_data():
    csolz = np.load('test/HICO/csolz.npy')
    csenz = np.load('test/HICO/csenz.npy')

    return csolz, csenz

@pytest.fixture(scope="session")
def read_HICO_data():
    l1_data = gas_corrections.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)

    return l1_data


@pytest.fixture(scope="session")
def read_OCSSW_lat_lon():
    with h5py.File('test/HICO/ozone/H2013010071847.L2.ozone.nc', 'r') as f:
        OCSSW_lat = np.array(f['/navigation_data/latitude'])
        OCSSW_lon = np.array(f['/navigation_data/longitude'])

    return OCSSW_lat, OCSSW_lon


@pytest.fixture(scope="session")
def read_OCSMART_lat_lon():
    OCSMART_lat = np.load('test/HICO/l1b_lat.npy')
    OCSMART_lon = np.load('test/HICO/l1b_lon.npy')

    return OCSMART_lat, OCSMART_lon


@pytest.fixture(scope="session")
def read_gas_transmittance_table():
    gas_transmittance_table = gas_corrections.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')

    return gas_transmittance_table


@pytest.fixture(scope="session")
def read_ozone_ancillary_data():
    koz = np.load('test/HICO/ozone/koz.npy')
    ozone_concentration = np.load('test/HICO/ozone/oz_concentration.npy')

    return koz, ozone_concentration


@pytest.fixture(scope="session")
def read_OCSMART_ozone_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/HICO/ozone/tg_sol_oz.npy')
    tg_sen_ocsmart = np.load('test/HICO/ozone/tg_sen_oz.npy')
    sensor_wavelengths = np.load('test/HICO/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture(scope="session")
def read_OCSSW_ozone_transmittance_benchmark_data():
    with h5py.File('test/HICO/ozone/H2013010071847.L2.ozone.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_no2_ancillary_data():
    kno2 = np.load('test/HICO/no2/k_no2.npy')
    no2_frac = np.load('test/HICO/no2/no2_frac.npy')
    no2_tropo = np.load('test/HICO/no2/no2_tropo.npy')
    no2_strat = np.load('test/HICO/no2/no2_strat.npy')

    return kno2, no2_frac, no2_tropo, no2_strat


@pytest.fixture(scope="session")
def read_OCSMART_no2_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/HICO/no2/tg_sol_no2.npy')
    tg_sen_ocsmart = np.load('test/HICO/no2/tg_sen_no2.npy')
    sensor_wavelengths = np.load('test/HICO/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture(scope="session")
def read_OCSSW_no2_transmittance_benchmark_data():
    with h5py.File('test/HICO/no2/H2013010071847.L2.no2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_co2_transmittance_benchmark_data():
    with h5py.File('test/HICO/co2/H2013010071847.L2.co2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_co_transmittance_benchmark_data():
    with h5py.File('test/HICO/co/H2013010071847.L2.co.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_ch4_transmittance_benchmark_data():
    with h5py.File('test/HICO/ch4/H2013010071847.L2.ch4.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_n2o_transmittance_benchmark_data():
    with h5py.File('test/HICO/n2o/H2013010071847.L2.n2o.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_o2_opt2_transmittance_benchmark_data():
    with h5py.File('test/HICO/o2/H2013010071847.L2.o2.oxaband_opt_2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_o2_opt3_transmittance_benchmark_data():
    with h5py.File('test/HICO/o2/H2013010071847.L2.o2.oxaband_opt_3.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_h2o_transmittance_benchmark_data():
    with h5py.File('test/HICO/h2o/H2013010071847.L2.h2o.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


def save_gas_transmittances(gas_transmittances, gas_subfolder):
    output_path = 'test/HICO/' + gas_subfolder + '/'
    np.save(output_path + 'solar_zenith.npy', gas_transmittances.solar_zenith)
    np.save(output_path + 'sensor_zenith.npy', gas_transmittances.sensor_zenith)
    np.save(output_path + 'total.npy', gas_transmittances.total)


def load_gas_transmittances(gas_subfolder):
    output_path = 'test/HICO/' + gas_subfolder + '/'
    solar_zenith_saved = np.load(output_path + 'solar_zenith.npy')
    sensor_zenith_saved = np.load(output_path + 'sensor_zenith.npy')
    total_saved = np.load(output_path + 'total.npy')

    return solar_zenith_saved, sensor_zenith_saved, total_saved


@pytest.mark.skip()
def test_ozone_OCSSW(read_ozone_ancillary_data, 
                     read_HICO_geometry_data, 
                     read_OCSMART_ozone_transmittance_benchmark_data, 
                     read_OCSSW_ozone_transmittance_benchmark_data):
    
    ancillary_data = gas_corrections.Ancillary_Data()
    ancillary_data.ozone_absorption_cross_section, ancillary_data.ozone_concentration = read_ozone_ancillary_data
    
    # TODO: Adapt approach from OCSMART ancillary.py to read from ancillary files directly instead of .npy files
    # TODO: This will require interpolation of ozone data to l1b grid
    l1_data = gas_corrections.L1_Data()
    l1_data.cos_solar_zenith, l1_data.cos_sensor_zenith = read_HICO_geometry_data
    l1_data.num_pixels = l1_data.cos_solar_zenith.shape[0] * l1_data.cos_solar_zenith.shape[1]
    l1_data.num_wavelengths = len(ancillary_data.ozone_absorption_cross_section)

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    gas_transmittances = gas_corrections.ozone_transmittance(l1_data=l1_data, ancillary_data=ancillary_data, use_gas_transmittance_lookup_table=False)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_ozone_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((1710, 1272, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((1710, 1272, 197))

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_ozone_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'ozone')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('ozone')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    # plt.plot(sensor_wavelengths, tg_sen_ocsmart[0, 0, :], '-g')
    # plt.plot(sensor_wavelengths, tg_sol_ocsmart[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths, tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/ozone/transmittance_comparison.png')


@pytest.mark.skip()
def test_co2_OCSSW(read_OCSSW_co2_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.co2_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table, use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co2_transmittance_benchmark_data

    sensor_wavelengths = gas_correction_manager.l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'co2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('co2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/co2/transmittance_comparison.png')


@pytest.mark.skip()
def test_co_OCSSW(read_OCSSW_co_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.co_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table, use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co_transmittance_benchmark_data

    sensor_wavelengths = gas_correction_manager.l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'co')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('co')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/co/transmittance_comparison.png')


@pytest.mark.skip()
def test_ch4_OCSSW(read_OCSSW_ch4_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.ch4_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table, use_gas_transmittance_lookup_table=False)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_ch4_transmittance_benchmark_data

    sensor_wavelengths = gas_transmittance_manager.l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'ch4')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('ch4')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/ch4/transmittance_comparison.png')


@pytest.mark.skip()
def test_n2o_OCSSW(read_OCSSW_n2o_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.n2o_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table, use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_n2o_transmittance_benchmark_data

    sensor_wavelengths = gas_transmittance_manager.l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'n2o')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('n2o')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/n2o/transmittance_comparison.png')


@pytest.mark.skip()
def test_no2_OCSSW(read_no2_ancillary_data, 
                   read_HICO_geometry_data, 
                   read_OCSMART_no2_transmittance_benchmark_data, 
                   read_OCSSW_no2_transmittance_benchmark_data):

    ancillary_data = gas_corrections.Ancillary_Data()
    ancillary_data.no2_absorption_cross_section, \
    ancillary_data.fraction_tropospheric_no2_above_200m, \
    ancillary_data.tropospheric_no2_concentration, \
    ancillary_data.stratospheric_no2_concentration = read_no2_ancillary_data

    # TODO: Adapt approach from OCSMART ancillary.py to read from ancillary files directly instead of .npy files
    # TODO: This will require interpolation of ozone data to l1b grid
    l1_data = gas_corrections.L1_Data()
    l1_data.cos_solar_zenith, l1_data.cos_sensor_zenith = read_HICO_geometry_data
    l1_data.num_pixels = l1_data.cos_solar_zenith.shape[0] * l1_data.cos_solar_zenith.shape[1]
    l1_data.num_wavelengths = len(ancillary_data.no2_absorption_cross_section)

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    gas_transmittances = gas_corrections.no2_transmittance(l1_data=l1_data, ancillary_data=ancillary_data, use_gas_transmittance_lookup_table=False)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_no2_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((1710, 1272, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((1710, 1272, 197))

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_no2_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'no2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('no2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths, tg_sen_ocsmart[0, 0, :], '-g')
    plt.plot(sensor_wavelengths, tg_sol_ocsmart[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths, tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith'])
    plt.savefig('test/HICO/no2/transmittance_comparison.png')


@pytest.mark.skip()
def test_o2_OCSSW_transmittance_table_option(read_OCSSW_o2_opt2_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, 
                                                          gas_transmittance_table=gas_transmittance_table, 
                                                          oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option().TRANSMITTANCE_TABLE, 
                                                          use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_o2_opt2_transmittance_benchmark_data
    
    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'o2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('o2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([750, 800])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/o2/transmittance_comparison_oxaband_opt_2.png')


@pytest.mark.skip()
def test_o2_OCSSW_surrounding_window_bands_option(read_OCSSW_o2_opt3_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, 
                                                          gas_transmittance_table=gas_transmittance_table, 
                                                          oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option().SURROUNDING_WINDOW_BANDS, 
                                                          use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_o2_opt3_transmittance_benchmark_data
    
    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'o2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('o2')

    # assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    # assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    # assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([750, 800])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/o2/transmittance_comparison_oxaband_opt_3.png')


@pytest.mark.skip()
def test_h2o_OCSSW(read_OCSSW_h2o_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_HICO_data('test/HICO/H2013010071847.L1B_ISS.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/HICO/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.h2o_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table, use_gas_transmittance_lookup_table=True)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_h2o_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'h2o')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('h2o')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d[:], tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d[:], tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_correction_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_correction_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([500, 900])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/HICO/h2o/transmittance_comparison.png')