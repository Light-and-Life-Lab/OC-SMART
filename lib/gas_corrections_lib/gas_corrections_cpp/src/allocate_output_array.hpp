#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

template<typename T>
pybind11::array_t<T> allocate_output_array(int n_rows, int n_cols)
{
    // Set up an empty array to populate, which we will return later
    pybind11::array_t<T> array = pybind11::array(
        pybind11::buffer_info(
            nullptr,                                // nullptr -> Ask numpy to allocate
            sizeof(T),                              // Size of the array elements
            pybind11::format_descriptor<T>::value,        // Format string
            2,                                      // Number of dimensions of this array
            {n_rows, n_cols},                       // Shape of this array
            {sizeof(T)*n_cols, sizeof(T)}           // Stride of this array
        )
    );

    return array;
}