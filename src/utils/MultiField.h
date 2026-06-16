/*
 * MultiField.h
 *
 *  Created on: Nov 25, 2023
 *      Author: lorenzo
 */

#ifndef SRC_UTILS_MULTIFIELD_H_
#define SRC_UTILS_MULTIFIELD_H_

#include <cassert>
#include <vector>

namespace ch {

template <class T>
class SpeciesView {
    const T *_data;
    int _bins;
    int _species;

public:
    SpeciesView(const T* data, int bins, int species)
        : _data(data), _bins(bins), _species(species) {}

    T operator[](int s) const {
        assert(s >= 0 && s < _species);
        return _data[s * _bins];
    }

    int size() const { return _species; }

    // Iterators for "for(auto v : view)" support
    struct iterator {
        const T* ptr;
        int bins;
        int s;

        T operator*() const { return ptr[s * bins]; }
        iterator& operator++() { ++s; return *this; }
        bool operator!=(const iterator& other) const { return s != other.s; }
    };

    iterator begin() const { return iterator{_data, _bins, 0}; }
    iterator end() const { return iterator{_data, _bins, _species}; }
};


template <class T>
class MultiField {
    std::vector<T> _data;
    int _bins;
    int _species;
    int _log2_N_per_dim;

public:
    MultiField() : _bins(0), _species(0) {}
    MultiField(int N_bins, int species) : _data(N_bins * species), _bins(N_bins), _species(species), _log2_N_per_dim(log2(N_bins)) {}

    T &operator()(int idx, int species) {
        assert(idx >= 0 && idx < _bins);
        assert(species >= 0 && species < _species);
        return _data[species * _bins + idx];
    }

    T operator()(int idx, int species) const {
        assert(idx >= 0 && idx < _bins);
        assert(species >= 0 && species < _species);
        return _data[species * _bins + idx];
    }

    void fill(const T &v) {
        std::fill(_data.begin(), _data.end(), v);
    }

    T *data() {
        return _data.data();
    }
    const T *data() const {
        return _data.data();
    }

    T field_sum(int idx) const {
        T tot = T{};
        for(int s = 0; s < _species; s++) {
            tot += operator()(idx, s);
        }

        return tot;
    }

    // zero-allocation species view
    SpeciesView<T> species_view(int idx) const {
        assert(idx >= 0 && idx < _bins);
        return SpeciesView<T>(_data.data() + idx, _bins, _species);
    }

    int bins() const {
        return _bins;
    }

    int species() const {
        return _species;
    }

    int size() const {
        return _data.size();
    }

    template<int dims>
    void fill_coords(int coords[dims], int idx) {
        for(int d = 0; d < dims; d++) {
            coords[d] = idx & (_bins - 1);
            idx >>= _log2_N_per_dim; // divide by N
        }
    }

    template<int dims>
    int cell_idx(int coords[dims]) {
        int idx = 0;
        int multiply_by = 1;
        for(int d = 0; d < dims; d++) {
            idx += coords[d] * multiply_by;
            multiply_by <<= _log2_N_per_dim; // multiply by N
        }
        return idx;
    }

    template<int dims>
    double cell_laplacian(int species, int idx, double dx) {
        if constexpr (dims == 1) {
            int idx_m = (idx - 1 + _bins) & (_bins - 1);
            int idx_p = (idx + 1) & (_bins - 1);

            return (this->operator()(idx_m, species)
                + this->operator()(idx_p, species)
                - 2.0 * this->operator()(idx, species)) / SQR(dx);
        } 
        else {
            int coords[dims];
            int coords_n[dims];
            fill_coords<dims>(coords, idx);
            memcpy(coords_n, coords, sizeof(coords));

            double sum = 0.0;

            for(int d = 0; d < dims; d++) {
                // minus direction
                coords_n[d] = (coords[d] - 1 + _bins) & (_bins - 1);
                sum += this->operator()(cell_idx<dims>(coords_n), species);

                // plus direction
                coords_n[d] = (coords[d] + 1) & (_bins - 1);
                sum += this->operator()(cell_idx<dims>(coords_n), species);
                coords_n[d] = coords[d];
            }

            return (sum - 2.0 * dims * this->operator()(idx, species)) / SQR(dx);
        }
    }
};

} /* namespace ch */

#endif /* SRC_UTILS_MULTIFIELD_H_ */
