/*
 * RhoMatrix.h.h
 *
 *  Created on: Nov 25, 2023
 *      Author: lorenzo
 */

#ifndef SRC_UTILS_RHOMATRIX_H_
#define SRC_UTILS_RHOMATRIX_H_

namespace ch {

template <class T>
class RhoMatrix {
    std::vector<T> _data;
    int _bins;
    int _species;

public:
    T &operator()(int idx, int species) {
        return _data[species * _bins + idx];
    }

    T operator()(int idx, int species) const {
        return _data[species * _bins + idx];
    }

    void fill(const T &v) {
        std::fill(_data.begin(), _data.end(), v);
    }

    T *data() {
        return _data.data();
    }

    T rho_tot(int idx) const {
        T tot = 0.;
        for(int s = 0; s < _species; s++) {
            tot += (*this)(idx, s);
        }

        return tot;
    }

    std::vector<T> rho_species(int idx) const {
        std::vector<T> res(_species);
        for(int s = 0; s < _species; s++) {
            res[s] = (*this)(idx, s);
        }

        return res;
    }

    int bins() const {
        return _bins;
    }

    RhoMatrix() {}
    RhoMatrix(int N_bins, int species) : _data(N_bins * species), _bins(N_bins), _species(species) {}
};

} /* namespace ch */

#endif /* SRC_UTILS_RHOMATRIX_H_ */
