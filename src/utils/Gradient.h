#pragma once
#include <array>

// Templated class to store a staggered gradient at a single point in D-dimensional space
template <std::size_t dims>
class Gradient {
public:
    std::array<double, dims> components;

    Gradient() {
        components.fill(0.0);
    }

    explicit Gradient(const std::array<double, dims>& values) : components(values) {}

    double& operator[](std::size_t index) {
        return components[index];
    }

    const double& operator[](std::size_t index) const {
        return components[index];
    }

    Gradient operator+(const Gradient& other) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++)
            result.components[i] = components[i] + other.components[i];
        return result;
    }

    Gradient operator-(const Gradient& other) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++) {
            result.components[i] = components[i] - other.components[i];
        }
        return result;
    }

    Gradient operator*(double scalar) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++) {
            result.components[i] = components[i] * scalar;
        }
        return result;
    }

    double operator*(Gradient &other) const {
        double result = 0;
        for (std::size_t i = 0; i < dims; i++) {
            result += components[i] * other[i];
        }
        return result;
    }

    Gradient& operator+=(const Gradient& other) {
        for (std::size_t i = 0; i < dims; ++i) {
            components[i] += other.components[i];
        }
        return *this;
    }

    Gradient& operator-=(const Gradient& other) {
        for (std::size_t i = 0; i < dims; ++i) {
            components[i] -= other.components[i];
        }
        return *this;
    }

    Gradient operator/(double scalar) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++) {
            result.components[i] = components[i] / scalar;
        }
        return result;
    }

    friend Gradient operator*(double scalar, const Gradient& grad) {
        return grad * scalar;
    }
};
