#include <array>

// Templated class to store a staggered gradient at a single point in D-dimensional space
template <std::size_t dims>
class Gradient {
public:
    std::array<double, dims> components;  // Stores the gradient components in each dimension

    // Constructor - Initializes to zero by default
    Gradient() {
        components.fill(0.0);
    }

    // Constructor - Initialize with an array of values
    explicit Gradient(const std::array<double, dims>& values) : components(values) {}

    double& operator[](std::size_t index) {
        return components[index];
    }

    const double& operator[](std::size_t index) const {
        return components[index];
    }

    // Addition
    Gradient operator+(const Gradient& other) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++)
            result.components[i] = components[i] + other.components[i];
        return result;
    }

    // Subtraction
    Gradient operator-(const Gradient& other) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++) {
            result.components[i] = components[i] - other.components[i];
        }
        return result;
    }

    // Scalar multiplication
    Gradient operator*(double scalar) const {
        Gradient result;
        for (std::size_t i = 0; i < dims; i++) {
            result.components[i] = components[i] * scalar;
        }
        return result;
    }

    // Dot product
    double operator*(Gradient &other) const {
        double result = 0;
        for (std::size_t i = 0; i < dims; i++) {
            result += components[i] * other[i];
        }
        return result;
    }

    // In-place addition
    Gradient& operator+=(const Gradient& other) {
        for (std::size_t i = 0; i < dims; ++i) {
            components[i] += other.components[i];
        }
        return *this;
    }

    // In-place subtraction
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

    // Print function for debugging
    void print() const {
        std::cout << "Gradient: (";
        for (std::size_t i = 0; i < dims; i++) {
            std::cout << components[i];
            if (i < dims - 1) std::cout << ", ";
        }
        std::cout << ")\n";
    }
};
