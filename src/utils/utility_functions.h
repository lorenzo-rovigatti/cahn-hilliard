#include <cstdint>

// copied from https://stackoverflow.com/a/67539748/5140209
template <class T2, class T1>
T2 cpp11_bit_cast(T1 t1) {
    static_assert(sizeof(T1) == sizeof(T2), "Types must match sizes");
    static_assert(std::is_pod<T1>::value, "Requires POD input");
    static_assert(std::is_pod<T2>::value, "Requires POD output");

    T2 t2;
    std::memcpy(std::addressof(t2), std::addressof(t1), sizeof(T1));
    return t2;
}

// copied from https://codereview.stackexchange.com/a/272170
inline bool safe_isnan(double val) noexcept {
    const auto x = cpp11_bit_cast<std::uint64_t>(val);
    return (x & 0x7FFFFFFFFFFFFFFFu) >= 0x7FF0000000000001u;
}

template<int dims>
constexpr double pow_dims(double dx) {
    double result = 1.0;
    for (int i = 0; i < dims; ++i) result *= dx;
    return result;
}
