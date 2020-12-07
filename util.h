#if !defined(UTIL_H)
#define UTIL_H

#include <iostream>
#include <string>

#define PANIC(...) std::cout << "Panicked at " << __FILE__ << "/" << __LINE__ << std::endl, panic(__VA_ARGS__);

using String = std::string;

inline void nl()
{
    std::cout << std::endl;
}

void panic()
{
    std::cout << std::endl;
    abort();
}

template <typename T, typename... Ts>
void panic(T var1, Ts... var2)
{
    std::cout << var1;
    panic(var2...);
}

// Adapted from iFreilicht on Stack Overflow
template <typename... Args>
std::string string_format(const std::string &format, Args... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    if (size <= 0)
    {
        std::cout << "Error during formatting." << std::endl;
        abort();
    }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

/*Template hack to convert smart pointers to raw pointers
Mainly intended for use with std::is_ptr - because is_ptr is false for smart pointers
Example:
    std::is_pointer<double*>::value; // true
    std::is_pointer<std::shared_ptr<double>>::value; // false

    std::is_pointer<to_raw_pointer<double*>::Raw>::value; // true
    std::is_pointer<to_raw_pointer<std::shared_ptr<double>>::Raw>::value; // true
*/
template <typename Ptr = void, typename Elem = void>
struct to_raw_pointer
{
    using Element = Elem;
    using Raw = Elem;
};

template <typename Elem>
struct to_raw_pointer<std::shared_ptr<Elem>>
{
    using Element = Elem;
    using Raw = Elem *;
};

template <typename Elem>
struct to_raw_pointer<std::unique_ptr<Elem>>
{
    using Element = Elem;
    using Raw = Elem *;
};

template <typename Elem>
struct to_raw_pointer<std::weak_ptr<Elem>>
{
    using Element = Elem;
    using Raw = Elem *;
};

template <typename Elem>
struct to_raw_pointer<Elem *>
{
    using Element = Elem;
    using Raw = Elem *;
};
#endif // UTIL_H
