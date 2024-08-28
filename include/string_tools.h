#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <tuple>


/**
 * @brief Checks if a string starts with a given substring.
 *
 * @param s The string to check.
 * @param beg The substring to check at the beginning of the string.
 * @return true If the string starts with the substring.
 * @return false If the string does not start with the substring.
 */
inline bool startswith(const std::string &s, const char *beg) {
    if (s.size()==0 && *beg != '\0') return false;
    for (size_t i=0; i<s.size(); ++i) {
        if (beg[i] == '\0') return true;
        if (s[i] != beg[i]) return false;
    }
    return true;
}


/**
 * @brief Overloads the output stream operator to print vectors.
 *
 * @tparam T The type of elements in the vector.
 * @param o The output stream.
 * @param x The vector to print.
 * @return std::ostream& The modified output stream.
 */
template<class T>
std::ostream& operator<<(std::ostream& o, std::vector<T> x) {
    o << "{";
    if (x.size()) o << x.front();
    for (size_t i=1; i<x.size(); ++i)  o << ", " << x[i];
    return o << "}\n";
}


/**
 * @brief Reads the next line from a file stream.
 *
 * @param fid The input file stream.
 * @return std::string The line read from the file.
 */
inline std::string read_line(std::istream &fid) {
    std::string line;
    getline(fid, line);
    return line;
}


/**
 * @brief Converts a string to a specified type.
 *
 * @tparam T The type to convert the string to.
 * @param s The string to convert.
 * @return T The converted value.
 */
template <typename T>
inline T from_string(const std::string &s) {
    std::istringstream ss(s);
    T t;  ss >> t;
    return t;
}


/**
 * @brief Extracts a tuple of arbitrary types from a vector of strings at a given index.
 *
 * @tparam A The types of the tuple elements.
 * @param args The vector of strings.
 * @param index The index from where to start extracting.
 * @return std::tuple<A...> The extracted tuple.
 */
template<typename... A>
std::tuple<A...> args_to_tuple(std::vector<std::string> args, int index) {
    auto end = args.begin() + index + sizeof...(A);
    return std::make_tuple(from_string<A>(*--end)...);
}


/**
 * @brief Converts a string to an unsigned integer.
 *
 * @param s The string to convert.
 * @return unsigned The converted unsigned integer.
 */
inline unsigned str2u32(const std::string &s) {
    return from_string<unsigned>(s);
}


/**
 * @brief Converts a string to a double.
 *
 * @param s The string to convert.
 * @return double The converted double.
 */
inline double str2dbl(const std::string &s) {
    return from_string<double>(s);
}


/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * @param s The string to trim.
 * @param ws The characters considered as whitespace (default: " \t\n\r").
 * @return std::string The trimmed string.
 */
inline std::string trim(std::string s, std::string ws=" \t\n\r") {
    if (s.empty()) return s;
    size_t a=s.find_first_not_of(ws);
    size_t b = s.find_last_not_of(ws)+1;
    return s.substr(a, b-a);
}


/**
 * @brief Returns the part of a string before the first occurrence of any
 * character in a given set.
 *
 * @param str The string to process.
 * @param spl The set of characters to split on.
 * @return std::string The substring before the first occurrence of any character in spl.
 */
inline std::string before(std::string str, std::string spl) {
    auto pos = str.find_first_of(spl);
    if (pos == std::string::npos) return str;
    return str.substr(0, pos);
}


/**
 * @brief Splits a string into a vector of substrings, similar to Python's split function.
 *
 * @param s The string to split.
 * @param delims The delimiters to split on (default: " \t\n\r").
 * @return std::vector<std::string> The vector of substrings.
 */
inline std::vector<std::string> split(std::string s, std::string delims=" \t\n\r")
{
    if (s.empty()) return {};
    std::vector<std::string> pieces;
    size_t begin=0, end=0;
    while (end != std::string::npos) {
        begin = s.find_first_not_of(delims, end);
        end   = s.find_first_of(delims, begin);
        if (begin != std::string::npos)
            pieces.push_back(s.substr(begin, end-begin));
    }
    return pieces;
}


/**
 * @brief Joins a vector of strings into a single string, with a space separator.
 *
 * @param v The vector of strings to join.
 * @return std::string The joined string.
 */
inline std::string join(const std::vector<std::string> &v)
{
    if (v.empty()) return "";
    std::string s(v.front());
    for (auto viter =v.begin()+1; viter != v.end(); ++viter) s += " " + *viter;
    return s;
}


/**
 * @brief Joins strings together from iterators to strings, with a space separator.
 *
 * @tparam _Iter The iterator type.
 * @param beg The beginning iterator.
 * @param end The ending iterator.
 * @return std::string The joined string.
 */
template <class _Iter>
inline std::string join(_Iter beg, _Iter end) {
    if (beg == end) return "";
    std::string s(*(beg++));
    for (; beg!=end; ++beg) s += " " + *beg;
    return s;
}


/**
 * @brief Converts a value to a string with optional width and precision.
 *
 * @tparam T The type of the value.
 * @param x The value to convert to a string.
 * @param width The width of the string (default: -1, no specific width).
 * @param precision The precision of the string (default: -1, no specific precision).
 * @return std::string The string representation of the value.
 */
template<typename T>
inline std::string make_string(const T& x, int width=-1, int precision=-1) {
    std::ostringstream o;
    if (width>=0)     o << std::setw(width) << std::setfill('0');
    if (precision>=0) o << std::fixed << std::setprecision(precision);
    o << x;
    return o.str();
}


/**
 * @brief Converts a string to lowercase.
 *
 * @param s The string to convert.
 * @return std::string The lowercase string.
 */
inline std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
            static_cast<int(*)(int)> (tolower));
    return s;
}


/**
 * @brief Finds the matching closing brace in a string.
 *
 * @param s The string to search.
 * @param i The index of the opening brace.
 * @return size_t The index of the matching closing brace, or std::string::npos if not found.
 */
inline size_t find_matching_brace(const std::string &s, size_t i) {
    char open = s[i];
    char close = ')';
    if      (open == '[') close = ']';
    else if (open == '{') close = '}';
    int p = 1;  // Open bracket counter.
    for (i=i+1; i<s.size(); ++i) {
        if      (s[i] == close) --p;
        else if (s[i] == open)  ++p;
        if (p == 0) return i;
    }
    return std::string::npos;
}

