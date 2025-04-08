/*
 * strings.h
 *
 *  Created on: Feb 28, 2020
 *      Author: lorenzo
 */

#ifndef UTILS_STRINGS_H_
#define UTILS_STRINGS_H_

#include <vector>
#include <sstream>

namespace ch {

namespace utils {

/**
 * Split a string into tokens.
 *
 * @param line
 * @param separators list of separators that will be used to split the string
 * @return a vector of strings containing the tokens
 */
std::vector<std::string> split(const std::string &str, const std::string &delims=" \t");

// taken from https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c/2072890
bool starts_with(const std::string &value, std::string beginning);
bool ends_with(const std::string &value, std::string ending);
bool contains(const std::string &value, std::string query);

// trim from start (in place)
void ltrim(std::string &s);

// trim from end (in place)
void rtrim(std::string &s);

// trim from both ends (in place)
void trim(std::string &s);

/**
 * Return a copy of the given string, without trailing and leading spaces
 *
 * @param source
 * @return a new string
 */
std::string trim_copy(std::string source);

}

}

#endif /* UTILS_STRINGS_H_ */
