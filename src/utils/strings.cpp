/*
 * strings.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: lorenzo
 */

#include "strings.h"

#include <bitset>
#include <cstring>
#include <algorithm>

namespace ch {

namespace utils {

std::vector<std::string> split(const std::string &str, const std::string &delims) {
	std::vector<std::string> output;
	output.reserve(15);

	const char *ptr = str.c_str();
	while(ptr) {
		auto base = ptr;
		ptr = std::strpbrk(ptr, delims.c_str());
		if(ptr) {
			// this check makes sure that no empty strings are added to the output
			if(ptr - base) {
				output.emplace_back(base, ptr - base);
			}
			ptr++;
		}
		else {
			output.emplace_back(base);
		}
	}

	if(output.size() == 0) {
		output.push_back("");
	}

	return output;
}

bool starts_with(const std::string &value, std::string beginning) {
	if(beginning.size() > value.size()) {
		return false;
	}
	return std::equal(beginning.begin(), beginning.end(), value.begin());
}

bool ends_with(const std::string &value, std::string ending) {
	if(ending.size() > value.size()) {
		return false;
	}
	return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool contains(const std::string &value, std::string query) {
	return value.find(query) != std::string::npos;
}

void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
		return !std::isspace(ch);
	}));
}

void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}

void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}

std::string trim_copy(std::string source) {
	trim(source);
	return source;
}

}

}
