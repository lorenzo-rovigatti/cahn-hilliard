/*
 * Object.h
 *
 *  Created on: Jul 20, 2023
 *      Author: lorenzo
 */

#ifndef SRC_OBJECT_H_
#define SRC_OBJECT_H_

#define GET_NAME(NAME) \
  std::string_view name() const override { return #NAME; }

#include "defs.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/fmt/fmt.h>
#include <type_traits>

namespace ch {

class Object {
public:
	Object();
	virtual ~Object();
	Object(const Object &other) = default;
	Object(Object &&other) = default;

	toml::node_view<toml::node> _config_node_view(toml::table &tbl, const std::string &path, bool mandatory) const;

	template<typename T>
	std::vector<T> _config_array_values(toml::table &tbl, const std::string &path, size_t output_size) const {
		static_assert(std::is_integral_v<T> || std::is_floating_point_v<T> || std::is_same_v<T, std::string>,
				"_config_array_values only accepts std::string, integral or floating point types");

		// this horror converts C++ types to toml++ supported types (integral types to int64_t, floating point types
		// to double and assumes that if T is not integral nor floating point then it is a string
		using toml_T = std::conditional_t<
				std::is_integral<T>::value,
				int64_t,
				std::conditional_t<
				std::is_floating_point<T>::value,
				double,
				std::string>
		>;

		std::vector<T> output;

		auto nv = _config_node_view(tbl, path, true);
		// if there is a single value then we assume the user want to use it for all the output_size cells of the output
		if(nv.is<toml_T>()) {
			T value = nv.value<toml_T>().value();
			output = std::vector<T>(output_size, value);
		}
		else if(nv.is_array()) {
			auto array = nv.as_array();
			if(!array->is_homogeneous()) {
				critical("The array specified in option '{}' does not contain values of the same type (i.e. homogeneous values)", path);
			}

			for(auto &&elem : *array) {
				if(!elem.template is<toml_T>()) {
					// we explicitly convert elem.source() to a std::string in order to print it
					std::stringstream ss;
					ss << elem.source();
					critical("Option '{}': cannot convert element at {}", path, ss.str());
				}
				T value = elem.template value<toml_T>().value();
				output.push_back(value);
			}
		}
		if(output.size() != output_size) {
			critical("Option '{}' should be a single value or an array of size {}", path, output_size);
		}

		return output;
	}

	template<typename T>
	T _config_value(toml::table &tbl, const std::string &path) const {
		auto nv = _config_node_view(tbl, path, true);
		return nv.value<T>().value();
	}

	template<typename T>
	T _config_optional_value(toml::table &tbl, const std::string &path, T default_value) const {
		auto nv = _config_node_view(tbl, path, false);
		return nv.value<T>().value_or(default_value);
	}

	template<typename FormatString, typename... Args>
	void info(const FormatString &fmt, Args&&...args) const {
		auto format_final = fmt::format("{} (source: {{}})", fmt);
		spdlog::info(format_final, std::forward<Args>(args)..., name());
	}

	template<typename FormatString, typename... Args>
	void critical(const FormatString &fmt, Args&&...args) const {
		auto format_final = fmt::format("{} (error source: {{}})", fmt);
		spdlog::critical(format_final, std::forward<Args>(args)..., name());
		exit(1);
	}

	virtual std::string_view name() const = 0;
};

} /* namespace ch */

#endif /* SRC_OBJECT_H_ */
