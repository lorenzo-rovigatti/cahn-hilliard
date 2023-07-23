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

namespace ch {

class Object {
public:
	Object();
	virtual ~Object();
	Object(const Object &other) = default;
	Object(Object &&other) = default;

	toml::node_view<toml::node> _config_node_view(toml::table &tbl, const std::string &path, bool mandatory) const;

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
