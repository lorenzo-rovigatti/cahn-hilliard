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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/fmt/fmt.h>
#include <string>

namespace ch {

class Object {
public:
	Object();
	virtual ~Object();
	Object(const Object &other) = default;
	Object(Object &&other) = default;

	template<typename FormatString, typename... Args>
	void critical(const FormatString &fmt, Args&&...args) {
		auto format_final = fmt::format("{} ({{}})", fmt);
		spdlog::critical(format_final, std::forward<Args>(args)..., name());
		exit(1);
	}

	virtual std::string_view name() const = 0;
};

} /* namespace ch */

#endif /* SRC_OBJECT_H_ */
