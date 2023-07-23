/*
 * Delta.h
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#ifndef SRC_UTILS_DELTA_H_
#define SRC_UTILS_DELTA_H_

#include "../Object.h"

namespace ch {

class Delta : public Object {
public:
	Delta(toml::table &table, std::string path);
	virtual ~Delta();
	Delta(const Delta &other) = default;
	Delta(Delta &&other) = default;

	operator double() const {
		return _delta;
	}

	GET_NAME(Delta utility class)

private:
	double _delta = 0.0;
};

} /* namespace ch */

#endif /* SRC_UTILS_DELTA_H_ */
