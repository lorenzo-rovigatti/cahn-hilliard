/*
 * defs.h
 *
 *  Created on: Jul 14, 2023
 *      Author: lorenzo
 */

#ifndef SRC_DEFS_H_
#define SRC_DEFS_H_

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

#define TOML_EXCEPTIONS 0
#define TOML_ENABLE_FORMATTERS 0
#include <toml++/toml.hpp>

#include <iostream>

#include "defs_CUDA.h"

#endif /* SRC_DEFS_H_ */
