/*
 * Delta.cpp
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#include "Delta.h"

#include "../defs.h"

namespace ch {

Delta::Delta(toml::table &table, std::string path) {
	auto nv = _config_node_view(table, path, true);
	if(nv.is_number()) {
		_delta = nv.value<double>().value();
	}
	else if(nv.is_table()){
		const toml::table &delta_table = *nv.as_table();
		double T = _config_value<double>(delta_table, "T");
		double salt = _config_optional_value<double>(delta_table, "salt", 1.0);
		int L_DNA = _config_optional_value<int>(delta_table, "sticky_size", 6);
		double delta_H = _config_value<double>(delta_table, "deltaH");
		double delta_S = _config_value<double>(delta_table, "deltaS");

		double delta_S_salt = 0.368 * (L_DNA - 1.0) * std::log(salt);
		double delta_G = delta_H - T * (delta_S + delta_S_salt);

		const double k_B = 1.9872036;
		_delta = 1.6606 * std::exp(-delta_G / (k_B * T));
	}
	else {
		critical("The {} option should be either a number or a table of values", path);
	}
}

Delta::~Delta() {

}

} /* namespace ch */
