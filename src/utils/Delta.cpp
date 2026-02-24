/*
 * Delta.cpp
 *
 *  Created on: Jul 23, 2023
 *      Author: lorenzo
 */

#include "Delta.h"

#include "../defs.h"
#include <spdlog/fmt/ostr.h>

namespace ch {

Delta::Delta(toml::node_view<const toml::node> nv) {
	if(nv.is_number()) {
		_delta = nv.value<double>().value();
	}
	else if(nv.is_table()){
		const toml::table &delta_table = *nv.as_table();

		_delta = _config_optional_value<double>(delta_table, "value", 0.0);
		if(_delta == 0.0) {
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
	}
	else {
		std::ostringstream oss;
    	oss << nv.node();
		critical("A delta must be specified either as a single number or as a table of values. This is the offending part of the input file that couldn't be parsed: {}", oss.str());
	}
}

Delta::Delta(const toml::table &table, std::string path) : Delta(_config_node_view(table, path, true)) {
	
}

Delta::~Delta() {

}

} /* namespace ch */
