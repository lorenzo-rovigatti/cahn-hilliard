/*
 * Object.cpp
 *
 *  Created on: Jul 20, 2023
 *      Author: lorenzo
 */

#include "Object.h"

namespace ch {

Object::Object() {

}

Object::~Object() {

}

toml::node_view<toml::node> Object::_config_node_view(toml::table &tbl, const std::string &path, bool mandatory) const {
	toml::path t_path(path);
	auto nv = tbl[t_path];

	if(mandatory && !nv) {
		critical("Mandatory option '{}' not found", path);
	}

	return nv;
}

} /* namespace ch */
