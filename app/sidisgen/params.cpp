#include "params.hpp"

#include <initializer_list>
#include <istream>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>

#include <TDirectory.h>

namespace {

// Remove whitespace from beginning and end of string.
std::string trim(std::string str) {
	auto not_whitespace = [](char c) {
		return !std::isspace(c);
	};
	str.erase(str.begin(), std::find_if(str.begin(), str.end(), not_whitespace));
	str.erase(std::find_if(str.rbegin(), str.rend(), not_whitespace).base(), str.end());
	return str;
}

// Search for comment character '#' and remove anything following it.
std::string trim_comment(std::string str) {
	str.erase(std::find(str.begin(), str.end(), '#'), str.end());
	return str;
}

}

std::string Type::to_string(Value const& value) const {
	std::ostringstream os;
	write_stream(os, value);
	if (!os) {
		throw std::runtime_error("Could not convert value to string.");
	}
	return os.str();
}

Params::Param::Param(
		Value const* default_value,
		std::initializer_list<char const*> tags,
		char const* usage,
		char const* brief,
		char const* doc) :
		type(default_value->type()),
		value(nullptr),
		default_value(default_value),
		tags(),
		used(false),
		usage(usage),
		brief(brief),
		doc(doc) {
	for (char const* tag : tags) {
		this->tags.insert(tag);
	}
}
Params::Param::Param(
		Type const& type,
		std::initializer_list<char const*> tags,
		char const* usage,
		char const* brief,
		char const* doc) :
		type(type),
		value(nullptr),
		default_value(nullptr),
		tags(),
		used(false),
		usage(usage),
		brief(brief),
		doc(doc) {
	for (char const* tag : tags) {
		this->tags.insert(tag);
	}
}

Value const& Params::get(char const* name) {
	Param& param = _params.at(name);
	if (param.value == nullptr && param.default_value == nullptr) {
		throw std::runtime_error(
			std::string("Required parameter '") + name + "' not provided.");
	} else if (param.value == nullptr) {
		return *param.default_value.get();
	} else {
		param.used = true;
		return *param.value.get();
	}
}

bool Params::set(char const* name, Value const* value) {
	Param& param = _params.at(name);
	bool old = (param.value != nullptr);
	if (param.type != value->type()) {
		throw std::runtime_error(
			std::string("Parameter '") + name + "' has different type.");
	}
	param.used = false;
	param.value.reset(value);
	return old;
}

void Params::check_format(Params const& other) const {
	for (auto const& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param = pair.second;
		auto it = other._params.find(name);
		if (it == other._params.end()) {
			throw std::runtime_error(
				std::string("Couldn't find parameter '") + name + "'.");
		}
		Param const& other_param = it->second;
		if (param.type != other_param.type) {
			throw std::runtime_error(
				std::string("Parameter '") + name + "' has different type.");
		}
		if (param.tags != other_param.tags) {
			throw std::runtime_error(
				std::string("Parameter '") + name + "' has different tags.");
		}
	}
}

void Params::check_complete() const {
	for (auto const& pair : _params) {
		if (pair.second.value != nullptr && !pair.second.used) {
			throw std::runtime_error(
				"Provided parameter '" + pair.first + "' was not used.");
		}
	}
}

void Params::check_can_provide(Params const& consumer) const {
	check_format(consumer);
	for (auto const& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param_prod = pair.second;
		Param const& param_cons = consumer._params.find(name)->second;
		if (!param_prod.type.can_provide(*param_prod.value, *param_cons.value)) {
			throw std::runtime_error(
				std::string("Parameter '") + name + "' is incompatible between "
				"producer "
				+ "(value " + param_prod.type.to_string(*param_prod.value) + ") and "
				+ "consumer "
				+ "(value " + param_prod.type.to_string(*param_cons.value) + ").");
		}
	}
}

void Params::check_equivalent(Params const& other) const {
	check_format(other);
	for (auto const& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param = pair.second;
		Param const& other_param = other._params.find(name)->second;
		if (!param.type.equivalent(*param.value, *other_param.value)) {
			throw std::runtime_error(
				std::string("Parameter '") + name + "' is not equivalent between source "
				+ "(value " + param.type.to_string(*param.value) + ") and dest "
				+ "(value " + param.type.to_string(*other_param.value) + ").");
		}
	}
}

void Params::clear_unused() {
	for (auto& pair : _params) {
		if (!pair.second.used) {
			pair.second.value.reset();
		}
	}
}

void Params::read_root(TDirectory& dir) {
	for (auto& pair : _params) {
		char const* name = pair.first.c_str();
		Param& param = pair.second;
		param.used = false;
		try {
			param.value = param.type.read_root(dir, name);
		} catch (...) {
			// TODO: Don't catch every exception here.
			throw std::runtime_error(
				std::string("Couldn't read parameter '")
				+ name + "' from ROOT directory.");
		}
	}
}

void Params::write_root(TDirectory& dir) const {
	for (auto& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param = pair.second;
		// TODO: Force writing the version parameter.
		if (param.value != nullptr) {
			try {
				param.type.write_root(dir, name, *param.value);
			} catch (...) {
				// TODO: Don't catch every exception here.
				throw std::runtime_error(
					std::string("Couldn't write parameter '")
					+ name + "' to ROOT directory.");
			}
		}
	}
}

void Params::read_stream(std::istream& is) {
	// Create a map matching parameter names to strings from the stream.
	std::map<std::string, std::string> map;
	while (is) {
		std::string line;
		std::getline(is, line);
		std::stringstream ss(trim(trim_comment(line)));
		std::string key;
		std::string value;
		ss >> key;
		std::getline(ss, value);
		value = trim(value);
		if (!key.empty()) {
			if (map.find(key) != map.end()) {
				throw std::runtime_error("Duplicate parameter '" + key + "'.");
			}
			map[key] = value;
		}
	}

	// Try to read each parameter in turn from the map.
	for (auto& pair : _params) {
		char const* name = pair.first.c_str();
		Param& param = pair.second;
		Type const& type = param.type;
		// Remove the parameter from the map once it's been read.
		auto map_it = map.find(name);
		if (map_it != map.end()) {
			// Strip all trailing whitespace or comments from param.
			std::istringstream ss(map_it->second);
			map.erase(map_it);
			try {
				param.value = type.read_stream(ss);
				if (!ss) {
					throw std::runtime_error(
						std::string("Could not read parameter '") + name
						+ "' from stream.");
				}
				param.used = false;
				std::string rem;
				std::getline(ss, rem);
				if (!rem.empty()) {
					throw std::runtime_error("");
				}
			} catch (std::exception const& e) {
				throw std::runtime_error(
					std::string("Failed to parse parameter '")
					+ name + "' from '" + ss.str() + "'.");
			}
		}
	}

	// Check for any leftover parameters.
	if (!map.empty()) {
		std::ostringstream ss_err;
		ss_err << "Unrecognized parameters";
		while (!map.empty()) {
			auto map_it = map.begin();
			ss_err << " '" << map_it->first << "'";
			map.erase(map_it);
		}
		ss_err << ".";
		throw std::runtime_error(ss_err.str());
	}
}

void Params::write_stream(std::ostream& os) const {
	for (auto& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param = pair.second;
		// TODO: Force writing the version parameter.
		if (param.value != nullptr) {
			os << name << ' ';
			param.type.write_stream(os, *param.value);
			os << std::endl;
			if (!os) {
				throw std::runtime_error(
					std::string("Could not write parameter '")
					+ name + "' to stream.");
			}
		}
	}
}

Params Params::filter(std::initializer_list<char const*> matches) {
	Params params;
	std::vector<std::regex> regex_list;
	for (char const* match : matches) {
		regex_list.emplace_back(match, std::regex_constants::ECMAScript);
	}
	for (auto& pair : _params) {
		char const* name = pair.first.c_str();
		Param const& param = pair.second;
		bool match = false;
		for (std::regex const& regex : regex_list) {
			if (std::regex_match(name, regex)) {
				match = true;
				goto found_match;
			}
			for (std::string const& tag : param.tags) {
				if (std::regex_match(tag.c_str(), regex)) {
					match = true;
					goto found_match;
				}
			}
		}
found_match:
		if (match) {
			params._params.emplace(name, param);
		}
	}
	return params;
}

