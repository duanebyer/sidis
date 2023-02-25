#include "params.hpp"

#include <algorithm>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

Filter const Filter::REJECT = Filter();
Filter const Filter::ACCEPT = Filter(Term{});

Filter& Filter::operator|=(Filter const& rhs) {
	_condition.insert(rhs._condition.begin(), rhs._condition.end());
	return *this;
}

Filter& Filter::operator&=(Filter const& rhs) {
	Condition condition_next;
	for (Term const& term_1 : _condition) {
		for (Term const& term_2 : rhs._condition) {
			Term term_union;
			std::set_union(
				term_1.begin(), term_1.end(),
				term_2.begin(), term_2.end(),
				std::inserter(term_union, term_union.begin()));
			condition_next.insert(term_union);
		}
	}
	_condition = condition_next;
	return *this;
}

Params::Param::Param(
		Value const* default_value,
		std::vector<std::string> tags,
		std::string usage,
		std::string brief,
		std::string doc) :
		type(default_value->type()),
		value(nullptr),
		default_value(default_value),
		tags(),
		used(false),
		usage(usage),
		brief(brief),
		doc(doc) {
	for (std::string const& tag : tags) {
		this->tags.insert(tag);
	}
}
Params::Param::Param(
		Type const& type,
		std::vector<std::string> tags,
		std::string usage,
		std::string brief,
		std::string doc) :
		type(type),
		value(nullptr),
		default_value(nullptr),
		tags(),
		used(false),
		usage(usage),
		brief(brief),
		doc(doc) {
	for (std::string const& tag : tags) {
		this->tags.insert(tag);
	}
}

Value const& Params::get(std::string const& name) {
	Param& param = _params.at(name);
	if (param.value == nullptr && param.default_value == nullptr) {
		throw std::runtime_error(
			"Required parameter '" + name + "' not provided.");
	} else if (param.value == nullptr) {
		return *param.default_value.get();
	} else {
		param.used = true;
		return *param.value.get();
	}
}

Value const& Params::get_soft(std::string const& name) const {
	Param const& param = _params.at(name);
	if (param.value == nullptr && param.default_value == nullptr) {
		throw std::runtime_error(
			"Required parameter '" + name + "' not provided.");
	} else if (param.value == nullptr) {
		return *param.default_value.get();
	} else {
		return *param.value.get();
	}
}

bool Params::set(std::string const& name, Value const* value) {
	Param& param = _params.at(name);
	bool old = (param.value != nullptr);
	if (value != nullptr && param.type != value->type()) {
		throw std::runtime_error(
			"Parameter '" + name + "' has different type.");
	}
	param.used = false;
	param.value.reset(value);
	return old;
}

bool Params::set_from(Params const& other, std::string const& name) {
	Param& param = _params.at(name);
	Param const& other_param = other._params.at(name);
	bool old = (param.value != nullptr);
	if (param.type != other_param.type) {
		throw std::runtime_error(
			"Parameter '" + name + "' has different type.");
	}
	if (param.tags != other_param.tags) {
		throw std::runtime_error(
			"Parameter '" + name + "' has different tags.");
	}
	if (param.default_value != other_param.default_value) {
		throw std::runtime_error(
			"Parameter '" + name + "' has different default.");
	}
	param.used = false;
	param.value = other_param.value;
	return old;
}

bool Params::set_from(Params const& other) {
	// TODO: This implementation could be made much more efficient by iterating
	// over the parameter map directly.
	bool old = false;
	for (std::string const& name : other.names()) {
		old |= this->set_from(other, name);
	}
	return old;
}

std::set<std::string> Params::names() const {
	std::set<std::string> result;
	for (auto const& pair : _params) {
		result.insert(pair.first);
	}
	return result;
}

void Params::check_format(Params const& other) const {
	for (auto const& pair : _params) {
		std::string const& name = pair.first;
		Param const& param = pair.second;
		auto it = other._params.find(name);
		if (it == other._params.end()) {
			throw std::runtime_error(
				"Could not find parameter '" + name + "'.");
		}
		Param const& other_param = it->second;
		if (param.type != other_param.type) {
			throw std::runtime_error(
				"Parameter '" + name + "' has different type.");
		}
		if (param.tags != other_param.tags) {
			throw std::runtime_error(
				"Parameter '" + name + "' has different tags.");
		}
		if (param.default_value != other_param.default_value) {
			throw std::runtime_error(
				"Parameter '" + name + "' has different default.");
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

void Params::check_equivalent(Params const& other) const {
	check_format(other);
	for (auto const& pair : _params) {
		std::string const& name = pair.first;
		Param const& param = pair.second;
		Param const& other_param = other._params.find(name)->second;
		if (((param.value == nullptr) ^ (other_param.value == nullptr))
				|| (param.value != nullptr && other_param.value != nullptr
					&& !param.type.equivalent(*param.value, *other_param.value))) {
			std::string param_str = "no value given";
			std::string other_param_str = "no value given";
			if (param.value != nullptr) {
				param_str = "value '" + param.value->to_string() + "'";
			}
			if (other_param.value != nullptr) {
				other_param_str = "value '" + other_param.value->to_string() + "'";
			}
			throw std::runtime_error(
				"Parameter '" + name + "' is not equivalent between source "
				"(" + param_str + ") and dest (" + other_param_str + ").");
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
	TDirectory* params_dir = dir.GetDirectory("params");
	if (params_dir == nullptr) {
		throw std::runtime_error("Could not open directory 'params'.");
	}
	for (auto& pair : _params) {
		std::string const& name = pair.first;
		Param& param = pair.second;
		param.used = false;
		try {
			param.value = param.type.read_root(*params_dir, name);
		} catch (...) {
			// TODO: Don't catch every exception here.
			throw std::runtime_error(
				"Could not read parameter '" + name + "' from ROOT directory.");
		}
	}
}

void Params::write_root(TDirectory& dir) const {
	TDirectory* params_dir = dir.mkdir("params", "params");
	if (params_dir == nullptr) {
		throw std::runtime_error("Could not create directory 'stats'.");
	}
	for (auto& pair : _params) {
		std::string const& name = pair.first;
		Param const& param = pair.second;
		// TODO: Force writing the version parameter.
		if (param.value != nullptr) {
			try {
				param.type.write_root(*params_dir, name, *param.value);
			} catch (...) {
				// TODO: Don't catch every exception here.
				throw std::runtime_error(
					"Could not write parameter '" + name + "' to ROOT "
					"directory.");
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
		std::string const& name = pair.first;
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
						"Could not read parameter '" + name + "' from stream.");
				}
				param.used = false;
				std::string rem;
				std::getline(ss, rem);
				if (!rem.empty()) {
					throw std::runtime_error("");
				}
			} catch (std::exception const& e) {
				throw std::runtime_error(
					"Failed to parse parameter '" + name + "' from '" + ss.str()
					+ "'.");
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
		std::string const& name = pair.first;
		Param const& param = pair.second;
		// TODO: Force writing the version parameter.
		if (param.value != nullptr) {
			os << name << ' ';
			param.type.write_stream(os, *param.value);
			os << std::endl;
			if (!os) {
				throw std::runtime_error(
					"Could not write parameter '" + name + "' to stream.");
			}
		}
	}
}

Params Params::filter(Filter const& filter) {
	Params params;
	for (auto& pair : _params) {
		std::string const& name = pair.first;
		Param const& param = pair.second;
		if (filter.check(param.tags.begin(), param.tags.end())) {
			params._params.emplace(name, param);
		}
	}
	return params;
}

