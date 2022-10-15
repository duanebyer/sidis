#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <map>
#include <set>

// Thoughts:
// * move physics checks out into main, and update how they work a bit.
//   * e.g. k_0_bar_cut (with min above soft threshold) should be a
//     pre-requisite for any other radiative cut, to reduce mistakes.
//   * k_0_bar_cut minimum should be either 0 or above soft threshold, if above
//     soft threshold, no generating nrad events
//   * Option to apply radiative cuts exclusively to radiative events above
//     `k_0_bar` threshold (still let non-radiative events and events below
//     threshold through), or to all events (remove non-radiative events
//     entirely).
// * Hashes for FOAMs. Parameter has an 'any' option, that is used to ignore the
//   hash by default. After generating, it gets filled in though, just like the
//   random seeds.
// * TBH specifying a seed is so useless for initialization phase because of
//   multithreading, just remove the capability and replace it with a hash.
// * Cuts can be in gen or in init.
// * Don't use initializer_list. It's bad

class TDirectory;
class TObject;
class Value;

// Abstract base class for types of values that can be stored in parameters.
class Type {
protected:
	Type() { }
public:
	// Prevent copying.
	Type(Type const&) = delete;
	virtual ~Type() { }
	Type& operator=(Type const&) = delete;

	// Whether two values are equivalent with each other. It's only valid to
	// call this on two values that have this as their type.
	virtual bool equivalent(Value const& value_1, Value const& value_2) const = 0;
	// Whether an 'object' from a 'producing' parameter file can be used by an
	// 'object' from a 'consuming' parameter file. It's only valid to call this
	// on two values that have this as their type.
	virtual bool can_provide(Value const& prod, Value const& cons) const = 0;

	// Read/write to ROOT files.
	virtual std::unique_ptr<Value> read_root(TDirectory& dir, char const* name) const = 0;
	virtual void write_root(TDirectory& dir, char const* name, Value const& value) const = 0;
	// Read/write to streams.
	virtual std::unique_ptr<Value> read_stream(std::istream& is) const = 0;
	virtual void write_stream(std::ostream& os, Value const& value) const = 0;

	std::string to_string(Value const& value) const;

	// Singleton equality.
	bool operator==(Type const& other) const {
		return this == &other;
	}
	bool operator!=(Type const& other) const {
		return !(*this == other);
	}
};

// Abstract base class for values that can be stored in parameters.
class Value {
private:
	Type const& _type;

protected:
	Value(Type const& type) : _type(type) { }

public:
	virtual ~Value() { }

	Type const& type() const {
		return _type;
	}

	// Convenience methods for quickly casting to a subtype.
	template<typename T>
	T const& as() const {
		return dynamic_cast<T const&>(*this);
	}
	template<typename T>
	T& as() {
		return dynamic_cast<T&>(*this);
	}
};

// Stores a collection of parameters that can be read from.
class Params final {
	struct Param final {
		// Type of the parameter.
		Type const& type;
		// Value provided by the parameter.
		std::shared_ptr<Value const> value;
		// Default value provided by the parameter when none other is available.
		std::shared_ptr<Value const> default_value;

		// Tags for filtering and selecting parameters.
		std::set<std::string> tags;

		// Whether the parameter has been read from.
		bool used;

		// Parameter metadata.
		std::string usage;
		std::string brief;
		std::string doc;

		Param(
			Value const* default_value,
			std::initializer_list<char const*> tags,
			char const* usage,
			char const* brief,
			char const* doc);
		Param(
			Type const& type,
			std::initializer_list<char const*> tags,
			char const* usage,
			char const* brief,
			char const* doc);
	};
	std::map<std::string, Param> _params;

public:

	// Get parameter type.
	Type const& type(char const* name) const {
		return _params.at(name).type;
	}
	// Get parameter metadata.
	char const* usage(char const* name) const {
		return _params.at(name).usage.c_str();
	}
	char const* brief(char const* name) const {
		return _params.at(name).brief.c_str();
	}
	char const* doc(char const* name) const {
		return _params.at(name).doc.c_str();
	}
	// Gets the value provided by the parameter. If parameter is empty, errors.
	Value const& get(char const* name);
	template<typename T>
	T const& get(char const* name) {
		return get(name).as<T>();
	}
	// Sets a value provided by the parameter. Overwrites any existing
	// parameter, returning whether it did so.
	bool set(char const* name, Value const* value);
	// Checks whether the parameter has been used.
	bool used(char const* name) const {
		return _params.at(name).used;
	}

	// Add new parameter.
	void add_param(
			char const* name,
			Value const* default_value,
			std::initializer_list<char const*> tags,
			char const* usage,
			char const* brief,
			char const* doc) {
		if (!_params.emplace(name, Param(default_value, tags, usage, brief, doc)).second) {
			throw std::runtime_error(
				std::string("Tried to add existing parameter '") + name + "'.");
		}
	}
	void add_param(
			char const* name, 
			Type const& type,
			std::initializer_list<char const*> tags,
			char const* usage,
			char const* brief,
			char const* doc) {
		if (!_params.emplace(name, Param(type, tags, usage, brief, doc)).second) {
			throw std::runtime_error(
				std::string("Tried to add existing parameter '") + name + "'.");
		}
	}

	// Gets a set of all parameter names.
	std::set<std::string> names() const;

	// Checks whether two sets of parameters have the same names + types + tags.
	void check_format(Params const& other) const;
	// Checks whether all provided parameters have been used.
	void check_complete() const;
	// Checks whether this set of parameters can provide objects for another set
	// of parameters to consume.
	void check_can_provide(Params const& consumer) const;
	// Checks whether this set of parameters is equal to another set of
	// parameters.
	void check_equivalent(Params const& other) const;

	// Clears unused parameters.
	void clear_unused();

	// Read/write parameters from ROOT file.
	void read_root(TDirectory& dir);
	void write_root(TDirectory& dir) const;
	// Read/write parameters from stream.
	void read_stream(std::istream& is);
	void write_stream(std::ostream& os) const;

	// Returns only the parameters matching the provided matcher. The matcher
	// can be a regex applied to both names and tags of parameters.
	Params filter(std::initializer_list<char const*> matches);
};

#endif

