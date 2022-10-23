#ifndef SIDISGEN_PARAMS_HPP
#define SIDISGEN_PARAMS_HPP

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <map>
#include <set>
#include <typeinfo>
#include <vector>

// Thoughts:
// * Cuts can be in gen or in init.

// TODO:
// Factor this out into its own lib. Remove the ROOT and stream stuff from Type.
// Instead, allow things like `write_stream` to be implemented in user code.
// * Make `Params` provide an iterator over (Name, Type, Value) tuples. This can
//   be transformed to get names(), values().

class TDirectory;
class TObject;
class Value;

// Type that can be constructed from any value and implicitly converted to any
// value.
class Any final {
	std::shared_ptr<void> _val;
	std::type_info const& _type;

public:
	Any() : _val(), _type(typeid(void)) { }
	template<typename T>
	Any(T const& val) : _val(new T(val)), _type(typeid(T)) { }

	template<typename T>
	T as() const {
		std::type_info const& other_type = typeid(
			typename std::remove_const<
				typename std::remove_reference<
					typename std::remove_const<T>::type >::type >::type);
		if (_type != other_type) {
			throw std::runtime_error(
				"Invalid type cast from " + std::string(_type.name()) + " to "
				+ std::string(other_type.name()) + ".");
		} else {
			return *static_cast<T const*>(_val.get());
		}
	}
	template<typename T>
	operator T() const {
		return this->as<T>();
	}
};

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

	// Read/write to ROOT files.
	virtual std::unique_ptr<Value> read_root(TDirectory& dir, std::string const& name) const = 0;
	virtual void write_root(TDirectory& dir, std::string const& name, Value const& value) const = 0;
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

	// Optional virtual function that can be overridden. The returned `Any` type
	// can be implicitly cast to anything with a dynamic type check, so this can
	// reduce boiler plate.
	virtual Any any() const {
		return Any();
	}

	// Convenience methods for easily calling methods from `Type`.
	std::string to_string() const {
		return _type.to_string(*this);
	}
	bool operator==(Value const& other) const {
		if (_type != other._type) {
			return false;
		} else {
			return _type.equivalent(*this, other);
		}
	}
	bool operator!=(Value const& other) const {
		return !(*this == other);
	}
};

class Filter final {
	// The tag is wrapped in this struct to make it easier to add negation in
	// the future, if desired.
	struct Factor {
		std::string tag;
		// For now, inverse conditions are not allowed.
		// bool invert;
		bool operator<(Factor const& rhs) const {
			return tag < rhs.tag;
		}
	};
	using Term = std::set<Factor>;
	using Condition = std::set<Term>;

	// Condition on tags, in the form:
	//     `(A1 && A2 && ...) || (B1 && B2 && ...) || ...`
	Condition _condition;

	Filter() : _condition{} { };
	explicit Filter(Factor const& factor) : _condition{ { factor } } { };
	explicit Filter(Term const& term) : _condition{ term } { };

public:
	static const Filter REJECT;
	static const Filter ACCEPT;

	explicit Filter(std::string const& tag) : Filter(Factor{ tag }) { };

	template<typename It>
	bool check(It begin, It end) const {
		bool accept_condition = false;
		for (Term const& term : _condition) {
			bool accept_term = true;
			for (Factor const& factor : term) {
				if (std::find(begin, end, factor.tag) == end) {
					accept_term = false;
					break;
				}
			}
			if (accept_term) {
				accept_condition = true;
				break;
			}
		}
		return accept_condition;
	}

	Filter& operator|=(Filter const& rhs);
	Filter& operator&=(Filter const& rhs);
	friend Filter operator|(Filter lhs, Filter const& rhs) {
		lhs |= rhs;
		return lhs;
	}
	friend Filter operator&(Filter lhs, Filter const& rhs) {
		lhs &= rhs;
		return lhs;
	}
};

inline Filter operator""_F(char const* tag, std::size_t len) {
	return Filter(std::string(tag, len));
}

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
			std::vector<std::string> tags,
			std::string usage,
			std::string brief,
			std::string doc);
		Param(
			Type const& type,
			std::vector<std::string> tags,
			std::string usage,
			std::string brief,
			std::string doc);
	};
	std::map<std::string, Param> _params;

public:

	// Get parameter type.
	Type const& type(std::string const& name) const {
		return _params.at(name).type;
	}
	// Get parameter metadata.
	std::string const& usage(std::string const& name) const {
		return _params.at(name).usage;
	}
	std::string const& brief(std::string const& name) const {
		return _params.at(name).brief;
	}
	std::string const& doc(std::string const& name) const {
		return _params.at(name).doc;
	}
	// Gets the value provided by the parameter. If parameter is empty, errors.
	Value const& get(std::string const& name);
	template<typename T>
	T const& get(std::string const& name) {
		return get(name).as<T>();
	}
	Value const& operator[](std::string const& name) {
		return get(name);
	}
	// Gets a value without actually using it. Useful for checks on a parameter
	// value that don't use the parameter to do anything useful.
	Value const& get_soft(std::string const& name) const;
	template<typename T>
	T const& get_soft(std::string const& name) const {
		return get_soft(name).as<T>();
	}
	// Sets a value provided by the parameter. Overwrites any existing
	// parameter, returning whether it did so.
	bool set(std::string const& name, Value const* value);
	// Sets with a value from another set of parameters. Checks tags, types, and
	// default values before assigning. Overwrites any existing parameter,
	// returning whether it did so.
	bool set_from(Params const& other, std::string const& name);
	// Sets with all values from another set of parameters. Checks tags, types,
	// and default values before assigning. Overwrites any existing parameters,
	// returning whether it did so.
	bool set_from(Params const& other);
	// Check if parameter has been set.
	bool is_set(std::string const& name) const {
		return _params.at(name).value != nullptr;
	}
	// Checks whether the parameter has been used.
	bool used(std::string const& name) const {
		return _params.at(name).used;
	}

	// Add new parameter.
	void add_param(
			std::string name,
			Value const* default_value,
			std::vector<std::string> tags,
			std::string usage,
			std::string brief,
			std::string doc) {
		if (!_params.emplace(name, Param(default_value, tags, usage, brief, doc)).second) {
			throw std::runtime_error(
				"Tried to add existing parameter '" + name + "'.");
		}
	}
	void add_param(
			std::string name, 
			Type const& type,
			std::vector<std::string> tags,
			std::string usage,
			std::string brief,
			std::string doc) {
		if (!_params.emplace(name, Param(type, tags, usage, brief, doc)).second) {
			throw std::runtime_error(
				"Tried to add existing parameter '" + name + "'.");
		}
	}

	// Gets a set of all parameter names.
	// TODO: Replace this with a proper iterator over `ParamView`, when the
	// parameter system is moved into its own package.
	std::set<std::string> names() const;

	// Checks whether two sets of parameters have the same names, types, tags,
	// and default values. Meta-data is not considered.
	void check_format(Params const& other) const;
	// Checks whether all provided parameters have been used.
	void check_complete() const;
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
	Params filter(Filter const& filter);
};

#endif

