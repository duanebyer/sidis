#ifndef SIDISGEN_EVENT_TYPE_HPP
#define SIDISGEN_EVENT_TYPE_HPP

// Enumeration for all types of supported events.
enum class EventType {
	// Non-radiative (without photon emission).
	NRAD,
	// Radiative
	RAD,
	// Exclusive
	EXCL,
};

int const NUM_EVENT_TYPES = 3;

// Short identifying names for each type of event.
inline char const* event_type_short_name(EventType type) {
	switch(type) {
	case EventType::NRAD:
		return "nrad";
	case EventType::RAD:
		return "rad";
	case EventType::EXCL:
		return "excl";
	default:
		return "<error>";
	}
}

// Longer identifying names (ex. for error message) for each type of event.
inline char const* event_type_name(EventType type) {
	switch (type) {
	case EventType::NRAD:
		return "non-radiative";
	case EventType::RAD:
		return "radiative";
	case EventType::EXCL:
		return "exclusive";
	default:
		return "<error>";
	}
}

#endif

