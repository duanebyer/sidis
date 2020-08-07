#include "sidis/sf_model/ww.hpp"

#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <utility>

#include <wstp.h>

#if WINDOWS_WSTP
#define WS_YIELD_CALLING_CONVENTION WINAPI
#else
#define WS_YIELD_CALLING_CONVENTION
#endif

#define SF_MODEL_DIR "./sidis/sf_model/"
#define WWSIDIS_FILE "./wwsidis/wwsidis.m"

using namespace sidis;
using namespace sidis::sf;
using namespace sidis::sf::model;

namespace {

struct TimeoutData {
	std::chrono::steady_clock::time_point timeout_begin;
	unsigned timeout;
};

// Timeout checking through the WSTP "yield function" feature.
int WS_YIELD_CALLING_CONVENTION link_yield(
		WSLINK link,
		WSYieldParameters param) {
	WSUserFunction f;
	TimeoutData* data = reinterpret_cast<TimeoutData*>(WSUserData(link, &f));
	unsigned duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - data->timeout_begin).count();
	if (duration > data->timeout) {
		return 1;
	} else {
		return 0;
	}
}

// Mandatory cleanup function for WSTP user data.
void link_data_cleanup(WSLINK link) {
}

}

struct WW::Impl {
	WSENV env = 0;
	WSLINK link = 0;
	TimeoutData timeout_data;
	unsigned timeout_init;
	unsigned timeout_packet;

	// These functions are called before any operations that wait on the WSTP
	// link, so that the timeout checks can be initialized correctly.
	void start_timeout_init() {
		timeout_data.timeout_begin = std::chrono::steady_clock::now();
		timeout_data.timeout = timeout_init;
	}
	void start_timeout_packet() {
		timeout_data.timeout_begin = std::chrono::steady_clock::now();
		timeout_data.timeout = timeout_packet;
	}

	// Make calls through the WSTP link.
	void load_library(char const* library_path);
	double evaluate_sf(
			char const* sf_name,
			double x,
			double z,
			double Q_sq,
			double ph_t);
};

WW::WW(WW&& other) noexcept : _impl(std::exchange(other._impl, nullptr)) { }
WW& WW::operator=(WW&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

WW::WW(unsigned timeout_packet, unsigned timeout_init) {
	// Initialize the WSTP environment and open a link to the kernel.
	int error;
	int packet;
	_impl = new Impl();
	_impl->timeout_packet = timeout_packet;
	_impl->timeout_init = timeout_init;
	_impl->env = WSInitialize((WSEnvironmentParameter) 0);
	if (_impl->env == (WSENV) 0) {
		throw WW::EnvironmentInitException();
	}
	_impl->link = WSOpenString(
		_impl->env,
		"-linkmode launch -linkname 'math -wstp'",
		&error);
	if (_impl->link == (WSLINK) 0 || error != WSEOK) {
		_impl->link = 0;
		throw WW::LinkInitException(error);
	}

	// Attach user data to the link, and create the timeout callback.
	WSSetUserData(
		_impl->link,
		reinterpret_cast<void*>(&(_impl->timeout_data)),
		link_data_cleanup);
	if (!WSSetYieldFunction(_impl->link, link_yield)) {
		error = WSError(_impl->link);
		throw WW::TimeoutInitException(error);
	}

	// Read the first input packet which is automatically generated.
	_impl->start_timeout_packet();
	packet = WSNextPacket(_impl->link);
	if ((error = WSError(_impl->link)) != WSEOK) {
		throw WW::ReceivePacketException(error);
	}
	if (packet != INPUTNAMEPKT) {
		throw WW::UnexpectedPacketException(packet);
	}
	if (!WSNewPacket(_impl->link)) {
		error = WSError(_impl->link);
		throw WW::NextPacketFailureException(error);
	}

	// Load the library. Try several possible paths.
	char const* library_path_1 = DATADIR SF_MODEL_DIR WWSIDIS_FILE;
	char const* library_path_2 = "../share/" SF_MODEL_DIR WWSIDIS_FILE;
	char const* library_path_4 = SF_MODEL_DIR WWSIDIS_FILE;
	char const* library_path_3 = WWSIDIS_FILE;
	char const* library_path_all =
		DATADIR SF_MODEL_DIR WWSIDIS_FILE ":"
		"../share/" SF_MODEL_DIR WWSIDIS_FILE ":"
		SF_MODEL_DIR WWSIDIS_FILE ":"
		WWSIDIS_FILE;
	try {
		_impl->load_library(library_path_1);
	} catch (WW::LibraryLoadException const&) {
		try {
			_impl->load_library(library_path_2);
		} catch (WW::LibraryLoadException const&) {
			try {
				_impl->load_library(library_path_3);
			} catch(WW::LibraryLoadException const&) {
				try {
					_impl->load_library(library_path_4);
				} catch(WW::LibraryLoadException const&) {
					throw WW::LibraryLoadException(library_path_all);
				}
			}
		}
	}
}

WW::~WW() {
	if (_impl->link != 0) {
		WSClose(_impl->link);
	}
	if (_impl->env != 0) {
		WSDeinitialize(_impl->env);
	}
	if (_impl != nullptr) {
		delete _impl;
	}
}

SfUU WW::sf_uu(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		0.,
		_impl->evaluate_sf(
			"wwsidis`Private`FUU",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUUcosphi",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUUcos2phi",
			x, z, Q_sq, ph_t),
	};
}

SfUL WW::sf_ul(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		_impl->evaluate_sf(
			"wwsidis`Private`FULsinphi",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FULsin2phi",
			x, z, Q_sq, ph_t),
	};
}

SfUT WW::sf_ut(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		0.,
		_impl->evaluate_sf(
			"wwsidis`Private`FUTf1tperp",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUTsin2phi",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUTh1tp",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUTsinphiS",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FUTh1",
			x, z, Q_sq, ph_t),
	};
}

SfLU WW::sf_lu(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		0.,
	};
}

SfLL WW::sf_ll(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		_impl->evaluate_sf(
			"wwsidis`Private`FLL",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FLLcosphi",
			x, z, Q_sq, ph_t),
	};
}

SfLT WW::sf_lt(Real x, Real z, Real Q_sq, Real ph_t) const {
	return {
		_impl->evaluate_sf(
			"wwsidis`Private`FLT",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FLTcos2phi",
			x, z, Q_sq, ph_t),
		_impl->evaluate_sf(
			"wwsidis`Private`FLTcosPhi",
			x, z, Q_sq, ph_t),
	};
}

void WW::Impl::load_library(char const* library_path) {
	int error;
	int packet;
	start_timeout_init();
	// Load the library into the kernel.
	WSPutFunction(link, "EvaluatePacket", 1);
	WSPutFunction(link, "Quiet", 1);
	WSPutFunction(link, "Block", 2);
	WSPutFunction(link, "List", 2);
	WSPutSymbol(link, "$ContextPath");
	WSPutSymbol(link, "Print");
	WSPutFunction(link, "Needs", 2);
	WSPutString(link, "wwsidis`");
	WSPutString(link, library_path);
	WSEndPacket(link);
	if ((error = WSError(link)) != WSEOK) {
		throw WW::SendPacketException(error);
	}

	// Receive the response to determine whether the library was loaded
	// correctly.
	packet = WSNextPacket(link);
	if ((error = WSError(link)) != WSEOK) {
		throw WW::ReceivePacketException(error);
	}
	if (packet != RETURNPKT) {
		throw WW::UnexpectedPacketException(packet);
	}
	char const* symbol;
	if (!WSGetSymbol(link, &symbol)) {
		error = WSError(link);
		throw WW::UnexpectedPacketContentsException(error);
	}
	if (std::strcmp(symbol, "$Failed") == 0) {
		WSReleaseSymbol(link, symbol);
		throw WW::LibraryLoadException(library_path);
	}
	if (std::strcmp(symbol, "Null") != 0) {
		WSReleaseSymbol(link, symbol);
		throw WW::UnexpectedPacketContentsException(0);
	}
	WSReleaseSymbol(link, symbol);
	if (!WSNewPacket(link)) {
		error = WSError(link);
		throw WW::NextPacketFailureException(error);
	}
}

double WW::Impl::evaluate_sf(
		char const* sf_name,
		double x,
		double z,
		double Q_sq,
		double ph_t) {
	int error;
	int packet;
	start_timeout_packet();
	// Check that we are in bounds.
	if (!std::isfinite(x)
			|| !std::isfinite(z)
			|| !std::isfinite(Q_sq)
			|| !std::isfinite(ph_t)) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(x > 0. && x <= 1.)
			|| !(z > 0. && z <= 1.)
			|| !(Q_sq > 0.)
			|| !(ph_t > 0.)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	// Send packet.
	WSPutFunction(link, "EvaluatePacket", 1);
	WSPutFunction(link, "Quiet", 1);
	WSPutFunction(link, sf_name, 5);
	WSPutString(link, "pi+");
	WSPutDouble(link, x);
	WSPutDouble(link, z);
	WSPutDouble(link, Q_sq);
	WSPutDouble(link, ph_t);
	WSEndPacket(link);
	if ((error = WSError(link)) != WSEOK) {
		throw WW::SendPacketException(error);
	}

	// Receive response.
	double value;
	packet = WSNextPacket(link);
	if ((error = WSError(link)) != WSEOK) {
		throw WW::ReceivePacketException(error);
	}
	if (packet != RETURNPKT) {
		throw WW::UnexpectedPacketException(packet);
	}
	if (!WSGetReal(link, &value)) {
		error = WSError(link);
		throw WW::UnexpectedPacketContentsException(error);
	}
	if (!WSNewPacket(link)) {
		error = WSError(link);
		throw WW::NextPacketFailureException(error);
	}

	return value;
}

WW::EnvironmentInitException::EnvironmentInitException() :
	std::runtime_error("Failed to initialize WSTP environment") { }
WW::LinkInitException::LinkInitException(int code) :
	std::runtime_error(
		std::string("Failed to initialize WSTP link (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::TimeoutInitException::TimeoutInitException(int code) :
	std::runtime_error(
		std::string("Failed to initialize WSTP timeout (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::LibraryLoadException::LibraryLoadException(char const* filepath) :
	std::runtime_error(
		std::string("Failed to load the WW-SIDIS library from \"")
		+ std::string(filepath)
		+ std::string("\"")),
	filepath(filepath) { }
WW::SendPacketException::SendPacketException(int code) :
	std::runtime_error(
		std::string("Failed to send packet to WSTP link (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::ReceivePacketException::ReceivePacketException(int code) :
	std::runtime_error(
		std::string("Failed to receive packet from WSTP link (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::UnexpectedPacketException::UnexpectedPacketException(int code) :
	std::runtime_error(
		std::string("Received unexpected packet from WSTP link (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::UnexpectedPacketContentsException::UnexpectedPacketContentsException(int code) :
	std::runtime_error(
		std::string("Found unexpected contents in packet (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }
WW::NextPacketFailureException::NextPacketFailureException(int code) :
	std::runtime_error(
		std::string("Failed to read the next packet (error ")
		+ std::to_string(code)
		+ std::string(")")),
	code(code) { }

