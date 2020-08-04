#ifndef __SIDIS_STRUCTURE_FUNCTION_HPP__
#define __SIDIS_STRUCTURE_FUNCTION_HPP__

#include "sidis/numeric.hpp"

namespace sidis {
namespace sf {

struct SfUU {
	Real F_UUL;
	Real F_UUT;
	Real F_UU_cos_phih;
	Real F_UU_cos_2phih;
};

struct SfUL {
	Real F_UL_sin_phih;
	Real F_UL_sin_2phih;
};

struct SfUT {
	Real F_UTL_sin_phih_m_phis;
	Real F_UTT_sin_phih_m_phis;
	Real F_UT_sin_2phih_m_phis;
	Real F_UT_sin_3phih_m_phis;
	Real F_UT_sin_phis;
	Real F_UT_sin_phih_p_phis;
};

struct SfLU {
	Real F_LU_sin_phih;
};

struct SfLL {
	Real F_LL;
	Real F_LL_cos_phih;
};

struct SfLT {
	Real F_LT_cos_phih_m_phis;
	Real F_LT_cos_2phih_m_phis;
	Real F_LT_cos_phis;
};

struct SfXU {
	Real F_UUL;
	Real F_UUT;
	Real F_UU_cos_phih;
	Real F_UU_cos_2phih;
	Real F_LU_sin_phih;

	SfXU(SfUU sf_uu, SfLU sf_lu) :
		F_UUL(sf_uu.F_UUL),
		F_UUT(sf_uu.F_UUT),
		F_UU_cos_phih(sf_uu.F_UU_cos_phih),
		F_UU_cos_2phih(sf_uu.F_UU_cos_2phih),
		F_LU_sin_phih(sf_lu.F_LU_sin_phih) { }

	operator SfUU() const {
		return {
			F_UUL,
			F_UUT,
			F_UU_cos_phih,
			F_UU_cos_2phih,
		};
	}
	operator SfLU() const {
		return {
			F_LU_sin_phih,
		};
	}
};

struct SfXL {
	Real F_UL_sin_phih;
	Real F_UL_sin_2phih;
	Real F_LL;
	Real F_LL_cos_phih;

	SfXL(SfUL sf_ul, SfLL sf_ll) :
		F_UL_sin_phih(sf_ul.F_UL_sin_phih),
		F_UL_sin_2phih(sf_ul.F_UL_sin_2phih),
		F_LL(sf_ll.F_LL),
		F_LL_cos_phih(sf_ll.F_LL_cos_phih) { }

	operator SfUL() const {
		return {
			F_UL_sin_phih,
			F_UL_sin_2phih,
		};
	}
	operator SfLL() const {
		return {
			F_LL,
			F_LL_cos_phih,
		};
	}
};

struct SfXT {
	Real F_UTL_sin_phih_m_phis;
	Real F_UTT_sin_phih_m_phis;
	Real F_UT_sin_2phih_m_phis;
	Real F_UT_sin_3phih_m_phis;
	Real F_UT_sin_phis;
	Real F_UT_sin_phih_p_phis;
	Real F_LT_cos_phih_m_phis;
	Real F_LT_cos_2phih_m_phis;
	Real F_LT_cos_phis;

	SfXT(SfUT sf_ut, SfLT sf_lt) :
		F_UTL_sin_phih_m_phis(sf_ut.F_UTL_sin_phih_m_phis),
		F_UTT_sin_phih_m_phis(sf_ut.F_UTT_sin_phih_m_phis),
		F_UT_sin_2phih_m_phis(sf_ut.F_UT_sin_2phih_m_phis),
		F_UT_sin_3phih_m_phis(sf_ut.F_UT_sin_3phih_m_phis),
		F_UT_sin_phis(sf_ut.F_UT_sin_phis),
		F_UT_sin_phih_p_phis(sf_ut.F_UT_sin_phih_p_phis),
		F_LT_cos_phih_m_phis(sf_lt.F_LT_cos_phih_m_phis),
		F_LT_cos_2phih_m_phis(sf_lt.F_LT_cos_2phih_m_phis),
		F_LT_cos_phis(sf_lt.F_LT_cos_phis) { }

	operator SfUT() const {
		return {
			F_UTL_sin_phih_m_phis,
			F_UTT_sin_phih_m_phis,
			F_UT_sin_2phih_m_phis,
			F_UT_sin_3phih_m_phis,
			F_UT_sin_phis,
			F_UT_sin_phih_p_phis,
		};
	}
	operator SfLT() const {
		return {
			F_LT_cos_phih_m_phis,
			F_LT_cos_2phih_m_phis,
			F_LT_cos_phis,
		};
	}
};

struct SfUX {
	Real F_UUL;
	Real F_UUT;
	Real F_UU_cos_phih;
	Real F_UU_cos_2phih;
	Real F_UL_sin_phih;
	Real F_UL_sin_2phih;
	Real F_UTL_sin_phih_m_phis;
	Real F_UTT_sin_phih_m_phis;
	Real F_UT_sin_2phih_m_phis;
	Real F_UT_sin_3phih_m_phis;
	Real F_UT_sin_phis;
	Real F_UT_sin_phih_p_phis;

	SfUX(SfUU sf_uu, SfUL sf_ul, SfUT sf_ut) :
		F_UUL(sf_uu.F_UUL),
		F_UUT(sf_uu.F_UUT),
		F_UU_cos_phih(sf_uu.F_UU_cos_phih),
		F_UU_cos_2phih(sf_uu.F_UU_cos_2phih),
		F_UL_sin_phih(sf_ul.F_UL_sin_phih),
		F_UL_sin_2phih(sf_ul.F_UL_sin_2phih),
		F_UTL_sin_phih_m_phis(sf_ut.F_UTL_sin_phih_m_phis),
		F_UTT_sin_phih_m_phis(sf_ut.F_UTT_sin_phih_m_phis),
		F_UT_sin_2phih_m_phis(sf_ut.F_UT_sin_2phih_m_phis),
		F_UT_sin_3phih_m_phis(sf_ut.F_UT_sin_3phih_m_phis),
		F_UT_sin_phis(sf_ut.F_UT_sin_phis),
		F_UT_sin_phih_p_phis(sf_ut.F_UT_sin_phih_p_phis) { }

	operator SfUU() const {
		return {
			F_UUL,
			F_UUT,
			F_UU_cos_phih,
			F_UU_cos_2phih,
		};
	}
	operator SfUL() const {
		return {
			F_UL_sin_phih,
			F_UL_sin_2phih,
		};
	}
	operator SfUT() const {
		return {
			F_UTL_sin_phih_m_phis,
			F_UTT_sin_phih_m_phis,
			F_UT_sin_2phih_m_phis,
			F_UT_sin_3phih_m_phis,
			F_UT_sin_phis,
			F_UT_sin_phih_p_phis,
		};
	}
};

struct SfLX {
	Real F_LU_sin_phih;
	Real F_LL;
	Real F_LL_cos_phih;
	Real F_LT_cos_phih_m_phis;
	Real F_LT_cos_2phih_m_phis;
	Real F_LT_cos_phis;

	SfLX(SfLU sf_lu, SfLL sf_ll, SfLT sf_lt) :
		F_LU_sin_phih(sf_lu.F_LU_sin_phih),
		F_LL(sf_ll.F_LL),
		F_LL_cos_phih(sf_ll.F_LL_cos_phih),
		F_LT_cos_phih_m_phis(sf_lt.F_LT_cos_phih_m_phis),
		F_LT_cos_2phih_m_phis(sf_lt.F_LT_cos_2phih_m_phis),
		F_LT_cos_phis(sf_lt.F_LT_cos_phis) { }

	operator SfLU() const {
		return {
			F_LU_sin_phih,
		};
	}
	operator SfLL() const {
		return {
			F_LL,
			F_LL_cos_phih,
		};
	}
	operator SfLT() const {
		return {
			F_LT_cos_phih_m_phis,
			F_LT_cos_2phih_m_phis,
			F_LT_cos_phis,
		};
	}
};

struct Sf {
	Real F_UUL;
	Real F_UUT;
	Real F_UU_cos_phih;
	Real F_UU_cos_2phih;
	Real F_UL_sin_phih;
	Real F_UL_sin_2phih;
	Real F_UTL_sin_phih_m_phis;
	Real F_UTT_sin_phih_m_phis;
	Real F_UT_sin_2phih_m_phis;
	Real F_UT_sin_3phih_m_phis;
	Real F_UT_sin_phis;
	Real F_UT_sin_phih_p_phis;
	Real F_LU_sin_phih;
	Real F_LL;
	Real F_LL_cos_phih;
	Real F_LT_cos_phih_m_phis;
	Real F_LT_cos_2phih_m_phis;
	Real F_LT_cos_phis;

	Sf(SfUU sf_uu, SfUL sf_ul, SfUT sf_ut, SfLU sf_lu, SfLL sf_ll, SfLT sf_lt) :
		F_UUL(sf_uu.F_UUL),
		F_UUT(sf_uu.F_UUT),
		F_UU_cos_phih(sf_uu.F_UU_cos_phih),
		F_UU_cos_2phih(sf_uu.F_UU_cos_2phih),
		F_UL_sin_phih(sf_ul.F_UL_sin_phih),
		F_UL_sin_2phih(sf_ul.F_UL_sin_2phih),
		F_UTL_sin_phih_m_phis(sf_ut.F_UTL_sin_phih_m_phis),
		F_UTT_sin_phih_m_phis(sf_ut.F_UTT_sin_phih_m_phis),
		F_UT_sin_2phih_m_phis(sf_ut.F_UT_sin_2phih_m_phis),
		F_UT_sin_3phih_m_phis(sf_ut.F_UT_sin_3phih_m_phis),
		F_UT_sin_phis(sf_ut.F_UT_sin_phis),
		F_UT_sin_phih_p_phis(sf_ut.F_UT_sin_phih_p_phis),
		F_LU_sin_phih(sf_lu.F_LU_sin_phih),
		F_LL(sf_ll.F_LL),
		F_LL_cos_phih(sf_ll.F_LL_cos_phih),
		F_LT_cos_phih_m_phis(sf_lt.F_LT_cos_phih_m_phis),
		F_LT_cos_2phih_m_phis(sf_lt.F_LT_cos_2phih_m_phis),
		F_LT_cos_phis(sf_lt.F_LT_cos_phis) { }

	Sf(SfXU sf_xu, SfXL sf_xl, SfXT sf_xt) :
		F_UUL(sf_xu.F_UUL),
		F_UUT(sf_xu.F_UUT),
		F_UU_cos_phih(sf_xu.F_UU_cos_phih),
		F_UU_cos_2phih(sf_xu.F_UU_cos_2phih),
		F_UL_sin_phih(sf_xl.F_UL_sin_phih),
		F_UL_sin_2phih(sf_xl.F_UL_sin_2phih),
		F_UTL_sin_phih_m_phis(sf_xt.F_UTL_sin_phih_m_phis),
		F_UTT_sin_phih_m_phis(sf_xt.F_UTT_sin_phih_m_phis),
		F_UT_sin_2phih_m_phis(sf_xt.F_UT_sin_2phih_m_phis),
		F_UT_sin_3phih_m_phis(sf_xt.F_UT_sin_3phih_m_phis),
		F_UT_sin_phis(sf_xt.F_UT_sin_phis),
		F_UT_sin_phih_p_phis(sf_xt.F_UT_sin_phih_p_phis),
		F_LU_sin_phih(sf_xu.F_LU_sin_phih),
		F_LL(sf_xl.F_LL),
		F_LL_cos_phih(sf_xl.F_LL_cos_phih),
		F_LT_cos_phih_m_phis(sf_xt.F_LT_cos_phih_m_phis),
		F_LT_cos_2phih_m_phis(sf_xt.F_LT_cos_2phih_m_phis),
		F_LT_cos_phis(sf_xt.F_LT_cos_phis) { }

	Sf(SfUX sf_ux, SfLX sf_lx) :
		F_UUL(sf_ux.F_UUL),
		F_UUT(sf_ux.F_UUT),
		F_UU_cos_phih(sf_ux.F_UU_cos_phih),
		F_UU_cos_2phih(sf_ux.F_UU_cos_2phih),
		F_UL_sin_phih(sf_ux.F_UL_sin_phih),
		F_UL_sin_2phih(sf_ux.F_UL_sin_2phih),
		F_UTL_sin_phih_m_phis(sf_ux.F_UTL_sin_phih_m_phis),
		F_UTT_sin_phih_m_phis(sf_ux.F_UTT_sin_phih_m_phis),
		F_UT_sin_2phih_m_phis(sf_ux.F_UT_sin_2phih_m_phis),
		F_UT_sin_3phih_m_phis(sf_ux.F_UT_sin_3phih_m_phis),
		F_UT_sin_phis(sf_ux.F_UT_sin_phis),
		F_UT_sin_phih_p_phis(sf_ux.F_UT_sin_phih_p_phis),
		F_LU_sin_phih(sf_lx.F_LU_sin_phih),
		F_LL(sf_lx.F_LL),
		F_LL_cos_phih(sf_lx.F_LL_cos_phih),
		F_LT_cos_phih_m_phis(sf_lx.F_LT_cos_phih_m_phis),
		F_LT_cos_2phih_m_phis(sf_lx.F_LT_cos_2phih_m_phis),
		F_LT_cos_phis(sf_lx.F_LT_cos_phis) { }

	operator SfXU() const {
		return {
			{
				F_UUL,
				F_UUT,
				F_UU_cos_phih,
				F_UU_cos_2phih,
			},
			{
				F_LU_sin_phih,
			},
		};
	}
	operator SfXL() const {
		return {
			{
				F_UL_sin_phih,
				F_UL_sin_2phih,
			},
			{
				F_LL,
				F_LL_cos_phih
			},
		};
	}
	operator SfXT() const {
		return {
			{
				F_UTL_sin_phih_m_phis,
				F_UTT_sin_phih_m_phis,
				F_UT_sin_2phih_m_phis,
				F_UT_sin_3phih_m_phis,
				F_UT_sin_phis,
				F_UT_sin_phih_p_phis,
			},
			{
				F_LT_cos_phih_m_phis,
				F_LT_cos_2phih_m_phis,
				F_LT_cos_phis,
			},
		};
	}
	operator SfUX() const {
		return {
			{
				F_UUL,
				F_UUT,
				F_UU_cos_phih,
				F_UU_cos_2phih,
			},
			{
				F_UL_sin_phih,
				F_UL_sin_2phih,
			},
			{
				F_UTL_sin_phih_m_phis,
				F_UTT_sin_phih_m_phis,
				F_UT_sin_2phih_m_phis,
				F_UT_sin_3phih_m_phis,
				F_UT_sin_phis,
				F_UT_sin_phih_p_phis,
			},
		};
	}
	operator SfLX() const {
		return {
			{
				F_LU_sin_phih,
			},
			{
				F_LL,
				F_LL_cos_phih,
			},
			{
				F_LT_cos_phih_m_phis,
				F_LT_cos_2phih_m_phis,
				F_LT_cos_phis,
			},
		};
	}
	operator SfUU() const {
		return {
			F_UUL,
			F_UUT,
			F_UU_cos_phih,
			F_UU_cos_2phih,
		};
	}
	operator SfUL() const {
		return {
			F_UL_sin_phih,
			F_UL_sin_2phih,
		};
	}
	operator SfUT() const {
		return {
			F_UTL_sin_phih_m_phis,
			F_UTT_sin_phih_m_phis,
			F_UT_sin_2phih_m_phis,
			F_UT_sin_3phih_m_phis,
			F_UT_sin_phis,
			F_UT_sin_phih_p_phis,
		};
	}
	operator SfLU() const {
		return {
			F_LU_sin_phih,
		};
	}
	operator SfLL() const {
		return {
			F_LL,
			F_LL_cos_phih,
		};
	}
	operator SfLT() const {
		return {
			F_LT_cos_phih_m_phis,
			F_LT_cos_2phih_m_phis,
			F_LT_cos_phis,
		};
	}
};

}
}

#endif

