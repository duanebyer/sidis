#include <catch2/catch.hpp>

#include <cmath>
#include <limits>
#include <sstream>
#include <tuple>

#include <sidis/math.hpp>

#include "rel_matcher.hpp"

using namespace sidis;

TEMPLATE_TEST_CASE(
		"√(1+x)-1 test values", "[math]",
		float, double, long double) {
	using Pair = std::tuple<TestType, TestType>;
	auto target = GENERATE(
		Pair {    0.L      ,   0L                                             },
		Pair {   -1.L      ,  -1.000000000000000000000000000000000000000L     },
		Pair {    1.L      ,   0.4142135623730950488016887242096980785697L    },
		Pair {   -0.5L     ,  -0.2928932188134524755991556378951509607152L    },
		Pair {    0.5L     ,   0.2247448713915890490986420373529456959830L    },
		Pair {   -1.25e-10L,  -6.250000000195312500012207031250953674316e-11L },
		Pair {    0.89e-10L,   4.449999999900987500004406056249754913121e-11L },
		Pair {    5.L      ,   1.449489742783178098197284074705891391966L     },
		Pair {    6.5L     ,   1.738612787525830567284848914004010669764L     },
		Pair {   10.L      ,   2.316624790355399849114932736670686683927L     },
		Pair {   25.L      ,   4.099019513592784830028224109022781989564L     },
		Pair {  100.L      ,   9.049875621120890270219264912759576186945L     },
		Pair { 2000.L      ,  43.73253849269008341597431497806218017572L      });

	TestType x = std::get<0>(target);
	TestType y = std::get<1>(target);
	TestType y_test = math::sqrt1p_1m(x);
	std::stringstream ss;
	ss << "x = " << x;
	INFO(ss.str());
	CHECK_THAT(
		y_test,
		RelMatcher<TestType>(y, 10. * std::numeric_limits<TestType>::epsilon()));
}

TEMPLATE_TEST_CASE(
		"Dilog test values", "[math]",
		float, double, long double) {
	using Pair = std::tuple<TestType, TestType>;
	auto target = GENERATE(
		Pair {   0.00L,  0L                                          },
		Pair {  -0.20L, -0.1908001377775356190369131537660839924181L },
		Pair {   0.25L,  0.2676526390827326069191838284878115758199L },
		Pair {  -0.35L, -0.3233703936260739552366847546742059221798L },
		Pair {   0.40L,  0.4492829744712816644647334023763193844553L },
		Pair {  -0.50L, -0.4484142069236462024430644059157743208343L },
		Pair {   0.80L,  1.074794600008248359395451922853995601055L  },
		Pair {  -0.90L, -0.7521631792172616203726927134268144689605L },
		Pair {   0.99L,  1.588625448076375327031229473980552467945L  },
		Pair {   1.00L,  1.644934066848226436472415166646025189219L  },
		Pair {   1.01L,  1.700732144324037076288456396604084281626L  },
		Pair {  -1.50L, -1.147380660375570754079976633862792129215L  },
		Pair {   2.L  ,  2.467401100272339654708622749969037783828L  },
		Pair {  -4.L  , -2.369939796998365831985537425350323048751L  },
		Pair {   5.L  ,  1.783719161266630627743559734721650413496L  },
		Pair { -16.L  , -5.427008531068195691294900694545378583535L  },
		Pair { 100.L  , -7.323953199000482242723970498308855784205L  });

	TestType x = std::get<0>(target);
	TestType y = std::get<1>(target);
	TestType y_test = math::dilog(x);
	std::stringstream ss;
	ss << "x = " << x;
	INFO(ss.str());
	CHECK_THAT(
		y_test,
		RelMatcher<TestType>(y, 10. * std::numeric_limits<TestType>::epsilon()));
}

