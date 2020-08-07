#include <catch2/catch.hpp>

#include <cmath>
#include <limits>
#include <tuple>

#include <sidis/math.hpp>

using namespace sidis;

TEMPLATE_TEST_CASE(
		"Dilog test values", "[math]",
		float, double, long double) {
	using Pair = std::tuple<TestType, TestType>;
	auto target = GENERATE(
		Pair {   0.00,  0                                          },
		Pair {  -0.20, -0.1908001377775356190369131537660839924181 },
		Pair {   0.25,  0.2676526390827326069191838284878115758199 },
		Pair {  -0.35, -0.3233703936260739552366847546742059221798 },
		Pair {   0.40,  0.4492829744712816644647334023763193844553 },
		Pair {  -0.50, -0.4484142069236462024430644059157743208343 },
		Pair {   0.80,  1.074794600008248359395451922853995601055  },
		Pair {  -0.90, -0.7521631792172616203726927134268144689605 },
		Pair {   0.99,  1.588625448076375327031229473980552467945  },
		Pair {   1.00,  1.644934066848226436472415166646025189219  },
		Pair {   1.01,  1.700732144324037076288456396604084281626  },
		Pair {  -1.50, -1.147380660375570754079976633862792129215  },
		Pair {   2.00,  2.467401100272339654708622749969037783828  },
		Pair {  -4.00, -2.369939796998365831985537425350323048751  },
		Pair {   5.00,  1.783719161266630627743559734721650413496  },
		Pair { -16.00, -5.427008531068195691294900694545378583535  },
		Pair { 100.00, -7.323953199000482242723970498308855784205  });

	TestType x = std::get<0>(target);
	TestType y = std::get<1>(target);
	TestType y_test = math::dilog(x);
	INFO("x = " + std::to_string(x));
	CHECK(y_test == Approx(y).epsilon(10. * std::numeric_limits<TestType>::epsilon()));
}

