#include <iostream>

#include <sidis/sidis.hpp>

// Compute the QED coupling constant at a certain value of Q^2.
int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr << "Usage: alpha_qed <Q² (GeV²)>" << std::endl;
		return 1;
	}
	sidis::Real Q_sq = std::stold(argv[1]);
	std::cout << "α = " << sidis::ph::alpha_qed(Q_sq) << std::endl;
	return 0;
}

