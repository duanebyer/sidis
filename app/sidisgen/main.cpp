#include <iostream>

#include <TCanvas.h>
#include <TF1.h>

#include <sidis/sidis.hpp>

using namespace sidis;

int main(int argc, char** argv) {
	char const* func_def = "sin(x)/x";
	std::cout << "Plotting a function " << func_def << std::endl << std::endl;
	TF1 f1("f1", func_def, -10., 10.);
	TCanvas c;
	f1.Draw();
	c.Print("f1.png");
	return 0;
}

