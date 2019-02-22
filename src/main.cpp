#include <iostream>
#include <getopt.h>

int main(int argc, char *argv[]) {
	int c;
	int index;
	char *cvalue = nullptr;

	while ((c = getopt(argc, argv, "")) != -1) {
		switch(c) {
		default:
			std::cout << "Unknown parameter: " << (char)c << std::endl;
			std::terminate();
		}
	}
	return 0;
}
