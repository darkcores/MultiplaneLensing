#include <iostream>
#include <getopt.h>

void print_help() {
	std::cout << "Program arguments" << std::endl;
	std::cout << "\t-i\tLens backprojecten" << std::endl;
	std::cout << "\t-d\tSelect cuda device" << std::endl;
	std::cout << "\t-l\tLens input file" << std::endl;
	std::cout << "\t-s\tSources input file" << std::endl;
}

int main(int argc, char *argv[]) {
	int c;
	int index;
	char *cvalue = nullptr;

	while ((c = getopt(argc, argv, "id:hl::s::")) != -1) {
		switch(c) {
		case 'i':
			std::cout << "Using inversion lens" << std::endl;
			break;
		case 'd':
			std::cout << "Cuda device: TODO" << std::endl;
			break;
		case 'l':
			std::cout << "Load lens: TODO" << std::endl;
			break;
		case 's':
			std::cout << "Load sources: TODO" << std::endl;
			break;
		case 'h':
			print_help();
			break;
		default:
			std::cout << "Unknown parameter: " << (char)c << std::endl;
			std::terminate();
		}
	}
	return 0;
}
