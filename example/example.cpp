#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>

// Multiplane CUDA API header
#include <context.h>

/**
 * Print help & options.
 */
void print_help();
/**
 * Parse program arguments.
 */
void parse_args(int argc, char *argv[]);
/**
 * Load data from input file.
 */
void load_file();
/**
 * Write calculated positions to file.
 */
void write_file(MultiPlaneContext *);

// Variables used in this file.
// To keep this nice and readable.

double angularUnit = 0;
// Lens redshifts.
std::vector<float> lens_z;
// Lens parameters.
std::vector<std::vector<PlummerParams>> lens_params;
// Source plane redshifts.
std::vector<float> source_z;
// Thetas
std::vector<Vector2D<float>> thetas;
// Mass params
std::vector<std::vector<float>> masses;

// Filenames
char *infile = nullptr;
char *outfile = nullptr;

void print_help() {
    std::cout << "Example program for calculating beta vectors"
              << "with multiple lens/source planes." << std::endl
              << std::endl;
    std::cout << "Usage: ./example -i data.txt -o points.txt" << std::endl
              << std::endl;
    std::cout << "Program arguments: " << std::endl;
    std::cout << "\t-d\tSelect cuda device" << std::endl;
    std::cout << "\t-i\tInput file" << std::endl;
    std::cout << "\t-o\tOutput file" << std::endl;
}

void parse_args(int argc, char *argv[]) {
    int c;

    while ((c = getopt(argc, argv, "d:hi:o:")) != -1) {
        switch (c) {
        case 'd': // Device select
            std::cout << "Cuda device: TODO" << std::endl;
            break;
        case 'i': // Input
            std::cout << "Load file: " << optarg << std::endl;
            infile = optarg;
            break;
        case 'o': // Output
            std::cout << "Output file: " << optarg << std::endl;
            outfile = optarg;
            break;
        case 'h': // Help
            print_help();
            std::exit(0);
            break;
        default:
            std::cout << "Unknown parameter: " << (char)c << std::endl;
            std::exit(1);
        }
    }

    if (!infile || !outfile) {
        std::cout << "No input or output file given, exiting." << std::endl;
        std::exit(1);
    }
}

void write_file(MultiPlaneContext *ctx) {
    printf("Writing %s\n", outfile);

	// Writing binary is a lot faster but not readable.  But this is
	// just an example, so it's not really important.
    std::ofstream out(outfile, std::ios::binary);
    out.precision(6);

    for (size_t i = 0; i < source_z.size(); i++) {
		size_t s = thetas.size();
		out.write((char *)&s, sizeof(size_t));
        auto points = ctx->getSourcePositions(0);
        for (auto &p : points) {
			out.write((char *)&p, sizeof(Vector2D<float>));
        }
    }
}

void load_file() {
    printf("Parsing %s\n", infile);

    std::ifstream in(infile);

    // Read angular unit
    in >> angularUnit;

    // Read lenses
    int numlenses = 0;
    in >> numlenses;
    std::cout << "Lenses: " << numlenses << std::endl;
    for (int i = 0; i < numlenses; i++) {
        double z = 0;
        in >> z;
        lens_z.push_back(z);
    }

    lens_params.resize(numlenses);
    for (int i = 0; i < numlenses; i++) {
        int lenses = 0;
        in >> lenses;
        std::cout << "Sublenses: " << lenses << std::endl;
        for (int j = 0; j < lenses; j++) {
            float x, y, angle;
            in >> x >> y >> angle;
            PlummerParams pp = {Vector2D<float>(x, y), angle, 1};
            lens_params[i].push_back(pp);
        }
    }

    // Read source planes
    int numsources = 0;
    in >> numsources;
    std::cout << "Source planes: " << numsources << std::endl;
    for (int i = 0; i < numsources; i++) {
        double z = 0;
        in >> z;
        source_z.push_back(z);
    }

    // Read thetas
    int numthetas = 0;
    in >> numthetas;
    std::cout << "Loading thetas: " << numthetas << std::endl;
    for (int i = 0; i < numthetas; i++) {
        float x, y;
        in >> x >> y;
        thetas.push_back(Vector2D<float>(x, y));
    }

    // Load masses
    masses.resize(numlenses);
    std::cout << "Loading masses" << std::endl;
    for (int i = 0; i < numlenses; i++) {
        int nmass = 0;
        in >> nmass;
        for (int j = 0; j < nmass; j++) {
            double m = 0;
            in >> m;
            masses[i].push_back(m);
        }
    }

    in.close();
}

int main(int argc, char *argv[]) {
    int error = 0;

    parse_args(argc, argv);
    load_file();

    // Setup
    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    MultiPlaneContext ctx(angularUnit, cosm);
    error = ctx.init(lens_z, lens_params, source_z);
    if (error)
        return error;

	std::vector<std::vector<Vector2D<float>>> th;
	th.push_back(thetas);
	th.push_back(thetas);
	th.push_back(thetas);
    error = ctx.setThetas(th);
    if (error)
        return error;

    // Calculate betas
    error = ctx.calculatePositions(masses);
    if (error)
        return error;

    write_file(&ctx);

    return 0;
}
