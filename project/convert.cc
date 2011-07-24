#include "config/config.h"
#include "include/time.h"
#include "modules/mesh/mesh.h"
#include "main.h"

#include <getopt.h>

int main (int argc, char * argv[]) {
    char input[500], output[500];
    input[0] = 0; output[0] = 0;
    bool update_c = false;
    float c = 1;

    while (1) {
	int option_index = 0;
	int ch = getopt(argc, argv, "i:o:c:h");

	if (ch == -1)
	    break;

	switch (ch) {
	    case 'o': strcpy(output, optarg); std::cout << output << std::endl; break;
	    case 'i': strcpy(input, optarg); break;
	    case 'c': update_c = true; c = atof(optarg); break;
	    case 'h': std::cout << "Usage: ./convert [-i <input>] [-o <output>] [-c <value>]" << std::endl;
		      exit(0);
	    case '?':
	    default:
		      abort();
	}
    }
    if (!strlen(input)) {
	std::cout << "No input file specified" << std::endl;
	return 1;
    }

    if (!strlen(output)) {
	std::cout << "No output file is specified" << std::endl;
	return 1;
	// char * inname = strtok(input, '.');
	// char * inext  = strtok(inname, '.');
	// strcpy(output, inname);
	// bool input_is_binary = false;
    }

    /*
     * Check extensions of the files:
     * crs - ascii
     * dat - binary
     */
    bool input_is_ascii  = (strstr(input, ".crs") ? true : false);
    bool output_is_ascii = (strstr(output, ".crs") ? true : false);

    if (input_is_ascii == output_is_ascii && update_c == false) {
	std::cout << "Both files are binary | ascii and no Laplace transformation" << std::endl;
	return 1;
    }

    std::cout << "Input (ascii = " << input_is_ascii << "): " << input << std::endl;
    std::cout << "Output (ascii = " << output_is_ascii << "): " << output << std::endl;
    SkylineMatrix A;
    A.load(input, input_is_ascii ? ASCII : BINARY, false);

    if (update_c) {
	const uint n = A.size();
	Vector x(n, 1), y(n);
	multiply(A, x, y);
	for (uint i = 0; i < n; i++)
	    A(i,i) = c - (y[i] - A(i,i));
    }

#if 0
    /* Add diagonal to matrix */
    for (uint i = 0; i < A.rows(); i++)
	A(i,i) += 100;
#endif

    dump(output, A, output_is_ascii ? ASCII : BINARY);

    return 0;
}
