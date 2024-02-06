#include <iostream>
#include <filesystem>
#include <string>

#include <cli.hpp>
#include <io.hpp>
#include "pex.cpp"   

int main(int argc, char** argv) {
    auto const opt = cli::parse_options(argc, argv);
    
    std::cout << "---{ welcome to floxer }---\n\n"
        << "- reference path: " << opt.reference_genome.c_str() << '\n'
        << "- query path: " << opt.queries.c_str() << '\n'
        << "- output path: " << opt.output_file.c_str() << "\n\n";

    std::cout << "---{ Reading input files..." << std::endl; 
    
    auto const input = io::read_inputs(opt.reference_genome, opt.queries);

    std::cout << "                              ...done }---\n\n"
        << "---{ Building FM-indices..." << std::endl;

    // TODO

    std::cout << "                              ...done }---\n\n"
        << "---{ Creating PEX search trees..." << std::endl;

    // TODO
    pex_tree p(12, 3, 2);
    p.debug_print();

    std::cout << "                              ...done }---\n\n"
        << "---{ Searching seeds..." << std::endl;

    // TODO

    std::cout << "                              ...done }---\n\n"
        << "---{ Verifiying hits..." << std::endl;

    // TODO

    std::cout << "                              ...done }---\n\n";

    return 0;
}
