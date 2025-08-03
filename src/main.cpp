#include<iostream>
#include<fstream>
#include<print>
#include"shh_code_generator.hpp"

constexpr size_t L = 4;

int main( int argc, char** argv )
{
    shh_code_generator< L > generator( "shhessian", false, true );
    return 0;
}
