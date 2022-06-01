#pragma once
#include <exception>
#include <iostream>
//CUSTOM EXCEPTIONS

struct arccos_out_of_range : public std::exception
{
    const char* what() const throw ()
    {
        return "Arccos's domain is [-1,1]";
    }
};

struct invalid_types_of_complex : arccos_out_of_range
{
    const char* what() const throw ()
    {
        return "Arccos's domain is [-1,1], therefore there are invalid complex roots";
    }
};

struct division_by_zero : public std::exception
{
    const char* what() const throw ()
    {
        return "You can't divide by zero";
    }
};

struct sqrt_of_negative_number : public std::exception
{
    const char* what() const throw ()
    {
        return "Square root's domain is [0, +inf)";
    }
};

struct unexcpected_number_of_roots : public std::exception
{
    const char* what() const throw ()
    {
        return "Received unexpected number of roots";
    }
};

struct nan_value : public std::exception
{
    const char* what() const throw ()
    {
        return "NaN value";
    }
};