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
