#pragma once
#include <string>

struct Point_d
{
	double x, y;
	double w;
};

Point_d* read_points(std::string path, int num);
