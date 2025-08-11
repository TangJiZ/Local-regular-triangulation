#include <fstream>
#include <sstream>

#include "Point_d.h"

Point_d* read_points(std::string path, int num)
{
	Point_d* points = new Point_d[num];

	std::ifstream file(path);
	std::string line;
	std::istringstream iss;
	int i = 0;

	for (int i = 0; i < num; i++)
	{
		getline(file, line);
		iss.clear();
		iss.str(line);
		iss >> points[i].x;
		iss >> points[i].y;
	}

	file.close();
	return points;
}