#include <iostream>

#include "Local_Regular.h"
#include "regular_common.h"

int main()
{
	Point_d a{ 0, 0 };
	Point_d b{ 1, 0 };
	Point_d c{ 0, 1 };
	Point_d d{ 2, 1 };
	Point_d e{ 5, -3 };
	Point_d f{ 2, -2 };
	Point_d g{ 4, -1 };
	Point_d points[7] = { a, b, c, d, e, f, g };
	Local_Regular lr(points, 7);

	return 0;
}