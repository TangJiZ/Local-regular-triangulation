#include <iostream>


#include "Local_Regular.h"
#include "regular_common.h"

int main()
{
	//Point_d a{ 0, 0 };
	//Point_d b{ 1, 0 };
	//Point_d c{ 0, 1 };
	//Point_d d{ 2, 1 };
	//Point_d e{ 5, -3 };
	//Point_d f{ 3, -1.5 };
	//Point_d g{ 4, -1 };
	//Point_d points[7] = { a, b, c, d, e, f, g };
	//Local_Regular lr(points, 7);

	int num = 20;
	int prop = 100;
	std::string path = "D:/Microsoft Visual Studio/code/Triangulation/data";
	std::string pathIn = path + "/" + std::to_string(num) + "-" + std::to_string(prop) + "%/" + std::to_string(num) + "_in.txt";
	Point_d* points = read_points(pathIn, num);
	Local_Regular lr(points, num);

	std::cout << "valid point num: " << lr.valid_num << std::endl;
	std::cout << "invalid point num: " << lr.invalid_num << std::endl;
	std::cout << "regular triangle num: " << lr.rtri.size() << std::endl;
	//for (int i = 0; i < lr.rtri.size(); i++)
	//{
	//	std::cout << "t" << i << " : " << "Point1: " << lr.rtri[i].a.x << " " << lr.rtri[i].a.y << "  Point2:" << lr.rtri[i].b.x << " " << lr.rtri[i].b.y << "  Point3: " << lr.rtri[i].c.x << " " << lr.rtri[i].c.y << std::endl;
	//}

	std::string outPath = "C:/Users/ÌÆ¼ÌÖÝ/Desktop/data";
	lr.saveAndPlot("lr", outPath, true);
	std::string path1 = path + "/" + std::to_string(num) + "-" + std::to_string(prop) + "%/" + std::to_string(num) + "_out.txt";
	Local_Regular lr2;
	lr2.read_triangle(path1);
	lr2.saveAndPlot("lr2", outPath, true);
	CompareResult cr = compare_RT(&lr, &lr2);
	compareRTSaveAndPlot(cr, "compareTest", outPath, true);

	return 0;
}