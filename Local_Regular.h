#pragma once
#include <vector>

#include "Point_d.h"

/// @brief powerͼ�б߽��洢�ṹ������Ԫ�ṹ
struct PCell
{
	Point_d q;				// �߽罻�㣬power���ԲԲ�ģ�p1,p2�д���infinate��ʱ��������������
	double r;				// �뾶
	Point_d* p1, * p2;		// �γɴ˽�������ھӽڵ�
	PCell* next;				// power��Ԫ�߽���ʱ�������һ����
};

/// @brief �ֲ���������regular�����ʷ���
class Local_Regular
{
/// @brief value
private:
	/// @brief �����
	Point_d* infinate;

	/// @brief ������㼯�͵㼯�нڵ�����
	Point_d* pts;	// ������㼯
	int n;			// �ڵ�����

	/// @brief �㼯��xy�����Сֵ��w�����ֵ
	double xmin, xmax, ymin, ymax, wmax;

	/// @brief �㼯������ز���
	int xnum, ynum;				// x��y�����������
	double xcell, ycell;		// ������Ԫ��С
	int* bucket;				// �洢ÿ�������е�һ���ڵ�
	int* next;					// �洢ÿ���ڵ��ڷ����е���һ���ڵ㣬��ڵ�ͬ����ʱ��������

public:


/// @brief function
private:
	/// @brief ��ȡ������㼯xy�����Сֵ��w���ֵ
	void getMinMax_ptr();

	/// @brief ������㼯���򻮷�
	void partition();

	/// @brief ����չ�����p���ڽ��ڵ㣬��֤�ڵ���������2������������Ϊ��Χ8��
	/// @param p �����ҽڵ�
	/// @param vp ��������p���ڽ��ڵ�
	void near_find(Point_d* p, std::vector<int>& vp);

	/// @brief power diagram��ʼ���죬Ҫ�󽻵���ʱ�������У�ֻ����ǰ�����ڵ�
	/// @param nf �ڽ��ڵ�
	/// @param first power diagram��������ͷ
	/// @param last power diagram��������β
	void Power_diagram_build_initial(std::vector<int>& nf, PCell*& first, PCell*& last);

public:
	/// @brief ���캯��ͬʱ�ڲ�ִ���ʷ�
	/// @param p ������㼯
	/// @param num �㼯�нڵ�����
	Local_Regular(Point_d* p, int num);
};

