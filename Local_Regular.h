#pragma once
#include <vector>

#include "Point_d.h"

/// @brief ��PCell��չȨֵԲ�ڵĽڵ㣬����洢
struct inPCell
{
	Point_d* p;
	inPCell* next;
};

/// @brief powerͼ�б߽��洢�ṹ������Ԫ�ṹ
struct PCell
{
	Point_d q;				// �߽罻�㣬power���ԲԲ�ģ�p1,p2�д���infinate��ʱ��������������
	Point_d* p1, * p2;		// �γɴ˽�������ھӽڵ�
	PCell* next;				// power��Ԫ�߽���ʱ�������һ����
	inPCell* ipc = nullptr;
	double r;				// ����չ���Ȩֵ��İ뾶
};

/// @brief Բ�������Χ��
struct Box
{
	int fy_min, fy_max, fx_min, fx_max;
};

/// @brief �ڵ���regular triangulaton�е��������ӽڵ�����
struct Nei
{
	Point_d* p;
	Nei* next;
};

/// @brief regular triangle�洢
struct RTriangle
{
	Point_d a, b, c;

	/// @brief ����==�����
	/// @param other ���Ƚ϶���
	/// @return �Ƿ���ͬ
	bool operator==(const RTriangle& other) const;
};

/// @brief RTriangleHash�ṹ�壬���ڰ���RTriangle��unorder_set�ȼ������ʹ���ʱ��hashģ������
struct RTriangleHash
{
	size_t operator()(const RTriangle& rt) const
	{
		double all = 0;
		all += rt.a.x + rt.a.y;
		all += rt.b.x + rt.b.y;
		all += rt.c.x + rt.c.y;
		/// �ض�double��С�����15λ��double�����ᵼ��hash�ļ��㲻ͬ���Ӷ�ʹ�ñȽ�ʱ��������
		all = floor(all * 1e15) / 1e15;
		return std::hash<double>{}(all);
	}
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

	std::vector<Nei*> nei;		// �洢ÿ���ڵ������ɺ�������ھӽڵ�

public:

	std::vector<RTriangle> rtri;		// �洢���������μ���
	int valid_num = 0;					// ��Ч�ڵ�����
	int invalid_num = 0;				// ��Ч�ڵ�����


/// @brief function
private:
	/// @brief ��ȡ������㼯xy�����Сֵ��w���ֵ
	void getMinMax_ptr();

	/// @brief ������㼯���򻮷�
	void partition();

	/// @brief ����չ�����p���ڽ��ڵ㣬��֤�ڵ���������2������������Ϊ��Χ8��ͬʱ��ȡ�����ľ�������Χ
	/// @param p �����ҽڵ�
	/// @param vp ��������p���ڽ��ڵ�
	/// @param box ���������ľ��������Χ�� 
	void near_find(Point_d* p, std::vector<int>& vp, Box& box);

	/// @brief power diagram��ʼ���죬Ҫ�󽻵���ʱ�������У�ֻ����ǰ�����ڵ�
	/// @param nf �ڽ��ڵ�
	/// @param c ���Ľڵ�
	/// @param first power diagram��������ͷ
	/// @param last power diagram��������β
	void power_diagram_build_initial(std::vector<int>& nf, Point_d c, PCell*& first, PCell*& last);

	/// @brief ��ʾpower diagram��
	/// @param iter power diagram����ʼָ��
	void show_pcell(PCell* iter);

	/// @brief ��power diagramͼ�в����½ڵ㣬������ͼ�ṹ���ж����ĵ��Ƿ�Ϊ���࣬�ɼ�������
	/// @param p ����ڵ�
	/// @param c ���Ľڵ�
	/// @param first power diagram��������ͷ��������ɺ����
	/// @param last power diagram��������β��������ɺ����
	/// @return ���ĵ�Ϊ����ڵ��޺�������false�����ĵ�Ƕ����������true 
	bool power_diagram_insert(Point_d* p, Point_d c, PCell* &first, PCell* &last);

	/// @brief ����������PCell�Ƿ���ֱ���жϵ�j����ж�
	/// @param cell ����������PCell
	/// @param A ֱ��ͳһ���̲���
	/// @param B 
	/// @param C 
	/// @param j �жϵ�j���ж�cell�Ƿ��j��ֱ�ߵ�ͬһ��
	/// @param c ���ĵ�c
	/// @return ��ͬһ��true������ͬһ��false
	bool infinate_same_side_judge(PCell* cell, double A, double B, double C, Point_d j, Point_d c);

	/// @brief ��ʾ������cell��չ���Բ���ڼнǷ�Χ�еĽڵ�
	/// @param iter power diagram����
	void show_inPCell(PCell* iter);

/// @brief ����power diagram���н�����չȨֵԲ�İ�����Χ����Ϊ�����γɰ��������Ͳ��������������
/// ------�����������------
	
	/// @brief ����pcell�İ�Χ��
	/// @param pcell ������power diagram����
	/// @param box ���ΰ�Χ��
	void pcell_box(PCell* pcell, Box& box);

	/// @brief �жϵ�p�Ƿ���cell��չ���Բ�ڣ�ͬʱ�����ĵ�c���������ɵ��γɵļн���
	/// @param p �жϵ�
	/// @param cell power diagram����
	/// @param c ���ĵ�
	/// @return ������չԲ�ͼн���true������false
	bool in_circle_and_angle(Point_d p, PCell* cell, Point_d c);

	/// @brief ��ȡ��ͬʱcell��չ���Բ�ͼн��еĽڵ㣬�洢��cell�е�ipc�У�����
	/// @param cell power diagram����
	/// @param exbox ���ĵ��ʼ����ʱ�Ѳ�ѯ�İ�Χ�з�Χ
	/// @param c ���ĵ�
	void findPoints_in_PCell(PCell* cell, Box exbox, Point_d c);

/// ------���������------
	
	/// @brief x�����ϣ�ֱ�ߺͳ������ཻ���֣�ÿ��xid��Ӧ��yid��Χ
	struct Bound
	{
		int x_min, x_max;	// ͨ��ֻʹ��x_min����ֱ�ߴ�ֱ��x��ʱ������Bound�߿ɱ�ʾ���з�Χ
		int y_min, y_max;
	};

	/// @brief ��������չֱ�ߺͽڵ��ⳤ���ο���Χ�ɵ�����Χ������Ϊ������Ӧ�ķ�������
	/// @param cell ������power diagram���㣬����infinate��
	/// @param c ���ĵ㣬ֱ��Ϊ���ĵ�ͷ�infinate���γ�
	/// @param vb Χ�ɵ�����Χ��x����������洢
	void pcell_bound_infinate(PCell* cell, Point_d c, std::vector<Bound>& vb);

	/// @brief ��ȡ������չ���ڵĽڵ㣬�洢��cell�е�ipc�У�����
	/// @param cell power diagram���㣬����infinate�ڵ�
	/// @param exbox ���ĵ��ʼ����ʱ�Ѳ�ѯ�İ�Χ�з�Χ
	/// @param c ���ĵ�
	void findPoints_in_PCell_infinate(PCell* cell, Box exbox, Point_d c);

/// ------******------

	/// @brief ����ڵ�󣬸�����Ķ������Ӧ��inPCell��Բ�ڽڵ㣩��������Ľ���Բ�ڽڵ�Ҫô���½���Բ�⣬Ҫô���½���Բ��
	/// @param start �����������ͷ
	/// @param end �����������β
	/// @param c ���Ľڵ�
	/// @param new1 �²�������1
	/// @param new2 �²�������2
	void update_inPCell(PCell* start, PCell* end, Point_d c, PCell* new1, PCell* new2);

	/// @brief �������������cell����ֱ�ߣ�������һ��������ֱ���ڲ�����Ľڵ�
	/// @param cell	power diagram���㣬����infinate
	/// @param c ���Ľڵ�
	/// @param A ����ֱ��ͨ�ò���
	/// @param B 
	/// @param C 
	/// @param p ������ֱ���ڲ�����һ��
	void compute_ext_line_and_point(PCell* cell, Point_d c, double& A, double& B, double& C, Point_d& p);

	/// @brief ��ȡ���������μ��ϣ�ÿ���ڵ�ֻ�ͱ�ű������Ľڵ��γ�������
	void get_triangle();

public:
	/// @brief ���캯��ͬʱ�ڲ�ִ���ʷ�
	/// @param p ������㼯
	/// @param num �㼯�нڵ�����
	Local_Regular(Point_d* p, int num);

	/// @brief Ĭ�Ϲ��캯��
	Local_Regular() {}

	/// @brief ��ȡ������
	/// @param path ��ȡ�ļ�·��
	void read_triangle(std::string path);

	/// @brief ���棬����������
	/// @param name �ļ�����������׺�������ļ��Զ�Ϊ�ļ���+txt��ͼ���ļ��Զ�Ϊ�ļ���+jpg
	/// @param path �ļ�����·��
	/// @param plot ѡ���Ƿ����ͼ��
	void saveAndPlot(std::string name, std::string path, bool plot);
};

