#ifndef Matrix_m
#define Matrix_m

#include<iostream>
#include<string>

using namespace std;

class CMatrix
{
private:
	int nrows, ncols;
	double** pp_rows;
public:
	static int  print;
	CMatrix(int r, int c);
	CMatrix();
	void destroy_matrix();
	void copy_matrix(const CMatrix &m);
	CMatrix(const CMatrix &m);
	~CMatrix();
	CMatrix& operator=(const CMatrix& m);
	void set_element(int a, int b, double value);
	CMatrix  operator+(CMatrix &m);
	CMatrix operator-(CMatrix &m);
	CMatrix  operator*(CMatrix &m);
	CMatrix  transpose();
	CMatrix get_cofactor(int r,int c);
	double get_determinant();
	CMatrix operator/(CMatrix &m);
	CMatrix operator/(double d);
	CMatrix operator*(double d);
	CMatrix operator+(double d);
	CMatrix operator-(double d);
	CMatrix operator-();
	CMatrix num_sub_mat(double d);
	CMatrix num_div_mat(double d);
	void print_matrix(string name);

	bool check_singularity();
	CMatrix(int r,int c,string type);
	CMatrix inv();
	double get_determinant_LU();

	friend void create_matrix(string&s);
	friend int check(string);
	friend void dop(string&s);
	friend void  detect_instruction(string&s);
};

#endif