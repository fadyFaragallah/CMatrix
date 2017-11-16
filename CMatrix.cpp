#include"CMatrix.h"
#include<vector>
vector <CMatrix> v ;
vector <string> names;

int CMatrix::print;

CMatrix::CMatrix(int r, int c)     
	{
		nrows = r;
		ncols = c;
		pp_rows = new double*[r];
		for (int i = 0; i<r; i++) pp_rows[i] = new double[c];
		for (int i = 0; i<r; i++)
			for (int j = 0; j<c; j++)
				pp_rows[i][j] = 0;
	}

CMatrix::CMatrix()
	{
		nrows = 0;
		ncols = 0;
		pp_rows = NULL;
	}

void CMatrix::destroy_matrix()
	{
		if (pp_rows)
		{
			for (int i = 0; i<nrows; i++)
				delete[] pp_rows[i];
			delete[] pp_rows;
		}
		nrows = 0;
		ncols = 0;
		pp_rows = NULL;
	}

void CMatrix::copy_matrix(const CMatrix &m) 
	{
		this->destroy_matrix();
		this->nrows = m.nrows;
		this->ncols = m.ncols;
		if ((nrows*ncols) == 0) { pp_rows = NULL; return; }
		pp_rows = new double*[nrows];
		for (int i = 0; i<nrows; i++) pp_rows[i] = new double[ncols];
		for (int i = 0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				pp_rows[i][j] = m.pp_rows[i][j];
	}

CMatrix::CMatrix(const CMatrix &m)
	{
		nrows = 0; ncols = 0; pp_rows = NULL;
		copy_matrix(m);
	}

CMatrix::~CMatrix()
	{
		destroy_matrix();
	}

CMatrix& CMatrix::operator=(const CMatrix& m) 
	{
		copy_matrix(m);
		return *this;
	}

void CMatrix::set_element(int a, int b, double value)
	{
		pp_rows[a][b] = value;
	}    

CMatrix CMatrix::operator+(CMatrix &m) // the name of the resulting matrix should be set manually
	{
		CMatrix c(nrows, ncols);
		for (int i = 0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				c.set_element(i, j, pp_rows[i][j] + m.pp_rows[i][j]);
		return c;
	}

CMatrix CMatrix::operator-(CMatrix &m)
	{
		CMatrix c(nrows, ncols);
		for (int i = 0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				c.set_element(i, j, pp_rows[i][j] - m.pp_rows[i][j]);
		return c;
	}

CMatrix CMatrix::operator*(CMatrix &m)
	{
		CMatrix c(nrows, m.ncols);
		for (int i = 0; i<nrows; i++)
		{
			for (int j = 0; j<m.ncols; j++)
			{
				double el = 0;
				for (int k = 0; k<ncols; k++)
					el += pp_rows[i][k] * m.pp_rows[k][j];
				c.set_element(i, j, el);
			}
		}
		return c;
	}

CMatrix CMatrix::transpose() 
	{
		CMatrix m(this->ncols, this->nrows);
		for (int i = 0; i<this->nrows; i++)
		{
			for (int j = 0; j<this->ncols; j++)
				m.set_element(j, i, pp_rows[i][j]);
		}
		return m;
	}

CMatrix CMatrix::get_cofactor(int r,int c)
	{
		CMatrix m(nrows-1,ncols-1);
		for(int i=0;i<m.nrows;i++)
		{
			for(int j=0;j<m.ncols;j++)
			{
				int sR = (i<r)?i:i+1;
				int sC = (j<c)?j:j+1;
				m.pp_rows[i][j] = pp_rows[sR][sC];
			}
		}
		return m;
	}

double CMatrix::get_determinant()
	{
		if(nrows==1 && ncols==1) return pp_rows[0][0];
		double values=0; int m=1;
		for(int i=0;i<nrows;i++)
		{
			values+=m*pp_rows[0][i]*get_cofactor(0,i).get_determinant();
			m*=-1;
		}
		return values;
	}

CMatrix CMatrix::operator/(CMatrix &m)
	{
	//	if(m.nrows!=m.ncols || this->nrows!=this->ncols) {cout<<"non-square matrix"<<endl; return;}
		double a=m.get_determinant();
	//	if (a==0) {cout<<"the determinant of the matrix=0"<<endl; return;}
		CMatrix x(m.nrows,m.ncols); int sign=1;
		for(int i=0;i<x.nrows;i++)
		{
			for(int j=0;j<x.ncols;j++)
			{
				if(i%2 != j%2) sign=-1; else sign=1;
				x.set_element(i,j,sign*m.get_cofactor(i,j).get_determinant());
				//sign*=-1;
			}
		}
		x=x.transpose();
		return *this*x*(1.0/a);
	}

CMatrix CMatrix::operator/(double d)
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<r1.nrows; i++)
			for (int j = 0; j<r1.ncols; j++)
				r1.set_element(i, j, double(this->pp_rows[i][j]) / d);
		return r1;
	}

CMatrix CMatrix::operator*(double d)
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<r1.nrows; i++)
			for (int j = 0; j<r1.ncols; j++)
				r1.set_element(i, j, this->pp_rows[i][j] * d);
		return r1;
	}

CMatrix CMatrix::operator+(double d)
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<r1.nrows; i++)
			for (int j = 0; j<r1.ncols; j++)
				r1.set_element(i, j, this->pp_rows[i][j] + d);
		return r1;
	}

CMatrix CMatrix::operator-(double d)
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<r1.nrows; i++)
			for (int j = 0; j<r1.ncols; j++)
				r1.set_element(i, j, this->pp_rows[i][j] - d);
		return r1;
	}

CMatrix CMatrix::operator-()
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<r1.nrows; i++)
			for (int j = 0; j<r1.ncols; j++)
				r1.set_element(i, j, -1 * this->pp_rows[i][j]);
		return r1;
	}

CMatrix CMatrix::num_sub_mat(double d)
	{
		return -*this + d;
	}

CMatrix CMatrix::num_div_mat(double d)
	{
		CMatrix r1(this->nrows, this->ncols);
		for (int i = 0; i<this->nrows; i++)
			for (int j = 0; j<this->ncols; j++)
				r1.set_element(i, j, 1.0 / this->pp_rows[i][j]);
		return r1*d;
	}

void CMatrix::print_matrix(string name)
	{
		if(print==0) return;
		cout << name << " =" << endl;
		for (int i = 0; i<nrows; i++)
		{
			for (int j = 0; j<ncols; j++)
				cout << "\t" << pp_rows[i][j];
			cout << endl;
		}
	}

////////////////////////////for the loving memory of jordan and guass/////////////////

CMatrix::CMatrix(int r,int c,string type)
{
		nrows = r;
		ncols = c;
		pp_rows = new double*[r];
		for (int i = 0; i<r; i++) pp_rows[i] = new double[c];
		for (int i = 0; i<r; i++){
			for (int j = 0; j<c; j++){
				if(i==j) pp_rows[i][j]=1;
				else pp_rows[i][j] = 0;
			}
		}
}

CMatrix CMatrix::inv()
{
	CMatrix m=*this;
	CMatrix x(m.nrows,m.ncols,"unity");
	for(int i=0;i<m.nrows;i++)
	{
		for(int j=0;j<m.ncols;j++)
		{
			double a=m.pp_rows[i][i];
			m.pp_rows[i][j]/=a;
			x.pp_rows[i][j]/=a;
		}
		for(int k=0;k<m.nrows;k++)
		{
			if(k==i) continue;
			else 
			{
				for(int z=0;z<m.ncols;z++)
				{
					double b=m.pp_rows[k][i];
					m.pp_rows[k][z]+=-1*b*m.pp_rows[i][z];
					x.pp_rows[k][z]+=-1*b*x.pp_rows[i][z];
				}
			}
		}
	}
	for(int i=0;i<m.nrows;i++)
	{
		for(int j=0;j<m.ncols;j++)
			cout<<m.pp_rows[i][j]<<" \t";
		cout<<endl;
	}
	for(int i=0;i<x.nrows;i++)
	{
		for(int j=0;j<x.ncols;j++)
			cout<<x.pp_rows[i][j]<<" \t";
		cout<<endl;
	}
	return x;
}

bool CMatrix::check_singularity()
{
	int zero=0;
	for(int i=0;i<this->nrows-1;i++)
	{
		double a=pp_rows[i+1][0]/pp_rows[i][0];
		for(int j=1;j<this->ncols;j++)
		{
			double b=pp_rows[i+1][j]/pp_rows[i][j];
			if(a!=b) {break;}
			if(j==this->ncols-1){zero=1;}
		}
		if(zero) break;
	}
	if(zero) return true;
	else if(this->get_determinant()==0) return true;
	else return false;
}

//////////////////////////////////parsing////////////////////////////////////////
int check(string name)
{
	int s = 0; unsigned int i;
	for (i = 0; i < names.size(); i++)
	{
		if (names[i] == name)
		{
			s = 1;
			break;
		}
	}
	if (s == 1)
		return i;
	else
		return -1;
}




void  dop(string&s)
{
	// updated on 29_10_2017 // to support cascading equal signs in operation line and the "./" operation
	s += ";";
	int j = 0, k = 0, t = 0, no_variables_equated = 0;
	string my_operation[100];

	string *ptr_variables_equated; string q = "";

	while (s[j] != ';')  //counting no_variables_equated, variables before the last equal sign
	{
		if (s[j] == '=') { no_variables_equated++; }
		j++;
	}
	// creating an array of string with size = no_variables_equated +1
	//it will be ended by "00"

	ptr_variables_equated = new string[no_variables_equated + 1];

	j = 0; t = 0; q = "";

	//filling the array with the names of variables to be equated
	//putting all the variables before the last equal sign in the array pointed by ptr_variables_equated

	while ((t<no_variables_equated))
	{

		if ((s[j] >= 'a' && s[j] <= 'z') || (s[j] >= 'A' && s[j] <= 'Z') || (s[j] >= '0' && s[j] <= '9') || (s[j] == '_'))
		{
			q += s[j];
			if (!((s[j + 1] >= 'a' && s[j + 1] <= 'z') || (s[j + 1] >= 'A' && s[j + 1] <= 'Z') || (s[j + 1] >= '0' && s[j + 1] <= '9') || (s[j + 1] == '_')))
			{
				ptr_variables_equated[t] = q; q = "";


				t++;
			}
		}

		j++;
	}

	//putting "00" in the last element of the array of strings pointed by ptr_variables_equated
	ptr_variables_equated[t] = "00";


	// i created an array of strings named my_operation to hold all the elements after the last equal sign
	// filling this array
	while (s[j] != ';')
	{
		if ( (s[j] >= '0' && s[j] <= '9')  || ( (s[j]== '-') && (s[j + 1] >= '0' && s[j + 1] <= '9') )
			|| ((s[j] == '.') && (s[j + 1] >= '0' && s[j + 1] <= '9'))       )

		{
			q += s[j];
			if (!(s[j + 1] >= '0' && s[j + 1] <= '9')  )
			{

				my_operation[k] = q;

				q = ""; k++;

			}

		}

		else if ((s[j] >= 'a' && s[j] <= 'z') || (s[j] >= 'A' && s[j] <= 'Z') || (s[j] == '_') || (s[j] >= '0' && s[j] <= '9'))

		{
			
			q += s[j];
			if (!((s[j + 1] >= 'a' && s[j + 1] <= 'z') || (s[j + 1] >= 'A' && s[j] <= 'Z') || (s[j + 1] == '_') || (s[j + 1] >= '0' && s[j + 1] <= '9')))
			{

				my_operation[k] = q;
				q = ""; k++;

			}

		}
		else if ((s[j] == '+') || (s[j] == '-') || (s[j] == '*') || (s[j] == '/') || (s[j] == '\'') || (s[j] == '.'))
		{
			if (s[j] == '.')
			{
				int xyz = j + 1;
				while (1)
				{
					if (s[xyz] == '/' || s[xyz]=='*' || s[xyz] == '+' || s[xyz] == '-')
					{
						my_operation[k] = s[j]; my_operation[k] += s[xyz]; // cout<< my_operation[k];
						break;
					}
					xyz++;
				}
				j = xyz + 1;

			}
			else  my_operation[k] = s[j];

			k++;
		}


		j++;
	}
	// putting "00" in the last element of my_operation array
	my_operation[k] = "00";

	//here, mr fady farag should use the array named my_operation to do the operations on the matrices
	//then, put the result in all the elements that are in the array pointed by ptr_variables_equated
	//you should loop in all these elements till you find "00"


	int a1 = -1, b1 = -1, tran = 0, inv = 0, dot = 0, defined = 0, dot_operand = -1, matched = 1;
	CMatrix op_res;

	if (my_operation[1] == "'")
	{
		tran = 1;
	}
	else if (my_operation[1] == "./" || my_operation[1] == ".*" || my_operation[1] == ".+" || my_operation[1] == ".-") { dot = 1; }

	for (unsigned int counter = 0; counter<v.size(); counter++)
	{
		if (names[counter] == my_operation[0] && (a1 == -1) && (!dot)) { a1 = check(names[counter]); /*cout<<a1<<endl;*/ }
		if ((!tran) && (!dot))
		{
			if (names[counter] == my_operation[2] && (b1 == -1)) { b1 = check(names[counter]);/*cout<<b1<<endl;*/ }
		}
		if (dot) {
			if (((int(my_operation[0][0]) >= 48) && (int(my_operation[0][0]) <= 57)) && (names[counter] == my_operation[2]))
			{
				a1 = atof(my_operation[0].c_str()); dot_operand = 0; b1 = check(names[counter]);
			}

			else  if ((int(my_operation[2][0]) >= 48 && int(my_operation[2][0]) <= 57) && (names[counter] == my_operation[0]))
			{
				b1 = atof(my_operation[2].c_str()); dot_operand = 2; a1 = check(names[counter]);
			}
		}
	}
	// cout<<a1<<" "<<b1<<" "<<tran<<" "<<inv<<" "<<endl;

	if ((((tran == 1)) && (a1 != -1)) || ((dot == 1) && (b1 != -1 || a1 != -1)) || ((inv == 0) && (tran == 0) && (a1 != -1) && (b1 != -1)))
	{
		defined = 1;
		if (my_operation[1] == "+")
		{
			if ((v[a1].nrows == v[b1].nrows) && (v[a1].ncols == v[b1].ncols))
				op_res = v[a1] + v[b1];
			else { cout << "unmatched dimensions of the two matrices" << endl; matched = 0; }
		}

		else if (my_operation[1] == "-")
		{
			if ((v[a1].nrows == v[b1].nrows) && (v[a1].ncols == v[b1].ncols))
				op_res = v[a1] - v[b1];
			else { cout << "unmatched dimensions of the two matrices" << endl; matched = 0; }
		}

		if (my_operation[1] == "*")
		{
			if ((v[a1].ncols == v[b1].nrows))
				op_res = v[a1] * v[b1];
			else { cout << "unmatched dimensions of the two matrices" << endl; matched = 0; }
		}

		if (my_operation[1] == "/")
		{
			//if(v[b1].get_determinant()==0) cout<<"the second matrix has determinant of zero"<<endl;
			if(v[b1].check_singularity()) {cout<<"the matrix "<<names[b1]<<" is singular."<<endl; matched=0;}
			else if ((v[a1].nrows == v[a1].ncols) && (v[b1].nrows == v[b1].ncols) && (v[a1].nrows == v[b1].ncols))
			op_res = v[a1] / v[b1];
			//op_res = v[a1]* v[b1].inv();
			else { cout << "unmatched dimensions of the two matrices" << endl; matched = 0; }
		}

		if (my_operation[1] == "'") { op_res = v[a1].transpose(); }

		if (my_operation[1] == "./")
		{
			if (dot_operand == 0) { op_res = v[b1].num_div_mat(a1); }
			else if (dot_operand == 2) { op_res = v[a1] / b1; }
		}

		if (my_operation[1] == ".*")
		{
			if (dot_operand == 0) { op_res = v[b1] * a1; }
			else if (dot_operand == 2) { op_res = v[a1] * b1; }
		}

		if (my_operation[1] == ".-")
		{
			if (dot_operand == 0) { op_res = v[b1].num_sub_mat(a1); }
			else if (dot_operand == 2) { op_res = v[a1] - b1; }
		}

		if (my_operation[1] == ".+")
		{
			if (dot_operand == 0) { op_res = v[b1] + a1; }
			else if (dot_operand == 2) { op_res = v[a1] + b1; }
		}


		//	 op_res.print_matrix();
	}
	else cout << "undefined matrices" << endl;
	if (defined&&matched)
	{
		for (int i = 0; ptr_variables_equated[i] != "00"; i++)
		{
			int place = check(ptr_variables_equated[i]);
			if (place != -1) { v[place] = op_res; }
			else {
				CMatrix temp;
				names.push_back(ptr_variables_equated[i]);
				temp = op_res;
				v.push_back(temp);
			}

		}
		op_res.print_matrix(ptr_variables_equated[0]);
	}
}



void create_matrix(string &s)
{
	s += ";";
	string q = "";
	int n_rows = 1, n_cols = 0, k = 0, j = 0, t = 0, no_variables_to_be_created = 0;
	double **p; string *ptr_variables_to_be_created;

	//i should create an array that contains the names of the variables to be created, the array's ptr is ptr_variables_to_be_created

	while (s[j] != '[')  //counting no_variables_to_be_created
	{
		if (s[j] == '=') { no_variables_to_be_created++; }
		j++;
	}

	// creating an array of string with size = no_variables_to_be_created +1 , bec it will be ended by "00"

	ptr_variables_to_be_created = new string[no_variables_to_be_created + 1];

	j = 0;


	//filling the array with the names of variables to be created
	while (s[j] != '[')
	{
		if ((s[j] >= 'a' && s[j] <= 'z') || (s[j] >= 'A' && s[j] <= 'Z') || (s[j] >= '0' && s[j] <= '9') || (s[j] == '_'))
		{
			q += s[j];
			if (!((s[j + 1] >= 'a' && s[j + 1] <= 'z') || (s[j + 1] >= 'A' && s[j + 1] <= 'Z') || (s[j + 1] >= '0' && s[j + 1] <= '9') || (s[j + 1] == '_')))
			{
				ptr_variables_to_be_created[t] = q; q = ""; t++;
			}
		}

		j++;
	}
	ptr_variables_to_be_created[t] = "00";
	t = 0;

	string name = q; q = "";

	while (s[j] != ']')
	{
		if (s[j] == ';')
		{
			if (s[j + 1] == ']')
				n_rows--;
			n_rows++;
		}
		j++;
	}

	j = 0;
	while (s[j] != ';')
	{
		if (((s[j] >= '0' && s[j] <= '9') || (s[j] == '-') || (s[j] == '.')))
		{
			if (!((s[j + 1] >= '0' && s[j + 1] <= '9') || (s[j + 1] == '-') || (s[j + 1] == '.')))
			{

				n_cols++;
			}
		}
		j++;

	}
	// cout<<n_rows<<" "<<n_cols<<"  "<<name<<endl;
	p = new double*[n_rows];

	for (int i = 0; i<n_rows; i++)
		p[i] = new double[n_cols];

	k = 0; t = 0; j = 0; q = "";
	while (s[j] != ']')
	{

		if (((s[j] >= '0' && s[j] <= '9') || (s[j] == '-') || (s[j] == '.')))
		{
			q += s[j];
			if (!((s[j + 1] >= '0' && s[j + 1] <= '9') || (s[j + 1] == '-') || (s[j + 1] == '.')))
			{
				p[k][t] = atof(q.c_str());
				//cout<<p[k][t]<<endl;
				q = "";
				t++;
			}

		}
		if (s[j] == ';') { k++; t = 0; }

		j++;

	}

	unsigned int counter = 0;


	while (ptr_variables_to_be_created[counter] != "00")
	{
		CMatrix c;
		if (check(ptr_variables_to_be_created[counter]) == -1)
		{
			names.push_back(ptr_variables_to_be_created[counter]);
			c.nrows = n_rows;
			c.ncols = n_cols;
			c.pp_rows = new double *[c.nrows];
			for (int qw = 0; qw < c.nrows; qw++)
			{
				c.pp_rows[qw] = new double[c.ncols];
			}
			for (int a = 0; a < c.nrows; a++)
			{
				for (int s = 0; s < c.ncols; s++)
				{
					c.pp_rows[a][s] = p[a][s];
				}
			}
			v.push_back(c);
			if (counter == 0)
				c.print_matrix(ptr_variables_to_be_created[0]);
		}
		else
		{
			for (int del = 0; del < v[check(ptr_variables_to_be_created[counter])].nrows; del++)
			{
				delete[]v[check(ptr_variables_to_be_created[counter])].pp_rows[del];
			}
			delete[]v[check(ptr_variables_to_be_created[counter])].pp_rows;
			///////deletion of old existed is complete.

			///////creation of new matrix in the same place of the old one.

			v[check(ptr_variables_to_be_created[counter])].nrows = n_rows;
			v[check(ptr_variables_to_be_created[counter])].ncols = n_cols;
			v[check(ptr_variables_to_be_created[counter])].pp_rows = new double *[v[check(ptr_variables_to_be_created[counter])].nrows];

			for (int qw = 0; qw < v[check(ptr_variables_to_be_created[counter])].nrows; qw++)
			{
				v[check(ptr_variables_to_be_created[counter])].pp_rows[qw] = new double[v[check(ptr_variables_to_be_created[counter])].ncols];
			}
			for (int a = 0; a < v[check(ptr_variables_to_be_created[counter])].nrows; a++)
			{
				for (int s = 0; s < v[check(ptr_variables_to_be_created[counter])].ncols; s++)
				{
					v[check(ptr_variables_to_be_created[counter])].pp_rows[a][s] = p[a][s];
				}
			}
			if (counter == 0)
				v[check(ptr_variables_to_be_created[counter])].print_matrix(ptr_variables_to_be_created[0]);
		}
		counter++;
	}
}







//now we must create another copies of the matrix c and put them in the vector, using ptr_variables_to_be_created,
// these matrices have the same elements as matrix c only differing in sring name 
// u shoud start creating matrices starting from ptr_variables_to_be_created[1] till you find "00" ...

void detect_instruction(string&s)
{
	if (s.find("[", 0) != -1) create_matrix(s);
	else dop(s);
}
