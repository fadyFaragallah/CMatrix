#include"CMatrix.h"
#include <fstream>


 void create_matrix(string&s) ;
 int check(string);
 void dop(string&s);
 void  detect_instruction(string &s);


int main(int argc,char*argv[])
{
	//int CMatrix::print=0;

if (argc == 1)
	{
		string h;
		while (getline(cin,h))
		{
			int equal=h.find("=");
			int pos_open_square=h.find("[");
			int pos_close_square=h.find("]");
			if(equal== -1 ) {continue;}
			if(pos_open_square!=-1 && pos_close_square== -1)
			{
				string x;
				while(getline(cin,x))
				{
					h=h+x;
					if(x.find(']')!=-1) {break;}
				}
			}
			if(h[h.length()-1]==';')  {CMatrix::print=0;}
			else CMatrix::print=1;
			detect_instruction(h);
		}
	}
	else if (argc == 2)
	{
		string s;
		ifstream file(argv[1]);
		while (getline(file, s))
		{
			int equal=s.find("=");
			int pos_open_square=s.find("[");
			int pos_close_square=s.find("]");
			if(equal== -1 ) {/*cin.ignore()*/; continue;}
			if(pos_open_square!=-1 && pos_close_square== -1)
			{
				string x;
				while(getline(file,x))
				{

					s+=x;
					if(x.find(']')!=-1) {break;}
				}
			}
			if(s[s.length()-1]==';')  {CMatrix::print=0;}
			else CMatrix::print=1;
			detect_instruction(s);
		}
	}

	return 0;
}
