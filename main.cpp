#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stack>
#include <cstdio>
#include <sstream>
#include "Matrix2017.h"
enum MI { MI_ZEROS, MI_ONES, MI_EYE, MI_RAND, MI_VALUE };
using namespace std;

Matrix* input_matrices= new Matrix[100];// array to hold matrices values
string* variables_names=new string[100];// array to hold matrices variable character
int in_cntr=0;

bool isChar(string s)
{
	//if (s.size() >1) return false;

	switch (s[0]) {
	case '+': return true;
	case '-': return true;
	case '*': return true;
	case '/': return true;
	case '(': return true;
	case ')': return true;
	case '^': return true;
	case '.': if(s[1]=='+') return true;
	else if(s[1]=='-') return true;
    else if(s[1]=='*') return true;
    else if(s[1]=='/') return true;
    else if(s[1]=='^') return true;
	default: return false;
	}
}

bool isLetter(string s)
{
    if (s[0]>='A' && s[0]<='Z')
        return true;
    else
        return false;
}



int main(int argc,char* argv[])
{
string matrix_char;
ifstream infile(argv[1]);
string contents[100];
int numberoflines=0;
int i=0;
while(getline(infile,contents[i]))// getting file contents
    {

        int pos =contents[i].find("[");
        int nbracketso =0;
        int nbracketsc =0;
        for(int j=0;j<contents[i].length();j++)
        {
            if(contents[i][j]=='[')
                nbracketso++;
            else if(contents[i][j]==']')
                nbracketsc++;
        }
        if(pos!=string::npos && nbracketso != nbracketsc)
        {
             getline(infile,contents[i+1]);
             contents[i] = contents[i] + ";" + contents[i+1];
        }


        if(contents[i]=="\n"||contents[i] == "\r" || contents[i] == "\t"||contents[i].find_first_not_of(' ') == std::string::npos)
        {
            i--;
            numberoflines--;
        }

        numberoflines++;

        i++;
    }


	int pos1;
for (i=0; i<numberoflines;i++)
{
    try
    {
int openbrackets =0, closedbrackets=0;
for (int h=0; h<contents[i].length();h++)
        {
        if(contents[i][h]=='(')
            openbrackets++;
        else if(contents[i][h]==')')
            closedbrackets++;
        }
    if (openbrackets!=closedbrackets)
        {
           throw("Syntax Error");
        }
	if (contents[i].find('[', 0) == string::npos && contents[i].find("ones", 0) == string::npos && contents[i].find("zeros", 0) == string::npos && contents[i].find("eye", 0) == string::npos && contents[i].find("rand", 0) == string::npos) //Only math expression
	{
	    string m;
	    string e;
		matrix_char = contents[i].substr(0, 1);
		if (contents[i].find('=',0)!=string::npos)
		{if (contents[i][contents[i].length()-1]!=';')
			m = contents[i].substr(contents[i].find('=', 0)+1, contents[i].length() - contents[i].find('=', 0)-1);
		else
			m = contents[i].substr(contents[i].find('=', 0) + 1, contents[i].length() - contents[i].find('=', 0)-2);}
        else{
            m = contents[i].substr(0,1);
        }
        e = evaluate(m);
        Matrix X(e);
        variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find('[', 0) != string::npos)
	{
	while (contents[i].find_first_of("+-*/^", contents[i].find('[', 0) + 1)!= string::npos) //find and replace all math expressions if square bracket is found
	{
		string d;
		int begin, end;
		pos1 = contents[i].find_first_of("+-*/^", contents[i].find('[', 0) + 1);
		begin = contents[i].find_last_of(" [];", pos1);
		end = contents[i].find_first_of(" [];", pos1);
		d = contents[i].substr(begin+1, end-begin-1);
		string f = evaluate(d);
		if(f.find("-",0)!=string::npos)
            f.replace(f.find("-",0),1,"$");
		contents[i].replace(begin+1, end - begin-1, f);
	}
	while(contents[i].find('$',0)!=string::npos)
        contents[i].replace(contents[i].find('$',0),1,"-");
	matrix_char = contents[i].substr(0,1);
	for (int k=1;k<contents[i].length();k++)
    {
        if(isLetter(contents[i].substr(k,1)))
        {
            for (int j=0;j<100;j++)
            {
                if (contents[i].substr(k,1)==variables_names[j])
                {
                    int opened_bracket = contents[i].find('[',k);
                    int closed_bracket = contents[i].rfind(']',k);
                    string space = contents[i].substr(k+1,opened_bracket-k-1);
                    string space2 = contents[i].substr(closed_bracket+1,k-closed_bracket-1);
                    if(space.find_first_not_of(' ') == std::string::npos || space2.find_first_not_of(' ') == std::string::npos)
                    contents[i].replace(k,1,input_matrices[j].getAltString());
                    else
                     contents[i].replace(k,1,input_matrices[j].getAltStringNoB());
                    break;
                }

            }
        }
    }
	string x;
    x = contents[i].substr(contents[i].find('[',0),contents[i].rfind(']')-contents[i].find('[',0)+1);
    x = concatenate(x);
    Matrix X(x);
    input_matrices[in_cntr] = X;
    variables_names[in_cntr] = matrix_char;
    in_cntr++;
	}

	else if (contents[i].find("rand", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = stoi(contents[i].substr(contents[i].find('(', 0) + 1, 1));
		int nC = stoi(contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1));
		Matrix X(nR, nC, MI_RAND,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("eye", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = stoi(contents[i].substr(contents[i].find('(', 0) + 1, 1));
		int nC = stoi(contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1));
		Matrix X(nR, nC, MI_EYE,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("zeros", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = stoi(contents[i].substr(contents[i].find('(', 0) + 1, 1));
		int nC = stoi(contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1));
		Matrix X(nR, nC, MI_ZEROS,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("ones", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = stoi(contents[i].substr(contents[i].find('(', 0) + 1, 1));
		int nC = stoi(contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1));
		Matrix X(nR, nC, MI_ONES,0.0);

        variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}

	//Output all the results
if(contents[i][contents[i].length()-2]!=';') //might be length()-1 due to difference in compilers ('\n') is counted in some compilers
cout << variables_names[in_cntr-1] << " = " << endl << input_matrices[in_cntr-1] << endl;
    }
    catch(char const* s)
    {
        cout<< "Error: " << s << endl;
    }
	}


return 0;
}
