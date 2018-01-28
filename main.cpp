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
string evaluate(string s1);

int countcommas(string x){
        string s = x;
        int commas=0;
        for(int i=0; i<s.length();i++){
            if(s[i]==','){
                commas++;
            }

        }
        return commas;
}
void addcomma(string& z){
    for(int j=0; j<z.length()-1;j++){
        if(z[j]==']'){
            int open=z.find_first_of('[',j);
            string str=z.substr(j+1,open-j-1);
            if(str.find_first_not_of(' ') != std::string::npos){

            }
            else{
                z.replace(j,1,"],");
            }
        }
    }
}
string concatenate (string s){

        string z =s;
        addcomma(z);
        Matrix matarr[2];
        int commapos,closedbracket,openedbracket;
        string tempconstruct;
        while (z.find(',') != string::npos)
        {
        commapos = z.find(',', 0);
        openedbracket=z.find_last_of('[',commapos);
        closedbracket=z.find_last_of(']',commapos);
        tempconstruct=z.substr(openedbracket, closedbracket - openedbracket+1);
        z.replace(openedbracket, closedbracket - openedbracket+1," ");
        matarr[0]=Matrix(tempconstruct);
        commapos = z.find(',', 0);
        openedbracket=z.find_first_of('[',commapos);
        closedbracket=z.find_first_of(']',commapos);
        tempconstruct=z.substr(openedbracket, closedbracket - openedbracket+1);
        z.replace(openedbracket, closedbracket - openedbracket+1," ");
        matarr[1]=Matrix(tempconstruct);
        matarr[0].addColumn(matarr[1]);
        string y = matarr[0].getAltString();
        z.replace(commapos,1,y);
        }
        return z;
        }


string replaceSin(string s)
{
	if (s.find("sin")!=string::npos)
	{
	    string mat_name ="ST";
		int begin = s.find('(', s.find("sin"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfSin = x1.length()+5;
		int start = (s.find("sin"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::sin(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfSin+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfSin+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	return NULL;
}
string replaceCos(string s)
{
	if (s.find("cos") != string::npos)
	{
	    string mat_name ="CT";
		int begin = s.find('(', s.find("cos"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfCos = x1.length()+5;
		int start = (s.find("cos"));
		Matrix Y(evaluate(x1));
         Matrix x;
		x=Matrix::cos(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfCos+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfCos+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	 return NULL;
}

string replaceTan(string s)
{
	if (s.find("tan") != string::npos)
	{
	    string mat_name ="TT";
		int begin = s.find('(', s.find("tan"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfTan = x1.length()+5;
		int start = (s.find("tan"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::tan(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfTan+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfTan+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	 return NULL;
}
string replaceLog(string s)
{
	if (s.find("log10") != string::npos)
{
        string mat_name ="LT";
		int begin = s.find('(', s.find("log10"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfLog10 = x1.length()+7;
		int start = (s.find("log10"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::log10(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfLog10+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfLog10+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
}
	return NULL;
}


string replaceLn(string s)
{
	if (s.find("log") != string::npos)
    {
	    string mat_name ="LLT";
		int begin = s.find('(', s.find("log"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfLn = x1.length()+5;
		int start = (s.find("tan"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::log(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfLn+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfLn+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	return NULL;
}


string replaceExp(string s)
{
	if (s.find("exp") != string::npos)
	{
	    string mat_name ="ET";
		int begin = s.find('(', s.find("exp"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfExp = x1.length()+5;
		int start = (s.find("exp"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::exp(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfExp+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfExp+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	return NULL;
}

string replaceSqrt(string s)
{
	if (s.find("sqrt") != string::npos)
	{
	    string mat_name ="SST";
		int begin = s.find('(', s.find("sqrt"));
		int end = s.find(')', begin+1);
		string x1 = s.substr(begin+1 , end-begin-1);
		int lengthOfSqrt = x1.length()+6;
		int start = (s.find("sqrt"));
		Matrix Y(evaluate(x1));
		Matrix x;
		x=Matrix::sqrt(Y);
//		char buffer[100];
//		sprintf(buffer, "%f", x);
//		string xstr = buffer;
        if(x.getn()==1)
		s.replace(start, lengthOfSqrt+1, x.getAltString());
		else
        {
            variables_names[in_cntr] = mat_name;
            mat_name+="X";
            input_matrices[in_cntr] = x;
            s.replace(start, lengthOfSqrt+1, variables_names[in_cntr]);
            in_cntr++;
        }
        return s;
	}
	return NULL;
}


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


Matrix solve(Matrix op1,Matrix op2, string ch)
{
    Matrix X;
	switch (ch[0])
	{
	case '+':
	    X = op1+op2;
	     return X;
	case '-': X = op1-op2;
	     return X;
	case '*': X = op1*op2;
	     return X;
	case '/': if(op2.getn()!=1) X = op1.div(op2);
              else X = op1/op2.getdDeterminant();
         return X;
	case '^': if (op2.getn()==1 && op1.getn()==1)
        X = pow(op1.getDeterminant(),op2.getDeterminant());
        else if(op2.getn()==1)
	    X = Matrix::pow(op1,op2.getDeterminant());
	    return X;
    case '.': if(ch[1]=='+') X = op1+op2.getDeterminant();
	else if(ch[1]=='-') X = op1-op2.getDeterminant();
    else if(ch[1]=='*') X = op1*op2.getDeterminant();
    else if(ch[1]=='/')if(op2.getn()==1) X = op1/op2.getDeterminant();
                        else if(op1.getn()==1) X = op1.rdivide(op1.getDeterminant(),op2);
                        else X = op1.rdivide(op1,op2);
    else if(ch[1]=='^') if (op2.getn()==1)
        X = Matrix::rpow(op1,op2.getDeterminant());
    return X;

	}
}

int preference(string ch)
{
	switch (ch[0])
	{
	case '+': return 1;
	case '-': return 1;
	case '*': return 2;
	case '/': return 2;
	case '^': return 3;
	case '(': return -1;
	case '.': if(ch[1]=='+') return 1;
	else if(ch[1]=='-') return 1;
    else if(ch[1]=='*') return 2;
    else if(ch[1]=='/') return 2;
    else if(ch[1]=='^') return 3;
	}
}

string evaluate(string s1)
{
	stack<string> operators;
	stack<Matrix> operands;
	string result;
	string s;
for (int i=0; i<s1.length();i++)
    if (s1[i]==' ')
    s1.erase(i,1);
	if(s1[0]=='-')
        s1.insert(0,"0");
    for (int i=0; i<s1.length();i++)
    if(s1[i]==' '&& s1[i+1] == '-')
        s1.insert(i+1,"0");
    else if(s1[i]=='('&& s1[i+1]=='-')
        s1.insert(i+1,"0");

	int pos = 0;
        while (pos != string::npos)
	{
		pos = s1.find_first_of("()+-/*^", pos);
		if (pos == -1) break;
		if (s1[pos - 1] == '.')
		{
			s1.insert(pos-1, " ");
			s1.insert(pos + 2, " ");
		}
		else
		{
            s1.insert(pos, " ");
            s1.insert(pos + 2, " ");
		}
		pos = s1.find_first_of("()+-/*^", pos) + 1;
	}

	if (s1.find("sin") != string::npos )
    {

			int pos = 0;
			while (pos != string::npos)
			{
				pos = s1.find("sin");
				if (pos == -1) break;
				pos += 3;
				s1 = replaceSin(s1);
			}
		}

    if (s1.find("cos") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("cos");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceCos(s1);
		}
	}

	if (s1.find("tan") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("tan");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceTan(s1);
		}
	}

	if (s1.find("log10") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("log10");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceLog(s1);
		}
	}


	if (s1.find("log") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("log");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceLn(s1);
		}
	}

	if (s1.find("exp") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("exp");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceExp(s1);
		}
	}

    if (s1.find("sqrt") != string::npos)
	{

		int pos = 0;
		while (pos != string::npos)
		{
			pos = s1.find("sqrt");
			if (pos == -1) break;
			pos += 3;
			s1 = replaceSqrt(s1);
		}
	}

	if (s1.find("Math error") != string::npos)
		return "Math error";



	istringstream is(s1);
	while (is >> s)
	{
		if (isChar(s))
		{
			if (s[0] == '(') operators.push("(");
			else if (s[0] == ')')
			{
				while (operators.top() != "(")
				{
					string ch = operators.top(); operators.pop();
					Matrix op2 = operands.top(); operands.pop();
					Matrix op1 = operands.top(); operands.pop();
					Matrix H;
                    H=solve(op1, op2, ch);
					operands.push(H);
				}
				operators.pop();
			}
			else
			{
				while (!operators.empty() && preference(s) <= preference(operators.top()))
				{

					string ch = operators.top(); operators.pop();
					Matrix op2 = operands.top(); operands.pop();
					Matrix op1 = operands.top(); operands.pop();

					operands.push(solve(op1, op2, ch));

				}

				operators.push(s);
			}

		}
		else if(isLetter(s))
        {
            for (int i=0;i<100;i++)
            {
                if(variables_names[i]==s)
                {
                    Matrix H(input_matrices[i]);
                    operands.push(H);
                    break;
                }
            }
        }
		else
		{
			Matrix op(s);
			operands.push(op);
		}

	}

	while (!operators.empty())
	{
		string ch = operators.top(); operators.pop();
		Matrix op2 = operands.top(); operands.pop();
		Matrix op1 = operands.top(); operands.pop();

		operands.push(solve(op1, op2, ch));
	}

//	char buffer[1000];
//	sprintf(buffer, "%f", operands.top().getDeterminate());
if(operands.top().getn()!=1)
	result = operands.top().getAltString();
	else
        result = operands.top().getAltStringNoB();

	return result;
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
		int nR = atof((contents[i].substr(contents[i].find('(', 0) + 1, 1).c_str()));
		int nC = atof((contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1)).c_str());
		Matrix X(nR, nC, MI_RAND,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("eye", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = atof((contents[i].substr(contents[i].find('(', 0) + 1, 1).c_str()));
		int nC = atof((contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1)).c_str());
		Matrix X(nR, nC, MI_EYE,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("zeros", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = atof((contents[i].substr(contents[i].find('(', 0) + 1, 1).c_str()));
		int nC = atof((contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1)).c_str());
		Matrix X(nR, nC, MI_ZEROS,0.0);
		variables_names[in_cntr] = matrix_char;
        input_matrices[in_cntr] = X;
        in_cntr++;
	}
	else if (contents[i].find("ones", 0) != string::npos)
	{
	    matrix_char = contents[i].substr(0,1);
		int nR = atof((contents[i].substr(contents[i].find('(', 0) + 1, 1).c_str()));
		int nC = atof((contents[i].substr(contents[i].find(',', 0) + 1, contents[i].find(')',0)-contents[i].find(',', 0) + 1)).c_str());
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
