#include "Matrix2017.h"
#include "iostream"
#include "string.h"
using namespace std;
int main()
{
	char** commands;
	ifstream infile("example.txt");
	Matrix* inputs= new Matrix[10];// array to hold matrices values
	string* variables=new string[60];// array to hold matrices variable charachter
	int in_cntr=0;// common index between 2 arrays
	Matrix m;
	string contents[20];//array contains file contents
	string line;
	int numberoflines=0;//index contains number of lines in file
	int i=0;//index for contents array
	int max_number_of_characters_in_a_command=0;
	int number_of_commands=0;
    while(getline(infile,contents[i]))// getting file contents
    {
    i++;
    numberoflines++;
    }
	commands = new char*[numberoflines];
for (int j=0,jcopy=0;j<numberoflines+1;j++)
{
       int pos =contents[j].find("[");
       if(pos!=-1)// "[" is found
       {   variables[in_cntr]=contents[j].substr(0,1);// save the matrix variable char
           line=contents[j].substr(pos);// getting only string with the matrix value
           inputs[in_cntr]=line;// saving matrix value
           in_cntr++;

       }
	else if (contents[j].find("=")!=-1){ int index =0;
	for(int i=0;i<contents[j].length();i++){
	if (contents[j].substr(i,1)==" " || contents[j].substr(i,1)==";") continue;
	index ++;

		}
	commands[jcopy] = new char[(index-1)];
	index =0;
	for(int i=0;i<contents[j].length();i++){
		if (contents[j].substr(i,1)==" " || contents[j].substr(i,1)==";") continue;
		commands[jcopy][index]=contents[j].at(i);
		index++;
		}
	jcopy++;
	}

number_of_commands=jcopy;
}

for(int i=0;i<number_of_commands;i++){
		for(int j=0;j<(sizeof commands[i] / sizeof(char));j++)
		cout<<commands[i][j]<<endl;
		cout<<"End of line "<<i<<endl;

}


	return 0;
}
