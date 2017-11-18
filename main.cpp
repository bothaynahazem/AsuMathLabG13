#include "Matrix2017.h"

#define INPUT_MATRICES_SIZE 100
#define VARIABLES_NAMES_SIZE 100

int main(int argc,char* argv[])
{
//    try
//    {


	char** commands;
	 //ifstream infile("phase_1_example.txt");
	 //ifstream infile("big_matrix_example.txt");
	 ifstream infile(argv[1]);

	Matrix* input_matrices= new Matrix[INPUT_MATRICES_SIZE];// array to hold matrices values
	string* variables_names=new string[VARIABLES_NAMES_SIZE];// array to hold matrices variable charachter
	//string variables_names=new string[60];// array to hold matrices variable charachter
//	char* variables_names=new char[60];// array to hold matrices variable charachter

	int in_cntr=0;// common index between 2 arrays

	Matrix m;
	string contents[50];//array contains file contents
	string line;

	int numberoflines=0;//index contains number of lines in file
	int i=0;//index for contents array
	//int max_number_of_characters_in_a_command=0;
	int number_of_commands=0;

    while(getline(infile,contents[i]))// getting file contents
    {
        int pos =contents[i].find("[");

        if(pos!=-1 && contents[i].find("]")==-1)
        {
			 getline(infile,contents[i+1],']');
			 contents[i]+=contents[i+1];
        }

        if(contents[i]==";")
        {
            i--;
            numberoflines--;
        }

        numberoflines++;

        i++;

    }

	commands = new char*[numberoflines];

    for (int j=0,jcopy=0;j<numberoflines+1;j++)
    {
       int pos =contents[j].find("[");
       if(pos!=-1)// "[" is found
       {
           variables_names[in_cntr]=contents[j].substr(0,1);// save the matrix variable char
           //make it capital
           line=contents[j].substr(pos);// getting only string with the matrix value
           input_matrices[in_cntr]=line;// saving matrix value
           in_cntr++;
       }

       else if (contents[j].find("=")!=-1)
       {
            int index =0;
            for(int i=0;i<contents[j].length();i++)
            {
                if (contents[j].substr(i,1)==" " || contents[j].substr(i,1)==";")
                    continue;
                index ++;
            }

            commands[jcopy] = new char[(index-1)];
            index =0;

            for(int i=0;i<contents[j].length();i++)
            {
                if (contents[j].substr(i,1)==" " || contents[j].substr(i,1)==";")
                continue;
                commands[jcopy][index]=contents[j].at(i);
                index++;
            }

            jcopy++;
        }

        number_of_commands=jcopy;
    }

    for (int i=0; i<number_of_commands; i++)
    {
        variables_names[i+2]=commands[i][0];  //i+2 --> because the first two are taken A & B
    }

if (variables_names[0]=="a")
    variables_names[0]="A";

Matrix first_matrix,second_matrix,result; //to carry out operations on

string current_var_name;
string first_operand; //for ex: A       /* why string ? because char didn't work as "variables_names is declared as string* and I don't know how to change it uptil now*/
string second_operand; //for ex: B

string test_var_name;
string test_cmd_name;

for (int a=0; a<number_of_commands; a++) //row
    {
        for (int b=0; b<(sizeof commands[i]/sizeof (char)); b++) //coloumn
        {
            if(commands[a][b]==commands[a][3]) // while iterating-->on finding the operation sign
            {
                switch(commands[a][3]) //switch on the operation (add+,subtract-,multiply*,divide/...etc)
                {
                    case '+': //add
                    {
                        for(int c=0; c<number_of_commands; c++)
                        {
                            current_var_name=variables_names[c];
                            first_operand=commands[a][2];
                            second_operand=commands[a][4];

                            if(current_var_name==first_operand)
                                first_matrix=input_matrices[c];

                            if (current_var_name==second_operand)
                                second_matrix=input_matrices[c];
                        }

                        result=first_matrix+second_matrix;

                        test_cmd_name=commands[a][0];
                        for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                        {
                            //on finding the matrix to be modified, modify it then break
                            test_var_name=variables_names[i];
                            if (test_var_name==test_cmd_name)
                            {
                                input_matrices[i]=result;
                                break;
                            }
                        }

                        cout<< commands[a][0] << "=" << "\n" <<result <<"\n\n";
                    }
                    break;

                    case '-': //subtract
                    {

                     for(int c=0; c<number_of_commands; c++)
                        {
                            current_var_name=variables_names[c];
                            first_operand=commands[a][2];
                            second_operand=commands[a][4];

                            if(current_var_name==first_operand)
                                first_matrix=input_matrices[c];

                            if (current_var_name==second_operand)
                                second_matrix=input_matrices[c];
                        }
                        result=first_matrix-second_matrix;

                        test_cmd_name=commands[a][0];
                        for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                        {
                            //on finding the matrix to be modified, modify it then break
                            test_var_name=variables_names[i];
                            if (test_var_name==test_cmd_name)
                            {
                                input_matrices[i]=result;
                                break;
                            }
                        }

                        cout<< commands[a][0] << "=" <<"\n" <<result <<"\n\n";
                    }
                    break;

                    case '*':
                    {
                        for(int c=0; c<number_of_commands; c++)
                        {
                            current_var_name=variables_names[c];
                            first_operand=commands[a][2];
                            second_operand=commands[a][4];

                            if(current_var_name==first_operand)
                                first_matrix=input_matrices[c];

                            if (current_var_name==second_operand)
                                second_matrix=input_matrices[c];
                        }
                        result=first_matrix*second_matrix;

                        test_cmd_name=commands[a][0];
                        for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                        {
                            //on finding the matrix to be modified, modify it then break
                            test_var_name=variables_names[i];
                            if (test_var_name==test_cmd_name)
                            {
                                input_matrices[i]=result;
                                break;
                            }
                        }

                        cout<< commands[a][0] << "=" << "\n" <<result <<"\n\n";
                    }
                    break;

                    case '/':
                    {

                    try{
                        for(int c=0; c<number_of_commands; c++)
                        {
                            current_var_name=variables_names[c];
                            first_operand=commands[a][2];
                            second_operand=commands[a][4];

                            if(current_var_name==first_operand)
                                first_matrix=input_matrices[c];

                            if (current_var_name==second_operand)
                                second_matrix=input_matrices[c];
                        }

                        result=first_matrix/second_matrix;


                        test_cmd_name=commands[a][0];
                        for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                        {
                            //on finding the matrix to be modified, modify it then break
                            test_var_name=variables_names[i];
                            if (test_var_name==test_cmd_name)
                            {
                                input_matrices[i]=result;
                                break;
                            }
                        }

                        cout<< commands[a][0] << "=" <<"\n" <<result <<"\n\n";
}
                    catch(char const* error){cout<<"Error: "<<error<<endl;}
                    }
                    break;

                    /*
                    this case could be split into another nested
                    switch (or if statements)case for other operations
                    starting with a period '.'
                    */
                    case '.':
                    {
                        if (commands[a][4] == '/') //rdivide operation
                        {
                            double inputDouble=0.0; //in case it's a double/Matrix
                            bool MatOverMat=false; //not used but we may use it later

                            for(int c=0; c<number_of_commands; c++)
                            {
                                current_var_name=variables_names[c];
                                string current_double;
                                current_double=commands[a][2];

                                inputDouble=atof(current_double.c_str());
                                first_operand=commands[a][5]; //the only operand

                                if(current_var_name==first_operand)
                                    first_matrix=input_matrices[c];
                            }
                             result=result.rdivide(inputDouble/*1 for example*/, first_matrix);
                        }

                            test_cmd_name=commands[a][0];
                            for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                            {
                                //on finding the matrix to be modified, modify it then break
                                test_var_name=variables_names[i];
                                if (test_var_name==test_cmd_name)
                                {
                                    input_matrices[i]=result;
                                    break;
                                }
                            }
                            cout<< commands[a][0] << "=" << "\n" <<result <<"\n\n";
                        }
                        break;

                    case '\'': //single quote --> Transpose
                    {
                        for(int c=0; c<number_of_commands; c++)
                        {
                            current_var_name=variables_names[c];
                            first_operand=commands[a][2];

                            if(current_var_name==first_operand)
                                first_matrix=input_matrices[c];
                        }
                        result=first_matrix.getTranspose();

                        test_cmd_name=commands[a][0];
                        for (int i=0; i<VARIABLES_NAMES_SIZE; i++)
                        {
                            //on finding the matrix to be modified, modify it then break
                            test_var_name=variables_names[i];
                            if (test_var_name==test_cmd_name)
                            {
                                input_matrices[i]=result;
                                break;
                            }
                        }
                        cout<< commands[a][0] << "=" << "\n" <<result <<"\n\n";
                    }
                    break;

                    default:
                        cout<< "ERROR: Invalid Operation, please try again.";
                    break;
                }
            }
        }
    }

//    }
//
//    catch (char* error){ cout<<"Error: "<<error<<endl; }
	return 0;
}


