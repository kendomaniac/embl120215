#ifndef __MAP_H__
#define __MAP_H__
#include <map>
#endif

#ifndef __FST_H__
#define __FST_H__
#include <fstream>
#endif 

#ifndef __IOS_H__
#define __IOS_H__
#include <iostream>
#endif 

#ifndef __STDL__
#define __STDL__
#include <stdlib.h>
#endif

#ifndef __STR_H__
#define __STR_H__
#include <string.h>
#endif

#ifndef __GNR_H__
	#define __GNR_H__
	#include "general.h"

#endif

#define MAXSIZE 2048


using namespace std;
using namespace general;



namespace general{

Parser::Parser(const std::string& delims, char* buffer)
{
this->delims=delims;
char* stoken=strtok(buffer,this->delims.c_str());
while (stoken!=NULL)
	{
	token.push_back(stoken);
	stoken=strtok(NULL,this->delims.c_str());
	}
}

Parser::Parser(const std::string& delims, const std::string& buffer)
{
this->delims=delims;
//per non modificare la stringa
std::string buffer1(buffer,0,buffer.size());
//per non modificare la stringa
char* stoken=strtok((char*)buffer1.c_str(),this->delims.c_str());
while (stoken!=NULL)
	{
	token.push_back(stoken);
	stoken=strtok(NULL,this->delims.c_str());
	}
}


void Parser::update(const std::string& delims, const std::string& buffer)
{
this->delims=delims;
token.clear();
std::string buffer1(buffer,0,buffer.size());
char* stoken=strtok((char*)buffer1.c_str(),this->delims.c_str());
while (stoken!=NULL)
	{
	token.push_back(stoken);
	stoken=strtok(NULL,this->delims.c_str());
	}
}

void Parser::update(const std::string& delims, char* buffer)
{
this->delims=delims;
token.clear();

char* stoken=strtok((char*)buffer,this->delims.c_str());
while (stoken!=NULL)
	{
	token.push_back(stoken);
	stoken=strtok(NULL,this->delims.c_str());
	}
}

void Parser::set(const int& i,std::string newstr)
{
token[i]=newstr;
}

int Parser::find(const std::string& test)
{
for (int i=0;i<(int)token.size();i++)
	{
	if (token[i]==test)
		{
		return i;
		}
	}
return -1;
}

int Parser::findD(const std::string& test,const std::string& test2)
{
for (int i=0;i<(int)token.size();i++)
	{
	if ((token[i]==test)||(token[i]==test2))
		{
		return i;
		}
	}
return -1;
}

int Parser::findS(const std::string& test)
{
int size= (int) test.size();
string comp="";
for (int i=0;i<(int)token.size();i++)
	{
	comp=token[i].substr(0,size);
	if (comp==test)
		{
		return i;
		}
	}
return -1;
}

std::string Parser::get(int i) 
{
if ((i>=0)&&(i<(int)token.size())) 
	return token[i]; 
else 
	return "";
}

ostream& operator<<(ostream& out,class Parser& data)	
		{
		out<<"-----------------------------------------";
		out<<"\nField  delims: "<<data.delims<<endl;
		for (int i=0;i<(int)data.token.size();i++)
			{
			out <<"Token["<<i<<"]: \""<<data.token[i]<<"\""<<endl;
			}
		out<<"\n-----------------------------------------\n";
		return out;
		}


std::string Parser::get_string()
{
std::string comp="";
for (int i=0;i<(int)token.size();i++)
	{
	comp+=token[i];
	}
return comp;
}

std::string Parser::get_string(char delim)
{
std::string comp="";
if (token.size()>0)  comp=token[0];

for (int i=1;i<(int)token.size();i++)
	{
	comp+=delim+token[i];
	}
return comp;
}
ParserB::ParserB(const std::string& delims, const std::string& buffer)
{

Parser();

this->delims=delims;
string temp="";
unsigned int i=0;
while (i<buffer.size())
	{
	if (buffer[i]!=delims[0])
		{
		temp+=buffer[i];
		}
	else
		{
		unsigned int k=1;
		bool fine=false;
		while ((k<delims.size())&&(!fine))
			{
			 if (buffer[i+k]!=delims[k])
				{
				fine=true;
				}
			k++;
			}
		if (fine==false)
			{
			token.push_back(temp);
			temp="";
			i=i+delims.size()-1;
			}
		else
			{
			temp+=buffer[i];
			}
		}
	i++;
	}
token.push_back(temp);
}

void ParserB::update(const std::string& delims, const std::string& buffer)
{
this->delims=delims;
string temp="";
unsigned int i=0;
while (i<buffer.size())
	{
	if (buffer[i]!=delims[0])
		{
		temp+=buffer[i];
		}
	else
		{
		unsigned int k=1;
		bool fine=false;
		while ((k<delims.size())&&(!fine))
			{
			 if (buffer[i+k]!=delims[k])
				{
				fine=true;
				}
			k++;
			}
		if (fine==false)
			{
			token.push_back(temp);
			temp="";
			i=i+delims.size()-1;
			}
		else
			{
			temp+=buffer[i];
			}
		}
	i++;
	}
token.push_back(temp);
}

}







extern "C" {
//It divides the reads in the input files in two files so that reads in the same pair-end are separated.
//It check the read name to discovered error in the pair-end. I.E. pair-end with more then 2 reads.
int StarParser(int argc, char **argv) {

//time_t time_1,time_4;
if (argc<3)
        {
        //std::cerr<<"\n\nUSE: Parser <input_file> <out_file> \n\n";
	return(-1);//error input parameters
        //exit(EXIT_FAILURE);
        }

//time(&time_1);

ifstream in(argv[1],ifstream::in);
if(!in)
        {
         // cout<<"*****Error opening the input file"<<argv[1] <<"*****\n\n";
        //exit(EXIT_FAILURE);
	  return(-2);//error opening input file
        }

ofstream out(argv[2],ofstream::out);
if(!out) 
        {
          //cout<<"*****Error opening the  output file *****\n\n"; 
	  return(-3);
          //exit(EXIT_FAILURE); 
        }



//cout<<"\n\nSTART EXECUTION..."<<endl;

char buffer[MAXSIZE];
unsigned int line=0;
char delimC[] = "\t :";
Parser parser;


map <string, pair<int,int> > fusion;

while (!in.eof()){

      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
            int num=in.gcount();
      if (buffer[num-1]!='\0'){
        buffer[num]='\0';
        num++;
      }
      if(buffer[0]!='\0'){
        line++;
//	if (line%1000001==0)
//	  cout<<"Processing line:"<<line<<endl;
	parser.update(delimC,buffer);
	if (parser.size()>6){
	  string gene1chr=parser.get(0);
	  string gene2chr=parser.get(3);

	  if  ((gene1chr[0]!='C')&&(gene1chr[0]!='c')){
	    gene1chr="chr"+gene1chr;
	    gene2chr="chr"+gene2chr;
	  }
	      
	  string id_fusion =  gene1chr+":"+parser.get(1)+":"+ parser.get(2)+":"+ gene2chr+":"+ parser.get(4)+":"+ parser.get(5);
#if DEBUG
	  cout<<fusion<<endl;
#endif
	  map <string, pair<int,int> >::iterator it;
	  if ((it=fusion.find(id_fusion))!=fusion.end()){
	    if (parser.get(6)=="-1")
	      it->second.first= it->second.first+1;
	    else
	     it->second.second= it->second.second+1;
	  }
	  else{
	    fusion[id_fusion]=make_pair<int,int>(parser.get(6)=="-1"?1:0,parser.get(6)=="-1"?0:1);
	  }
	}
	else{
//	  cerr<<"\n\t Line "<<line<<"has a wrong format\n";
	  line--;
	}
      }
}
in.close();

out<<"gene1.chr\tgene1.start\tgene1.strand\tgene2.chr\tgene2.start\tgene2.strand\tn.spanning\tn.encompassing"<<endl;
for (map <string, pair<int,int> >::iterator it=fusion.begin();it!=fusion.end();++it){
  parser.update(delimC,it->first);
  out<<parser.get(0)<<"\t"<<parser.get(1)<<"\t"<<parser.get(2)<<"\t"<<parser.get(3)<<"\t"<<parser.get(4)<<"\t"<<parser.get(5)<<"\t"<<it->second.second<<"\t"<<it->second.first<<"\n";
}

out.close();

//cout<<"correct reads. "<<line<<endl;
//time(&time_4);
//cout<<"\n\nEND EXECUTION"<<endl;
//cout<<"\nResults are saved in: "<<argv[2]<<endl;
//cout<<"\n=========================== TIME ===========================\n\n\t";
//cout<<"Total time required: "<<(time_4-time_1)<<"s."<<endl;
//cout<<"\n=========================== TIME ===========================\n\n";
return EXIT_SUCCESS;
}
}
