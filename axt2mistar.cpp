/*
 * axt2mistar
 * Date: Dec-11-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

   string line;
   string lineSp1;
   string lineSp2;

   igzstream myFile;
   string usage=string(""+string(argv[0])+" [chr name] [name sample]  [axt file]"+
			"\nThis program will parse an axt alignment and print a mistar file\n\n");
   
    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }
    
    




   string chrname    = string(argv[argc-3]);
   string nameSample = string(argv[argc-2]);
   string axtfile    = string(argv[argc-1]);

    string programLine;
    for(int i=0;i<(argc);i++){ 
	programLine+=(string(argv[i])+" ");
    }

   myFile.open(axtfile.c_str(), ios::in);
   cout<<"#MISTAR"<<endl;    	
   cout<<"#PG:"<<programLine<<endl;
   cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
   cout<<"#DATE: "<<getDateString()<<endl;	
   cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<nameSample<<endl;
   
   if (myFile.good()){
     while ( getline (myFile,line)){
	 if(strBeginsWith(line,"#"))
	     continue;
	 
	 vector<string> allToks = allTokens(line,' ');

	 if(allToks[1]  != "chr"+chrname){
	     cerr<<"The chromosome name does not match the one provided on the command line"<<line<<endl;
	     return 1;
	 }

	 unsigned int  startC = destringify<unsigned int>(allToks[2]);
	 unsigned int  endC   = destringify<unsigned int>(allToks[3]);
	 unsigned int  coordC =startC;
	 unsigned int  lastC = startC;
	 bool lastCharWasC=false;

	 if(allToks.size() == 9){
	     getline (myFile,lineSp1);
	     getline (myFile,lineSp2);

	     string  lastToPrintS="";
	     for(unsigned i=0;i<lineSp1.size();i++){
		 if(lineSp1[i] == '-'){

		 }else{
		     
		     if(lineSp2[i] == '-'){
			 
		     }else{
			 //cout<<coordC<<"\t"<<lineSp1[i]<<lineSp2[i]<<endl;
			 string toprintS=chrname+"\t"+stringify(coordC)+"\t";
			 
			 char cRef= char(toupper(lineSp1[i]));
			 char cAlt= char(toupper(lineSp2[i]));
			 string toprint;
			 if(cRef == cAlt){
			     toprintS+=stringify(cRef)+",N\t";
			     toprint = "1,0";
			 }else{
			     toprintS+=stringify(cRef)+","+stringify(cAlt)+"\t";
			     toprint = "0,1";
			 }
			 
			 
			 bool cpgFlag=false;
			 if( (lastC+1) == coordC){
			     if(lastCharWasC && (cRef == 'G' || cAlt == 'G' )){
				 cpgFlag=true;				 
			     }
			 }
			 
			 if(cpgFlag){
			     if(!lastToPrintS.empty())
				 cout<<lastToPrintS<<":1"<<endl;			     
			     toprintS+="0,0:0\t0,0:0\t"+toprint+":1";
			     cout<<toprintS<<endl;
			     lastToPrintS="";
			 }else{
			     if(!lastToPrintS.empty())
				 cout<<lastToPrintS<<":0"<<endl;
			     toprintS+="0,0:0\t0,0:0\t"+toprint;
			     lastToPrintS=toprintS;		       			 
			 }

			 //cout<<toPrintS<<endl;
			 
			 //cpg
			 if(cRef == 'C' || cAlt == 'C')
			     lastCharWasC = true;
			 else
			     lastCharWasC = false;

			 lastC = coordC;
		     }
		     coordC++;
		 }

	     }//end loop
	     if(!lastToPrintS.empty())
		 cout<<lastToPrintS<<":0"<<endl;

	     if((coordC-1) != endC){
		 cerr<<"The block did not end with the proper coordinate got:"<<coordC<<" expected "<<endC<<endl;
		 return 1;
	     }

	     // cout<<startC<<"\t"<<endC<<endl;
	     // return 1;
	     getline (myFile,line);	 
	 }else{
	     cerr<<"Wrong number of fields for line "<<line<<endl;
	     return 1;
	 }
     }
     myFile.close();
   }else{
       cerr << "Unable to open file "<<axtfile<<endl;
       return 1;
    }

   cerr << "Program axt2mistar finished gracefully "<<axtfile<<endl;

    return 0;
}

