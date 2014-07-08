/*
 * testMistar.cpp
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"
#include "MistarParser.h"

using namespace std;

int main (int argc, char *argv[]) {






    
    cout<<"opening "<<argv[1]<<endl;
    
    MistarParser mp (argv[1]);
    AlleleRecords * test;

    cout<<"Header\n\n\n\n "<<mp.getHeader()<<endl;
    cout<<"\n\n\n\n "<<endl;
    while(mp.hasData()){
	cout<<"data"<<endl;
	test = mp.getData();
	cout<<test->coordinate<<endl;
    }


    return 0;
}

