/*
 * tests
 * Date: Feb-23-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <queue>
#include <map>
#include <algorithm>
#include <unistd.h>

#include "utils.h"
#include "MistarParser.h"
#include "SumStatAvgCoa.h"
#include "SumStatD.h"

using namespace std;

//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread

*/     
#define checkResults(string, val) {             \
 if (val) {					\
     cerr<<"Failed with "<<val<<" at "<<string<<endl;   \
   exit(1);                                     \
 }                                              \
}

map<unsigned int, int>       threadID2Rank;

queue< vector<string>  * >  * queueFilesToprocess;

pthread_mutex_t  mutexTHREADID   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexQueue      = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter    = PTHREAD_MUTEX_INITIALIZER;
bool doneReading;

vector<string> * populationNames;
bool allowUndefined;

// template <typename STAT> //type 

// template <class > 
// class Results{  
// private:
//     T val; 
// public:
//   static int count;
//   Test()
//   {
//     count++;
//   }
//   // some other stuff in class
// };

// template <class STAT> //type 
// class Results{  
// private:
//     //T val; 
// public:
//     vector<STAT> results;
//     //static int count;
//     //Test(){
//     //count++;
//     //}
//     // some other stuff in class
// };

template <typename STAT> //type 
void *mainComputationThread(void * argc){

    vector<STAT *> * results =  static_cast<vector<STAT *> *>( argc ) ;

    int   rc;
    // int   stackIndex;
    string freqFileNameToUse;
    int rankThread=0;
    vector < string >  * dataToUse;

    rc = pthread_mutex_lock(&mutexTHREADID);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];
    
    cerr<<"Thread #"<<rankThread <<" is starting"<<endl;

    rc = pthread_mutex_unlock(&mutexTHREADID);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    // cerr<<"Thread #"<<rankThread <<" started and is requesting mutex"<<endl;

    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;

   
    // cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;
    // cerr<<"Thread #"<<(unsigned int)pthread_self() <<" started "<<endl;
    //cout<<"Thread # "<<rankThread<<"TRes "<<results<<endl;

    // cout<<"Thread "<<(unsigned int)pthread_self()<<" taking mutex queue "<<endl;
    if(!queueFilesToprocess->empty()){    
	//cout<<"Thread #"<<rankThread <<" is requesting data"<<endl;
	//cout<<"Thread "<<(unsigned int)pthread_self()<<" taking mutex queue "<<endl;
 	foundData=true;
	dataToUse = queueFilesToprocess->front();
 	queueFilesToprocess->pop();
 	//cerr<<"Thread #"<<rankThread<<" is <<endl;
    }

    
  

    if(!foundData){
 	if(doneReading){

	    rc = pthread_mutex_unlock(&mutexQueue);
	    checkResults("pthread_mutex_unlock()\n", rc);

	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
	    return NULL;	
 	}else{
	    //  cout<<"Queue is empty, thread #"<<rankThread<<" will sleep for 5 seconds and wait for data"<<endl;

	    rc = pthread_mutex_unlock(&mutexQueue);
	    checkResults("pthread_mutex_unlock()\n", rc);

 	    sleep(2);

 	    goto checkqueue;	   
 	}
    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }


    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    //cout<<"Thread #"<<rankThread<<" is starting computations"<<endl;
    //cout<<vectorToString(*populationNames)<<endl;
    // cout<<"Thread #"<<rankThread<<" is starting computations2"<<endl;
    unsigned int count=0;
    // for(unsigned int i=0;i<=dataToUse->size();i++){
    // 	//count++;
    // 	cout<<"mcp "<<i<<"\t"<<(dataToUse->at(i))<<endl;
    // }

    STAT * statComputer = new STAT(populationNames);
    // cout<<"Thread #"<<rankThread <<" addrt stat "<<statComputer<<endl;


    MistarParser mp                  (dataToUse,*populationNames);    
    AlleleRecords  * currentRecord;

    while(mp.hasData()){
	currentRecord = mp.getData() ;

	statComputer->computeStatSingle(currentRecord,allowUndefined);

	//cout<<test->print()<<endl;
    }
    //cout<<"Thread #"<<rankThread <<" read "<<count<<" records "<<results->size()<<endl;
    //cout<<"Thread #"<<rankThread <<" done "<<endl;
    // delete(dataToUse);
    //////////////////////////////////////////////////////////////
    //                  END COMPUTATION                         //
    //////////////////////////////////////////////////////////////
    
    delete dataToUse;
    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
    //cout<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    results->push_back(statComputer);

    //outputToPrint.push_back(toAddToLog);
    
    //cout<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);


    goto checkqueue;	   


    

    
    //cout<<"Thread "<<rankThread<<" ended "<<endl;
    return NULL;
    
}



template <class STAT> //type 
class parallelP{  

public:
    void launchThreads(string filename,int numberOfThreads,int sizeBins, bool produceEveryBoot );
};//end class parallelP

//template <class STAT> //type 

// template <typename STAT> 
// vector<STAT> results;

template <class STAT> //type 
void parallelP<STAT>::launchThreads(string filename,int numberOfThreads,int sizeBins , bool produceEveryBoot){

    doneReading=false;
    queueFilesToprocess = new queue< vector< string >  * >()  ;

    //MistarParser mp                  (filename);
    igzstream myFile;
    myFile.open(filename.c_str(), ios::in);
    vector<string> lines;
    // vector<string> populationNames;
    unsigned int numberPopulations;
   
    allowUndefined  = false;
    //populationNames = mp.getPopulationsNames();

    pthread_mutex_init(&mutexTHREADID,   NULL); 
    pthread_mutex_init(&mutexQueue,      NULL); 
    pthread_mutex_init(&mutexCounter,    NULL);

    pthread_t             thread[numberOfThreads];  
    int                   rc=0;   
    vector<STAT * > * results=new vector<STAT *>();
    //cout<<"res "<<results<<endl;
    //launchThreads(numberOfThreads,*thread);
    for(int i=0;i<numberOfThreads;i++){
        rc = pthread_create(&thread[i], NULL, mainComputationThread<STAT>, results); 
        checkResults("pthread_create()\n", rc); 
    }                                          

    //threads are running here  
    map<string,int> chr2index;
    int chr2indexCurrent=0;
    //    int indexBin;
    int lastBin=-1;
    int chrBin=-1;

    // AlleleRecords  * currentRecord;
    vector< string  > * vecForBin;



   
    if (myFile.good()){
	string line;
	while ( getline (myFile,line)){
	    //while(mp.hasData()){
	    // cout<<"line"<<line<<endl;

	    if(line[0] == '#'){
	       
		if(strBeginsWith(line, "#chr")){
		    vector<string> fields=allTokens(line,'\t');
		    if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
		    if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
		    if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
		    if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
		    if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }

		    for(unsigned int i=3;i<fields.size();i++){
			populationNames->push_back(fields[i]);
			numberPopulations++;
		    }

		    //break;
		}else{
		  
		}
		continue;
	    
	    }
	    

	    //first field chr
	    //second field coord
	    string chrS;
	    string coordS;
	    unsigned int  coordSUI;

	    bool  inChr   = true;
	    bool  inCoord = false;

	    for(unsigned int f=0;f<line.size();f++){
		// cout<<"line["<<f<<"] "<<line[f] <<"\t"<<inChr<<"\t"<<inCoord<<endl;
		if(line[f] == '\t'){

		    if(inChr){
			inChr  =false;
			inCoord=true;
		    }else{
			if(inCoord){
			    //done reading
			    break;
			}else{
			    cerr<<"Parsing error "<<line<<endl;
			    exit(1);
			}
		    }

		}else{
		    if(inChr){
			chrS+=line[f];
		    }else{
			if(inCoord){
			    coordS+=line[f];
			}else{
			    cerr<<"Parsing error "<<line<<endl;
			    exit(1);
			}
		    }
		}
	    }

	    // cout<<"#"<<chrS<<"\t#"<<coordS<<"#"<<endl;
	    coordSUI=destringify<unsigned int>(coordS);

	    // cout<<chrS<<"\t"<<coordS<<endl;
	    // currentRecord = mp.getData() ;
	    // cout<<currentRecord->chr<<"\t"<<coordSUI<<endl;

	    if( chr2index.find(chrS) == chr2index.end() ){//new chr
		chr2index[ chrS ]  = (chr2indexCurrent++);//adding
		// if(lastBin != -1)
		// 	currentBin = lastBin+1+currentBin;
		if(chrBin == -1){
		    chrBin = 0;
		}else{
		    chrBin = lastBin +1; //a step above the last bin
		}
		cerr<<"new chr found, processing chr #: "<<(chrS)<<endl;

	    }
	
	    int currentBin = chrBin+(coordSUI/sizeBins);	


	    if(lastBin != currentBin){
		//cout<<"2: "<<currentRecord->chr<<"\tc="<<coordSUI<<"\tbin="<<(currentBin)<<endl;
		//cout<<"new bin"<<endl;
		if( (currentBin%10)==0){
		    //cerr<<"processing : "<<chrBin<<":"<<coordSUI<<endl;
		    cerr<<"processing chr #: "<<(chrS)<<":"<<coordSUI<<endl;
		}
		if(lastBin == -1){ //first bin
		    //cout<<"new bin first"<<endl;
		    //re-init
		    //vecForBin =  new vector< AlleleRecords > ();		
		    vecForBin =  new vector< string > ();		

		    lastBin=currentBin;
		}else{

		    int rc = pthread_mutex_lock(&mutexQueue);
		    checkResults("pthread_mutex_lock()\n", rc);
		    //cout<<"new bin old\t"<<int(queueFilesToprocess->size())<<"\tresadr\t"<<results<<"\tressize\t"<<results->size()<<endl;
		
		    bool needToAskMutex=false;
		    //add old
		    while( int(queueFilesToprocess->size()) > numberOfThreads ){
			needToAskMutex=true;
			//unlock mutex
			rc = pthread_mutex_unlock(&mutexQueue);
			checkResults("pthread_mutex_unlock()\n", rc);

			//cout<<"Queue is full main threads will sleep for 10 seconds and wait for threads to finish"<<endl;
			sleep(4);
		    }

		    if(needToAskMutex){
			rc = pthread_mutex_lock(&mutexQueue);
			checkResults("pthread_mutex_lock()\n", rc);
		    }

		    
		    queueFilesToprocess->push(vecForBin);


		    rc = pthread_mutex_unlock(&mutexQueue);
		    checkResults("pthread_mutex_unlock()\n", rc);
		
		    //cout<<"pushing "<<vecForBin<<endl;
		    // for(unsigned int j=0;j<vecForBin->size();j++){
		    //     cout<<"v "<<(vecForBin->at(j))<<endl;
		    // }
		    //re-init
		    vecForBin =  new vector< string  > ();		
		    lastBin=currentBin;
		    //
		}
	    }else{
		//cout<<"old bin"<<endl;
	    }
	
	    //cout<<currentRecord->chr<<"\t"<<coordSUI<<endl;
	    //vecForBin->push_back(*currentRecord);
	    vecForBin->push_back(line);

	}//done reading

	myFile.close();

    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }
    doneReading=true;
    //waiting for threads to finish
    for (int i=0; i <numberOfThreads; ++i) {
        rc = pthread_join(thread[i], NULL); 
        checkResults("pthread_join()\n", rc);
    }
    //cout<<"ALL DONE1"<<endl;
    pthread_mutex_destroy(&mutexQueue);   
    pthread_mutex_destroy(&mutexCounter);
    pthread_mutex_destroy(&mutexTHREADID);
    //    cout<<"ALL DONE2"<<endl;



    //DO jacknifing
    if(results->size()>1){
	STAT * allResults=new STAT((*results->at(0)));
	// cout << "VEC "<<vectorToString( *((*results->at(0)).populationNames) )<<endl;
	// cout<<"done all1"<<endl;
	// cout<<allResults->print()<<endl;
	// cout<<"done all2"<<endl;
	vector<STAT * > * jacknife=new vector<STAT *>();

	for (unsigned int i=1; i<results->size() ; ++i) {    
	    //cout<<"add\t"<<i<<endl;	
	    (*allResults)+=(*results->at(i));
	    // //cout<<"done dd1"<<endl;
	    // cout<<allResults->print()<<endl;
	    // cout<<"done dd2"<<endl;
	}

	if(produceEveryBoot){
	    for (unsigned int i=0; i<results->size() ; ++i) {    
		// cout<<i<<"\n#####\n"<<endl;
		//	cout<<results->at(i)<<endl;
		cout<<"---------------------------"<<endl;
		cout<<results->at(i)->print()<<endl;
		//	cout<<i<<"\n#####\n"<<endl;
	    }
	}

	for (unsigned int i=0; i<results->size() ; ++i) {    
	    STAT * test =new STAT (*allResults); //creating a copy
	    *test-=(*results->at(i)); //removing ith block
	    jacknife->push_back(test); 
	    //cout<<"ji "<<i<<endl<<test->print()<<endl;
	}	    
	
	//cout<<"done all1"<<endl;
	//COMMENT allResults contains a matrix of AvgCoaResults
	//add a method for jacknife in allResults
	cout<<allResults->printWithBootstraps(jacknife);
	//cout<<allResults->print()<<endl;
	//cout<<"done all2"<<endl;
	
    }
    pthread_exit(NULL); 
    //cout<<"ALL DONE3"<<endl;
}

int main (int argc, char *argv[]) {

    // queueFilesToprocess->push("file1");
    // queueFilesToprocess->push("file2");
    // queueFilesToprocess->push("file3");
    // queueFilesToprocess->push("file4");
    int numberOfThreads =       1;
    int sizeBins        = 1000000;
    string program="none";
    bool produceEveryBoot=false;
    const string usage=string(string(argv[0])+
                              "\nThis program computers summary stats on mistar files\n\n"+
                              +" <options> [mistar file]"+"\n\n"+

                              "\t"+"-p [program]"  +"\t\t" +"Program to use: (Default: "+program+"\n"+
			      "\t"+"    paircoacompute"+"\tTo compute pairwise average coalescence\n"+
			      "\t"+"    dstat"+"\t\tTo compute triple-wise D-statistics\n\n"+

			      "\t"+"-b\t"  +"\t\t" +"Print every subregion separately (Default: "+booleanAsString(produceEveryBoot)+")\n"+
                              "\t"+"-t [threads]"  +"\t\t" +"Threads to use (Default: "+stringify(numberOfThreads)+")\n"+
                              "\t"+"-s [size bin]" +"\t\t" +"Size of bins (Default: "+stringify(sizeBins)+")\n"+
			      
			      "");

    populationNames = new 	vector<string>  ();
	
    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    for(int i=1;i<(argc-1);i++){ 
        if(string(argv[i]) == "-t" ){
            numberOfThreads=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-s" ){
	    sizeBins=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-p" ){
	    program=string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-b" ){
	    produceEveryBoot=true;
            continue;
        }

        cerr<<"Wrong option "<<argv[i]<<endl;
        return 1;

    }

    if(program == "paircoacompute"){
	parallelP<SumStatAvgCoa> pToRun;
	pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,produceEveryBoot);
    }else{
	if(program == "dstat"){
	    parallelP<SumStatD> pToRun;
	    pToRun.launchThreads(string(argv[argc-1]),numberOfThreads,sizeBins,produceEveryBoot);	
	}else{
	    cerr<<"Wrong program "<<program<<endl;
	    return 1;

	}
    }
    return 0;
}

