/*
 * mistar
 * Date: Jun-04-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


//#define COMPUTEDIV
#define DEBUG1
#define REMOVEANC

#include <iostream>
#include <fstream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <RInside.h>                    // for the embedded R via RInside


#include "VecAllPairDistanceResult.h"
#include "utils.h"
#include "MistarParser.h"
#include "MistarPairwiseDiff.h"
#include "GenomicRange.h"
#include "mistarOperations.h"
#include "mistarArmadillo.h"
// #include "ComputeGMM.h"

using namespace std;
// using namespace arma;



int main (int argc, char *argv[]) {
    // vector< vector<double> > testvec = vector< vector<double> >(4, vector<double>(4));
    
    // mat testmat= vectorDouble2vec(testvec);
    // cout<<testmat<<endl;

    // return 0;
    string individualid; //id of the individual with admixture to detect
    string bedFileRegions;
    string fastaIndex  =  "/mnt/454/Altaiensis/users/gabriel/faidx/index.hg19.fai"  ;
    int components =2;
    vector<chrinfo> chrFound;
    uint64_t genomeLength;
    string outprefix;
    vector<string> unadmixedInd;
    vector<string> allInd;

    string usageMandatory=string("")+
	"\t\t"+"--admx   [individual id]"   + "\t"+"Id of the individual with potential admixture (Default: none)\n"+
	"\t\t"+"--regions [sorted bed file]" + "\t"+"Genomic regions to iterate on (Default: none)\n"+
	"\t\t"+"--fai     [file]"            + "\t\t"+"Fasta index for genome (produced by \"samtools faidx\") (default: "+fastaIndex+")\n"+
	"\t\t"+"-o        [out prefix]"      + "\t\t"+"Output prefix (default: none)\n"+
	"\n";

    string clusteringOpt=string("")+
	"\t\t"+"-k      [comp]"            + "\t\t"+"Components for admixted individual (default: "+stringify(components)+")\n"+
	"\n";
	
	    
    string usage=string(""+string(argv[0])+" <mistar file (tabixed)> "+
                        "\n\nOptions:\n\n"+
			"\tMandatory options:\n"+
			usageMandatory+
			"\tClustering options:\n"+
			clusteringOpt);

    
    if(argc == 1 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
        cerr << "Usage  "<<usage<<endl;
        return 1;       
    }

    for(int i=1;i<(argc-1);i++){ 

	if(string(argv[i]) == "-k"){
	    components=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-o"){
	    outprefix=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--admx"){
	    individualid=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--regions"){
	    bedFileRegions=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--fai"){
	    fastaIndex=string(argv[i+1]);
	    i++;
	    continue;
	}

        cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
        return 1; 	

    }

    if(bedFileRegions.empty()){
	cerr<<"The --region option is mandatory"<<endl;
	return 1;
    }

    if(individualid.empty()){
	cerr<<"The --admx option is mandatory"<<endl;
	return 1;
    }

    if(outprefix.empty()){
	cerr<<"The -o option is mandatory"<<endl;
	return 1;
    }

    if(fastaIndex.empty()){
	cerr<<"The --fai option cannot be empty"<<endl;
	return 1;
    }




    















    // if(1){
    
    //    cerr<<"ok "<<endl;
    //opening mistar file
    /////////////////////////////////////////
    //      Initialize variables           //
    /////////////////////////////////////////
    //read fasta index
    readFastaIndex(fastaIndex,chrFound,genomeLength);
    //read bed regions
    map< string, vector<GenomicRange> * > * rangesToConsider = readBEDSortedfile(bedFileRegions);


    MistarParser mp   (string(argv[argc-1]),string(argv[argc-1])+".tbi","1",0,1);
    


    //detecting individual
    unsigned int individualIndex=0;
    bool found=false;

    //vector<AllPairDistanceResult * > vectorAPDR;


    for(unsigned i=0;i<(mp.getPopulationsNames()->size());i++){
	if(i != 1)//skipping anc field
	     allInd.push_back( mp.getPopulationsNames()->at(i) );

	if(mp.getPopulationsNames()->at(i) == individualid){
	    individualIndex=i;
	    found=true;
	    //break;
	}else{
	    if(i != 1)//skipping anc field
		unadmixedInd.push_back( mp.getPopulationsNames()->at(i) );
	}
    }

    if(!found){
	cerr<<"Cannot find individual "<<individualIndex<<endl;
	return 1;
    }

    cerr<<"ok1 "<<endl;




    unsigned int totalAmountRegions=0;
    for(unsigned int i=0;i<(chrFound.size());i++){
	//cout<<chrFound[i].name<<endl;
	if(rangesToConsider->find( chrFound[i].name ) == rangesToConsider->end())
	    continue;
	vector<GenomicRange> * touse=rangesToConsider->at( chrFound[i].name );
	for(unsigned int j=0;j<(touse->size());j++){
	    //cout<<touse->at(j)<<endl;
	    totalAmountRegions++;
	}
    }

    //mat matDistAdmx = zeros<mat>(mp.getPopulationsNames()->size()-2,totalAmountRegions);//minus the individual and the anc.

    cerr<<"ok2 "<<endl;


    // initialize vector of all pairwise distance results
    AllPairDistanceResult    * allAPDR ;     //sum of all distances
    VecAllPairDistanceResult * vectorAPDR = new VecAllPairDistanceResult();   //vector of all distances

    allAPDR = new AllPairDistanceResult ( int(mp.getPopulationsNames()->size()), 
					  *(mp.getPopulationsNames()) );
    

    ////////////////////////////////////////////
    //     end Initialize variables           //
    ////////////////////////////////////////////


















    ////////////////////////////////////////////
    // computing differences for all regions  //
    ////////////////////////////////////////////

#ifndef DEBUG1
    cerr<<"DEBUG1"<<endl;
    //    return 1;
    unsigned int currentRegion=0;
    //This matrix is to store the distance of the admixed individual to the remaining ones
    //matDistAdmx = zeros<mat>(mp.getPopulationsNames()->size()-2,totalAmountRegions);//minus the individual and the anc.
   
    cerr<<"Computing divergence over regions"<<endl;

    for(unsigned int i=0;i<(chrFound.size());i++){
	if(rangesToConsider->find( chrFound[i].name ) == rangesToConsider->end())
	    continue;
	vector<GenomicRange> * touse=rangesToConsider->at( chrFound[i].name );

	for(unsigned int j=0;j<(touse->size());j++){
	    //compute on those:
	    cout<<"region: "<<touse->at(j)<<endl;

	    //mp.repositionIterator(
	    AllPairDistanceResult *  apdr = pairwiseDifferences_calc(&mp,
								     touse->at(j).getChrName(),
								     touse->at(j).getStartCoord(),
								     touse->at(j).getEndCoord());
	    // cout<<"apdr1 #"<<endl<<*apdr<<"#"<<endl;
	    //BEGIN to remove
	    // stringstream ssapdr;
	    // ssapdr<<*apdr;
	    
	    // // cout<<"ss "<<ssapdr.str()<<endl;
	    // AllPairDistanceResult *  apdr2 = new AllPairDistanceResult(ssapdr.str())  ;
	    // cout<<"apdr2 "<<endl<<*apdr2<<endl;
	    //END to remove

	    // return 1;
	    vectorAPDR->push_back( apdr );


	}

    }
    cerr<<"done"<<endl;

    
    /////////////////////////////////////////////
    // END computing distance for all regions  //
    /////////////////////////////////////////////

        






    


    //saving distance results
    cerr<<"ok3 "<<endl;

    ofstream outputFiled;
    outputFiled.open( string(outprefix+".apdr.dat").c_str() );

    if(!outputFiled.good()){
	cerr<<"Cannot write to "+string(outprefix+".dist.dat")<<endl;
	return 1;
    }
	
    outputFiled <<*vectorAPDR<<endl;
    outputFiled.close();


#endif

    cerr<<"ok4"<<endl;



#ifdef DEBUG1
    

    //VecAllPairDistanceResult vectorAPDR;
    ifstream inputFiled;
    cerr<<"Opening distance rest"<<endl;
    inputFiled.open(string(outprefix+".apdr.dat").c_str() , ios::in);
    string apdrdat="";
    if (inputFiled.is_open()){
	string line;
	while ( getline (inputFiled,line)){
	    apdrdat+=line+"\n";
	}
	inputFiled.close();
   }else{
       cerr << "Unable to open file "<< string(outprefix+".apdr.dat") <<endl;
       return 1;
    }

    vectorAPDR = new VecAllPairDistanceResult (apdrdat) ;
    cerr<<"..done"<<endl;

    if(0)
    for(unsigned int i=0;i<(vectorAPDR->size());i++){
    //for(unsigned int i=7;i<(8);i++){
	
	cerr<<"Computing tree "<<i<<" of "<<(vectorAPDR->size())<<endl;
	ofstream outputFileSmall;
	outputFileSmall.open( string(outprefix+stringify(i)+".nw").c_str() );
	if(!outputFileSmall.good()){
	    cerr<<"Cannot write to "+string(outprefix+stringify(i)+".nw")<<endl;
	    return 1;
	}
	
	//cout<<(unadmixedInd.size() -1  )<<endl;
	// arma::mat allDistPrint = allDist;
	// cout.precision(11);
	// cout.setf(ios::fixed);
	// //allDistPrint.raw_print(cout);

	arma::mat allDist2 =  vectorDouble2mat( vectorAPDR->at(i)->computeDNADist("JC69") );
	
	//cout<<"all dist "<<1000*allDist<<endl;
	//cout<<vectorToString(unadmixedInd)<<endl;
	outputFileSmall <<*neighborJoinFromDist( allDist2,
						 allInd,
						 allInd.size())<<endl;
	outputFileSmall.close();

	
    }

    // cout<<*vectorAPDR<<endl;
    //  return 1;
#endif





     cerr<<"ok5"<<endl;


    ////////////////////////////////////////////////
    //      Computing divergence over data        //
    ////////////////////////////////////////////////


    vector< arma::mat > allDistanceMatrices;
    arma::mat matDistAdmx;
    arma::mat matDistAdmx_raw;

#ifdef REMOVEANC	

    matDistAdmx     = zeros<mat>(mp.getPopulationsNames()->size()-3,totalAmountRegions);//minus the individual and the anc.
    matDistAdmx_raw = zeros<mat>(mp.getPopulationsNames()->size()-2,totalAmountRegions);//minus the individual and the anc.

#else

    matDistAdmx     = zeros<mat>(mp.getPopulationsNames()->size()-2,totalAmountRegions);//minus the individual and the anc.
    matDistAdmx_raw = zeros<mat>(mp.getPopulationsNames()->size()-1,totalAmountRegions);//minus the individual and the anc.

    // matDistAdmx = zeros<arma::mat>(mp.getPopulationsNames()->size()-3,totalAmountRegions);//minus the individual and the anc.

#endif

    vector<double> normalizingFactor;

    for(unsigned int i=0;i<(vectorAPDR->size());i++){
	//cout<<i<<"\t"<<totalAmountRegions<<endl;
	const AllPairDistanceResult * apdr;
	apdr = vectorAPDR->at(i);
	vector< vector<double> > dist = apdr->computeDNADist("JC69");

	arma::mat matDist= vectorDouble2mat(dist);
	allDistanceMatrices.push_back(matDist);
	
	// cout<<"test1 "<<i<<endl;

	arma::vec distForInd=matDist.col(individualIndex-1);
	distForInd.shed_row(individualIndex-1);

	matDistAdmx_raw.col(i)=distForInd;
	
	//removing anc
#ifdef REMOVEANC	
	distForInd.shed_row(0);
#endif
	// cout<<"test2 "<<i<<endl;
	

	if(1){
	    //normalize
	    normalizingFactor.push_back(accu(distForInd));
	    distForInd=distForInd/accu(distForInd);

	}

	

	// cout<<"post "<<distForInd<<endl;


	//cout<<matDistAdmx.col(currentRegion)<<endl;
	// cout<<"1: "<<allAPDR->getNumberOfPopulations()<<endl;
	// cout<<"2: "<<apdr->getNumberOfPopulations()<<endl;
	
	//return 1;

	//cout<<1000.0*matDist<<endl;
	//cout<<matDistAdmx.col(currentRegion-1)<<endl;
	//delete(apdr);
	(*allAPDR)+=*apdr;

	//cout<<"allAPDR "<<i<<"\t"<<(*allAPDR)<<endl;

	matDistAdmx.col(i)=distForInd;
	//cout<<"test "<<i<<endl;

    }

  
    ////////////////////////////////////////////////////
    //      end Computing divergence over data        //
    ///////////////////////////////////////////////////










    cerr<<"ok5 "<<endl;




    

    // for(unsigned int i=0;i<normalizingFactor.size();i++)
    // 	cerr<<"i "<<normalizingFactor[i]<<endl;

    // return 1;
    
    ofstream outputFilet;
    outputFilet.open( string(outprefix+".dist.dat").c_str() );
    if(!outputFilet.good()){
	cerr<<"Cannot write to "+string(outprefix+".dist.dat")<<endl;
	return 1;
    }
	
    outputFilet <<"ARMA_MAT_TXT_FN008"<<endl;
    outputFilet <<matDistAdmx.n_cols<<" "<< matDistAdmx.n_rows<<endl;   
    
    // outputFilet.precision(11);
    // outputFilet.setf(ios::fixed);
    matDistAdmx.t().raw_print(outputFilet);
    //outputFilet << matDistAdmx.t();
    outputFilet.close();
    

    cout<<"done writing"<<endl;

    /////////////////////////////////////////
    //      Computing GMM over data        //
    /////////////////////////////////////////

    cerr<<"MAT IN "<<matDistAdmx.t()<<endl;
    cerr<<"Calling "<<endl;
    //    ComputeGMM * cgmm = new ComputeGMM (matDistAdmx.t(), components);
    // cerr<<"MU "<<endl<<cgmm->getMu()<<endl;

    arma::mat mu;
    arma::mat probGauss;
    arma::vec pi;

    RInside R(argc, argv);		// create an embedded R instance
    
    try {

	cerr<<"1"<<endl;
	arma::mat inputmat;
	inputmat.load( string(outprefix+".dist.dat") );
	string txt = "suppressMessages(library(EMCluster))";
	R.parseEvalQ(txt);              // load library, no return value
	cerr<<"2"<<endl;
	R["M"] = inputmat;
	cerr<<"3"<<endl;
	txt = "ret <- init.EM(M, nclass = "+stringify(components)+")";
	R.parseEval(txt);  // assign mat. M to NumericMatrix
	
	cout<<"done"<<endl;
	txt = "ret$Mu;";
	// cout<<R.parseEval (txt)<<endl;
	// // ar
	arma::mat mu_  = 	R.parseEval (txt);
	// cout<<"MU "<<mu<<endl ;

	txt = "ret$pi;";
	arma::vec pi_   = 	R.parseEval (txt);
	// cout<<"pi "<<pi<<endl ;

	txt = "e.step(as.matrix(M),emobj=ret)$Gamma;";
	arma::mat P_ =  R.parseEval(txt); 
	// cout<<P<<endl ;

	mu=mu_;
	pi=pi_;
	probGauss=P_;

    } catch(std::exception& ex) {
        std::cerr << "ComputeGMM::Exception caught : " << ex.what() << std::endl;
    } catch(...) {
        std::cerr << "ComputeGMM::Unknown exception caught" << std::endl;
    }

    // cerr<<mu<<endl;
    // cerr<<pi<<endl;
    cout.precision(11);
    cout.setf(ios::fixed);
    cout<<"mu"<<endl;
    mu.raw_print(cout);
    cout<<"pi"<<endl;
    pi.raw_print(cout);
    // inputmat.raw_print(cout);

    probGauss.raw_print(cout);
    //return 1;
    // mat dataamdxWindow;
    // // bool status =  data.load("40000.mat");
    // bool status =  dataamdxWindow.load( string(outprefix+".dist.dat") );
    
    // if(status == true){
    // 	cout << "loaded okay" << endl;
    // }else{
    // 	cout << "problem with loading" << endl;
    // }

    // //cout<<data<<endl;
    // // Compute the gaussian
    // dataamdxWindow=dataamdxWindow.t();
    // const int gaussians = components;
    // const size_t maxIterations = 250;
    // const double tolerance     = 1e-10;
    // const double perturbation  = 1e-30;
    // EMFit<> em(maxIterations, tolerance, perturbation);

    // cerr<<"before gmm"<<endl;
    // GMM<> gmm(size_t(gaussians), dataamdxWindow.n_rows, em);
    // cerr<<"after gmm"<<endl;
    
    
    // Timer::Start("em");
    // int trials=10;
    // double likelihood = gmm.Estimate(dataamdxWindow, trials);
    // Timer::Stop("em");
    // cerr<<"after em"<<endl;
    
    // cout << "Log-likelihood of estimate: " << likelihood << ".\n";
    // cout<<gmm.Weights()<<endl;
    // //cout<<vectorToString()<<endl;
    
    // // Save results.
    // gmm.Save(outprefix+".gmm.out");

    // for(int k=0;k<components;k++){
    // 	cout<<"k "<<gmm.Means()[k]<<endl;
    // }

    /////////////////////////////////////////////
    //      END computing GMM over data        //
    /////////////////////////////////////////////


    cerr<<"done gmm"<<endl;
     // return 0;


    arma::mat allDist =  vectorDouble2mat( allAPDR->computeDNADist("JC69") );
    
    
    ofstream outputFileRawTree;
    outputFileRawTree.open( string(outprefix+".raw.nw").c_str() );
    if(!outputFileRawTree.good()){
	cerr<<"Cannot write to "+string(outprefix+".raw.nw")<<endl;
	return 1;
    }
       
    //cout<<(unadmixedInd.size() -1  )<<endl;
    arma::mat allDistPrint = allDist;
    cout.precision(11);
    cout.setf(ios::fixed);
    allDistPrint.raw_print(cout);

    //cout<<"all dist "<<1000*allDist<<endl;
    cout<<"unad "<<vectorToString(unadmixedInd)<<endl;
    Tree *  unadmixedTree = neighborJoinFromDist( allDist,
						  allInd,
						  allInd.size());
    outputFileRawTree << *unadmixedTree <<endl;


    outputFileRawTree.close();



    cerr<<"before shed"<<endl;

    //DELETING COLUMN OF INDIVIDUAL
    cerr<<individualIndex<<endl;
    cerr<<allDist.n_rows<<endl;
    cerr<<allDist.n_cols<<endl;
    cout<<allDist<<endl;

    allDist.shed_col(individualIndex-1);
    allDist.shed_row(individualIndex-1);
    cout.precision(11);
    cout.setf(ios::fixed);
    allDist.raw_print(cout);

    cout<<"allDist after shed: "<<endl<<allDist<<endl;
    
    
    if(1){

	vector<double> scalingFact (components,0.0); 
	// vector<double> sumScalling (components,0.0);
	vector<double> sumProb;

	for(unsigned int i=0;i<(vectorAPDR->size());i++){ //for each window
	    for(int k=0;k<components;k++){
		sumProb.push_back(accu(probGauss.col(k)));
	    }
	    
	}

	for(unsigned int i=0;i<(vectorAPDR->size());i++){ //for each window
	    for(int k=0;k<components;k++){
		scalingFact[k]+=normalizingFactor[i] *(probGauss(i,k) / sumProb[k]);
	    }
	}
	
	for(int k=0;k<components;k++){
	    arma::vec unscalledMu = mu.row(k).t() * scalingFact[k];
	    //cout<<k<<"\t"<<mu<<"\t"<<unscalledMu<<endl;
	    //unscalledMu.raw_print(cout);
	    vector<double> distForNode;
	    vector<string> distNames;
	    //put in tree
	    unadmixedTree->addNewNodeUsingDist("comp",distForNode,distNames);
	    return 1;
	}


    }else{
  
    











    //ADDING the new names
    for(int k=0;k<components;k++){
	cerr<<"k1 "<<k<<endl;
	cerr<<pi<<endl;
	cerr<<pi[k]<<endl;
	//#FIX
	unadmixedInd.push_back(individualid+"#"+stringify(k+1)+"_"+stringify(100.0*pi[components-k-1])+"%" );	
    }
    
    //de-normalize
    // double denormFact= allDist.col( allDist.n_rows);
    cout<<"unadmixedInd "<<vectorToString(unadmixedInd)<<endl;

    vector<double> sumProb;
    vector<arma::vec> componentsToAdd;

    
    for(int k=0;k<components;k++){
	sumProb.push_back(accu(probGauss.col(k)));
	cout<<accu(probGauss.col(k))<<endl;
	//componentsToAdd
	//vec A = zeros<vec>( mp.getPopulationsNames()->size() - 3 );
	vec A = zeros<vec>( mp.getPopulationsNames()->size() - 2 );
	cout<<"k "<<k<<"\n"<<A<<endl;
	componentsToAdd.push_back(A);
    }
    

    //computing a weighted 
    for(unsigned int i=0;i<(vectorAPDR->size());i++){ //for each window

	for(int k=0;k<components;k++){
	    // cout<<"1: "<<i<<"-"<<k<<"\t"<<probGauss(i,k)<<"\t"<<sumProb[k]<<"\n"<<matDistAdmx_raw.col(i)<<endl;
	    // cout<<"2: "<<(individualIndex-1)<<endl;
	    arma::vec temp=matDistAdmx_raw.col(i);
	    //temp.shed_col(individualIndex-1);
	    cout<<"temp "<<temp<<endl;
	    componentsToAdd[k] = componentsToAdd[k] + ((probGauss(i,k) / sumProb[k])  * temp);
	    //cout<<"i "<<endl<<matDistAdmx_raw.col(i)<<endl;
	}
    }

    cout<<"componentsToAdd"<<endl;
    for(int k=0;k<components;k++){
	cout.precision(11);
	cout.setf(ios::fixed);
    	cout<<"k "<<k<<endl;
	componentsToAdd[k].raw_print();
    }


     // return 1;






    // //empty apdr
    // vector<AllPairDistanceResult *> vectorComponents;
    // for(int k=0;k<components;k++){
    // 	cerr<<"k2 "<<k<<endl;
    // 	AllPairDistanceResult * newapdr = new  AllPairDistanceResult(int(mp.getPopulationsNames()->size()), 
    // 								     *(  mp.getPopulationsNames()) );
    // 	cerr<<"k3 "<<newapdr->getNumberOfPopulations()<<endl;
    // 	cerr<<"k3 "<<newapdr<<endl;

    // 	vectorComponents.push_back(newapdr);
    // }

    //#FIX
    //arma::Col< size_t >  	labels;

    //    gmm.Classify( dataamdxWindowcout<<labels<<endl;

    // for(unsigned int i=0;i<vectorAPDR->size() ;i++){
    // 	// cerr<<"i1 "<<i<<endl;
    // 	int id=labels[i];
    // 	// cerr<<id<<endl;
    // 	// cerr<<vectorComponents[ id ]<<endl;
    // 	// cerr<<vectorAPDR[i]->getNumberOfPopulations()<<endl;
    // 	*(vectorComponents[ id ]) += *(vectorAPDR->at(i));
    // }
    


    //resize    
    cerr<<allDist<<endl;
    allDist.raw_print(cout);

    allDist.resize(allDist.n_rows+components,
		   allDist.n_cols+components);
    cerr<<"post resize"<<endl;
    cerr<<allDist<<endl;
    allDist.raw_print(cout);

    //adding components

    for(int k=0;k<components;k++){
	cerr<<"k3 "<<k<<endl;
    	//vec admixed=gmm.Means()[k];
	// arma::mat admixedMat = vectorDouble2mat ( vectorComponents[k]->computeDNADist("JC69") );
	// arma::vec toadd = admixedMat.col(individualIndex-1);
	// toadd.shed_row( individualIndex-1 );
	// cerr<<toadd<<endl;

	// // cerr<<sqrt(accu(pow(toadd,2)))<<endl;
	

	// cout<<mu<<endl;
	// cout<<mu.row(k).t()<<endl;
	// cout<<accu(componentsToAdd[k])*mu.row(k).t()<<endl;

	// cout<<componentsToAdd[k]<<endl;
	// return 1;
	arma::vec toadd2 = componentsToAdd[k];
	//arma::vec toadd2 = accu(componentsToAdd[k])*mu.row(k).t();
	cout.precision(11);
	cout.setf(ios::fixed);
	cout<<"toadd2"<<endl;
	toadd2.raw_print(cout);


	//	for(int k1=0;k1<components;k1++)
	toadd2.insert_rows(toadd2.n_rows,components);

	cerr<<"toadd2 "<<toadd2<<endl;
	//for(int k1=1;k1<=components;k1++)
	allDist.col(allDist.n_cols - (k+1) )  = toadd2;
	allDist.row(allDist.n_rows - (k+1) )  = toadd2.t();
    }


    //cerr<<allDist<<endl;
    allDist.raw_print(cout);






    //estimating distance between components

    vector<unsigned int> indexOfMinElemt;//index of the closest elem. for each component
    for(int k1=0;k1<components;k1++){

	//find non-null closest to k1
	vec distK1=allDist.col(allDist.n_cols - (k1 + 1) );
	double       minvalk1=DBL_MAX;
	unsigned int minindk1=0;

	for(unsigned int i=0;i<distK1.n_rows;i++){
	    if(distK1[i] < minvalk1 &&
	       distK1[i] != 0 ){
		minvalk1 = distK1[i];
		minindk1 = i;
	    }
	}
	cerr<<distK1<<endl;
	cerr<<minvalk1<<endl;
	cerr<<minindk1<<endl;

	indexOfMinElemt.push_back(minindk1);
	// for(int k2=0;k2<components;k2++){
	//     vec distK2=allDist.col(allDist.n_cols - (k1 + 1) );	    
	// }

    }

    for(int k1=0;k1<components;k1++){
	for(int k2=0;k2<components;k2++){
	    allDist( (allDist.n_rows - (k1+1) ),(allDist.n_rows - (k2+1) )) = allDist(indexOfMinElemt[k1],indexOfMinElemt[k2]);
	}	
    }

    //cerr<<"final "<<allDist<<endl;
    cout<<"FINAL"<<endl;
    cout.precision(11);
    cout.setf(ios::fixed);
    allDist.raw_print(cout);
    //return 1;

    ofstream outputFileTree;
    outputFileTree.open( string(outprefix+".nw").c_str() );
    if(!outputFileTree.good()){
	cerr<<"Cannot write to "+string(outprefix+".dist.dat")<<endl;
	return 1;
    }
	
    cout<<vectorToString(unadmixedInd)<<endl;
    //cout<<(unadmixedInd.size() -1  )<<endl;

    outputFileTree <<*neighborJoinFromDist( allDist,
					    unadmixedInd,
					    unadmixedInd.size() )<<endl;
    
    //return 1;
    outputFileTree.close();
    //}

    // // for(unsigned int cols=0;cols<dataamdxWindow.n_cols ;cols++){
    // // 	cout<<dataamdxWindow.col(cols).t()<<endl;

    // // 	cout<<gmm.Probability( dataamdxWindow.col(cols),0 )/gmm.Probability(gmm.Means()[0],0)<<endl;
    // // 	cout<<gmm.Probability( dataamdxWindow.col(cols),1 )/gmm.Probability(gmm.Means()[1],1)<<endl;

    // // }
    // MistarParser mp2   (string(argv[argc-1]),string(argv[argc-1])+".tbi","1",0,1);

    // AllPairDistanceResult * componentsDist  [components];
    // for(int i=0;i<(components);i++){
    // 	componentsDist[i] =      new AllPairDistanceResult    ( int(mp2.getPopulationsNames()->size()), *(mp2.getPopulationsNames()) );
    // }



    // unsigned int indexOfRegion=0;
    
    // for(unsigned int i=0;i<(chrFound.size());i++){

    // 	if(rangesToConsider->find( chrFound[i].name ) == rangesToConsider->end())
    // 	    continue;
    // 	vector<GenomicRange> * touse=rangesToConsider->at( chrFound[i].name );

    // 	for(unsigned int j=0;j<(touse->size());j++){
    // 	    //compute on those:
    // 	    cout<<touse->at(j)<<endl;

    // 	    //mp.repositionIterator(
    // 	    AllPairDistanceResult *  apdr = pairwiseDifferences_calc(&mp2,
    // 	    							     touse->at(j).getChrName(),
    // 	    							     touse->at(j).getStartCoord(),
    // 	    							     touse->at(j).getEndCoord());
    // 	    size_t indexd=labels[indexOfRegion];
    // 	    *(componentsDist[ indexd ]) += *apdr;

    // 	    indexOfRegion++;

    // 	}

    // }
    // cerr<<"done"<<endl;



    // for(int i=0;i<(components);i++){
    // 	//vector< vector<double> > dist = componentsDist[ i ]->computeDNADist("JC69");
    // 	mat matDist= vectorDouble2mat(  componentsDist[ i ]->computeDNADist("JC69") );
    // 	matDist=1000.0*matDist;
    // 	cout<<matDist<<endl;		
    // }

    // for(int i=0;i<(components);i++){
    // 	//vector< vector<double> > dist = componentsDist[ i ]->computeDNADist("JC69");
    // 	cout<<*neighborJoinFromDist(    componentsDist[ i ]->computeDNADist("JC69"),
    // 				        componentsDist[ i ]->getPopulationNamesNoAnc(),
    // 				        componentsDist[ i ]->getNumberOfPopulations()-1)<<endl;
	
    // }
    }

    return 0;



    // return 0;
}

