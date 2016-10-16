/*
 * bam2mistar
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */




#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>

#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "utils.h"
#include "ReadTabix.h"

using namespace BamTools;
using namespace std;

#define MAXCOV 250
#define MIN(a,b) (((a)<(b))?(a):(b))

bool hasLineToPrint;
bool hasLineToPrintCoverage;
string lineToPrint;
char   offsetQual=33;
static unsigned int limitToReOpenFP = 200; //if the coordinate is this far away, we will re-open the file pointer



class mistarVisitor : public PileupVisitor {
  
    public:
    mistarVisitor(const RefVector& references, Fasta * fastaReference,bool useQCFail,const int minBaseQual,const int minCovcutoff,const int maxCovcutoff)
            : PileupVisitor()
	    , m_references(references)
	    , m_fastaReference(fastaReference)
	    , m_useQCFail(useQCFail)
	    , m_minBaseQual(minBaseQual)
	    , m_minCovcutoff(minCovcutoff)
	    , m_maxCovcutoff(maxCovcutoff)
	      //, m_out(out)
        { 
	    previousPosAlign=0;
	    hasCpreviousPos = false;




	    firstLine=true;
	}
        ~mistarVisitor(void) { 
	    
	    if(hasLineToPrint){
		if(hasLineToPrintCoverage)
		    cout<<lineToPrint<<":0"<<endl;
	    }
	    hasLineToPrint=true;		

	}
  
    // PileupVisitor interface implementation




    public:
	// prints coverage results ( tab-delimited )
        void Visit(const PileupPosition& pileupData) {
	    char referenceBase = 'N';
	    char altBase       = 'N';

    	    unsigned int posAlign    = pileupData.Position+1;
    	    int          posAlignInt = int(pileupData.Position+1);

	    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
		cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
		return;
	    }
	    referenceBase = toupper(referenceBase);


	    //cout<<
	    if(referenceBase == 'N')
		return;


	    //m_coverageCounter->at(MIN(pileupData.PileupAlignments.size(),MAXCOV))++;
	    int countRef=0;
	    int countAlt=0;
	    string currentLine;
	    string chimpString;
	    string ancString;
	    int countTotal ;
	    bool coverageAccept;
		
	    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
		if(pileupData.PileupAlignments[i].Alignment.IsFailedQC() && 
		   !m_useQCFail){
		    continue;
		}
		//skip deletion in the reads/insertion in the reference
		if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		    pileupData.PileupAlignments[i].IsNextInsertion ){
		    continue;
		}
		char b =       pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		//char q = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		char q   = int(pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual);

		if(q<m_minBaseQual)//base on read does not pass quality cutoff
		    continue;

		if(b != referenceBase){
		    if(altBase == 'N'){
			altBase = b;
			countAlt++;
		    }else{
			if(b == altBase){
			    countAlt++;
			}else{ //tri-allelic site
			    goto skiptonextpos;
			}
		    }
		}else{
		    countRef++;
		}

	    }

	    countTotal     = countRef + countAlt ;
	    coverageAccept = true;
	    if(countTotal < m_minCovcutoff)
		coverageAccept=false;

	    if(countTotal > m_maxCovcutoff)
		coverageAccept=false;


	    // cerr<<posAlign<<"\talt1\t"<<altBase<<endl;


	    if(countRef != 0 || countAlt != 0 ){
		//		char alt=altBase;//(toprint->getAlt()=="."?'N':toprint->getAlt()[0]);
		

		//unresolved ancestral allele (A,C,G,T)
		//if(!isResolvedDNA(allel_chimp)){
		chimpString="0,0:0";					
		// }
		// //resolved ancestral allele
		// else{
		//     if(allel_chimp == referenceBase){//no diff between chimp and reference
		// 	chimpString="1,0:"+string(cpgEPO?"1":"0");
		//     }else{
		// 	if(altBase == 'N'){//no alt defined, the chimp becomes the alt			    
		// 	    altBase = allel_chimp;
		// 	    chimpString="0,1:"+string(cpgEPO?"1":"0");
		// 	}else{
		// 	    if(altBase == allel_chimp){//alt is chimp
		// 		chimpString="0,1:"+string(cpgEPO?"1":"0");
		// 	    }else{ //tri-allelic site, discard
		// 		//				continue;
		// 		goto skiptonextpos;
		// 	    }
		// 	}
		//     }
		// }
		
		// cerr<<posAlign<<"\talt2\t"<<altBase<<endl;

		//if(!isResolvedDNA(allel_anc)){
		    ancString="0,0:0";					
		// }
		// //resolved ancestral allele
		// else{
		//     if(allel_anc == referenceBase){//no diff between ancestor and reference
		// 	ancString="1,0:"+string(cpgEPO?"1":"0");
		//     }else{
		// 	if(altBase == 'N'){//no alt defined, the ancestor becomes the alt			    
		// 	    altBase = allel_anc;
		// 	    ancString="0,1:"+string(cpgEPO?"1":"0");			    
		// 	}else{
		// 	    if(altBase == allel_anc){//alt is ancestor
		// 		ancString="0,1:"+string(cpgEPO?"1":"0");
		// 	    }else{ //tri-allelic site, discard
		// 		//continue;
		// 		goto skiptonextpos;
		// 	    }
		// 	}
		//      }
		// }

		
		// cerr<<posAlign<<"\talt3\t"<<altBase<<endl;
		currentLine=m_references[pileupData.RefId].RefName + "\t" +stringify(posAlign)+"\t"+ stringify(referenceBase) + ","+stringify(altBase)+"\t"+chimpString+"\t"+ancString+"\t"+stringify(countRef)+","+stringify(countAlt);
		
		// cout<<toprint->getChr()<<"\t"<< toprint->getPosition()<<"\t"<<
		//     toprint->getRef()<<","<<
		//     alt<<"\t"<<
		//     chimpString<<"\t"<<
		//     ancString<<"\t"<<
		//     pairCount.first<<","<<pairCount.second<<":"<<(toprint->isCpg()?"1":"0")<<endl;	


		if( hasLineToPrint && //case of CpG
		    previousPosAlign == (posAlign-1) &&
		    hasCpreviousPos                  &&
		    ( referenceBase == 'G' || altBase  == 'G')){
		    //currentLine+=":1";
		    if(hasLineToPrintCoverage)
			cout<<lineToPrint<<":1"<<endl;
		    if(coverageAccept)
			cout<<currentLine<<":1"<<endl;
		
		    hasLineToPrint=false;
		    lineToPrint ="";
		}else{//no CpG
		    if(hasLineToPrint){
			if(hasLineToPrintCoverage)
			    cout<<lineToPrint<<":0"<<endl;
		    }
		    hasLineToPrint         = true;		
		    lineToPrint            = currentLine;		
		    hasLineToPrintCoverage = coverageAccept;
		}


	    }else{

		//no line to print
		if(hasLineToPrint){
		    if(hasLineToPrintCoverage)
			cout<<lineToPrint<<":0"<<endl;
		}
		hasLineToPrint=false;		

	    }



	skiptonextpos:
	    previousPosAlign = posAlign;
	    hasCpreviousPos  = ( referenceBase == 'C' || altBase  == 'C');

        }
        
    private:
    RefVector m_references;
    Fasta * m_fastaReference;
    unsigned int previousPosAlign;
    bool hasCpreviousPos;
    bool m_useQCFail;
    int m_minBaseQual;
    int m_minCovcutoff;
    int m_maxCovcutoff;
    

    //EPO stuff
    // string epoFile;
    // string epoFileidx;
    // string epoChr;
    // unsigned int epoCoord;

    // ReadTabix * rtEPO ;
    // string lineFromEPO;
    // bool lineLeftEPO;
    // bool cpgEPO;
    bool firstLine;
    char allel_chimp;
    char allel_anc;

};


int main (int argc, char *argv[]) {

    bool useQCFail = false;

    hasLineToPrint=false;
    lineToPrint ="";
    int minBaseQual=0;
    int minCovcutoff=0;
    int maxCovcutoff=1000;

    string usage=string(""+string(argv[0])+" <options> [name sample] [fasta file] [bam file]"+
			"\n\nThis program produces a  mistar matrix given a BAM file\n\n"+
			"\tOptions\n"+	     
			"\t\t"+"--qc"+"\t\t\t"+"Use even the reads that have failed quality control (Default: "+boolStringify(useQCFail)+" )\n"
			"\t\t"+"--qual"+"\t\t\t"+"Base quality cutoffs on a PHRED scale (Default: "+stringify(minBaseQual)+" )\n"
			"\t\t"+"--minCov [cov]" +"\t\t"+"Minimal coverage  (default: "+stringify(minCovcutoff)+")\n"+
			"\t\t"+"--maxCov [cov]" +"\t\t"+"Maximal coverage  (default: "+stringify(maxCovcutoff)+")\n"+

			// "\t\t"+"--het"+"\t"+"Produce two fasta files containing the alleles for het sites (Default: "+boolStringify(printRoot)+" )\n"		
			""
	);

    if(argc < 4 ||
       (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    //all but last arg
    for(int i=1;i<(argc-4);i++){ 

	if( string(argv[i]) == "--qc"){
	    useQCFail=true;
	    continue;
	}

	if( string(argv[i]) == "--qual"){
	    minBaseQual=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}
	
	if(strcmp(argv[i],"--minCov") == 0 ){
	    minCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--maxCov") == 0 ){
	    maxCovcutoff=destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	cerr<<"Unknown option "<<string(argv[i])<<endl;
	return 1;
    }


    string bamfiletopen = string(argv[argc-1]);
    string fastaFile    = string(argv[argc-2]);//fasta file
    string nameSample   = string(argv[argc-3]);//fasta file


    //    string epoFile  = string(argv[argc-1]);
    // string epoFileidx = epoFile+".tbi";
    // string epoChr;
    // unsigned int epoCoord;

    // ReadTabix * rtEPO ;
    // string lineFromEPO;
    // bool lineLeftEPO;
    // bool cpgEPO=false;
    // bool firstLine=true;
    // char allel_chimp;
    // char allel_anc;
    // string regionToUse  = string(argv[argc-3]);//chr:st-end


    // size_t dotdotpos = regionToUse.find("-");
    // if(dotdotpos != string::npos ) {
    // 	regionToUse = regionToUse.substr(0,dotdotpos)+".."+regionToUse.substr(dotdotpos+1);
    // }

    BamReader reader;

     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM file"<<bamfiletopen<< endl;
    	return 1;
     }



     // if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
     // 	 cerr << "Could not open input index BAM files." << endl;
     // 	 return 1;
     // }

    // retrieve reference data
     const RefVector  references = reader.GetReferenceData();

    
     Fasta fastaReference;
     if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){
	 
	 cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	 return false;

     }


     mistarVisitor* cv = new mistarVisitor(references,&fastaReference,useQCFail,minBaseQual, minCovcutoff, maxCovcutoff );
     PileupEngine pileup;
     pileup.AddVisitor(cv);


     string programLine;
     for(int i=0;i<(argc);i++){ 
	 programLine+=(string(argv[i])+" ");
     }
     cout<<"#MISTAR"<<endl;    	
     cout<<"#PG:"<<programLine<<endl;
     cout<<"#GITVERSION: "<<returnGitHubVersion(argv[0],"")<<endl;
     cout<<"#DATE: "<<getDateString()<<endl;
     cout<<"#BAM2MISTAR:"<<endl;
     //cout<<"#chr\tcoord\tREF,ALT"<<"\t"<<nameSample<<endl;
     cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t"<<nameSample<<endl;

     BamAlignment al;
     unsigned int numReads=0;
     while ( reader.GetNextAlignment(al) ) {

	 pileup.AddAlignment(al);
	 numReads++;
     }

     if( hasLineToPrint )
	 cout<<lineToPrint<<":0"<<endl;


     cerr<<"Program bam2mistar terminated gracefully, looked at "<<numReads<< " reads"<<endl;

     //clean up
     pileup.Flush();
     reader.Close();
     fastaReference.Close();

     delete cv;
    

    return 0;
}

