// $Id: input.C,v 1.1 2002/06/18 11:26:07 nlambros Exp nlambros $
//#include <libc.h>
#include <iostream>            // C++ 
#include <fstream>             // C++
#include <iomanip>             // C++
#include <cmath>               // C++
#include <assert.h>            // C++
#include <stdio.h>             // C++
#include <string>              // C++
#include <vector>              // C++ 
//FUN 1
using namespace std;

//**************************************
static string app(const int num){

	 static char temp[20];
	 sprintf(temp,"%d",num);
	 return string(temp);
}
//**************************************



class Generate{


private:


    double   rMax ;                   // maximum Radius
    string _gauge  ;

    // misc
    int maxPossibleL ;
    int    lTotalMax ;              //e.g. 4 for S--F, 5 for S--G  
    int    n2eMax ;

    static string    _dat ;
    static string    _inp ;
    static string    _inpDir ;
    static string    _out ;
    static string    _egv ;


    string    _bspeg  ;               // In-Out  filename
    string    _knotFile ;             // Write knots in _KnotFile
    string    _flnm ;                 // Coefficients files for atoms withcore
    int      splinesOrder ;           // Order of B-Splines
    double   zNuclear  ;              // Atomic Number 
    double   rMin ;
    int      noPoints ;                // # of grid points for integral evaluations
    int      gridBsplines ;            // 0 is sine-scale 1 is exp-scale for B-splines
                                       //   basis diagonalization   
    int      gridWaveF ;               // 0 is sine-scale 1 is exp-scale for WaveFunction
                                       //   representation   
    int      noS, noP, noD ;           // # of S, P, D orbitals in frosen core
    int      redmass;                  // 1 = hydrogenic, 2 = positronium
    int      useCoreWaveF ;            // 1 uses Hydrogenic,not 1 BSP-coef input
    int      lCoreMax ;                // lCoreMax+1 is largest core symmetry
    int      lMax ;                    // l(+1) is max l in H-orb used 
    int      lFrom ;                   // for aedmx is used 
    int      noIterationsMax ;         // for HF calculations
    vector<int>           zk ;         // Effective charges for core orbits
    vector<double>      alp1 ;         //  
    vector<double>       r01 ;         //
    vector<int>  noCollocationPoints ; //# of H-like orbitals
    vector<double>  firstKnotPoint ;   //
    double   criterion ;               // convergence for HF calculations 
    double   fractionInitial ;         //  HF
    double   fractionFinal ;           //  HF
    int     idBspeg ;                  //   ??
    double xInitial ;                  // grid parameter
    int     nRinExp ;                  // grid parameter
    string _bspwf   ;                  // wf2e.inp file
    string _cfg     ;                  // configuration file 
    string _hxbsp   ;                  // r12.inp   file
    string _diaglh  ;                  // d2e.inp   file
    string _diagout ;                  // wf2e.inp  file
    string _hdmx    ;                  // dmx1e.inp file
    string _osael   ;                  // dmx2e.inp file   

    vector<int>    Nvalence ;

public:

  Generate(const string & InputFile, const double & R, const string & gauge);
  
  void  writeBspeg(int & nIterations, int & Lmin) ; // write diag1e.inp        
  void  writeGrid() ;               // write grid.inp         
  void  writeBspwf(int & Lmin) ;    // write wf1e.inp           
  void  writeHxbsp() ;              // write r12.inp            
  void  writeDiagl() ;              // write diag2e.inp
  void  writeDiagout(int & Lmax) ;  // write wf2e.inp
  void  writeAedmx(int & Lmin ) ;   // write dmx1e.inp
  void  writeOsael(int & Lmax) ;    // write dmx2e.inp

private:

  void  writeFormated2(ofstream &outFile, int width,
		      int noPerLine, int total, 
		      string fileName, string _tail) ;

  void  writeFormated3(ofstream &outFile, int width,
		      int noPerLine, int total, 
		      string fileName, string _tail) ;

  void  writeFormated4(ofstream &outFile, int width,
		      int noPerLine, int total, 
		      string fileName, string _tail) ;

  void  writeFormated5(ofstream &outFile, int width,
		      int noPerLine, int total, 
		      string fileName, string _tail ) ;
} ;

// Class Definition
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// static variables

string Generate::   _dat = ".dat" ;
string Generate::   _inp = ".inp" ;
string Generate::_inpDir = "inp/" ;
string Generate::   _out = ".out" ;
string Generate::   _egv = "egvh" ;

void error(char *) ;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Generate::Generate(const string & InitFile, const double & R, const string & gauge)
    :rMax(R), _gauge(gauge), 
    splinesOrder(9),         //  Bsplines order
    noS(0),                  //  Max Principal number for the s-core electrons
    noP(0),                  //  Max Principal number for the p-valence electrons
    noD(0),                  //  Max Principal number for the d-valence electron
    useCoreWaveF(1),         //  (= 1) HF calculations, ( = 0) no HF
    lTotalMax(5),            //  ( lTotalMax - 1 ) == Max angular number symmetry
    lCoreMax(0),             //  (lCoreMax + 1) ==  Max angular number for core-electrons 
    noIterationsMax(1),      //  No of iterations for FCHF calculation
    criterion(0.1e-13),      //  HF criterion
    fractionInitial(0.9),    //  
    fractionFinal(0.7),
    idBspeg(6)
{
    int i ;
    ifstream fInput ;
  
    _flnm = "core"  ;
    fInput.open(InitFile.c_str(), ios::in) ;  
    clog << "  INPOUT FILE WAS  : " << InitFile <<"\n" ;

//xxx       B-splines parameter

    fInput >> zNuclear     >> redmass    >> maxPossibleL >> lTotalMax ;
    fInput >> rMin         >> noPoints   >> splinesOrder ;
    fInput >> gridBsplines >> xInitial   >> nRinExp ;

    int   tmp1 ;
    double tmp2 ;

  // number of the collocation points
  for ( i = 0 ; i < lTotalMax ; ++i) { 
       fInput >> tmp1 ;
       noCollocationPoints.push_back(tmp1) ;

  }

  // First non-zero of the grid knot points
  for (  i = 0 ; i < lTotalMax ; ++i){       
       fInput >> tmp2 ;
       firstKnotPoint.push_back(tmp2) ;
  }


  // Polarization Potential parameters


  // polarization potential

  for ( i = 0 ; i < lTotalMax ; ++i ){
       fInput >> tmp2 ;
       alp1.push_back(tmp2) ;
  }


  // r effective
  for ( i = 0 ; i < lTotalMax ; ++i ){

       fInput >> tmp2 ;
       r01.push_back(tmp2) ;
  }


  // z effective
  for ( i = 0 ; i < lTotalMax ; ++i ){
       fInput >> tmp1 ;
       zk.push_back(tmp1) ;
  }

  cout << " Reading input file :        " << InitFile  << "\n"
       << " Atomic Number           Z = " << zNuclear  << "\n"
       << " Box   Radius  (au).     R = " << rMax      << "\n"
       << " # Symmetries           nw = " << lTotalMax << "\n" ;

  // if no core, next does not read any
  // Input Core Files, used for  Mg, Ca, ..., ...,
  // 
  //  if( lCoreMax > 0 ) fInput >> _flnm ;


  // dmx1e programm

  fInput >> n2eMax ;


  // for the dmx1e programm. Principal quantum number depending the core

  for ( i = 0 ; i < maxPossibleL ; ++i )      Nvalence.push_back(1) ;

   // read name of input files. 

  // d1e, core, knot, wf1e, dmx1e 

  fInput >> _bspeg  >> _flnm  >> _knotFile  >> _bspwf >> _hdmx  ;

  //  cfg, r12,  d2e, wf2e, dmx2e 

  fInput >> _cfg    >> _hxbsp >> _diaglh  >> _diagout  >> _osael ;

  cout << " outpout file is " << _bspeg <<"\n" ;

  //##################################################################
  //
  //  Specialize to two-electron atomic systems: 
  //
  //       H-, He, Mg, Ca
  //


     if( zNuclear == 1 || zNuclear == 2 ) {         // He , H-

	  lCoreMax        = 0 ;
	  noS             = 0 ;
	  noP             = 0 ;
	  noD             = 0 ;

     }
     else if( zNuclear == 5 ) {                  // Boron - like
	  
	  lCoreMax        = 2 ;
	  noS             = 2 ;
	  noP             = 0 ;
	  noD             = 0 ;

     }
     else if( zNuclear == 12 ) {                  // Mg - like  
	  
	  lCoreMax        = 2 ;
	  noS             = 2 ;
	  noP             = 1 ;
	  noD             = 0 ;

     }
     else if( zNuclear == 20 ) {                 // Ca
	  
	  // To be improved

	  lCoreMax        = 2 ;
	  noS             = 3 ;
	  noP             = 2 ;
	  noD             = 0 ;

     }
     else {
	 
	 cerr << " Atom with Z = " << zNuclear 
	      << " Not supported yet. Modify input files by hand. \n" ;
	 exit(-1) ;
     }
	 
      Nvalence[0] = noS + 1 ;
      Nvalence[1] = noP + 1 ;
      Nvalence[2] = noD + 1 ;

  fInput.close() ;

}
/*******************************************************************/
//xx This method writes the input files for the bspbeg program
//xx The files bspeg.din
void Generate::writeBspeg(int & nIterations, int & Lmin) {


    cout << "        Lmin = " << Lmin << "\n" ;
    int i ;
    // Use Hartree-Fock or not ?

    if( nIterations > 1) {

	useCoreWaveF    = 1 ;	
	noIterationsMax = nIterations ;
	for (  int j = 0 ; j < lTotalMax ; ++j)   alp1[j] = 0.0 ;

    }
    else{

	useCoreWaveF    = 0 ;
	noIterationsMax = 1 ;	;

    }


     // File name
     string   _bspegInputFile = _inpDir + _bspeg  + _inp ;	

     // define input stream and configure

     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //   useCoreWaveF == 1 ::
     //
     //    lCoreMax == 0  (H-, He)
     //  
     //               Use Hydrogenic WaveFunctions 
     //                                  for the valence electron
     //    else 
     //               First  Generate Core HF waveFunctions      
     //
     //   and then             
     //
     //   useCoreWaveF == 0 :: 
     //
     //                Generate WaveFunctions for the valence electron
     //                using existing Core WaveFunctions calculated 
     //                previously                 
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     ofstream bspegInp ;
     bspegInp.open( _bspegInputFile.c_str() , ios::out) ;
     bspegInp.setf(ios::left, ios::adjustfield) ;


     // Now write in File
//     bspegInp << _bspeg + _out << "\t \t" << _knotFile + _dat << "\n" ;
//     writeFormated2(bspegInp, 16, 5, maxPossibleL, _bspeg, _dat) ;

     
     bspegInp << setw(12)  << zNuclear << setw(12) << rMin << setw(12) <<  rMax << '\n'
	      << setw(12)  <<noPoints  << setw(12) << redmass          << '\n'
	      << setw(12)  << noS      << setw(12) << noP   << setw(12) <<  noD << '\n'
	      << setw(12)  << useCoreWaveF << setw(12) << lTotalMax << '\n'
	      << setw(12)  << Lmin +1  << setw(12) << lCoreMax  << setw(12)  << Lmin + lTotalMax
	      << setw(12)  << noIterationsMax << '\n' ;
  
     for ( i = 0 ; i < lTotalMax ; ++i)   bspegInp << setw(12) << alp1[i] ;
     bspegInp << '\n' ;
     for (i = 0 ; i< lTotalMax; ++i)      bspegInp << setw(12) << r01[i]  ;
     bspegInp << '\n' ;

     //  int ntc = noS + noP + noD;
     //WARNING!!! Next needs reimplementation...
     for (i = 0 ; i< lTotalMax; ++i)  bspegInp << setw(12) << zk[i] ;
     bspegInp << '\n' ;  

     for (i = 0 ; i< lTotalMax; ++i) bspegInp << setw(12) << noCollocationPoints[i] ;
     bspegInp << splinesOrder << '\n';

     for (i = 0; i< lTotalMax; ++i) bspegInp << setw(12) << firstKnotPoint[i] ;
     bspegInp << '\n' ;

     bspegInp << setw(12) << criterion       
	      << setw(12) << fractionInitial 
	      << setw(12) << fractionFinal   
	      << setw(12) << idBspeg         << '\n' ;

//     writeFormated2(bspegInp, 16, 5, lCoreMax, _flnm, _dat);

//     if (lCoreMax == 0)   bspegInp << '\n' ;


     bspegInp.close() ;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeGrid(){


  //       gridBsplines      == 0   sine-like grid              (hard-coded)
  //                         == 1   exp -like grid
  //       xInitial          == 0   First Knot point            (hard-coded)
  //       nrp               == 1   Exponent for exp-like grid. (hard-coded)

     // File name
     string   _gridFile = _inpDir + "grid"  + _inp ;	

     // define input stream and configure

     ofstream rinInp ;
     rinInp.open( _gridFile.c_str() , ios::out) ;
     rinInp.setf(ios::left, ios::adjustfield) ;


//  rinInp.open("grid.inp") ;
/*
  rinInp << gridBsplines   << "\t"  << gridWaveF   << "\n"
         << xInitial       << "\t"  << nRinExp     << "\n"
         << rMin           << "\t"  << rMax        << "\t"
	 << noPoints       << "\n" ;
*/

   /*  fxd Boundary code */

  rinInp << gridBsplines   << "\n"
         << xInitial       << "\t"  << rMin        << "\t"  << rMax        << "\n"
	 << noPoints       << "\t"  << nRinExp     << "\n" ;

  rinInp.close() ;

} 
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeBspwf(int & Lmin){

     // File name
     string   _bspwfInputFile = _inpDir + _bspwf  + _inp ;	

     // define input stream and configure

     ofstream bspwfInp ;
     bspwfInp.open( _bspwfInputFile.c_str() , ios::out) ;
     bspwfInp.setf(ios::left, ios::adjustfield) ;

     // Now write in File

//     writeFormated2( bspwfInp, 16, 5, maxPossibleL, _bspeg, _dat) ;
//     writeFormated2( bspwfInp, 16, 5, maxPossibleL, _bspwf, _dat) ;

//     bspwfInp << _bspwf    +   _out         << "\t" 
//	      << _bspwf    + "-en"   + _dat << "\t"
//	      << _knotFile + "WaveF" + _dat << '\n' ;

     bspwfInp << setw(12) << 1
	      << setw(12) << Lmin + lTotalMax
	      << '\n' ;
     int i ;
     for ( i = 0 ; i < Lmin + lTotalMax ; ++i)   
	 bspwfInp << setw(12) << Nvalence[i] ;
     bspwfInp << '\n';
     bspwfInp << setw(5) << 1                  //to plot P(1,Lmin)
	      << setw(5) << Lmin 
	      << '\n' ;




     bspwfInp.width(0) ;
     bspwfInp.close() ;

  
}
/*******************************************************************/
void Generate::writeHxbsp()
{
     // File name
     string   _hxbspInputFile = _inpDir + _hxbsp  + _inp ;	
     
     // define input stream and configure

     ofstream hxbspInp ;
     hxbspInp.open( _hxbspInputFile.c_str() , ios::out) ;
     hxbspInp.setf(ios::left, ios::adjustfield) ;

     // Now write in File
 
//     hxbspInp << _bspwf  + "-en"   + _dat            << "\t"
//	      << _cfg    + _inp                      << "\t"
//	      << _hxbsp  + _out                      << '\n' ;

//     writeFormated2(hxbspInp, 16, 5, lTotalMax, _bspwf, _dat) ;
//     writeFormated4(hxbspInp, 16, 5, lTotalMax, _hxbsp, _dat) ;

     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     //
     //xxx   From here on only from hx_mix.f code is read 
     //xxx   note ::
     //
     //    noCollocationPoints[0] :: actually needs total number of B-splines
     //                            it is assumed  that 
     //      number of B-splines == noCollocationPoints[0]                       
     //      
     //     alp[0] :: it is assumed that alp == 0
     //    
     //    ** Reform this part in combination with hx_mix.f  07 Oct 99
     //
     //
     //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     hxbspInp << setw(12) << alp1[0] 
	      << setw(12) << redmass - 1  ;     // hard-coded  alp1
     hxbspInp << '\n' ;
     int i ;
     for ( i = 0 ; i < lTotalMax ; ++i) hxbspInp << setw(12) << r01[i] ;
     hxbspInp << '\n' ;


     hxbspInp << setw(12) << splinesOrder           
	      << setw(12) << noCollocationPoints[0] 
	      << setw(12) << rMax                   
	      << '\n' ;

     for ( i = 0 ; i< lTotalMax; ++i) 
	 hxbspInp << setw(12) << firstKnotPoint[i] ;

     hxbspInp << '\n';

     hxbspInp.close() ;

}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeDiagl()
{
     // File name
     string   _diaglInputFile = _inpDir + _diaglh  + _inp ;	

     // define input stream and configure
     ofstream diaglInp ;
     diaglInp.open( _diaglInputFile.c_str() , ios::out) ;
     diaglInp.setf(ios::left, ios::adjustfield) ;

     // Now write in File

//     diaglInp << _cfg  +  _inp << '\n' ;

//     writeFormated4(diaglInp, 16, 5, lTotalMax, _diaglh,".eng") ;
//     writeFormated4(diaglInp, 16, 5, lTotalMax, _diaglh, _dat) ;
//     writeFormated4(diaglInp, 16, 5, lTotalMax,  _hxbsp, _dat);

     diaglInp.width(0) ;
     diaglInp.close() ;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeDiagout(int & Lmax)
{


    int i ;

    for( i = 0 ; i <= Lmax + lTotalMax ; i++ ) {

     // File name
     string   _diagoutInputFile = _inpDir + _diagout + app(i) + _inp ;

     // define input stream and configure

     ofstream diagoutInp ;
     diagoutInp.open( _diagoutInputFile.c_str() , ios::out) ;	
     diagoutInp.setf(ios::left,ios::adjustfield);

     diagoutInp << i <<'\n' ;
     diagoutInp <<setw(12) << 1 << setw(12)<<n2eMax << '\n' ;  
     
     diagoutInp.close() ;

    }
}
/*******************************************************************/
void Generate::writeAedmx(int & Lmin)
{
     // File name
     string   _aedmxInputFile = _inpDir + _hdmx + _inp ;

     // define input stream and configure

     ofstream aedmxInp ;
     aedmxInp.open( _aedmxInputFile.c_str() , ios::out) ;
     aedmxInp.setf(ios::left,ios::adjustfield);


     aedmxInp << setw(12)  <<   Lmin  
              << setw(12)  << ( Lmin + lTotalMax - 1 ) 
	      << setw(12)  <<   _gauge 
	      << '\n' ;

     int i ;
     for ( i = 0 ; i < lTotalMax ; ++i)   
	 aedmxInp << setw(12) << Nvalence[Lmin + i] ;
     aedmxInp << '\n';
     

     for ( i = 0 ; i < lTotalMax  ; ++i) {

	 if(( noCollocationPoints[i] - 2 ) < 2) {
	     cerr << " Incosistent number of B-splines and order. Exiting. \n" 
		  << " N = " << noCollocationPoints[i]              << "\n" 
		  << " k = " << splinesOrder                        << "\n" ;
	     exit(-1) ;
	 }

	 aedmxInp << setw(12) << ( noCollocationPoints[i] - 2 ) ;

	 }
     aedmxInp << '\n';
     aedmxInp.close() ;
}
/*******************************************************************/
void Generate::writeOsael( int & Lmax)
{

    int i;
     // File name
    
    for( i = 0 ; i <= Lmax + lTotalMax ; i++ ) {

     string   _osaelInputFile = _inpDir + _osael + app(i) + _inp ;

     // define input stream and configure

     ofstream osaelInp ;

     osaelInp.open( _osaelInputFile.c_str() , ios::out) ;
     osaelInp.setf(ios::left,ios::adjustfield);


     osaelInp << setw(12)  <<   i 
	      << setw(12)  <<  _gauge 
	      << '\n' ;
	     
     //     osaelInp << _hdmx + _dat << "\t" << _osael + _out <<"\n" ;
     // writeFormated4(osaelInp, 16, 5, lTotalMax, _diagout, _dat) ;
     // writeFormated3(osaelInp, 16, 5, lTotalMax - 1 , _hdmx, _dat) ;
     // close
     osaelInp.close() ;

    }


}
/*******************************************************************/
//                 Formated Routines 
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeFormated2( ofstream &outFile, int width, 
			       int noPerLine,     int total,  
			       string _fileName, string _tail){

     string  _outFileName ;

     outFile.setf(ios::left, ios::adjustfield) ;	

     for ( int i = 0 ; i < total ; ) {
	  for (int j = 0 ; ( j < noPerLine && i < total) ; ++i, ++j){

	       _outFileName = _fileName + "-" + app(i) + _tail ;
	       _outFileName.resize(width, ' ') ;	

	       outFile << _outFileName ;
	  }
	       outFile << '\n' ;
     }
}	
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeFormated3( ofstream &outFile, int width, 
			       int noPerLine,     int total,  
			       string  _fileName, string _tail)
{
     string _outFileName ;
     outFile.setf(ios::left,ios::adjustfield) ;

     for ( int i = 0 ; i < total ; ) {
	  for (int j = 0 ; ( j < noPerLine && i < total) ; ++i, ++j) {

	       _outFileName = _fileName  + "-"
		                        + app(i) 
		                        + app(i+1) 
		                        + _tail ;

	       _outFileName.resize(width, ' ') ;

	       outFile << _outFileName ;
	  }

	  outFile <<'\n' ;
     }	

    
}	
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeFormated4( ofstream &outFile, int width, 
			       int noPerLine,     int total,  
			       string  _fileName, string _tail)
{

     string _outFileName ;
     outFile.setf(ios::left,ios::adjustfield) ;

     int count = 0 ;
     for ( int i = 0 ; i < total ; ++i) {

	  count++ ;

 	  if( count > noPerLine ){
	       count = 1 ;
	       outFile << '\n' ;
	  }

	  _outFileName =  _fileName + "-"
	                           + app(1) 
	                           + app(i) + _tail ;

	  _outFileName.resize(width, ' ') ;

	  outFile << _outFileName ;

	  count++ ;

 	  if( count > noPerLine ){
	       count = 1 ;
	       outFile << '\n' ;
	  }	


	  _outFileName =  _fileName + "-"
	                            + app(3) 
	                            + app(i) + _tail ;

	  _outFileName.resize(width, ' ') ;

	  outFile << _outFileName ;



     }
	  outFile << '\n' ;
}	
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void Generate::writeFormated5( ofstream &outFile, int width, 
			       int noPerLine,     int total,  
			       string  _fileName, string _tail)
{

     string  _outFileName ;
     outFile.setf(ios::left,ios::adjustfield) ;

     int count = 0 ;

     for ( int i = 0 ; i < total ; ++i){
	  
	  count++ ;	

	  if( count > noPerLine ){
	       count = 0 ;
	       outFile << '\n' ;
	  }
	  
	  _outFileName =  _fileName  + "-"
	                             + app(i) 
	                             + app(i+1) + _tail ;

	  _outFileName.resize(width, ' ') ;

	  outFile << _outFileName ;


	  count++ ;

	  if( count > noPerLine ){
	       count = 0 ;
	       outFile << '\n' ;
	  }



	  _outFileName = _fileName   + "-"
	                         + app(i+1) 
	                         + app(i)  + _tail ;

	  _outFileName.resize(width, ' ') ;
	  
	  
	  outFile << _outFileName ;

     }
     
     outFile << "\n" ;
}	
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
void error(char *message)
{
  cerr << message << '\n' ;
  exit(-1) ;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
void   Options(char *option) ;
void   Help() ;
void   getOptions(int argC, char** argV, string & FileName, double & R, int & nIterations, 
                  int & Lmin, string & gauge ) ;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// AS  :: (A)tomic (S)tructure
// DME :: (D)ipole (M)atrix (E)lements

int main(int argc,char **argv) {


     cout <<"\n" ;
     string InitFile ; //= "atom.inp" ;    //default value. 
     int nIterations, Lmin ;
     double R ;
     string gauge ;
  
     getOptions(argc, argv, InitFile, R, nIterations, Lmin, gauge) ;


     Generate Atom(InitFile, R, gauge) ;


     //------------------------- One-Electron Input-Output AS files (fxd-free)
  
     Atom.writeBspeg( nIterations, Lmin ) ;       // diag1e.inp
     Atom.writeGrid()  ;                          // grid.inp
     Atom.writeBspwf( Lmin) ;                     // wf1e.inp


     //------------------------- Two-Electron Input-Output AS files
     Atom.writeHxbsp() ;         // r12.inp      (fxd -free)
     //     Atom.writeDiagl() ;         // diag2e.inp
     //------------------------- Two-Electron Input-Output DME files 
     //     Atom.writeDiagout( Lmin ) ;     // wf2e.inp
     Atom.writeAedmx( Lmin )   ;     //  dmx1e.inp   ( fxd - free   )
     //     Atom.writeOsael( Lmin )   ;     //  dmx2e.inp

     return 0 ;  

}
//FUN 2
//*****************************************************************
void getOptions( int          argC, 
		 char       **argV, 
		 string & FileName, 
		 double        & R, 
		 int & nIterations,
		 int &        Lmin,
		 string    & gauge 
    ) {   

     
    FileName = "A.inp" ;    //default value. 
    nIterations = 1 ;
    Lmin = 0 ;
    R  = 100. ;
    gauge = "l" ;

    for (int i = 1; i < argC; i++) {

	if ( strncmp(argV[i],"-o", 2) == 0 ) {

	    Options(argV[0]);
	}
	else if ( strncmp(argV[i],"-h", 2 ) == 0) { 
	    Help();
	}
	else if( strncmp(argV[i],"-f", 2 ) == 0) {
	    if (++i >= argC )
		Options(argV[0]) ;

	    FileName = argV[i] ;

	}
	else if( strncmp(argV[i],"-i", 2 ) == 0) {
	    if (++i >= argC )
		Options(argV[0]) ;

	    nIterations = atoi(argV[i]) ;

	}
	else if( strncmp(argV[i],"-l", 2 ) == 0) {
	    if (++i >= argC )
		Options(argV[0]) ;

	    Lmin = atoi(argV[i]) ;

	}
	else if( strncmp(argV[i],"-r", 2 ) == 0) {
	    if (++i >= argC )
		Options(argV[0]) ;

	    R = atof(argV[i]) ;

	}
	else if( strncmp(argV[i],"-g", 2 ) == 0) {
	    if (++i >= argC )
		Options(argV[0]) ;

	    gauge = argV[i] ;

	}
	else{

	    cout<< " Look -o option for better results ...   \n" ;
	    exit(1);

	}
    }


}
// FUN 3
//*****************************************************************
void Options(char *option) {

  cout<<"Be carefull.Options are for the time being:\n"
      <<"1. -op     : Options Display\n"
      <<"2. -help   :Informations about the programm\n"
      <<"3. -f      :Initialization File [A.inp]\n"
      <<"4. -l      : ( Lmin - Lmin +lTotalMax) radial SE are solved [1]\n"
      <<"5. -r      :Box radius in a.u. [100]\n"
      <<"6. -g      :gauge selection '[length]','velocity'\n" ;
  exit(-11) ;
}
// FUN 4
//*****************************************************************
void Help(){

  cout<<"This programm creates the input  files for the Bsplines"<<endl;
  cout<<"programm.It needs initially a sample input file.   "<<endl;
  cout<<"For Negative Hydrogen :   InputH    file.          "<<endl;
  cout<<"For Helium            :   InputHe   file.          "<<endl;
  cout<<"Produces the bspeg.din............................."<<endl;
  cout<<"..................................................."<<endl;
  cout<<"Programms that are connected with this files are   "<<endl;
  cout<<"the following:bspeg.f,bspwf.f, hxbsp.f,diagl.f     "<<endl;
  cout<<"aedmx.f,diagout.f,osael.f,norm.f and their subrtines"<<endl;
  cout<<"More informations for their dependecies ora the    "<<endl;
  cout<<"corresponding makefile.Stop here.I kourastnka.     "<<endl;
  cout<<"                                July 96 MPQ days   "<<endl;
  cout<<"Look the Inputhelp file in the same directory for  "<<endl;
  cout<<" brief notes for some of the important parameters  "<<endl;
  cout<<endl;  cout<<endl;
  exit(1);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//EOF



