#include <cmath>
//#include <math.h> 
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_struct_ls.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <assert.h>

#ifdef M_PI
  #define PI M_PI
#else
  #define PI 3.14159265358979
#endif

//#define H5
//#define CLUSTER
//#define JUQUEEN
#ifndef NCOL
  #define NCOL 3
#endif
#define LOCAL

#ifdef H5
  #include "hdf5.h"
  #define RANK 2
  char buff[50];
#else
  int ncol = NCOL;
#endif

#define ZERO 0.0
#define TOL 1e-80
  
// 0.1177585 is the mean gap at load 0.01 for N=512
// 0.1177799 is the mean gap at load 0.01 for N=1024
//double slip = 0.014*0.1177585;
double slip = 0;

// #define VERBOSE

using namespace std;

/* Macro to evaluate a function F in the grid point (i,j) */
#define Eval(F,i,j) (F( (ilower[0]+(i))*h, (ilower[1]+(j))*h ))
#define bcEval(F,i,j) (F( (bc_ilower[0]+(i))*h, (bc_ilower[1]+(j))*h ))

#ifdef H5
  void readH5(string name, int pi, int pj, int num_procs, int myid, double * field);
  void writeH5(string name, int pi, int pj, int num_procs, int myid, double * field);
  bool DoesDatasetExist(hid_t, const std::string& rDatasetName);
#endif

int optionU0, optionF;

/* Boundary condition */
double U0(double x, double y)
{
   switch (optionU0)
   {
      case 3:
         if (x <= 0.0) { return 1.0;   }
         if (x >= 1.0) { return 0.0;   }
         else          { return 0.0;   }
      default:
	 if (x <= 0.0) { return 1.0;   }
         if (x >= 1.0) { return 0.0;   }
         else          { return 0.0;   }
   }
}

/* Boundary condition */
double P(double x, double y)
{
   return 1.0-x;
}

/* Right-hand side */
double F(double x, double y)
{
   switch (optionF)
   {
      case 1:
         return 0.0;
      default:
         return 0.0;
   }
}

#ifdef CLUSTER
  enum clusterMode {CONT,NONCONT};

  struct clusterInfo{
    int elem;
    int spanningX, spanningY;  
    int surface;
    int clusterNumber;
  };
  clusterInfo newCluster(int e, int x, int y, int s, int n) 
    { clusterInfo c; c.elem=e; c.spanningX=x; c.spanningY=y; c.surface=s; c.clusterNumber=n; return(c); }
  void addToCluster(vector <clusterInfo> &, int);
  bool clusterOrderNonCont (const clusterInfo & c1, const clusterInfo & c2);
  bool clusterOrderCont (const clusterInfo & c1, const clusterInfo & c2);  
  void clusterFind(vector <clusterInfo> & c, const clusterMode mode);

  int getNeighbor(int , int );
  void initVars();
  void cleanup(int);
  void dumpSurface(int,string);

  const int indCont = 1, indNonCont = 0;   // integers indicating contact / noncontact
  int clusterNumber;                        // running number of clusters
  int totNonCond;                           // how many conducting/nonconducting cells
  
  bool *cont;
  int *contact, *sitesToTest, *checkX, *checkY;
  vector <clusterInfo> ContClusters, nonContClusters;
  int periodicX = 0, periodicY = 1;  
#endif /* CLUSTER */


int n, N;

int main (int argc, char *argv[])
{
 #ifdef H5
   #ifdef CLUSTER
   cerr << "Cluster calculation only works on 1 processor; can't do H5 file, currently. Terminating...\n\n"; exit(0);
   #endif
 #endif
  
   int k;

   int myid, num_procs;

   int pi, pj;
   double h, h2;
   long long int ilower[2], iupper[2];//CM-Fix

   int solver_id;
   int n_pre, n_post;
   int rap, relax, skip;
   //int time_index;
   int iter, iterLocal, iterMax;
   string readIn;
   int alternate_solver, its_alt;

   long long int num_iterations, totalIterations, iter_output = 1024, iter_local = 4096;//CM-Fix
   double final_res_norm = 1e10*TOL;
   double currentTotLoc = 0.0, currentTotAll = 0, currentTotAllOld = 1e80;
   double curColumnRelDevMax, curColumnRelDevMean;
   vector <double> currentCol;
   
 #ifdef LOCAL
   double chiSquTotLoc, chiSquTotAll;
   double deltaPress;
   double chi_dummy;
   long iIterLocal;
   double chiSquNow, chiSquPlu, chiSquMin, *chi;
   double dp, dP;
   bool localDone = 0;
 #endif //LOCAL
   
   int vis;   
 
   HYPRE_StructGrid     grid;
   HYPRE_StructStencil  stencil;
   HYPRE_StructMatrix   A;
   HYPRE_StructVector   b;
   HYPRE_StructVector   x;
   HYPRE_StructSolver   solver;
   HYPRE_StructSolver   precond;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   MPI_Status status; 
   
   double startTimeProgram = MPI_Wtime();
   
   /* Set default parameters */
   n         = 1024;
   optionU0  = 3;   // b.c. == zero
   optionF   = 1;   // rhs  == zero
   solver_id = 11;
   n_pre     = 1;
   n_post    = 1;
   rap       = 0;
   relax     = 1;
   skip      = 0;
   readIn.clear();
   alternate_solver = 0;
   its_alt   = 0;
   iter      = 512;
   iterLocal = 64;
   iterMax   = 16536;
   
   vis       = 0;
   
   /* Parse command line */
   {
      int arg_index = 0;
      bool print_usage = 0;

      if (argc == 1) print_usage = 1;
      
      while ((arg_index < argc) && (!print_usage))
      {
         if ( strcmp(argv[arg_index], "-n") == 0 )
         {
            arg_index++;
            n = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-U0") == 0 )
         {
            arg_index++;
            optionU0 = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-readin") == 0 )
         {
            arg_index++;
            readIn = argv[arg_index++];
         }
         else if ( strcmp(argv[arg_index], "-iter") == 0 )
         {
            arg_index++;
            iter= atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-iterLocal") == 0 )
         {
            arg_index++;
            iterLocal= atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-iterMax") == 0 )
         {
            arg_index++;
            iterMax= atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-alternate") == 0 )
         {
            arg_index++;
            alternate_solver = atoi(argv[arg_index++]);
	    its_alt = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-F") == 0 )
         {
            arg_index++;
            optionF = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-solver") == 0 )
         {
            arg_index++;
            solver_id = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-v") == 0 )
         {
            arg_index++;
            n_pre = atoi(argv[arg_index++]);
            n_post = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-rap") == 0 )
         {
            arg_index++;
            rap = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-relax") == 0 )
         {
            arg_index++;
            relax = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-skip") == 0 )
         {
            arg_index++;
            skip = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-slip") == 0 )
         {
            arg_index++;
            slip = atof(argv[arg_index++]);
         }
         else if (( strcmp(argv[arg_index], "-help") == 0 ) || ( strcmp(argv[arg_index], "--help") == 0 ) || ( strcmp(argv[arg_index], "-?") == 0 ) || ( strcmp(argv[arg_index], "-h") == 0 ))
         {
            print_usage = 1;
            break;
         }
         else
         {
	    arg_index++;
         }
      }

      if ((print_usage) && (myid == 0))
      {
         printf("\n");
         printf("Usage: %s [<options>]\n", argv[0]);
         printf("\n");
         printf("  -n  <n>             : problem size per processor (default: 8)\n");
         printf("  -U0 <U0>            : choice for the boundary condition (default: 3)\n");
         printf("  -F  <F>             : choice for the right-hand side (default: 1) \n");
         printf("  -solver <ID>        : solver ID\n");
         printf("                        0  - SMG \n");
         printf("                        1  - PFMG\n");
         printf("                        10 - CG with SMG precond (default)\n");
         printf("                        11 - CG with PFMG precond\n");
         printf("                        17 - CG with 2-step Jacobi\n");
         printf("                        18 - CG with diagonal scaling\n");
         printf("                        19 - CG\n");
         printf("                        30 - GMRES with SMG precond\n");
         printf("                        31 - GMRES with PFMG precond\n");
         printf("                        37 - GMRES with 2-step Jacobi\n");
         printf("                        38 - GMRES with diagonal scaling\n");
         printf("                        39 - GMRES\n");
         printf("                        40 - BiCGSTAB with SMG precond\n");
         printf("                        41 - BiCGSTAB with PFMG precond\n");
         printf("                        47 - BiCGSTAB with 2-step Jacobi\n");
         printf("                        48 - BiCGSTAB with diagonal scaling\n");
         printf("                        49 - BiCGSTAB\n");
         printf("  -v <n_pre> <n_post> : number of pre and post relaxations\n");
         printf("  -rap <r>            : coarse grid operator type\n");
         printf("                        0 - Galerkin (default)\n");
         printf("                        1 - non-Galerkin ParFlow operators\n");
         printf("                        2 - Galerkin, general operators\n");
         printf("  -relax <r>          : relaxation type\n");
         printf("                        0 - Jacobi\n");
         printf("                        1 - Weighted Jacobi (default)\n");
         printf("                        2 - R/B Gauss-Seidel\n");
         printf("                        3 - R/B Gauss-Seidel (nonsymmetric)\n");
         printf("  -skip <s>           : skip levels in PFMG (0 or 1)\n");
         printf("  -iter <iter>        : proceed in sets of <iter> iterations (default: 512)\n");
	 printf("  -iterLocal <iterLocal>  : proceed in sets of <iterLocal> iterations of the local solver (default: 64)\n");
	 printf("  -iterMax <iterMax>  : end program after <iterMax> iterations, to restart (default: 16536)\n");
	 printf("  -readin <filename>  : read in initial condition for pressure\n");
	 printf("  -alternate <ID2> <its_alt> : insert a set of <iter> iterations of solver <ID2> every <its_alt> sets\n");
	 printf("  -slip <slip>        : work with a (negative) sliplength of <slip>, i.e., decrease the gap by <slip> everywhere\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }
   }
      
   /* Figure out the processor grid (N x N).  The local
      problem size is indicated by n (n x n). pi and pj
      indicate position in the processor grid. */
   
   /*  first check whether we have a quadratic number of processors */
   k = 0;
   while (k * k < num_procs) { ++k; }
   /*   k will be an integer that is not less than the square root of num_procs */
   if (k * k == num_procs) { N = k; } // N = (int) sqrt(num_procs)
   else                    
   { 
     if (myid == 0) 
     {
	cerr << " ##### Error: must use a square processor grid, therefore you can only use 1, 4, 16, 64, 256, 1024, 4096, ... cores.\n";
	cerr << "       Other square numbers (not powers of 2) are possible, but not advisable because each node has an even number of cores.\n\n";
	cerr << "       You specified " <<  num_procs << "...\n\n";
     }
     MPI_Finalize();
     exit(1);
   }
   //N  = (int) sqrt(num_procs);   
   
   h  = 1.0 / (N*n-1);
   h2 = h*h;
   pj = myid / N;
   pi = myid - pj*N;
   
   if (myid == 0) { cout << " # running on " << N << " x " << N << " = " 
                         << num_procs << " processors, each holding " 
			 << n << " x " << n << " points.\n\n"; }

   if (myid == 0) { cout << " # using a slip length of "<< slip << "\n\n"; }
			 
   /* Define the nodes owned by the current processor (each processor's
      piece of the global grid) */
   ilower[0] = pi*n;
   ilower[1] = pj*n;
   iupper[0] = ilower[0] + n-1;
   iupper[1] = ilower[1] + n-1;
   
#ifdef SINGLE   
   float *S = NULL, *loc_S = NULL;
   loc_S   = (float *) calloc((n*n), sizeof(float));
#else
   double *S = NULL, *loc_S = NULL;
   loc_S   = (double *) calloc((n*n), sizeof(double));
#endif   
//    double *topdata, *botdata;
//    topdata = (double *) calloc((N*n), sizeof(double));
//    botdata = (double *) calloc((N*n), sizeof(double));
      
   int solver_id_ini = solver_id;
   if (myid == 0)      
   {      
     cout << " ### Running with solver "<< solver_id_ini;
     if (its_alt) { cout << ", inserting a set of iterations with solver "<<alternate_solver<<" every "<<its_alt*iter<<" iterations."; }
     cout << endl << endl;
   }
   
   MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
   double t1 = MPI_Wtime();
   
#ifdef H5
//    string filen;
//    int filenSize = 0;
//    while (myid == 0) // let only rank == root try to open the file, so that there's not 2000 filehandles
//    {
//      std::ifstream readDummy;
//      sprintf(buff,"_%04d",N*n);     
//      filen.clear(); filen.append("gap"); filen.append(charbuffer); filen.append(".h5");
//      readDummy.open(filen.c_str()); 
//      if (readDummy.is_open()) // search for gap file 
//      { readDummy.close(); filenSize = filen.size(); break; } // file exists, now try to open it with HDF5 (first have to close it again)
//      else 
//      {
//       sprintf(charbuffer,"_%04d_S%d", N*n, 4957); // try another file name format
//       filen.clear(); filen.append("gap"); filen.append(charbuffer); filen.append(".h5");
//       readDummy.open(filen.c_str()); 
//       if (readDummy.is_open()) // search for dump file with coarser resolution for continuation run
//       { readDummy.close(); filenSize = filen.size(); break; } // file exists, now try to open it with HDF5 (first have to close it again)
//       cerr << " ##### Error: file not found! Terminating...\n\n";
//      }
//    }
//    MPI_Bcast(&filenSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    if (filenSize == 0) { MPI_Abort(MPI_COMM_WORLD, -1); }
//    else 
//    {
//      MPI_Bcast(&(filen[0]), filenSize, MPI_CHAR, 0, MPI_COMM_WORLD);
//    }
   
   string filen = "gap.cluster.h5";
   // string filen = "gap"; sprintf(buff,"._%04d.h5",N*n); filen.append(buff);  
   readH5(filen, pi, pj, num_procs, myid, loc_S);   
   
//    filen = "gap";
//    sprintf(buff,"_%04d.testbef.dat",N*n);
//    filen.append(buff);  
//    ofstream testout(filen.c_str());
//    /* transpose the read-in values */
//    for (long ix = 0; ix < n; ++ix) 
//    {
//      for (long iy = 0; iy < n; ++iy)
//        testout << 1.0*ix/(n-1) << " " << 1.0*iy/(n-1) << " " << loc_S[ix*n+iy] << endl;
//      testout << endl;
//    }
//    testout.close();
      
//    // change gap values to conductivity
//    for (long ix = 0; ix < n; ++ix) 
//      for (long iy = 0; iy < n; ++iy)
//        loc_S[ix*n+iy] = pow(loc_S[ix*n+iy],3)/12.0;
   
   double *dummy = (double *) calloc((n*n), sizeof(double));
   memcpy(dummy,loc_S,sizeof(double)*n*n);
      
   for (long ix = 0; ix < n; ++ix) {
     for (long iy = 0; iy < n; ++iy) {
       //loc_S[(n-1-iy)*n+ix] = pow(loc_S[ix*n+iy],3)/12.0;
       //loc_S[ix*n+iy] = pow(loc_S[iy*n+ix],3)/12.0;
       //loc_S[ix*n+iy] = pow(dummy[iy*n+ix],3)/12.0; /* transpose the read-in values */
       loc_S[ix*n+iy] = pow(max(dummy[ix*n+iy]-slip,0.0),3)/12.0;
   } }
   free(dummy);
     
//    filen = "gap";
//    sprintf(buff,"_%04d.testaft.dat",N*n);
//    filen.append(buff);  
//    testout.open(filen.c_str());
//    /* transpose the read-in values */
//    for (long ix = 0; ix < n; ++ix) 
//    {
//      for (long iy = 0; iy < n; ++iy)
//        testout << 1.0*ix/(n-1) << " " << 1.0*iy/(n-1) << " " << loc_S[ix*n+iy] << endl;
//      testout << endl;
//    }
//    testout.close();
   
     
#else //H5  --- not using HDF5, and thus only one processor reading in --- allows for cluster analysis
   
   if (myid == 0)
   {     
      double dummy;
      string restOfLine;
      string file = "gap";
      
      char buffer[50];
      string filename; filename.clear();
      filename.append(file);
      sprintf(buffer,"_%04d.dat",N*n);
      filename.append(buffer);  
      ifstream readDisp ( filename.c_str() ); readDisp.precision(8);    
      if(!readDisp){ cerr << "\n ##### Error opening file " << filename << ". Terminating...\n"; exit(1); }  
      else         { cout << "\n ### Reading column "<<ncol<<" of file " << filename << "...\n\n"; }  
          
      /* generate Matrix S */
    #ifdef SINGLE
      S = (float *) calloc((N*n)*(N*n), sizeof(float));
    #else
      S = (double *) calloc((N*n)*(N*n), sizeof(double));
    #endif      
      for (long ix = 0; ix < N*n; ++ix) 
      {
	for (long iy = 0; iy < N*n; ++iy)
	{
	  for (long icol = 1; icol <= ncol; ++icol) { readDisp >> dummy; }
	  getline(readDisp,restOfLine);
// 	#ifdef JUQUEEN
// 	  readDisp >> dummy; getline(readDisp,restOfLine);
// 	#else
// 	  readDisp >> x >> y >> dummy; getline(readDisp,restOfLine);
// 	#endif
	#ifndef CLUSTER
          S[ix*N*n+iy] = (dummy == 0) ? ZERO : pow(dummy,3)/12.0;
          //S[ix*N*n+iy] = (dummy == 0) ? ZERO : dummy; ///NOTE TODO
	#else          
	  S[ix*N*n+iy] = (dummy == 0) ? ZERO : dummy;
	#endif
	  //S[ix*N*n+iy] = pow(dummy,3)/12.0;
	  //S[ix*N*n+iy] = dummy;
	#ifdef VERBOSE
	  cout.precision(8);
	  cout.setf( std::ios::scientific ); //cout.setf(ios::fixed);
	  // cout.setf(ios::fixed,ios::floatfield); // cout.width(10);
	  cout << "x: " << x << "\t" << "y: " << y << "\t"	  
	       << "ix: " << ix << "\t" << "iy: " << iy << "\t" 
	       << "dummy: " << dummy << "\t" << "suscep: " << S[ix*N*n+iy] << endl;
	#endif // VERBOSE
	}
      }
      readDisp.close();
            
//       cout << "### Global matrix: \n"; 
//       for (long ii=0; ii<N*n; ii++) 
//       { 
// 	for (long jj=0; jj<N*n; jj++) 
// 	{ cout << "S["<<ii<<"]["<<jj<<"]: " << S[ii*(N*n)+jj] << " \n"; }
// 	cout << endl;
//       }
      
      #ifdef CLUSTER
	// ===================================================================== 
	// =================== cluster analysis ================================  
	// =====================================================================
	contact     = (int*)  calloc (N*n*N*n,sizeof(int));
	cont        = (bool*) calloc (N*n*N*n,sizeof(bool));  
	sitesToTest = (int*)  calloc (N*n*N*n,sizeof(int));
	checkX      = (int*)  calloc (N*n,sizeof(int));
	checkY      = (int*)  calloc (N*n,sizeof(int));
	
	for (long ix = 0; ix < N*n; ++ix) 
	{
	  for (long iy = 0; iy < N*n; ++iy)
	  {
	  #ifdef SINGLE
	    contact[ix*N*n+iy] = (S[ix*N*n+iy] < 1e-30) ? 1 : 0;
	  #else
	    contact[ix*N*n+iy] = (S[ix*N*n+iy] < 1e-40) ? 1 : 0;
	  #endif
	    cont[ix*N*n+iy] = contact[ix*N*n+iy];      
	  }
	}  
	dumpSurface(0,filename);
	clusterFind(ContClusters, CONT);  
	clusterFind(nonContClusters, NONCONT);
	
	bool problemCont = 0;
// 	for (uint iCluster = 0; iCluster < ContClusters.size(); ++iCluster) {
// 	  if (ContClusters[iCluster].spanningY) { problemCont = 1; }
// 	  if (ContClusters[iCluster].spanningX && problemCont) {
// 	    cout << "\n### Contact percolates in both directions. No current is possible. Exiting...\n\n"; 
// 	    cleanup(0); 
// 	    exit(0); 
// 	  }    
// 	}
	bool problemNonCont = 1;
	vector <int> eligibleClusters;
	for (uint iCluster = 0; iCluster < nonContClusters.size(); ++iCluster) 
	{
	  if (nonContClusters[iCluster].spanningX) 
	  { 
	    problemNonCont = 0; 
	    for (long i = 0; i < n*N*n*N; ++i) 
	    { 
	      if (contact[i] == nonContClusters[iCluster].clusterNumber) { contact[i] = 0; cont[i] = 0; }
	    }
	    eligibleClusters.push_back(nonContClusters[iCluster].clusterNumber);
	  }    
	}
	if (problemNonCont) 
	{
	  if (problemCont) 
	  {
	    cout << "\n ### Contact percolates in y-directions, and simultaneously," 
		 << "non-contact does not percolate in x-direction (i.e. no stripes, either). "
		 << "No current is possible. Exiting...\n\n"; 
	    cleanup(0); 
	    exit(0); 
	  } 
	  else 
	  {
	    cout << "\n ### There is no non-contact cluster percolating in x-direction. "
		 << "No current is possible. Exiting...\n\n"; 
	    cleanup(0); 
	    exit(0);
	  }      
	}    
	// call contact everything that does not belong to a percolating cluster of non-contact
	for (long ix = 0; ix<n*N; ++ix) 
	{
	  for (long iy = 0; iy<n*N; ++iy) 
	  {
	    if (contact[ix*n*N+iy] != 0) { cont[ix*n*N+iy] = 1; S[ix*N*n+iy] = ZERO; }
	  } 
	}
	dumpSurface(9,filename);
      #ifndef JUQUEEN
	///* dump the modified surface as gap_4096.cluster.dat for use with JUQUEEN (but not ON JUQUEEN) *///
	file = "gap";      
	filename.clear();
	filename.append(file);
	sprintf(buffer,"_%04d.cluster.dat",N*n);
	filename.append(buffer);  
	ofstream writeDisp ( filename.c_str() ); writeDisp.precision(8);    
	if(!writeDisp){ cerr << "\n ##### Error opening file " << filename << ". Terminating...\n"; exit(1); }  
	else          { cout << "\n ### Writing to file " << filename << "...\n\n"; }
	for (long ix = 0; ix < N*n; ++ix) 
	{
	  for (long iy = 0; iy < N*n; ++iy)
	  {
	    writeDisp << S[ix*N*n+iy] << endl;
	  }
	  writeDisp << endl;
	}
	writeDisp.close();
	cleanup(0);
	
	#ifdef SINGLE
	cout << "\n ### ...done.\n\n";
	exit(0);
	#endif
      #endif      
	// now calculate conductivity from height profile
	for (long ix = 0; ix<n*N; ++ix) 
	{
	  for (long iy = 0; iy<n*N; ++iy) 
	  {
	    if (contact[ix*n*N+iy] == 0) { S[ix*N*n+iy] = pow(S[ix*N*n+iy],3)/12.0; }
	  } 
	}
	cleanup(0);
	// ===================================================================== 
	// =================== end cluster analysis ============================
	// =====================================================================      
      #endif /* CLUSTER */
      
   } // end if (myid == 0)  

   int COLS = N*n;
   int ROWS = N*n;
   
   int NPROWS = N; /* number of horizontal blocks */
   int NPCOLS = N; /* number of vertical blocks   */
   
   int BLOCKROWS = n; /* number of rows in each _block_ */
   int BLOCKCOLS = n; /* number of cols in each _block_ */
      
   MPI_Datatype blocktype;
   MPI_Datatype blocktype2;   
      
   MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_DOUBLE, &blocktype2);
   MPI_Type_create_resized( blocktype2, 0, sizeof(double), &blocktype);
   MPI_Type_commit(&blocktype);
   MPI_Type_commit(&blocktype2);

   int disps[NPROWS*NPCOLS];  // vector of gaps, one for each block
   int counts[NPROWS*NPCOLS]; // vector of number or elements, one for each block
   for (long ii=0; ii<NPROWS; ii++) 
   {
     for (long jj=0; jj<NPCOLS; jj++) 
     {       
       disps[ii*NPCOLS+jj] = ii*BLOCKROWS+jj*ROWS*BLOCKCOLS;    /// NOTE proc 1 located on the RIGHT of proc 0       
       //disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;  /// NOTE proc 1 located ABOVE proc 0
       counts [ii*NPCOLS+jj] = 1;
       //if (myid == 0)
       //cout << " ii: " <<  ii << " jj: " << jj << "; disps["<<ii*NPCOLS+jj<<"]: " << disps[ii*NPCOLS+jj] 
       //                                       << "; counts["<<ii*NPCOLS+jj<<"]: " << counts[ii*NPCOLS+jj] << endl;
     }
   }
   
   MPI_Scatterv(S, counts, disps, blocktype, loc_S, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
   
#endif //H5      
   double t2 = MPI_Wtime()-t1;
 
//    /* each proc prints its "loc_S" out, in order */
//    for (long proc = 0; proc < num_procs; proc++) 
//    {
//      if (proc == myid) 
//      {
//        cout << "Proc # " << myid << "; local Matrix:\n";
//        for (long ii=0; ii<n; ii++) 
//        {
//          for (long jj=0; jj<n; jj++)
// 	 //{ cout << loc_S[ii*n+jj] << " "; }
// 	 { cout << "id: " << myid << "; pi: " << pi << ", pj: " << pj << "; loc_S["<<ii<<"]["<<jj<<"]: " << loc_S[ii*n+jj] << " \n"; }	 
//          cout << endl;
//        }
//        printf("\n");
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//     }
//    MPI_Finalize();
//    exit(1);
   
   vector < vector <double> > s;
   s.resize(n+2);
   for (long ii = 0; ii < n+2; ++ii) { s[ii].resize(n+2,0.0); }   

   /* interior points */
   for (long ii = 1; ii < n+1; ++ii)
   { for (long jj = 1; jj < n+1; ++jj)
     { s[ii][jj] = loc_S[(ii-1)*n+(jj-1)]; } }
         
   /* set up neighbors */
   int topNeigh = 0, botNeigh = 0, leftNeigh = 0, rightNeigh = 0;
   if (N > 1)
   {
     ///* NOTE proc 1 located to the RIGHT of proc 0 */ 
     leftNeigh  = (pi > 0)   ? myid - 1 : MPI_PROC_NULL; // if on boundary, has no left neighbor
     rightNeigh = (pi < N-1) ? myid + 1 : MPI_PROC_NULL; // if on boundary  has no right neighbor
     botNeigh   = (myid - N + N*N)%(N*N); // which proc is below me? periodicity!
     topNeigh   = (myid + N)      %(N*N); // which proc is above me? periodicity!
     ///* NOTE proc 1 located ABOVE proc 0 */
     //leftNeigh  = (pj > 0)   ? myid - N : MPI_PROC_NULL; // if on boundary, has no left neighbor
     //rightNeigh = (pj < N-1) ? myid + N : MPI_PROC_NULL; // if on boundary  has no right neighbor
     //botNeigh   = (pi > 0)   ? myid - 1 : myid + (N-1); // which proc is below me? periodicity!
     //topNeigh   = (pi < N-1) ? myid + 1 : myid - (N-1); // which proc is above me? periodicity!     
   }   
      
   //cout << " ### system value for MPI_PROC_NULL: " << MPI_PROC_NULL << endl;
      
//    /* check if neighbors are set up properly */      
//    for (long i = 0; i < num_procs; ++i)  
//    {
//       if (myid == i)
//       cout << "myid: " << myid << " " << "\t" 
//            << "pi: " << pi << "\t" << "pj: " <<  pj << "\t"
//            << "loc_S[0]: " << loc_S[0] << "\t"
//            << "loc_S[1]: " << loc_S[1] << "\t"
//            << "top: " << topNeigh << "\t" << "bot: " << botNeigh << "\t"
//            << "left: " << leftNeigh << "\t" << "right: " << rightNeigh << "\t"
//            << endl;
//       MPI_Barrier(MPI_COMM_WORLD);
//    }
//    MPI_Finalize();
//    exit(1);
      
   double *sendbuf = (double *) calloc(n, sizeof(double));
   double *recvbuf = (double *) calloc(n, sizeof(double));
   double *recvbuflarge = NULL;
   if (myid == 0) { recvbuflarge = (double *) calloc(N*N*n, sizeof(double)); } // to hold the current in each column, for each proc
   
   /* send top row to top neighbor's bottom row */
   for (long ii = 0; ii < n; ++ii) { sendbuf[ii] = loc_S[ii*n+(n-1)]; recvbuf[ii] = 0.0; }
   MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, topNeigh, 4, 
                 recvbuf, n, MPI_DOUBLE, botNeigh, 4, MPI_COMM_WORLD, &status );
     
   for (long ii = 1; ii < n+1; ++ii) { s[ii][0] = recvbuf[ii-1]; }
     
   /* send bottom row to bottom neighbor's top row */
   for (long ii = 0; ii < n; ++ii) { sendbuf[ii] = loc_S[ii*n]; recvbuf[ii] = 0.0; }
   MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, botNeigh, 3, 
                 recvbuf, n, MPI_DOUBLE, topNeigh, 3, MPI_COMM_WORLD, &status );  
     
   for (long ii = 1; ii < n+1; ++ii) { s[ii][n+1] = recvbuf[ii-1]; }
      
   /* choose cells right/left of the right/left boundaries equal to 
      their neighbors so that harmonic mean will give correct value */
   
   /* send left row to left neighbor's right row */
   for (long jj = 0; jj < n; ++jj) { sendbuf[jj] = loc_S[0*n+jj]; recvbuf[jj] = 0.0; }
//    if (myid == 2) { for (long jj = 0; jj < n; ++jj) { cout << sendbuf[jj] << endl; } }
   MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, leftNeigh, 1, 
                 recvbuf, n, MPI_DOUBLE, rightNeigh, 1, MPI_COMM_WORLD, &status );  
      
   if (pi == N-1) /* right boundary, x = 1 */   ///* NOTE proc 1 located to the RIGHT of proc 0 */ 
   //if (pj == N-1) /* right boundary, x = 1 */ ///* NOTE proc 1 located ABOVE proc 0 */
     for (long jj = 0; jj < n; ++jj) { recvbuf[jj] = loc_S[(n-1)*n+jj]; }
   for (long jj = 1; jj < n+1; ++jj) { s[n+1][jj] = recvbuf[jj-1]; }
      
   /* send right row to right neighbor's left row */
   for (long jj = 0; jj < n; ++jj) { sendbuf[jj] = loc_S[(n-1)*n+jj]; recvbuf[jj] = 0.0; }
   MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, rightNeigh, 2, 
                 recvbuf, n, MPI_DOUBLE, leftNeigh, 2, MPI_COMM_WORLD, &status );  
   
   if (pi == 0) /* left boundary, x = 0 */   ///* NOTE proc 1 located to the RIGHT of proc 0 */ 
   //if (pj == 0) /* left boundary, x = 0 */ ///* NOTE proc 1 located ABOVE proc 0 */
     for (long jj = 0; jj < n; ++jj) { recvbuf[jj] = loc_S[0*n+jj]; }
   for (long jj = 1; jj < n+1; ++jj) { s[0][jj] = recvbuf[jj-1]; }
   
//    /* check if recv'd and set up properly */
//    for (long i = 0; i < num_procs; ++i)  
//    {
//       if (myid == i)      
//       for (long ii = 0; ii < n+2; ++ii)
//       {
// 	for (long jj = 0; jj < n+2; ++jj)
// 	{
// 	  cout << "id: " << myid << "; pi: " << pi << ", pj: " << pj << "; s["<<ii<<"]["<<jj<<"]: " << s[ii][jj] << endl;
// 	}
// 	cout << endl;
//       }
//       MPI_Barrier(MPI_COMM_WORLD);
//    }   
//      
//    free(S);
//    free(loc_S);
//    
//    MPI_Finalize();
//    exit(1);
   
   free(S);
   free(loc_S);
   
   MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
   
   
// #ifdef H5
//    // ===============================================================
//    // ======================= HDF5 test =============================
//    // ===============================================================
//    
//    /*
//     * HDF5 APIs definitions
//     */ 	
//    hid_t       file_id, dset_id;         /* file and dataset identifiers */
//    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
//    hsize_t     dimsf[RANK];                 /* dataset dimensions */
//    hsize_t     chunk_dims[RANK];            /* chunk dimensions */
//    double      *data;                    /* pointer to data buffer to write */
//    hsize_t	count[RANK];	          /* hyperslab selection parameters */
//    hsize_t	stride[RANK];
//    hsize_t	block[RANK];
//    hsize_t	offset[RANK];
//    hid_t	plist_id;                 /* property list identifier */
//    herr_t	status_h5;
//     
//    MPI_Info info  = MPI_INFO_NULL;
//    
//    /* 
//     * Set up file access property list with parallel I/O access
//     */
//    plist_id = H5Pcreate(H5P_FILE_ACCESS);
//    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
// 
// #define H5FILE_NAME "gap.h5"
// #define DATASETNAME "gap"
// #define RANK 2
//    
//    /*
//     * Create a new file collectively and release property list identifier.
//     */
//    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//    H5Pclose(plist_id);
// 
//    /*
//     * Create the dataspace for the dataset.
//     */
//    dimsf[0] = N*n;
//    dimsf[1] = N*n;
//    chunk_dims[0] = n;   
//    chunk_dims[1] = n;   
//    filespace = H5Screate_simple(RANK, dimsf, NULL); 
//    memspace  = H5Screate_simple(RANK, chunk_dims, NULL);    
//         
//    /*
//     * Create chunked dataset.
//     */
//    plist_id = H5Pcreate(H5P_DATASET_CREATE);
//    status_h5 = H5Pset_layout( plist_id, H5D_CHUNKED );
//    assert(status_h5 >= 0);
//    status_h5 = H5Pset_chunk(plist_id, RANK, chunk_dims);
//    assert(status_h5 >= 0);
//    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
//    H5Pclose(plist_id);
//    H5Sclose(filespace);
// 
//    /* 
//     * Each process defines dataset in memory and writes it to the hyperslab
//     * in the file.
//     */
//    count[0] = 1;               // one block in x-direction
//    count[1] = 1 ;              // one block in y-direction
//    stride[0] = 1;              // the blocks have no gaps between them in x-direction
//    stride[1] = 1;              // the blocks have no gaps between them in y-direction
//    block[0] = chunk_dims[0];   // contiguous x-size of block is chunk_dims[0]
//    block[1] = chunk_dims[1];   // contiguous x-size of block is chunk_dims[1]
//    offset[0] = pi*n;           // x-index [in data file] of (0,0) entry [in memory] of block is pi*n
//    offset[1] = pj*n;           // y-index [in data file] of (0,0) entry [in memory] of block is pj*n
//    // // if reordering should be required: ((N-1)-pi)*n;
//    
//    /*
//     * Select hyperslab in the file.
//     */
//    filespace = H5Dget_space(dset_id);
//    status_h5 = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
//    assert(status_h5 >= 0);
//    
//    /*
//     * Create property list for collective dataset write.
//     */
//    plist_id = H5Pcreate(H5P_DATASET_XFER);
//    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//  
//    /*
//     * Initialize data buffer 
//     */
//    data = (double *) calloc(chunk_dims[0]*chunk_dims[1],sizeof(double));
//    for (long i = 0; i < num_procs; ++i)  
//    {
//      if (myid == i)
//      {
//        for (ulong j=0; j < (ulong)chunk_dims[0]*chunk_dims[1]; j++) 
//        {
// 	 ulong ii = j/chunk_dims[0]+1; // the +1 comes from the fact that s[][] is padded
// 	 ulong jj = j%chunk_dims[1]+1;
// 	 data[j] = s[ii][jj];
// 	 //cout << myid << " " << j << " " << data[j] << " " << ii << " " << jj << " " << s[ii][jj] << " " << endl;
//        }
//        //cout << endl;
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
//    
//    /*
//     * Write to file
//     */   
//    status_h5 = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
//    assert(status_h5 >= 0);
//    free(data);
// 
//    /*
//     * Close/release resources.
//     */
//    H5Dclose(dset_id);
//    H5Sclose(filespace);
//    H5Sclose(memspace);
//    H5Pclose(plist_id);
//    H5Fclose(file_id);
//    
//    

//    /*
//     * Now we begin the read section of this example.
//     */
// 
//    /*
//     * Open file_id and dataset using the default properties.
//     */
//    plist_id = H5Pcreate(H5P_FILE_ACCESS);
//    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
//    file_id = H5Fopen (H5FILE_NAME, H5F_ACC_RDONLY, plist_id);
//    H5Pclose(plist_id);
// 
//    dset_id = H5Dopen (file_id, DATASETNAME, H5P_DEFAULT);
//    filespace = H5Dget_space (dset_id);
//     
//    plist_id = H5Pcreate(H5P_DATASET_XFER);
//    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//           
//    double *rdata = (double *) calloc(chunk_dims[0]*chunk_dims[1],sizeof(double));
//     
//    /*
//     * Define and select the hyperslab to use for reading.
//     */
//    memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 
//     
//    count[0] = 1;               // one block in x-direction
//    count[1] = 1 ;              // one block in y-direction
//    stride[0] = 1;              // the blocks have no gaps between them in x-direction
//    stride[1] = 1;              // the blocks have no gaps between them in y-direction
//    block[0] = chunk_dims[0];   // contiguous x-size of block is chunk_dims[0]
//    block[1] = chunk_dims[1];   // contiguous x-size of block is chunk_dims[1]
//    offset[0] = pi*n;           // x-index [in data file] of (0,0) entry [in memory] of block is pi*n
//    offset[1] = pj*n;           // y-index [in data file] of (0,0) entry [in memory] of block is pj*n
//    // // if reordering should be required: ((N-1)-pi)*n;
//   
//    status_h5 = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, stride, count, block);
//    assert(status_h5 >= 0);
// 
//    /*
//     * Read the data using the previously defined hyperslab.
//     */
//    status_h5 = H5Dread (dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, rdata);
//    assert(status_h5 >= 0);
//  
//    for (long j=0; j<num_procs; ++j)
//    {
//      if (myid == j)
//      {
//        /*
// 	* Output the data to the screen.
// 	*/
//        printf ("\nData as written to disk by hyberslabs:\n");
//        for (ulong i=0; i < (ulong)chunk_dims[0]*chunk_dims[1]; i++)
//        {
// 	 cout << myid << " " << i << " " << rdata[i] << " " << endl;
//        }
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
//  
//    /*
//     * Close and release resources.
//     */
//    status_h5 = H5Pclose (plist_id);  assert(status_h5 >= 0);
//    status_h5 = H5Dclose (dset_id);   assert(status_h5 >= 0);
//    status_h5 = H5Fclose (file_id);   assert(status_h5 >= 0);
//    status_h5 = H5Sclose (filespace); assert(status_h5 >= 0);
//    status_h5 = H5Sclose (memspace);  assert(status_h5 >= 0);
//     
//    free(rdata);
//
//    if (myid == 0) { cout << " ### end HDF5 test.\n\n"; }
//    MPI_Finalize();
//    exit(0);
// 
//    // ===============================================================
//    // ======================= end HDF5 test =========================
//    // ===============================================================
// #endif //H5
   
   double t3 = MPI_Wtime()-t2 -t1;
   if (myid == 0) { cout << " # time for read-in and distribution was " << t2 << " and " << t3 << " seconds.\n\n"; }
         
   /* 1. Set up a grid */
   {
      /* Create an empty 2D grid object */
      HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

      /* Add a new box to the grid */
      HYPRE_StructGridSetExtents(grid, ilower, iupper);
      
      /* Make grid periodic in the second dimension, with period N*n */
      long long int periodic[2] = {0,N*n};//CM-Fix
      HYPRE_StructGridSetPeriodic (grid, periodic); 

      /* This is a collective call finalizing the grid assembly.
         The grid is now ``ready to be used'' */
      HYPRE_StructGridAssemble(grid);
   }

   /* 2. Define the discretization stencil */
   /* Define the geometry of the stencil */
   long long int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};//CM-Fix

   /* Create an empty 2D, 5-pt stencil object */
   HYPRE_StructStencilCreate(2, 5, &stencil);

   /* Assign stencil entries */
   for (long i = 0; i < 5; i++)
      HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
   
   /* 3. Set up Struct Vectors for b and x */
   {
      double *values = (double *) calloc((n*n), sizeof(double));

      /* Create an empty vector object */
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
      
      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      
      k = 0;
      /* Set the values of b in left-to-right, bottom-to-top order */      
      for (long j = 0; j < n; j++)
         for (long i = 0; i < n; i++, k++)
            values[k] = h2 * Eval(F,i,j);
      HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

      free(values);

      /* Assembling is postponed since the vectors will be further modified */
   }

   /* 3.5 Set up Struct Vector for x */      
   {
      double *values = (double *) calloc((n*n), sizeof(double));
      
      /* Create an empty vector object */
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(x);
         
      MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
      t1 = MPI_Wtime();
      if (!readIn.empty())
      {

#ifdef H5
	ifstream readPress ( readIn.c_str() );
	if(!readPress)
	{ 
	  if (myid == 0) { cerr << "\n ##### Error opening file " << readIn << ". Terminating...\n\n"; }
	  free(values); 
	  MPI_Finalize();
	  exit(1);
	}
	else { readPress.close(); } // file exists, read it with HDF5
// 	readH5(readIn.c_str(), pi, pj, num_procs, myid, values);
	
	double *pressReadin = (double *) calloc((n*n), sizeof(double));
	readH5(readIn.c_str(), pi, pj, num_procs, myid, pressReadin);
	
	/* transpose the read-in values */
	for (long ii = 0; ii < n; ++ii)
	  for (long jj = 0; jj < n; ++jj)
	    values[jj*n+ii] = pressReadin[ii*n+jj];
	free(pressReadin);
	
#else // H5
	
// 	///* everyone reads its own file --- only possible if the restart run has the same # of procs *///
// 	{
// 	  double x, y, dummy;
// 	  string restOfLine;
// 	  string file = "test2";
// 	  
// 	  char buffer[50];
// 	  string filename; filename.clear();
// 	  filename.append(file);
// 	  sprintf(buffer,"_%05d.press.%06d.dat",N*n,myid);
// 	  filename.append(buffer);  
// 	  ifstream readPress ( filename.c_str() ); readPress.precision(8);    
// 	  if(!readPress){ cerr << "\n ##### Error opening file " << filename << ". Terminating...\n"; exit(1); }  
// 	  else          { cout << "\n ### Proc " << myid << " reading file " << filename << "...\n\n"; }  
// 	
// 	  k = 0;
// 	  /* Set x to old pressure */
// 	  for (long j = 0; j < n; ++j)
// 	  {
// 	    for (long i = 0; i < n; ++i, k++)
// 	    {
// 	      readPress >> x >> y >> dummy; getline(readPress,restOfLine);
// 	      values[k] = dummy;
// 	    #ifdef VERBOSE
// 	      cout.precision(8);
// 	      cout.setf( std::ios::scientific ); //cout.setf(ios::fixed);
// 	      // cout.setf(ios::fixed,ios::floatfield); // cout.width(10);
// 	      cout << " x: " << x << " \t" << "y: " << y << " \t"  
// 		    << " i: " << i << " \t" << "j: " << j << " \t"
// 		    << " dummy: " << dummy << " \t" << "pressReadIn: " << values[k] << endl;
// 	    #endif // VERBOSE
// 	    }
// 	  }
// 	  readPress.close();
// 	}
	///* proc 0 reads its global file and then scatters the data across all other procs *///
	{
	  double *pressReadin = NULL;
          if (myid == 0) 
	  { 
	    pressReadin = (double *) calloc((N*n*N*n), sizeof(double));
          
	    double dummy;
	    string restOfLine;
	    ifstream readPress ( readIn.c_str() ); readPress.precision(8);    
	    if(!readPress)
	    { 
	      cerr << "\n ##### Error opening file " << readIn << ". Terminating...\n";
	      free(pressReadin); free(values); 
	      MPI_Finalize();
	      exit(1);
	    }  
	    else          
	      cout << "\n ### reading file " << readIn << "...\n\n";
	  
	    /* Read in old pressure */
	    for (long i = 0; i < N*n; ++i)
	      for (long j = 0; j < N*n; ++j)
	      {
		k = j*N*n+i;
		for (long icol = 1; icol <= ncol; ++icol) { readPress >> dummy; }
		getline(readPress,restOfLine);
// 	      #ifdef JUQUEEN
// 		readPress >> dummy; getline(readPress,restOfLine);
// 	      #else
// 		readPress >> x >> y >> dummy; getline(readPress,restOfLine);
// 	      #endif
		pressReadin[k] = dummy;
	      }
	    readPress.close();
	  }
	  for (long ii=0; ii<NPROWS; ii++)
	    for (long jj=0; jj<NPCOLS; jj++) 
	    {       
	      /// NOTE for some reason needs the procs ordering needs to be transposed
	      disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;  /// NOTE proc 1 located ABOVE proc 0
	      counts [ii*NPCOLS+jj] = 1;
	    }
	  MPI_Scatterv(pressReadin, counts, disps, blocktype, values, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  free(pressReadin);
	}
#endif // H5
      }
      else
      {      
	/* Set x = 0 */
	for (long i = 0; i < (n*n); i++) {
	  values[i] = 0.0; // Eval(P,i%n,0); /// TODO insert linear initial guess
         //left side :Anle 
        /*
	 if (i%n == 0) {
	   values[i] = 1.;
	 }
	 else if ((i+1)%n == 0) {
	   values[i] = 0.;
	 } 
	 else{
	   //get the index of ix
	   int ix = i/n;
	   values[i] = (n-ix-1)*1./(n-1);
	}
	*/
       }	 
      }   
      HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

      free(values); 
      MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
      t2 = MPI_Wtime()-t1;
      if (myid == 0) { cout << " # time for setting up pressure was " << t2 << " seconds.\n\n"; }
   }

   /* 4. Set up a Struct Matrix */
   {
      /* Create an empty matrix object */
      HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
      
      /* Use symmetric storage? */
      HYPRE_StructMatrixSetSymmetric(A, 0);
      
      /* Indicate that the matrix coefficients are ready to be set */
      HYPRE_StructMatrixInitialize(A);

      /* Set the stencil values in the interior. Here we set the values
         at every node. We will modify the boundary nodes later. */
      long long int stencil_indices[5] = {0, 1, 2, 3, 4}; /* labels correspond //CM-Fix
                                                          to the offsets */
      double *values = (double *) calloc(5*(n*n), sizeof(double));
     
      k = 0;
      /* The order is left-to-right, bottom-to-top */
      for (long j = 0; j < n; j++)
      {
        for (long i = 0; i < n; i++, k+=5)
	{
	    int ix = i+1; // ilower[0]+i+1;
	    int iy = j+1; // ilower[1]+j+1;
            
	    //modified by Anle:
 	   // values[k+1] = 2. / (1./s[ix][iy] + 1./s[ix-1][iy]); // p[i-1][j]
 	   // values[k+2] = 2. / (1./s[ix][iy] + 1./s[ix+1][iy]); // p[i+1][j]
 	   // values[k+3] = 2. / (1./s[ix][iy] + 1./s[ix][iy-1]); // p[i][j-1]
 	   // values[k+4] = 2. / (1./s[ix][iy] + 1./s[ix][iy+1]); // p[i][j+1]
 	   // values[k+0] = - (values[k+1] + values[k+2] + values[k+3] + values[k+4]);
	    
	    if (s[ix][iy] > 0.0) // point open
	    { // make sure that 1/(1/0) does not give undefined result. Also, treat points surrounded by zero conductivity
	      values[k+1] = (s[ix-1][iy] > 0.0) ? 2. / (1./s[ix][iy] + 1./s[ix-1][iy]) : 0.0; // p[i-1][j]
	      values[k+2] = (s[ix+1][iy] > 0.0) ? 2. / (1./s[ix][iy] + 1./s[ix+1][iy]) : 0.0; // p[i+1][j]
	      values[k+3] = (s[ix][iy-1] > 0.0) ? 2. / (1./s[ix][iy] + 1./s[ix][iy-1]) : 0.0; // p[i][j-1]
	      values[k+4] = (s[ix][iy+1] > 0.0) ? 2. / (1./s[ix][iy] + 1./s[ix][iy+1]) : 0.0; // p[i][j+1]
	      if ((values[k+1] == 0) && (values[k+2] == 0) && (values[k+3] == 0) && (values[k+4] == 0)) { values[k+0] = 1.0; } 
	      else { values[k+0] = - (values[k+1] + values[k+2] + values[k+3] + values[k+4]); }
	    }
	    else                 // point closed --- identity matrix, TODO eliminate the adjacent points' reaching in
	    {
	      values[k+1] = values[k+2] = values[k+3] = values[k+4] = 0.0;
	      values[k+0] = 1;
	    }
	}
      }
      
      HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 5,
                                     stencil_indices, values);
      free(values);      
   }    

   /* 5. Set the boundary conditions, while eliminating the coefficients
         reaching ouside of the domain boundary. We must modify the matrix
         stencil and the corresponding rhs entries. */
   {
      long long int bc_ilower[2];//CM-Fix
      long long int bc_iupper[2];//CM-Fix

      long long int stencil_indices[5] = {0, 1, 2, 3, 4};//CM-Fix
      double *values, *bvalues;

      long long int nentries;//CM-Fix
      nentries = 5;
      
      values  = (double *) calloc(nentries*n, sizeof(double));
      bvalues = (double *) calloc(n, sizeof(double));

      /* The stencil at the boundary nodes is 1-0-0-0-0. Because
         we have I x_b = u_0; */
      for (long i = 0; i < nentries*n; i += nentries)
      {
         values[i] = 1.0;
         for (long j = 1; j < nentries; j++)
            values[i+j] = 0.0;
      }

      /* Processors at x = 0 */
      if (pi == 0)
      {
         bc_ilower[0] = pi*n;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (long j = 0; j < n; j++) 
           bvalues[j] = bcEval(U0,0,j);
         
         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 1 */
      if (pi == N-1)
      {
         bc_ilower[0] = pi*n + n-1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         /* Modify the matrix */
         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, values);

         /* Put the boundary conditions in b */
         for (long j = 0; j < n; j++) 
	   bvalues[j] = bcEval(U0,0,j);
         
         HYPRE_StructVectorSetBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Recall that the system we are solving is:
         [A_ii 0; 0 I] [x_i ; x_b] = [b_i - A_ib u_0; u_0].
         This requires removing the connections between the interior
         and boundary nodes that we have set up when we set the
         5pt stencil at each node. We adjust for removing
         these connections by appropriately modifying the rhs.
         For the symm ordering scheme, just do the top and right
         boundary */

      /* Processors at x = 0, neighbors of boundary nodes */
      if (pi == 0)
      {
         bc_ilower[0] = pi*n + 1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         stencil_indices[0] = 1;

         /* Modify the matrix */
         for (long j = 0; j < n; j++)
            bvalues[j] = 0.0;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                        stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (long j = 0; j < n; j++)
	 {
	   int ix =   2;
	   int iy = j+1;
	   bvalues[j] = bcEval(U0,-1,j) * (-2. / (1./s[ix][iy] + 1./s[ix-1][iy])); 
	 }

         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      /* Processors at x = 1, neighbors of boundary nodes */
      if (pi == N-1)
      {
         bc_ilower[0] = pi*n + (n-1) - 1;
         bc_ilower[1] = pj*n;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + n-1;

         stencil_indices[0] = 2;
         
         /* Modify the matrix */
         for (long j = 0; j < n; j++)
            bvalues[j] = 0.0;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, 1,
                                        stencil_indices, bvalues);

         /* Eliminate the boundary conditions in b */
         for (long j = 0; j < n; j++)
	 {
	    int ix = n-1;
	    int iy = j+1;
	    bvalues[j] = bcEval(U0,1,j) * (-2. / (1./s[ix][iy] + 1./s[ix+1][iy]));
	 }

         HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, bvalues);
      }

      free(values);
      free(bvalues);
   }

   /* Finalize the vector and matrix assembly */
   HYPRE_StructMatrixAssemble(A);
   HYPRE_StructVectorAssemble(b);
   HYPRE_StructVectorAssemble(x);
    
   //change here: 
//    char filename[255];
//    sprintf(filename, "%s.%06d.%s", "_matrixA", myid, "txt");
//    HYPRE_StructMatrixPrint (filename, A, 1);
//    sprintf(filename, "%s.%06d.%s", "_vectorX", myid, "txt");
//    HYPRE_StructVectorPrint (filename, x, 1);
//    sprintf(filename, "%s.%06d.%s", "_vectorB", myid, "txt");
//    HYPRE_StructVectorPrint (filename, b, 1);

   num_iterations = iter;
   totalIterations = 0;
   
   double maxTimeIter = 0, maxTimeOutput = 0, maxTimeLocal = 0;
      
   while (1)
   {
     double startTimeIter = MPI_Wtime();
     if (iter) // did the user specify a nonzero number of iterations?
     {
      if (its_alt)
      {
	/// each set has "iter" iterations.
	/// every "its_alt" sets of iterations with "solver_id_ini", insert one set with "alternate_solver" 
	if ((solver_id == solver_id_ini) && ((totalIterations%(its_alt*iter))==0) && totalIterations) 
	{ solver_id = alternate_solver; }
	else                 
	{ solver_id = solver_id_ini; }
	
	if (myid == 0) { cout << " ### running with solver " << solver_id << "...\n\n"; }
      }

      /* 6. Set up and use a solver */
      if (solver_id == 0) /* SMG */
      {
	  /* Options and setup */
	  HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
	  HYPRE_StructSMGSetMemoryUse(solver, 0);
	  HYPRE_StructSMGSetMaxIter(solver, iter);
	  HYPRE_StructSMGSetTol(solver, TOL ); // NOTE was 1.0e-06
	  HYPRE_StructSMGSetRelChange(solver, 0);
	  HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
	  HYPRE_StructSMGSetNumPostRelax(solver, n_post);
	  HYPRE_StructSMGSetPrintLevel(solver, 1);
	  HYPRE_StructSMGSetLogging(solver, 1);
	  HYPRE_StructSMGSetup(solver, A, b, x);

	  /* Solve */
	  HYPRE_StructSMGSolve(solver, A, b, x);

	  /* Get info and release memory */
	  HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
	  HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
	  HYPRE_StructSMGDestroy(solver);
      }

      if (solver_id == 1) /* PFMG */
      {
	  /* Options and setup */
	  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
	  HYPRE_StructPFMGSetMaxIter(solver, iter);
	  HYPRE_StructPFMGSetTol(solver, TOL ); // NOTE was 1.0e-06
	  HYPRE_StructPFMGSetRelChange(solver, 0);
	  HYPRE_StructPFMGSetRAPType(solver, rap);
	  HYPRE_StructPFMGSetRelaxType(solver, relax);
	  HYPRE_StructPFMGSetNumPreRelax(solver, n_pre);
	  HYPRE_StructPFMGSetNumPostRelax(solver, n_post);
	  HYPRE_StructPFMGSetSkipRelax(solver, skip);
	  HYPRE_StructPFMGSetPrintLevel(solver, 1);
	  HYPRE_StructPFMGSetLogging(solver, 1);
	  HYPRE_StructPFMGSetup(solver, A, b, x);

	  /* Solve */
	  HYPRE_StructPFMGSolve(solver, A, b, x);

	  /* Get info and release memory */
	  HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
	  HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
	  HYPRE_StructPFMGDestroy(solver);
      }

      /* Preconditioned CG */
      if ((solver_id > 9) && (solver_id < 20))
      {
	  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
	  HYPRE_StructPCGSetMaxIter(solver, iter );   // NOTE was 200
	  HYPRE_StructPCGSetTol(solver, TOL );   // NOTE was 1.0e-06
	  HYPRE_StructPCGSetTwoNorm(solver, 1 );
	  HYPRE_StructPCGSetRelChange(solver, 0 );
	  HYPRE_StructPCGSetPrintLevel(solver, 2 );

	  if (solver_id == 10)
	  {
	    /* use symmetric SMG as preconditioner */
	    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructSMGSetMemoryUse(precond, 0);
	    HYPRE_StructSMGSetMaxIter(precond, 1);
	    HYPRE_StructSMGSetTol(precond, 0.0);
	    HYPRE_StructSMGSetZeroGuess(precond);
	    HYPRE_StructSMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructSMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructSMGSetPrintLevel(precond, 0);
	    HYPRE_StructSMGSetLogging(precond, 0);
	    HYPRE_StructPCGSetPrecond(solver,
				      HYPRE_StructSMGSolve,
				      HYPRE_StructSMGSetup,
				      precond);
	  }

	  else if (solver_id == 11)
	  {
	    /* use symmetric PFMG as preconditioner */
	    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructPFMGSetMaxIter(precond, 1);
	    HYPRE_StructPFMGSetTol(precond, 0.0);
	    HYPRE_StructPFMGSetZeroGuess(precond);
	    HYPRE_StructPFMGSetRAPType(precond, rap);
	    HYPRE_StructPFMGSetRelaxType(precond, relax);
	    HYPRE_StructPFMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructPFMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructPFMGSetSkipRelax(precond, skip);
	    HYPRE_StructPFMGSetPrintLevel(precond, 0);
	   HYPRE_StructPFMGSetLogging(precond, 0);
	    HYPRE_StructPCGSetPrecond(solver,
				      HYPRE_StructPFMGSolve,
				      HYPRE_StructPFMGSetup,
				      precond);
	  }

	  else if (solver_id == 17)
	  {
	    /* use two-step Jacobi as preconditioner */
	    HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructJacobiSetMaxIter(precond, 2);
	    HYPRE_StructJacobiSetTol(precond, 0.0);
	    HYPRE_StructJacobiSetZeroGuess(precond);
	    HYPRE_StructPCGSetPrecond( solver,
					HYPRE_StructJacobiSolve,
					HYPRE_StructJacobiSetup,
					precond);
	  }

	  else if (solver_id == 18)
	  {
	    /* use diagonal scaling as preconditioner */
	    precond = NULL;
	    HYPRE_StructPCGSetPrecond(solver,
				      HYPRE_StructDiagScale,
				      HYPRE_StructDiagScaleSetup,
				      precond);
	  }

	  /* PCG Setup */
	  HYPRE_StructPCGSetup(solver, A, b, x );

	  /* PCG Solve */
	  HYPRE_StructPCGSolve(solver, A, b, x);

	  /* Get info and release memory */
	  HYPRE_StructPCGGetNumIterations( solver, &num_iterations );
	  HYPRE_StructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );
	  HYPRE_StructPCGDestroy(solver);

	  if (solver_id == 10)
	  {
	    HYPRE_StructSMGDestroy(precond);
	  }
	  else if (solver_id == 11 )
	  {
	    HYPRE_StructPFMGDestroy(precond);
	  }
	  else if (solver_id == 17)
	  {
	    HYPRE_StructJacobiDestroy(precond);
	  }
      }

      /* Preconditioned GMRES */
      if ((solver_id > 29) && (solver_id < 40))
      {
	  HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);

	  /* Note that GMRES can be used with all the interfaces - not
	    just the struct.  So here we demonstrate the
	    more generic GMRES interface functions. Since we have chosen
	    a struct solver then we must type cast to the more generic
	    HYPRE_Solver when setting options with these generic functions.
	    Note that one could declare the solver to be
	    type HYPRE_Solver, and then the casting would not be necessary.*/
	  
	  HYPRE_GMRESSetMaxIter((HYPRE_Solver) solver, iter ); // NOTE was 500
	  HYPRE_GMRESSetKDim((HYPRE_Solver) solver,30);
	  HYPRE_GMRESSetTol((HYPRE_Solver) solver, TOL );   // NOTE was 1.0e-06
	  HYPRE_GMRESSetPrintLevel((HYPRE_Solver) solver, 2 );
	  HYPRE_GMRESSetLogging((HYPRE_Solver) solver, 1 );

	  if (solver_id == 30)
	  {
	    /* use symmetric SMG as preconditioner */
	    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructSMGSetMemoryUse(precond, 0);
	    HYPRE_StructSMGSetMaxIter(precond, 1);
	    HYPRE_StructSMGSetTol(precond, 0.0);
	    HYPRE_StructSMGSetZeroGuess(precond);
	    HYPRE_StructSMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructSMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructSMGSetPrintLevel(precond, 0);
	    HYPRE_StructSMGSetLogging(precond, 0);
	    HYPRE_StructGMRESSetPrecond(solver,
					HYPRE_StructSMGSolve,
					HYPRE_StructSMGSetup,
					precond);
	  }

	  else if (solver_id == 31)
	  {
	    /* use symmetric PFMG as preconditioner */
	    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructPFMGSetMaxIter(precond, 1);
	    HYPRE_StructPFMGSetTol(precond, 0.0);
	    HYPRE_StructPFMGSetZeroGuess(precond);
	    HYPRE_StructPFMGSetRAPType(precond, rap);
	    HYPRE_StructPFMGSetRelaxType(precond, relax);
	    HYPRE_StructPFMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructPFMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructPFMGSetSkipRelax(precond, skip);
	    HYPRE_StructPFMGSetPrintLevel(precond, 0);
	    HYPRE_StructPFMGSetLogging(precond, 0);
	    HYPRE_StructGMRESSetPrecond( solver,
					  HYPRE_StructPFMGSolve,
					  HYPRE_StructPFMGSetup,
					  precond);
	  }

	  else if (solver_id == 37)
	  {
	    /* use two-step Jacobi as preconditioner */
	    HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructJacobiSetMaxIter(precond, 2);
	    HYPRE_StructJacobiSetTol(precond, 0.0);
	    HYPRE_StructJacobiSetZeroGuess(precond);
	    HYPRE_StructGMRESSetPrecond( solver,
					  HYPRE_StructJacobiSolve,
					  HYPRE_StructJacobiSetup,
					  precond);
	  }

	  else if (solver_id == 38)
	  {
	    /* use diagonal scaling as preconditioner */
	    precond = NULL;
	    HYPRE_StructGMRESSetPrecond( solver,
					HYPRE_StructDiagScale,
					HYPRE_StructDiagScaleSetup,
					precond);
	  }

	  /* GMRES Setup */
	  HYPRE_StructGMRESSetup(solver, A, b, x );

	  /* GMRES Solve */
	  HYPRE_StructGMRESSolve(solver, A, b, x);

	  /* Get info and release memory */
	  HYPRE_StructGMRESGetNumIterations(solver, &num_iterations);
	  HYPRE_StructGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
	  HYPRE_StructGMRESDestroy(solver);

	  if (solver_id == 30)
	  {
	    HYPRE_StructSMGDestroy(precond);
	  }
	  else if (solver_id == 31)
	  {
	    HYPRE_StructPFMGDestroy(precond);
	  }
	  else if (solver_id == 37)
	  {
	    HYPRE_StructJacobiDestroy(precond);
	  }
      }
      
      /* preconditioned BiCGSTAB */
      if ((solver_id > 39) && (solver_id < 50))
      { 
	  HYPRE_StructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
	  HYPRE_StructBiCGSTABSetMaxIter(solver, iter );   // NOTE was 200
	  HYPRE_StructBiCGSTABSetTol(solver, TOL );   // NOTE was 1.0e-06
	  HYPRE_StructBiCGSTABSetPrintLevel(solver, 2 );

	  if (solver_id == 40)
	  {
	    /* use symmetric SMG as preconditioner */
	    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructSMGSetMemoryUse(precond, 0);
	    HYPRE_StructSMGSetMaxIter(precond, 1);
	    HYPRE_StructSMGSetTol(precond, 0.0);
	    HYPRE_StructSMGSetZeroGuess(precond);
	    HYPRE_StructSMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructSMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructSMGSetPrintLevel(precond, 0);
	    HYPRE_StructSMGSetLogging(precond, 0);
	    HYPRE_StructBiCGSTABSetPrecond(solver,
				      HYPRE_StructSMGSolve,
				      HYPRE_StructSMGSetup,
				      precond);
	  }

	  else if (solver_id == 41)
	  {
	    /* use symmetric PFMG as preconditioner */
	    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructPFMGSetMaxIter(precond, 1);
	    HYPRE_StructPFMGSetTol(precond, 0.0);
	    HYPRE_StructPFMGSetZeroGuess(precond);
	    HYPRE_StructPFMGSetRAPType(precond, rap);
	    HYPRE_StructPFMGSetRelaxType(precond, relax);
	    HYPRE_StructPFMGSetNumPreRelax(precond, n_pre);
	    HYPRE_StructPFMGSetNumPostRelax(precond, n_post);
	    HYPRE_StructPFMGSetSkipRelax(precond, skip);
	    HYPRE_StructPFMGSetPrintLevel(precond, 0);
	    HYPRE_StructPFMGSetLogging(precond, 0);
	    HYPRE_StructBiCGSTABSetPrecond(solver,
				      HYPRE_StructPFMGSolve,
				      HYPRE_StructPFMGSetup,
				      precond);
	  }

	  else if (solver_id == 47)
	  {
	    /* use two-step Jacobi as preconditioner */
	    HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &precond);
	    HYPRE_StructJacobiSetMaxIter(precond, 2);
	    HYPRE_StructJacobiSetTol(precond, 0.0);
	    HYPRE_StructJacobiSetZeroGuess(precond);
	    HYPRE_StructBiCGSTABSetPrecond( solver,
					HYPRE_StructJacobiSolve,
					HYPRE_StructJacobiSetup,
					precond);
	  }

	  else if (solver_id == 48)
	  {
	    /* use diagonal scaling as preconditioner */
	    precond = NULL;
	    HYPRE_StructBiCGSTABSetPrecond(solver,
				      HYPRE_StructDiagScale,
				      HYPRE_StructDiagScaleSetup,
				      precond);
	  }

	  /* PCG Setup */
	  HYPRE_StructBiCGSTABSetup(solver, A, b, x );

	  /* PCG Solve */
	  HYPRE_StructBiCGSTABSolve(solver, A, b, x);

	  /* Get info and release memory */
	  HYPRE_StructBiCGSTABGetNumIterations( solver, &num_iterations );
	  HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm( solver, &final_res_norm );
	  HYPRE_StructBiCGSTABDestroy(solver);
		
	  if (solver_id == 40)
	  {
	    HYPRE_StructSMGDestroy(precond);
	  }
	  else if (solver_id == 41 )
	  {
	    HYPRE_StructPFMGDestroy(precond);
	  }
	  else if (solver_id == 47)
	  {
	    HYPRE_StructJacobiDestroy(precond);
	  }     
      }

      if ((solver_id != 0)  && // SMG
	  (solver_id != 1)  && // PFMG
	  (solver_id != 10) && // CG + SMG preconditioner
	  (solver_id != 11) && // CG + PFMG preconditioner
	  (solver_id != 17) && // CG + Jacobi preconditioner
	  (solver_id != 18) && // CG + diagonal scaling
	  (solver_id != 19) && // CG
	  (solver_id != 30) && // GMRES + SMG preconditioner
	  (solver_id != 31) && // GMRES + PFMG preconditioner
	  (solver_id != 37) && // GMRES + Jacobi preconditioner
	  (solver_id != 38) && // GMRES + diagonal scaling
	  (solver_id != 39) && // GMRES
	  (solver_id != 40) && // BiCGSTAB + SMG preconditioner
	  (solver_id != 41) && // BiCGSTAB + PFMG preconditioner
	  (solver_id != 47) && // BiCGSTAB + Jacobi preconditioner
	  (solver_id != 48) && // BiCGSTAB + diagonal scaling
	  (solver_id != 49))   // BiCGSTAB
      { cerr << "\n\n##### ERROR: no valid solver selected! Terminating...\n\n"; MPI_Finalize(); exit(1);}
      
//       /* calculate the flow for the current configuration */
/// 	// need one additional layer in both directions
      
      totalIterations += num_iterations;      
      //if ((totalIterations % (iter)) == 0)                    { vis = 1; } // output after every set
    #ifndef JUQUEEN
      if (totalIterations > iter_output) { iter_output +=1024; vis = 1; } // output every so often, but not at first, and not on JUQUEEN
    #endif
      //if ((final_res_norm < TOL*3) && (num_iterations != iter) ) { vis = 2; } // output when tolerance reached
      //if ((final_res_norm < TOL*3) && (num_iterations != iter) && (final_res_norm > 0.0)) { vis = 2; } // output when tolerance reached
       
     } // end if (iter) // check if nonzero iterations were requested
      
     //modified here:

      /* 7. postprocessing, communication and output */     
      {
	double *values = (double *) calloc(n*n, sizeof(double));

	vector < vector <double> > current_x;
	current_x.resize(n+1);
	for (long ii = 0; ii < n+1; ++ii) { current_x[ii].resize(n,0.0); }
	
	vector < vector <double> > current_y;
	current_y.resize(n);
	for (long ii = 0; ii < n; ++ii)   { current_y[ii].resize(n+1,0.0); }

	vector < vector <double> > pressLoc;
	pressLoc.resize(n+2);
	for (long ii = 0; ii < n+2; ++ii) { pressLoc[ii].resize(n+2,0.0); }
	
      #ifdef LOCAL
        double startTimeLocal = MPI_Wtime();
        iIterLocal = 0;
	while (1)
	{
      #endif //LOCAL
	
	/* get the local solution */
	HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);

	k = 0;
	for (long j = 0; j < n; j++) // NOTE pressure needs to be output transposed
	  for (long i = 0; i < n; i++) 
	    pressLoc[i+1][j+1] = max(0.0,min(values[k++],1.0)); /// NOTE safeguard: make sure that the pressure is <= 1
	
	long long int bc_ilower[2];//CM-Fix
	long long int bc_iupper[2];//CM-Fix
	
	///* send top row to top neighbor's bottom row */
	bc_ilower[0] = pi*n;         bc_iupper[0] = bc_ilower[0] + n-1;
	bc_ilower[1] = pj*n + (n-1); bc_iupper[1] = bc_ilower[1];
	HYPRE_StructVectorGetBoxValues(x, bc_ilower, bc_iupper, sendbuf);
	for (long ii = 0; ii < n; ++ii) { recvbuf[ii] = 0.0; }
	MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, topNeigh, 4, 
		      recvbuf, n, MPI_DOUBLE, botNeigh, 4, MPI_COMM_WORLD, &status );
	for (long ii = 1; ii < n+1; ++ii) { pressLoc[ii][0] = recvbuf[ii-1]; }
	  
	///* send bottom row to bottom neighbor's top row */
	bc_ilower[0] = pi*n;            bc_iupper[0] = bc_ilower[0] + n-1;
	bc_ilower[1] = pj*n;            bc_iupper[1] = bc_ilower[1];
	HYPRE_StructVectorGetBoxValues(x, bc_ilower, bc_iupper, sendbuf);
	for (long ii = 0; ii < n; ++ii) { recvbuf[ii] = 0.0; }
	MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, botNeigh, 3, 
		      recvbuf, n, MPI_DOUBLE, topNeigh, 3, MPI_COMM_WORLD, &status );  
	for (long ii = 1; ii < n+1; ++ii) { pressLoc[ii][n+1] = recvbuf[ii-1]; }
	    
	///* send left row to left neighbor's right row */
	bc_ilower[0] = pi*n;            bc_iupper[0] = bc_ilower[0];
	bc_ilower[1] = pj*n;            bc_iupper[1] = bc_ilower[1] + n-1;
	HYPRE_StructVectorGetBoxValues(x, bc_ilower, bc_iupper, sendbuf);
	for (long jj = 0; jj < n; ++jj) { recvbuf[jj] = 0.0; }
	MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, leftNeigh, 1, 
		      recvbuf, n, MPI_DOUBLE, rightNeigh, 1, MPI_COMM_WORLD, &status );  
	/* right boundary, x = 1 */
	for (long jj = 1; jj < n+1; ++jj) { pressLoc[n+1][jj] = recvbuf[jj-1]; }
	if (pi == N-1) { for (long jj = 0; jj < n; ++jj) { pressLoc[n+1][jj+1] = pressLoc[n][jj+1]; } }
	    
	///* send right row to right neighbor's left row */
	bc_ilower[0] = pi*n + n-1;      bc_iupper[0] = bc_ilower[0];
	bc_ilower[1] = pj*n;            bc_iupper[1] = bc_ilower[1] + n-1;
	HYPRE_StructVectorGetBoxValues(x, bc_ilower, bc_iupper, sendbuf);
	for (long jj = 0; jj < n; ++jj) { recvbuf[jj] = 0.0; }
	MPI_Sendrecv( sendbuf, n, MPI_DOUBLE, rightNeigh, 2, 
		      recvbuf, n, MPI_DOUBLE, leftNeigh, 2, MPI_COMM_WORLD, &status );  
	/* left boundary, x = 0 */
	for (long jj = 1; jj < n+1; ++jj) { pressLoc[0][jj] = recvbuf[jj-1]; }
	if (pi == 0) { for (long jj = 0; jj < n; ++jj) { pressLoc[0][jj+1] = pressLoc[1][jj+1]; } }
	
	/* calculate current */
	for (long j = 0; j < n; j++)
	{
	  for (long i = 0; i < n+1; i++)
	  {
	    int ix = i+1; /// s is the _extended_ susceptibility field, the edges padded by boundary condition entries
	    int iy = j+1; /// pressLoc is the _extended_ solution == pressure field, edges also padded
	    
	    current_x[i][j] = (pressLoc[ix][iy]-pressLoc[ix-1][iy]) * 2. / (1./s[ix][iy] + 1./s[ix-1][iy]);
// 	    cout << " i: " << i << " j: " << j 
// 		<< " p["<<ix<<"]["<<iy<<"]: " << pressLoc[ix][iy] << " " << pressLoc[ix-1][iy]
// 		<< " c: " << current_x[i][j]
// 		<< endl;
	  }
	  if (pi == 0)   { current_x[0][j] = current_x[1][j]; }   // here, p=const==1, and thus c_y == 0 by definition, and c_x == const
	  if (pi == N-1) { current_x[n][j] = current_x[n-1][j]; } // here, p=const==0, and thus c_y == 0 by definition, and c_x == const
	}
	  
	cout.precision(6); cout.setf( std::ios::scientific );
	double KahanCorr = 0.0, corrected_next_term = 0.0, new_sum = 0.0, sum = 0.0;
	sum = (pi == 0) ? -current_x[1][0] : -current_x[0][0];
	
	// calculate total current
	currentTotLoc = 0.0;
	for (long i = 0; i < n; ++i)
	{
	  if ((pi == 0) && (i == 0)) { continue; } // cx[0][j] is == cx[1][j]
	  for (long j = 0; j < n; ++j)  
	  {
// 	      currentTotLoc -= current_x[i][j];        // calculate overall total
	    if ( ((i == 1) && (j == 0)) && (pi == 0) ) { continue; } // first term is already in "sum"
	    if ( ((i == 0) && (j == 0)) && (pi != 0) ) { continue; } // dito, for more than one proc
	    corrected_next_term = -current_x[i][j] - KahanCorr;
	    new_sum             = sum + corrected_next_term;
	    KahanCorr           = (new_sum - sum) - corrected_next_term;
	    sum                 = new_sum;           // calculate overall total using Kahan Summation
	  }
	}
// 	  cout << " ### proc " << myid << " 1currentTotLoc = " << currentTotLoc << " sum = " << sum << endl;
	currentTotLoc = sum;
	
	/* add total current from all processors */
	//MPI_Reduce(&currentTotLoc, &currentTotAll, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Allreduce (&currentTotLoc, &currentTotAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	currentTotAll /= (N*n-1);
	
	for (long j = 0; j < n+1; j++)
	{
	  for (long i = 0; i < n; i++)
	  {
	    int ix = i+1; /// s is the _extended_ susceptibility field, the edges padded by boundary condition entries
	    int iy = j+1; /// pressLoc is the _extended_ solution == pressure field, edges also padded
	    
	    current_y[i][j] = (pressLoc[ix][iy]-pressLoc[ix][iy-1]) * 2. / (1./s[ix][iy] + 1./s[ix][iy-1]);
	    if (pi == 0)   { current_y[0][j]   = 0.0; } // here, p=const==1, and thus c_y == 0 by definition
	    if (pi == N-1) { current_y[n-1][j] = 0.0; } // here, p=const==0, and thus c_y == 0 by definition
	  }
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	chiSquTotLoc = 0.0;
	for (long i = 0; i < n; ++i)
	  for (long j = 0; j < n; ++j)
	  {
	    chi_dummy = current_x[i][j]-current_x[i+1][j]
		      + current_y[i][j]-current_y[i][j+1];
	    chiSquTotLoc += chi_dummy*chi_dummy;
	  }
	MPI_Allreduce (&chiSquTotLoc, &chiSquTotAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	chiSquNow = chiSquTotAll;

// 	cout << " i j   pressLoc        s           s_x        current_x     current_y      chi " << endl;
// 	cout << "===============================================================================" << endl;
// 	for (long i = 0; i < n; ++i)
// 	  for (long j = 0; j < n; ++j)
// 	    cout << " " << i << " " << j << " " 
// 	         << pressLoc[i+1][j+1] << " "
// 		 << s[i+1][j+1] << " "
// 		 << 2./(1./s[i+1][j+1] + 1./s[i+1-1][j+1]) << " " 
// 		 << -current_x[i][j] << " " << -current_y[i][j] << " "  
// 		 //<< -current_x[i+1][j] << " " << current_y[i][j] << " " << -current_y[i][j+1] << " " 
// 		 << -current_x[i][j]+current_x[i+1][j] - current_y[i][j]+current_y[i][j+1] << " "
// 		 << endl;
	//MPI_Finalize(); exit(0);
	
      #ifdef LOCAL
	  //if (1) { break; } // never go into local solver
	  //if (0) { break; } // always go into local solver
	  if ( (totalIterations < iter_local) || (localDone) ) { break; } // go into local solver every so often 
	  else // LOCAL SOLVER
	  {
	    deltaPress = 1.0/(N*n);
	    
	    ++iIterLocal;
	    if (iIterLocal == iterLocal) { localDone = 1; iter_local += 4096; }
	    
	    if (myid == 0)
	    cout << "  " 
	         << iIterLocal << "  "
		 << chiSquNow << "  " << currentTotAll << "  " 
		 << chiSquNow/(currentTotAll*currentTotAll) << "  "
		 << "_LOCAL_" << "  "
		 << endl;
		
	    //updateCurrentGlobal(); // done above through MPI communication, through array "values"
	    //updateChiSquTot();     // done above
	    
	    for (long i = 0; i < n; ++i)
	    {
	      for (long j = 0; j < n; ++j)
	      {
		int ix = i+1; /// s is the _extended_ susceptibility field, the edges padded by boundary condition entries
		int iy = j+1; /// pressLoc is the _extended_ solution == pressure field, edges also padded
		
		if ((pi == 0)   && (i == 0))   { pressLoc[1][iy] = 1.0; continue; } // here, p=const==1
		if ((pi == N-1) && (i == n-1)) { pressLoc[n][iy] = 0.0; continue; } // here, p=const==0
		
	      #ifdef CLUSTER
		if (s[ix][iy] == 0) { pressLoc[ix][iy] = 0.0; continue; } // NOTE: Neumann BC
	      #endif
	      
		for (long m = 0; m < 3; ++m) // three "function calls" without having to communicate data
		{
		  switch(m) 
		  {
		    case 0: { chi = &chiSquNow; dP = 0.0;         break; }
		    case 1: { chi = &chiSquPlu; dP = deltaPress;  break; }
		    case 2: { chi = &chiSquMin; dP = -deltaPress; break; }
		    default: { if (myid == 0) cerr << " ##### Something went wrong! Terminating...\n\n"; MPI_Finalize(); exit(0); }
		  }
		  dp = pressLoc[ix][iy] + dP;
		  *chi = ( pressLoc[ix-1][iy] - dp )*2./(1./s[ix][iy] + 1./s[ix-1][iy])  //current_x[i][j]
		       - ( dp - pressLoc[ix+1][iy] )*2./(1./s[ix][iy] + 1./s[ix+1][iy])  //current_x[i+1][j]      
		       + ( pressLoc[ix][iy-1] - dp )*2./(1./s[ix][iy] + 1./s[ix][iy-1])  //current_y[i][j] 
		       - ( dp - pressLoc[ix][iy+1] )*2./(1./s[ix][iy] + 1./s[ix][iy+1]); //current_y[i][j+1]   
		  *chi *= *chi;
		}
		//chiSquNow = updateChiSquLocal(ix, iy, 0.0);
		//chiSquPlu = updateChiSquLocal(ix, iy, deltaPress);
		//chiSquMin = updateChiSquLocal(ix, iy, -deltaPress);
		  
		double curva = chiSquPlu + chiSquMin - 2.0*chiSquNow;
		double slope = (chiSquPlu - chiSquMin) / 2.0;
		
		assert(curva >= 0.0); 
		if (curva == 0.0) { continue; }
		pressLoc[ix][iy] -= deltaPress * slope / curva;
	      } // end loop over i
	    }   // end loop over j  
// 	  MPI_Finalize(); exit(0);
	    
	    k = 0;
	    for (long j = 0; j < n; j++) 
	      for (long i = 0; i < n; i++) 
		values[k++] = max(0.0,min(pressLoc[i+1][j+1], 1.0));  /// NOTE safeguard: make sure that the pressure is 0<=p<= 1
	      
	    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
	    HYPRE_StructVectorAssemble(x);  
	    
	    //updateCurrentGlobal(); // done above through MPI communication, through array "values"
	    //updateChiSquTot();     // done above
	    
	  } // end else local solver
	  
	} // end while iIterLocal    
	if (myid == 0)
	cout << "  " 
	     << iIterLocal+1 << "  "
	     << chiSquNow << "  " << currentTotAll << "  " 
	     << chiSquNow/(currentTotAll*currentTotAll) << "  "
	     << "_LOCAL_" << "  "
	     << endl << endl;
	maxTimeLocal = max(maxTimeLocal, MPI_Wtime()-startTimeLocal);
	
	if (localDone) { vis = 1; localDone = 0; } // output after local solver has done its job
      #endif //LOCAL
	
	// calculate current column totals
	for (long i = 0; i < n; ++i)
	{
	  sendbuf[i] = 0.0; 
	  if ((pi == 0) && (i == 0)) { continue; } // cx[0][j] == cx[1][j] and does not need to be calculated twice
	  for (long j = 0; j < n; ++j) { sendbuf[i] -= current_x[i][j]; }        // calculate column totals
	}
	MPI_Gather(sendbuf, n, MPI_DOUBLE, recvbuflarge, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (myid == 0) 
	{
	  double dummy;
	  int counter;
	  
	  vector <double> curColAve (3,0.0), curColRelDiff (3,1e80), curColRelDevMean (3,0.0), curColRelDevMax (3,0.0);
	  // index 0 is all data
	  // index 1 is the data with outliers automatically removed (using quartile analysis)
	  // index 2 is the data without every n'th data point because those seem to be outliers.
	  
	  currentCol.resize(N*n, 0.0);
	  for (long ii = 0; ii < num_procs; ++ii) 
	    for (long jj = 0; jj < n; ++jj) 
	      currentCol[(ii % N)*n+jj] += recvbuflarge[ii*n+jj]; // put data from all procs together
	  
  	  for (long ii = 1; ii < N*n; ++ii) 
	    curColAve[0] += currentCol[ii];
	  curColAve[0] /= (N*n-1);
	  
	  // now get a robust average, removing outliers automatically. 
	  // WARNING doing things automatically may have unintended consequences.
	  // first sort a copy of the dataset
	  vector <double> currentColSorted = currentCol; sort(currentColSorted.begin(), currentColSorted.end());
	  // next divide the dataset into quartiles [ http://en.wikipedia.org/wiki/Quartile ]
	  // median is the mean between N*n/2-1 and N*n/2 ( e.g., for N*n=64 it's the mean of [31] and[32] ), 
	  // similarly the 1st and 3rd quartile; need their difference (the "interquartile range")
	  double Q1  = (currentColSorted[1*N*n/4-1] + currentColSorted[1*N*n/4])/2;
	  double Q3  = (currentColSorted[3*N*n/4-1] + currentColSorted[3*N*n/4])/2;
	  double IQR = Q3-Q1;
	  //cout << Q1 << " " << Q3 << " " << IQR << endl;
	  //cout << currentColSorted[1*N*n/4-1] << " " << currentColSorted[1*N*n/4] << " " << 1*N*n/4-1 << " " <<  1*N*n/4 << endl;
	  //cout << currentColSorted[3*N*n/4-1] << " " << currentColSorted[3*N*n/4] << " " << 3*N*n/4-1 << " " <<  3*N*n/4 << endl;
	  counter = 0;
	  for (long ii = 1; ii < N*n; ++ii) 
	  { 
	    dummy = currentColSorted[ii]; // only do the expensive vector dereferencing once
	    // ignore points outside the "fence"
	    if ((dummy < Q1-3*IQR) && (dummy > Q3+3*IQR)) { continue; }
	    curColAve[1] += dummy;
	    ++counter;
	  }
	  curColAve[1] /= counter;
	  
	  // get another average, ignoring certain problematic points
	  counter = 0;
	  for (long ii = 1; ii < N*n; ++ii) 
	  { // ignore points 511, 1023, 1535, ... (for n=512); this seems to be where the outliers live
	    if ((ii+1) % (n)) { continue; }
	    curColAve[2] += currentCol[ii];
	      ++counter;
	  }
	  curColAve[2] /= counter;
	  
	  //Anle
	  cout.precision(8); cout.setf( std::ios::scientific );
	  cout  << " ### currentTot = " << currentTotAll;
	  ofstream outCurrentTot("totCurrent.dat", ios_base::app);
	  outCurrentTot.unsetf(ios::fixed); outCurrentTot.setf( std::ios::scientific );
	  outCurrentTot << currentTotAll << endl; 
	  outCurrentTot.close();
	  
	  for (long jj=0; jj<(int)curColAve.size(); ++jj)
	  {
	    curColRelDiff[jj] = 0.0; // check that currentCol[?] is not zero, else all will diverge
	    curColRelDevMean[jj] = curColRelDevMax[jj] = (currentCol[1]) ? abs(currentCol[1]-curColAve[jj])/currentCol[1] : curColAve[jj]; 
	    for (long ii = 2; ii < N*n; ++ii) 
	    {
	      dummy = (currentCol[ii]) ? abs(currentCol[ii]-curColAve[jj])/currentCol[ii] : curColAve[jj]; // only do the vector dereferencing and computations once
	      //curColRelDiff[jj]  += abs(currentCol[ii]-currentCol[ii-1])/currentCol[ii]; // not used
	      curColRelDevMean[jj] += dummy;
	      curColRelDevMax[jj]   = max(curColRelDevMax[jj], dummy);
	    }
	    curColRelDevMean[jj] /= (N*n-1); // curColRelDiff[jj] /= (N*n-2); // not used
	      
	    cout  << "; mean/max dev: " << curColRelDevMean[jj] << "/" << curColRelDevMax[jj];
	  }
	  cout << endl << endl;
	  cout.precision(6); cout.unsetf( std::ios::scientific );
	  
	  curColumnRelDevMax  = curColRelDevMax[2];  // use the mean with some points deliberately (NOT automatically) excluded
	  curColumnRelDevMean = curColRelDevMean[2]; // dito
	}
	MPI_Bcast(&curColumnRelDevMean, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // not used
	MPI_Bcast(&curColumnRelDevMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if ( (abs((currentTotAll-currentTotAllOld)/currentTotAll) < 1e-5) && // current is converged
	      (curColumnRelDevMax < 1e-4) )                                  // maximum deviation < 1e-4 in each column
	{ vis = 2; } 
      #ifdef JUQUEEN
	if ( 1800.0-(MPI_Wtime()-startTimeProgram) < 2*(maxTimeIter+maxTimeLocal+max(maxTimeOutput,10.0)) ) { vis = 2; } // make sure to output before the 30mins are up
      #endif  
      
      #ifdef LOCAL
	if (!localDone) // don't include the local solver in the timing
	{ 
      #endif
	  maxTimeIter = max(maxTimeIter, MPI_Wtime()-startTimeIter); 
      #ifdef LOCAL
	}
      #endif
      
	if (vis)
	{
	  double startTimeOutput = MPI_Wtime();
	  if (myid == 0)
	  {
	    cout << " # outputting solution...";
	    string file;
	    char buffer[50];
	    sprintf(buffer,"currentRow.%05d.dat",N*n);
	    file.append(buffer);
	    ofstream curCol( file.c_str()); curCol.precision(8); 
	    
	   // test from Anle

	    // output current averaged over the number of rows
	    for (long ii = 1; ii < N*n; ++ii) 
	    { 
	      curCol.unsetf(ios::fixed); curCol.setf( std::ios::scientific );
	      curCol << currentCol[ii]/(N*n) << endl;
	    } 
	    curCol.close();
	  }
#ifdef H5
          string filen;
	  char buffer[50];
	  	  
	  filen = "pressure";
	  sprintf(buffer,"_%04d.h5",N*n);
	  filen.append(buffer);
	  k = 0;
	  for (long j = 0; j < n; j++)
	    for (long i = 0; i < n; i++)
	      values[k] = max(0.0, min(values[k],1.0)); ++k;
// 	  writeH5(filen, pi, pj, num_procs, myid, values);
	  
	  double *pressOutput = (double *) calloc((n*n), sizeof(double));
	  /* transpose the values */
	  for (long ii = 0; ii < n; ++ii)
	    for (long jj = 0; jj < n; ++jj)
	      pressOutput[ii*n+jj] = values[jj*n+ii];
	  writeH5(filen, pi, pj, num_procs, myid, pressOutput);  
	  free(pressOutput);
	  	  
	  filen = "current";
	  sprintf(buff,"_%04d.h5",N*n);
	  filen.append(buffer);
	  
	  double * curLoc = (double *) calloc((n+2)*(n+2),sizeof(double));
	  for (long i = 0; i < n; i++)
	  {
	    for (long j = 0; j < n; j++)
	    {
		/* NOTE all that flows in must flow out (incompressibility). 
			Therefore adding the inflows is the same as 1/2 sum |j|, 
			which is the flow through a node by its 4 neighbors */
		curLoc[i*n+j]  = abs(current_x[i][j]) + abs(current_x[i+1][j]);  // inflow in x
		curLoc[i*n+j] += abs(current_y[i][j]) + abs(current_y[i][j+1]);  // inflow in y
		curLoc[i*n+j] /= 2;
	    }
	  }
	  writeH5(filen, pi, pj, num_procs, myid, curLoc);
	  free(curLoc);
#else // H5
	  FILE *file;
	  char filename[255];

	  /* save solution (pressure) with global unknown numbers */
	  sprintf(filename, "%s_%05d.%s.%06d.%s", argv[0] , N*n,  "press", myid, "dat");
	  if ((file = fopen(filename, "w")) == NULL)
	  {
	    printf("Error: can't open output file %s\n", filename);
	    free(values);
	    MPI_Finalize();
	    exit(1);
	  }
	  
	  k = 0;
	  for (long j = 0; j < n; j++)
	  {
	    for (long i = 0; i < n; i++)
	    {
	      /// NOTE the added shift is to match the coordinates. 
	      ///      i.e., for N procs per row, the first covers the n entries in [0:(n-1)/(N*n-1)], the next [(n)/(N*n-1):1]. 
	      ///      e.g., for N*n = 128, and N = 2 (np = 4, n = 64), [0:63/127] and [64/127:1]
	      ///      e.g., for N*n = 128, and N = 4 (np = 16, n = 32), [0:31/127], [32/127:63/127], [64/127:95/127], and [96/127:1]
	      fprintf(file, "%.8f %.8f %.16e \n", 1.0*i/(N*n-1)+1.0*pi*n/(N*n-1), 1.0*j/(N*n-1)+1.0*pj*n/(N*n-1), values[k]);
	      //fprintf(file, "%.8f %.8f %.14e \n", 1.0*i/(N*n-1), 1.0*j/(N*n-1), values[k]);
	      //fprintf(file, "%.8f %.8f %.14e %06d \n", 1.0*i/(N*n-1), 1.0*j/(N*n-1), values[k], pj*N*n*n+pi*n+j*N*n+i);
	      ++k;
	    }
	    fprintf(file, "\n");
	  }
	  fflush(file);
	  fclose(file);
	  
	  /* save solution (current) with global unknown numbers */
	  sprintf(filename, "%s_%05d.%s.%06d.%s", argv[0] , N*n,  "current", myid, "dat");
	  if ((file = fopen(filename, "w")) == NULL)
	  {
	    printf("Error: can't open output file %s\n", filename);
	    free(values);
	    MPI_Finalize();
	    exit(1);
	  }
	  
	  double flowThru;
	  for (long i = 0; i < n; i++)
	  {
	    for (long j = 0; j < n; j++)
	    {
		/* NOTE all that flows in must flow out (incompressibility). 
			Therefore adding the inflows is the same as 1/2 sum |j|, 
			which is the flow through a node by its 4 neighbors */
		flowThru  = abs(current_x[i][j]) + abs(current_x[i+1][j]);  // inflow in x
		flowThru += abs(current_y[i][j]) + abs(current_y[i][j+1]);  // inflow in y
		fprintf(file, "%.8f %.8f %.14e \n", 1.0*i/(N*n-1)+1.0*pi*n/(N*n-1), 1.0*j/(N*n-1)+1.0*pj*n/(N*n-1), flowThru/2);
		//if ((pi == 0) && (i == 0)) { continue; }
		//else { currentTotLoc += current_x[i][j]; }
	    }
	    fprintf(file, "\n");
	  }
	  //cout << "# proc " << myid << " 2currentTotLoc = " << currentTotLoc << endl;
	  fflush(file);
	  fclose(file);
#endif // H5
	#ifndef JUQUEEN
	  vis = 0; // on JUQUEEN: only output at convergence or the end of the time the program can run (30 mins)  /// NOTE
	#endif
	  if (myid == 0) { cout << "\n #...done.\n\n"; }
	  maxTimeOutput = max(maxTimeOutput, MPI_Wtime() - startTimeOutput);
	}
	free(values); 
	if (totalIterations > iterMax) { vis = 3; } // code needs to be restarted once in a while. Don't go on indefinitely
      }
          
      if (myid == 0)
      {
	  printf("\n");
	  printf(" Iterations = %d\n", totalIterations); //num_iterations
	  printf(" Final Relative Residual Norm = %e\n", final_res_norm);
	  printf("\n");
      }      
//       //if ( (final_res_norm < TOL*3) && (num_iterations != iter) && (final_res_norm > 0.0) ) // if converged, exit loop
//       if ( (final_res_norm < TOL*3) && (num_iterations != iter) ) // if converged, exit loop
//       { 
// 	if (myid == 0) { cout << " # hypre convergence reached. \n # All done...\n\n"; }
// 	break;
//       } 
      if ( (abs((currentTotAll-currentTotAllOld)/currentTotAll) < 1e-5) &&      // current is converged
	        (curColumnRelDevMax < 1e-4) )                                      // maximum deviation < 1e-4 in each column
      { 
	if (myid == 0) 
	{ 
	  cout  << " # change in total current = " << abs((currentTotAll-currentTotAllOld)/currentTotAll) 
		<< " ==> convergence in current achieved." << endl
		<< " # average deviation from mean in each column = " << curColumnRelDevMean << "." << endl
		<< " # maximum deviation from mean = " << curColumnRelDevMax << "." << endl;
	  cout  << " # All done...\n\n";
	}
	break;
      } 
      else { currentTotAllOld = currentTotAll; }
      cout.flush();
      cerr.flush();
      
      if (myid == 0)
      {
        cout << " # maxTimeIter = " << maxTimeIter << " ("<<iter<<");"
             << " maxTimeOutput = " << maxTimeOutput << ";"
	     << " maxTimeLocal = " << maxTimeLocal << " ("<<iterLocal<<");"
             << " total time program = " << (MPI_Wtime()-startTimeProgram) << " \t"
             << endl << endl;
      }

    //#ifdef JUQUEEN
      if ( vis >= 2) { break; } /// NOTE
    //#endif
      
//       sprintf(filename, "%s.%06d.%s", "_vectorX_after", myid, "txt");
//       HYPRE_StructVectorPrint (filename, x, 1);
   
   } // end while (1)
   
   free(sendbuf);
   free(recvbuf);
   free(recvbuflarge);

#ifndef H5
   MPI_Type_free(&blocktype);
   MPI_Type_free(&blocktype2);
#endif

   /* Free memory */
   HYPRE_StructGridDestroy(grid);
   HYPRE_StructStencilDestroy(stencil);
   HYPRE_StructMatrixDestroy(A);
   HYPRE_StructVectorDestroy(b);
   HYPRE_StructVectorDestroy(x);
   
   /* Finalize MPI */
   MPI_Finalize();

   if (vis == 3) { system ("touch ___restart"); } // indicate to the launching script that restart is required
   return (0);
}

#ifdef CLUSTER

int getNeighbor(int s, int j) 
{  // identify neighbor elements, non-periodic boundaries
  switch(j) {
    case 0: if ((s % (n*N)) == 0)     { return(((periodicY) ? (s-1+n*N):-1)); } else { return(s-1);  } // top neighbor
    case 1: if ((s % (n*N)) == n*N-1) { return(((periodicY) ? (s+1-n*N):-1)); } else { return(s+1);  } // bottom
    case 2: if ((s / (n*N)) == 0)     { return(((periodicX) ? (s-n*N)+n*N*n*N:-1)); } else { return(s-n*N); } // left
    case 3: if ((s / (n*N)) == n*N-1) { return(((periodicX) ? (s+n*N)-n*N*n*N:-1)); } else { return(s+n*N); } // right

//     case 0: if ((s % (n*N)) == 0)     { return(-1); } else { return(s-1);  } // top neighbor
//     case 1: if ((s % (n*N)) == n*N-1) { return(-1); } else { return(s+1);  } // bottom
//     case 2: if ((s / (n*N)) == 0)     { return(-1); } else { return(s-n*N); } // left
//     case 3: if ((s / (n*N)) == n*N-1) { return(-1); } else { return(s+n*N); } // right
    
    default: return(-1);
  }      
}

void addToCluster(vector <clusterInfo> & c, int s) 
{ // add element to current cluster, update the extrema
  (c.back()).elem++; 
  checkX[ s / (n*N) ] = 1; 
  checkY[ s % (n*N) ] = 1;
}

void clusterFind(vector< clusterInfo >& c, clusterMode mode) 
{  
  int comp;
  if (mode == CONT) {
    clusterNumber =  10;
    cout << "\n # Finding contact clusters..." << endl;    
    comp = indCont; 
  } else {
    clusterNumber =  -10;
    cout << "\n # Finding non-contact clusters..." << endl;
    comp = indNonCont; 
  }
  int numSitesToTest = 0;
  int initialSite = 0; 
  
  while (1) // not all done
  {
    while (initialSite < n*N*n*N) // find cell to start next cluster
    {
      if (contact[initialSite] == comp) // starting site identified
      {
	if (mode == CONT) { ++clusterNumber; } else { --clusterNumber; }
	contact[initialSite] = clusterNumber; 
	break; 
      } else { ++initialSite; }            // try next cell
    } if (initialSite == n*N*n*N) { break; } // if all cells have been treated, exit
    //cout << "# starting site for cluster "<<clusterNumber<<": " << initialSite << endl;
  
    // --- reset spanning check
    for (long ix = 0; ix < n*N; ++ix) { checkX[ix] = 0; }
    for (long iy = 0; iy < n*N; ++iy) { checkY[iy] = 0; }
    
    c.push_back(newCluster(0,0,0,0,clusterNumber)); 
    sitesToTest[numSitesToTest++] = initialSite; 
    addToCluster(c, initialSite);
    while (numSitesToTest > 0) // not all members of current cluster identified
    {
      int site = sitesToTest[--numSitesToTest];
      int neighborNum = 0;
      for (long j = 0; j < 4; ++j ) // loop over neighbors
      {
	int neighbor = getNeighbor(site, j); 
	if ((contact[neighbor] == comp) || (contact[neighbor] == clusterNumber)) // neighbor part of same cluster
	   { ++neighborNum; } 
	if ((neighbor >= 0) && (contact[neighbor] == comp)) // neighbor not on boundary and fits test
	{ 
	  contact[neighbor] = clusterNumber;
	  sitesToTest[numSitesToTest++] = neighbor;           // check if neighbor has further neighbors
	  addToCluster(c, neighbor);                          // add current neighbor to cluster counter
	} // end if another member of current cluster identified	
      }   // end for neighbors checked
      if (neighborNum < 4) { c.back().surface += 1; }       // check if current site is a boundary site
    }     // end while all cluster members found
    
    // --- check if cluster is spanning
    for (long ix = 0; ix < n*N; ++ix) { (c.back()).spanningX += checkX[ix]; } 
    for (long iy = 0; iy < n*N; ++iy) { (c.back()).spanningY += checkY[iy]; }
    (c.back()).spanningX /= (n*N); // will be zero if not all indices are covered
    (c.back()).spanningY /= (n*N); // dito    
        
  }       // end while all done
  
  if (mode == CONT) { sort(c.begin(),c.end(),clusterOrderCont); }
  else              { sort(c.begin(),c.end(),clusterOrderNonCont); }
}

bool clusterOrderNonCont (const clusterInfo & c1, const clusterInfo & c2) 
{ // sorts clusters
//   // --- sort by size
//   if (c1.elem > c2.elem) {return true; } 
//   else                   {return false; }

  // --- sort contact clusters by "whether-spanning in x-direction" (if no, no percolation)
  if (c1.spanningX > c2.spanningX) { return true;  }
  else                             { return false; }  
}

bool clusterOrderCont (const clusterInfo & c1, const clusterInfo & c2) 
{ // sorts clusters
  // --- sort contact clusters by "whether-spanning in y-direction" (if yes, no percolation)
  if (c1.spanningY > c2.spanningY) { return true;  }
  else                             { return false; }  
}

void dumpSurface(int num, string file) 
{ 
  file.erase(file.length()-3,3);
  char buffer[50];
  sprintf(buffer,"clusterPic.%d.out",num);
  file.append(buffer);
  ofstream dumpPic( file.c_str()); 
  
  dumpPic << "P1" << endl << n*N << " " << n*N << endl;
  // NOTE: ix must vary faster than ix, otherwise the output will be transposed, 
  //       and iy must run run backwards, otherwise it will be inverted 
  ///      (all relative to gnuplot!!!)
  for (long iy = n*N-1; iy >= 0; --iy) {
    for (long ix = 0; ix < n*N; ++ix) {
      dumpPic << cont[ix*n*N + iy];
      if ((num != 0) && (num != 9)) { dumpPic << " "; }
    } dumpPic << endl;
  }   dumpPic << endl;
  dumpPic.close();
  
//   string command; 
//   command.append("pnmtopng ");
//   command.append(file);
//   command.append(" > ");
//   command.append(file.replace(file.length()-4,4,".png"));
//   cout << "# executing " << command << endl << endl;
//   system (command.c_str());  
}

void cleanup(int error) {
  free (contact);
  free (cont);
  free (sitesToTest);
  free (checkX);
  free (checkY);
  if (error > 0) { exit(error); }  
}
#endif /* CLUSTER */     

#ifdef H5
// from: http://stackoverflow.com/questions/1686869/searching-a-hdf5-dataset
bool DoesDatasetExist(hid_t mFileId, const std::string& rDatasetName)
{
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
    // This is a nice method for testing existence, introduced in HDF5 1.8.0
    htri_t dataset_status = H5Lexists(mFileId, rDatasetName.c_str(), H5P_DEFAULT);
    return (dataset_status>0);
#else
    bool result=false;
    // This is not a nice way of doing it because the error stack produces a load of 'HDF failed' output.
    // The "TRY" macros are a convenient way to temporarily turn the error stack off.
    H5E_BEGIN_TRY
    {
	hid_t dataset_id = H5Dopen(mFileId, rDatasetName.c_str());
	if (dataset_id>0)
	{
	    H5Dclose(dataset_id);
	    result = true;
	}
    }
    H5E_END_TRY;
    return result;
#endif
}

void readH5(string name, int pi, int pj, int num_procs, int myid, double * field)
{
   // ===============================================================
   // =================== read input from HDF5 ======================
   // ===============================================================
   
   /*
    * HDF5 APIs definitions
    */ 	
   hid_t       file_id, dset_id;         /* file and dataset identifiers */
   hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
   hsize_t     chunk_dims[RANK];            /* chunk dimensions */
   hsize_t	count[RANK];	          /* hyperslab selection parameters */
   hsize_t	stride[RANK];
   hsize_t	block[RANK];
   hsize_t	offset[RANK];
   hid_t	plist_id;                 /* property list identifier */
   herr_t	status_h5;
    
   MPI_Info info  = MPI_INFO_NULL;

   if (myid == 0) { cout << "\n ### Reading file " << name << "...\n\n"; }
   
   /*
    * Open file_id and dataset using the default properties.
    */
   plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
   file_id = H5Fopen (name.c_str(), H5F_ACC_RDONLY, plist_id); // name.c_str() is the file name
   H5Pclose(plist_id);

   std::string fname = name;
   std::size_t pos = fname.find(".cluster.h5"); // position for first occurence of string
   std::string dataset = fname.substr(0,pos);   // everything up to first occurence
   
   // check whether the dataset exists
   string dset_name = name;
   if (DoesDatasetExist(file_id, dataset))
   { dset_id = H5Dopen (file_id, dataset.c_str(), H5P_DEFAULT); } // dataset.c_str() is the dataset name
   else if (DoesDatasetExist(file_id, dset_name))
   { dset_id = H5Dopen (file_id, dset_name.c_str(), H5P_DEFAULT); } // name.c_str() is the dataset name
   else if (DoesDatasetExist(file_id, dset_name.replace(dset_name.end()-3,dset_name.end(),"_S4957.h5"))) // append _S4957
   { dset_id = H5Dopen (file_id, dset_name.c_str(), H5P_DEFAULT); } // dset_name.c_str() is the dataset name
   else
   { std::cerr << " ##### Error: dataset " << dataset << " or " << name << " or " << dset_name << " not found! Terminating...\n\n"; exit(-1); }
   
   filespace = H5Dget_space (dset_id);
    
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
   /*
    * Define and select the hyperslab to use for reading.
    */
   chunk_dims[0] = n;   
   chunk_dims[1] = n;   
   memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 
    
   count[0] = 1;               // one block in x-direction
   count[1] = 1 ;              // one block in y-direction
   stride[0] = 1;              // the blocks have no gaps between them in x-direction
   stride[1] = 1;              // the blocks have no gaps between them in y-direction
   block[0] = chunk_dims[0];   // contiguous x-size of block is chunk_dims[0]
   block[1] = chunk_dims[1];   // contiguous x-size of block is chunk_dims[1]
   offset[0] = pi*n;           // x-index [in data file] of (0,0) entry [in memory] of block is pi*n
   offset[1] = pj*n;           // y-index [in data file] of (0,0) entry [in memory] of block is pj*n
   // // if reordering should be required: ((N-1)-pi)*n;
   
   status_h5 = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, stride, count, block);
   assert(status_h5 >= 0);

   /*
    * Read the data using the previously defined hyperslab.
    */
   status_h5 = H5Dread (dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, field);
   assert(status_h5 >= 0);
   
//    for (long j=0; j<num_procs; ++j)
//    {
//      if (myid == j)
//      {
//        /*
// 	* Output the data to the screen.
// 	*/
//        printf ("\nData as written to disk by hyberslabs:\n");
//        for (ulong i=0; i < (ulong)chunk_dims[0]*chunk_dims[1]; i++)
//        {
// 	 cout << myid << " " << i << " " << field[i] << " " << endl;
//        }
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
 
   /*
    * Close and release resources.
    */
   status_h5 = H5Pclose (plist_id);  assert(status_h5 >= 0);
   status_h5 = H5Dclose (dset_id);   assert(status_h5 >= 0);
   status_h5 = H5Fclose (file_id);   assert(status_h5 >= 0);
   status_h5 = H5Sclose (filespace); assert(status_h5 >= 0);
   status_h5 = H5Sclose (memspace);  assert(status_h5 >= 0);
}

void writeH5( string name, int pi, int pj, int num_procs, int myid,double * field)
{   
   // ===============================================================
   // =================== write output to HDF5 ======================
   // ===============================================================
   
   if (myid == 0) { cout << "\n # ...writing to file " << name << "..."; }
  
   /*
    * HDF5 APIs definitions
    */ 	
   hid_t       file_id, dset_id;         /* file and dataset identifiers */
   hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
   hsize_t     dimsf[RANK];                 /* dataset dimensions */
   hsize_t     chunk_dims[RANK];            /* chunk dimensions */
   hsize_t	count[RANK];	          /* hyperslab selection parameters */
   hsize_t	stride[RANK];
   hsize_t	block[RANK];
   hsize_t	offset[RANK];
   hid_t	plist_id;                 /* property list identifier */
   herr_t	status_h5;
   double * data = NULL;
    
   MPI_Info info  = MPI_INFO_NULL;
   
   /* 
    * Set up file access property list with parallel I/O access
    */
   plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

   /*
    * Create a new file collectively and release property list identifier.
    */
   file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); // name.c_str() is the file name
   H5Pclose(plist_id);

   /*
    * Create the dataspace for the dataset.
    */
   dimsf[0] = N*n;
   dimsf[1] = N*n;
   chunk_dims[0] = n;   
   chunk_dims[1] = n;   
   filespace = H5Screate_simple(RANK, dimsf, NULL); 
   memspace  = H5Screate_simple(RANK, chunk_dims, NULL);    
        
   /*
    * Create chunked dataset.
    */
   plist_id = H5Pcreate(H5P_DATASET_CREATE);
   status_h5 = H5Pset_layout( plist_id, H5D_CHUNKED );
   assert(status_h5 >= 0);
   status_h5 = H5Pset_chunk(plist_id, RANK, chunk_dims);
   assert(status_h5 >= 0);
   dset_id = H5Dcreate(file_id, name.c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT); // name.c_str() is the data set name
   H5Pclose(plist_id);
   H5Sclose(filespace);

   /* 
    * Each process defines dataset in memory and writes it to the hyperslab
    * in the file.
    */
   count[0] = 1;               // one block in x-direction
   count[1] = 1 ;              // one block in y-direction
   stride[0] = 1;              // the blocks have no gaps between them in x-direction
   stride[1] = 1;              // the blocks have no gaps between them in y-direction
   block[0] = chunk_dims[0];   // contiguous x-size of block is chunk_dims[0]
   block[1] = chunk_dims[1];   // contiguous x-size of block is chunk_dims[1]
   offset[0] = pi*n;           // x-index [in data file] of (0,0) entry [in memory] of block is pi*n
   offset[1] = pj*n;           // y-index [in data file] of (0,0) entry [in memory] of block is pj*n
   // // if reordering should be required: ((N-1)-pi)*n;
   
   /*
    * Select hyperslab in the file.
    */
   filespace = H5Dget_space(dset_id);
   status_h5 = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
   assert(status_h5 >= 0);
   
   /*
    * Create property list for collective dataset write.
    */
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
 
//    /*
//     * Initialize data buffer 
//     */
//    data = (double *) calloc(chunk_dims[0]*chunk_dims[1],sizeof(double));
//
//    for (long i = 0; i < num_procs; ++i)  // for 2d field to be outputted:
//    {
//      if (myid == i)
//      {
//        for (ulong j=0; j < (ulong)chunk_dims[0]*chunk_dims[1]; j++) 
//        {
// 	 ulong ii = j/chunk_dims[0]+1; // the +1 comes from the fact that s[][] is padded
// 	 ulong jj = j%chunk_dims[1]+1;
// 	 data[j] = field[ii][jj];
// 	 //cout << myid << " " << j << " " << data[j] << " " << ii << " " << jj << " " << s[ii][jj] << " " << endl;
//        }
//        //cout << endl;
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }

//    for (long k = 0; k < num_procs; ++k)  // for 1d field to be outputted:
//    {
//      if (myid == k)
//      {
//        for (long i=0; i < n; i++) 
//        for (long j=0; j < n; j++) 
//        {
// 	 cout << myid << " " << "field["<<i*n+j<<"]: " << field[i*n+j] << " " << endl;
//        }
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
   
   /*
    * Write to file
    */   
   status_h5 = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, field);
   assert(status_h5 >= 0);
   free(data);

   /*
    * Close/release resources.
    */
   H5Dclose(dset_id);
   H5Sclose(filespace);
   H5Sclose(memspace);
   H5Pclose(plist_id);
   H5Fclose(file_id);
   
}

#endif // H5
