#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void printHelp(string);
int nx = 1024, ny = 1024;
bool auto_n = true;
string filename = "gap0.dat";


int main(const int argc, const char *argv[]){
  // process command line inputs
  if (argc > 1) {
    //cout << "checking args\n";//DEBUG
    vector<string> argument(argc);
    for (int i = 0; i<argc; ++i) argument[i] = argv[i];
    
    for (int i = 1; i<argc; ++i) {
      if (argument[i]=="-h" || argument[i]=="--help") printHelp(argument[0]);
    }
    
    for (int i = 1; i<argc-1; ++i) {
      //cout << argument[i] << "\n";//DEBUG
      if (argument[i]=="-i") filename = argument[i+1];
      else if (argument[i]=="-nx") {
	nx = stoi(argument[i+1]);
	auto_n = false;
      }
      else if (argument[i]=="-ny") {
	ny = stoi(argument[i+1]);
	auto_n = false;
      }
    }
  }
  if (!auto_n) cout << "reading " + filename << " with manual resolution (" << nx << "," << ny << ")\n";

  //read in the file and output the file one more time
  vector<double> gapOld, gapNew; 
  gapNew.resize(nx*ny);
  
  
  // read the old part  
    std::ifstream configIn(filename);
    if (!configIn.is_open()) {
      cerr << "ERROR: file not found!\n";
      cerr << "use argument --help for more information.\n";
      exit(0); // check if the file exists
    }

    double dummy, gapRead;
    std::string str;
    char dummyChar; // used to hold the "#" character in gap file

    if (auto_n) {
      configIn >> dummyChar >> nx >> ny;
      cout << "reading " + filename << " with auto resolution (" << nx << "," << ny << ")\n";
    }
    else {configIn >> dummyChar >> dummy >> dummy;}

    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            for (int i = 1; i < 3; i++)
                configIn >> dummy; // skip other column
            configIn >> gapRead;
            getline(configIn, str);
            int ii = ix * ny + iy;
            //locCont[ii] = (gapRead < 1.0e-9);
            if (gapRead < 0) gapRead = 0.;
	    gapNew[ii] = gapRead;
	}
        getline(configIn, str); // skip the extra line
    }
    configIn.close();

    string num = to_string(ny);
    if (ny<1000) num = "0"+num;
    if (ny<100) num = "0"+num;
    if (ny<10) num = "0"+num;
    ofstream output("gap_"+num+".dat");
    for (int ix = 0; ix < nx; ++ix){
      for (int iy = 0; iy < ny; ++iy){
        int ii = ix*ny + iy;
        output << ix*1./nx << "\t" << iy*1./ny << "\t" << gapNew[ii] << endl;
      }
      output << endl;
    }
    output.close();
}

void printHelp(string name) {
  //cout << "printHelp()\n";//DEBUG
  cout << "Usage:\n";
  cout << name + " [options]\n";
  cout << "\nOptions:\n";
  cout << "-i  : path to the gap file to be converted (default: gap0.dat)\n";
  cout << "-nx : number of grid points in x direction (default: 1024)\n";
  cout << "-ny : number of grid points in y direction (default: 1024)\n";
  cout << "\nExample:\n";
  cout << name << " -i gap0.dat -nx 1024 -ny 1024\n";
  exit(0);
}
