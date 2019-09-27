//#include <mpi.h>
#include <string>
#include <iostream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include "parse.h"
#include "clus_str.h"
#include "getvalue.h"
#include "ctype.h"
#include "version.h"
#include "plugin.h"
#include <random>
#include <vector> 
#include <typeinfo> //used to check the type of a variable

using namespace std;

/*
//define constants
const int binNum=30; //66; //int(argv[2])-int(argv[1]); //46;
const int NumEle=4;
typedef vector< vector<double> > matrix;*/

class wlce
{
public:
//define constants
const int binNum=30; //66; //int(argv[2])-int(argv[1]); //46;
const int NumEle=4;
typedef vector< vector<double> > matrix;

//define new functions
void readLog(double gamma,vector<double> dos, vector<int> hist, vector<int> count, matrix sro, char *logfile);
void writeLog(double gamma,vector<double> dos, vector<int> hist, vector<int> count, matrix sro, char *logfile);
//void writeLog(double gamma,array<double,binNum> dos, array<int,binNum>hist, array<int,binNum> count, double sro[binNum][NumEle], char *logfile);
int isFlat(vector<int> hist,int maskedN,int N_step,int maskedThreshold);
int bin(double E, double lb, double ub, int binSize);
float etot(const char* cmd);
//void swap_write_struct(int atom_a,int atom_b, Structure &str, Array<AutoString> &atom_label, rMatrix3d &axes, char *file);
Structure swapped_str(int atom_a,int atom_b, Structure &str);
float cal_etot(Structure ideal_str,LinkedList<MultiCluster> clusterlist,SpaceGroup spacegroupi, Array<Arrayint> labellookup, Array<Real> eci);
//vector<double> short_range_order_bcc(Structure str, double latt, int ith_NN, int numEle);
vector<vector<double>> short_range_order(Structure str, double latt, int ith_NN, int numEle, int cryStru);
void resizeVec( vector<vector<double> > &vec , const unsigned short ROWS , const unsigned short COLUMNS );
int countSameElement(Arrayint labels,int elementLabel);
vector<int> random_pick_two_numbers(vector<int> components);

private:
int genRandom(int range_from, int range_to);
float maxArray(vector<float>rest, float maskedValue);
Vector3d<double> NNNbcc(double latt,int ithNN,int jth_ithNN);
Vector3d<double> NNNfcc(double latt,int ithNN,int jth_ithNN);
void moveInBox(Vector3d<double> &Atom, double Box[3]);
int matchVector(Vector3d<double> atom1, Vector3d<double> atom2, double Box[3]);
};

