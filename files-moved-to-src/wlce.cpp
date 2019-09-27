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

#include "wlce-class.cpp"

using namespace std;

//define constants
const int binNum=30; 
typedef vector< vector<double> > matrix;


int main(int argc, char **argv) {
wlce wfun;
double lb, ub, gamma_start,gamma_end;
int cryStruct,numEle,binSize,N_skip;
int show_detail_step, cal_short_range_order, restart;
char *log="log.wlce";
char *restartLog="restart.log.wlce";
const char * split = ",";
char * p;
char ch[1001];
int *ptr_numEle;
ptr_numEle = (int*)( &numEle );
//-------------read input data-----------------
 fstream fin("incar", ios::in);
 int row=0;
  if (!fin) {cout<<"Sorry bad file"; exit(0);}
  while (!fin.eof())
   {
     row=row+1;
     fin.getline(ch,100);
     if (row==2) { p=strtok(ch,split);  cryStruct=atoi(p);}
     if (row==3) { p=strtok(ch,split);  numEle=atoi(p);}//*ptr_NumEle=atoi(p);}
     if (row==4) { p=strtok(ch,split);  lb=atof(p); }
     if (row==5) { p=strtok(ch,split);  ub=atof(p); }
     if (row==6) { p=strtok(ch,split);  gamma_start=atof(p); }
     if (row==7) { p=strtok(ch,split);  gamma_end=atof(p); }
     if (row==8) { p=strtok(ch,split);  show_detail_step=atoi(p);}
     if (row==9) { p=strtok(ch,split);  cal_short_range_order=atoi(p);}
     if (row==10) { p=strtok(ch,split);  restart=atoi(p);}
     if (row==11) { p=strtok(ch,split);  binSize=atoi(p);}
     if (row==12) { p=strtok(ch,split);  N_skip=atoi(p);}
     if (row>1) {cout <<"I finished reading input info"<<row<<" "<<ch<<'\n';}
  }
  fin.close();
//---------------------------------------------
vector<double> dos;
vector<int> hist; 
vector<int> count; 
vector<int> component;
dos.resize(binSize);
hist.resize(binSize);
count.resize(binSize);
component.resize(numEle);
matrix sro; 
wfun.resizeVec(sro,binSize,numEle*numEle);
float E0=0;
float E1=0;
char *strfilename;
char *restart_str_file="restart.struct";
if (restart)
  {strfilename=restart_str_file;}
else
  {strfilename="init.struct";}
char *latticefilename="lat.in";
char *ecifile="eci.out";
//-----------------read lattice parameter---------------------
double latt;
fstream flatt(latticefilename, ios::in);
int line=0;
if (!flatt) {cout<<"Sorry bad file"; exit(0);}
while (!flatt.eof()){
  line=line+1;
  flatt.getline(ch,100);
  if (line==1) { p=strtok(ch,split);  latt=atof(p);}
  cout <<"lattice para: "<<latt<<'\n';
  break;
}
flatt.close();
//-----------------read lattice file and get atom labels and space groups-------------------
int readocc=0;
  Structure lattice;                                                             
  Array<Arrayint> labellookup;                                                   
  Array<AutoString> label;                                                       
  rMatrix3d axes;                                                                
  {                                                                              
    ifstream latticefile(latticefilename);                                       
    if (!latticefile) ERRORQUIT("Unable to open lattice file");                  
    if (readocc) {                                                               
      Array<Array<Real> > occ;                                                   
      parse_rndstr_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &occ, &labellookup, &label, latticefile, &axes);
    }                                                                            
    else {                                                                       
      parse_lattice_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &labellookup, &label, latticefile, &axes);                                      
    }                                                                            
    wrap_inside_cell(&lattice.atom_pos,lattice.atom_pos,lattice.cell);           
  } 
  SpaceGroup spacegroup;                                                                 
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
//--------------------read initial structure files------------------------
  Structure str;
  ifstream strfile(strfilename);
  if (!strfile) ERRORQUIT("Unable to open structure file");
  int strnum=1;
  while (!strfile.eof()) {
    parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,label,strfile);
    skip_to_next_structure(strfile);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    rMatrix3d supercell=(!lattice.cell)*str.cell;
    rMatrix3d rounded_supercell=to_real(to_int(supercell));
    rMatrix3d transfo=lattice.cell*rounded_supercell*(!str.cell);
    str.cell=transfo*str.cell;
    //cout << typeid(str.atom_type).name() << endl;
   // cout<<str.cell(0,0)<<endl;
   // cout<<str.atom_type(100)<<" "<<str.atom_type.get_size()<<" "<<str.atom_type(1)<<" "<<str.atom_type[1]<<endl;
   // cout<<label<<endl;
   // cout<<str.atom_pos<<endl;
  //  cout<<labellookup<<endl;
/*    for (int i=0; i<str.atom_pos.get_size(); i++) {
      cout<<label(i) <<endl;
     // cout <<str.atom_pos.get_size()<<'\n';
    }*/
  }
for (int i=0;i<numEle;i++){
  component[i]=wfun.countSameElement(str.atom_type,i);
}
// Arrayint str.atom_type;
//   cout<<component[0]<<" "<<component[1]<<endl;
//--------------------read clusterlist from cluster.out and eci from eci.out-------------------------------
  LinkedList<MultiCluster> clusterlist;
  ifstream clusterfile("clusters.out"); 
  if (!clusterfile) ERRORQUIT("Unable to open clusters.out file");
  read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);

  LinkedListIterator<MultiCluster> icluster(clusterlist);
  Array<Real> eci;
  if (strlen(ecifile)>0) {
    ifstream file(ecifile);
    if (!file) ERRORQUIT("Unable to open ECI file.");
    LinkedList<Real> ecilist;
    while (skip_delim(file)) {
      Real e; 
      file >> e;
      ecilist << new Real(e);
    }         
    LinkedList_to_Array(&eci,ecilist);
    if (clusterlist.get_size()!=eci.get_size()) ERRORQUIT("Number of ECI does not match number of clusters."); 
  }
//--------------write str.out------------------
//char *outfile="str.out";
//------------------------------------------------------------------------
double gamma=gamma_start;
double gamma_final=gamma_end;
int column;
//array<double,binNum> dos; 
////array<double,binNum,NumEle> sro; //define short range parameters
//double sro[binNum][NumEle];
//array<int,binNum> count,hist;
if (restart) {
  fstream logf(restartLog,ios::in);     
  int row=0;         
  char *strings[11];
  if (!logf) {cout<<"Sorry bad file"; exit(0);}
  while (!logf.eof()){                    
     row=row+1;         
     logf.getline(ch,1000);         
     //cout <<row<<" "<<ch<<endl;  
     if (row==1) { p=strtok(ch,split);  gamma=atof(p); }
     if (row>2 and row<binSize+3) { p=strtok(ch,split);  dos[row-3]=atof(p); }
     if (row>binSize+3 and row<2*binSize+4) { p=strtok(ch,split);  hist[row-binSize-4]=atof(p); }
     if (row>2*binSize+4 and row<3*binSize+5) {
        column=0;
        strings[column]=strtok(ch,split);        
        while (strings[column] != NULL){ 
          sro[row-2*binSize-5][column]=stod(strings[column]);
          if (column==numEle*numEle) {count[row-2*binSize-5]=atoi(strings[column]);} //cout<<count[row-2*binSize-5]<<endl;}
          strings[++column] = strtok(NULL, split);} //cout <<count[row-2*binNum-5]<<endl; }
        for (int i=0;i<numEle*numEle;i++){
          sro[row-2*binSize-5][i]=sro[row-2*binSize-5][i]*count[row-2*binSize-5];}
     }
   } 
 // readLog(&gamma,dos,hist,count,sro,restartLog);
 // cout << gamma <<endl;
 // for (int i=0;i<60;i++){
 // cout << sro[i]<<endl;}
  }
else {
  for (int i = 0; i < binSize; i++) // ...initialize it
  {
    dos[i]  = 0;
    count[i]= 0;
    hist[i] = 1;
    for (int j1 = 0; j1 < numEle; j1++){
    for (int j2 = 0; j2 < numEle; j2++){
      sro[i][j1*numEle+j2]  = 0;}
   // count[i]<<","<<""
  }
}}
cout <<"I finished preparing the calculations, but do not start it yet!"<<'\n'; 
E0=wfun.cal_etot(str, clusterlist,spacegroup, labellookup, eci);
cout<<"Hello!"
cout <<"first E0 calcualated!"<<E0<<endl;
Structure initial_struct;
vector<int> two_atom_ids(2);
int atom_a, atom_b;
//int num_species=str.atom_pos.get_size()/numEle;
while (E0 <lb or E0>ub){
  two_atom_ids=wfun.random_pick_two_numbers(component);
  atom_a=two_atom_ids[0];
  atom_b=two_atom_ids[1];

// Structure initial_struct;
 initial_struct=wfun.swapped_str(atom_a,atom_b,str);                             
 E1=wfun.cal_etot(initial_struct, clusterlist,spacegroup, labellookup, eci); 
 if (((E0 < lb) and (E1>E0)) or ((E0 > lb) and (E1<E0))){
   str=initial_struct;
   E0=E1;}
  cout <<E0<<endl;
  cout <<"I'm finding the E0 between El and Eu! \n";
 //cout <<E0<<" "<<bin(E0)<<'\n';
}
int bin_now;
int N=0; //record iteration step
Structure ideal_str;
while (gamma>gamma_final){
  N +=1;
  two_atom_ids=wfun.random_pick_two_numbers(component);
  atom_a=two_atom_ids[0]; 
  atom_b=two_atom_ids[1]; 
  ideal_str=wfun.swapped_str(atom_a,atom_b,str);
  E1=wfun.cal_etot(ideal_str, clusterlist,spacegroup, labellookup, eci);
  if (E1>lb and E1<ub){
    srand(time(NULL));
    double r = ((double) rand() / (RAND_MAX));  //(rand() % 10000) / 10000.0;// ((double) rand() / (RAND_MAX)) + 1;
    if (r<exp(dos[wfun.bin(E0,lb,ub,binSize)]-dos[wfun.bin(E1,lb,ub,binSize)])){
      str=ideal_str;
      E0=E1;}
  } //}
  bin_now=wfun.bin(E0,lb,ub,binSize);
  dos[bin_now] = dos[bin_now]+gamma;
  hist[bin_now] += 1;
  if (cal_short_range_order){
    count[bin_now] += 1;
    for (int j1=0;j1<numEle;j1++){
    for (int j2=0;j2<numEle;j2++){
      sro[bin_now][j1*numEle+j2]   = sro[bin_now][j1*numEle+j2]+wfun.short_range_order(str,latt,1,numEle,cryStruct)[j1][j2];}}
  }
  if (show_detail_step) {
    cout <<E0<<" "<<bin_now<<" ";
    for (int i1=0;i1<numEle;i1++){
    for (int i2=0;i2<numEle;i2++){
      cout<<sro[bin_now][i1*numEle+i2]/count[bin_now]<<" ";}}
    cout<<count[bin_now]<<'\n';}
    //cout<<N<<'\n';}
  if (wfun.isFlat(hist,N,N_skip,2)){
    gamma=gamma/2.0;
    for (int i = 0; i < binSize; i++){ hist[i] = 1;}
    N=0;
  //}
  //print logs when gamma changes
  wfun.writeLog(gamma,dos,hist,count,sro, log);
  ofstream wfile(restart_str_file);
  write_structure(str, label, axes, wfile);
  }
}
 if ((N%N_skip)==0){
  wfun.writeLog(gamma,dos,hist,count,sro, log);
  ofstream wfile(restart_str_file);
  write_structure(str, label, axes, wfile);
  }
return 0;
}

