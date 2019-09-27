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

#include "wlce-class.hpp"

using namespace std;
/*
//define constants
const int binNum=30; //66; //int(argv[2])-int(argv[1]); //46;
const int NumEle=4;
typedef vector< vector<double> > matrix;
*/

//----------------define functions here-----------------------
Vector3d<double> wlce::NNNbcc(double latt,int ithNN,int jth_ithNN){         
Vector3d<double> offset;                                              
switch(ithNN){ 
  case 1:
    switch(jth_ithNN){
      case 0: //offset(0.5*latt,0.5*latt,0.5*latt);      
        offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;
      case 1:         
        offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;
      case 2:         
        offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;
      case 3:         
        offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;
      case 4:         
        offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;
      case 5:         
        offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;
      case 6:         
        offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;
      case 7:         
        offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;
    }   
  case 2:
    exit;
}
return offset;
}

Vector3d<double> wlce::NNNfcc(double latt,int ithNN,int jth_ithNN){
Vector3d<double> offset;
switch(ithNN){                                          
  case 1:                  
    switch (jth_ithNN){                                   
     case 0: 
       offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.0*latt; break;
     case 1:                                       
       offset[0]= 0.5*latt; offset[1]= 0.0*latt; offset[2]= 0.5*latt; break;
     case 2:                                       
       offset[0]= 0.0*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;
     case 3:                                       
       offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.0*latt; break;
     case 4:                                       
       offset[0]= 0.5*latt; offset[1]= 0.0*latt; offset[2]=-0.5*latt; break;
     case 5:                                       
       offset[0]= 0.0*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;
     case 6:                                       
       offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.0*latt; break;
     case 7:                                       
       offset[0]=-0.5*latt; offset[1]= 0.0*latt; offset[2]= 0.5*latt; break;
     case 8:                                       
       offset[0]= 0.0*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;
     case 9:                                       
       offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.0*latt; break;
     case 10:                                      
       offset[0]=-0.5*latt; offset[1]= 0.0*latt; offset[2]=-0.5*latt; break;
     case 11:                                      
       offset[0]= 0.0*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;
    }
  case 2:
    exit;
}
return offset;
}

vector<vector<double>> wlce::short_range_order(Structure str, double latt, int ith_NN, int numEle, int cryStru){
//int occ_i;
//vector<double> loc_avg_i_j(numEle), glo_avg_i_j(numEle);
vector<vector<double>> loc_avg_i_j, glo_avg_i_j;
loc_avg_i_j.resize(numEle,vector<double>(numEle,0));
glo_avg_i_j.resize(numEle,vector<double>(numEle,0));
Vector3d<double> offset,neighbour;
double e1,e2;
//double latt=3.57;
//cout <<"SRO_fcc, I reade latt="<<latt<<'\n';
double Box[3]={str.cell(0,0),str.cell(1,1),str.cell(2,2)};
//array<int,12> occ_j;
vector<int> occ_i(numEle,0);
vector<vector<int>> occ_j;
//occ_i.resize(numEle);
int N_neighbor = 12;
if (cryStru==1){N_neighbor=8; occ_j.resize(8,vector<int>(numEle,0));}
else if (cryStru==2){N_neighbor=12; occ_j.resize(12,vector<int>(numEle,0));}
for (int i=0; i<str.atom_pos.get_size(); i++) {
  e1=(str.atom_type[i]+1)%numEle;
  for (int i1=0;i1<numEle;i1++){occ_i[i1]=0;}
  occ_i[e1]=1;
  for (int j1=0;j1<N_neighbor;j1++){
  for (int i1=0;i1<numEle;i1++){occ_j[j1][i1]=0;}}
  for (int j=0; j<N_neighbor; j++){     
    if (cryStru==1) {offset=NNNbcc(latt,ith_NN,j);}
    else if (cryStru==2){offset=NNNfcc(latt,ith_NN,j);}
    for (int i2=0;i2<3;i2++){neighbour[i2] =str.atom_pos[i][i2]+offset[i2];}
    moveInBox(neighbour,Box);    
    int k=0;     
    while(k<str.atom_pos.get_size() and !matchVector(str.atom_pos[k],neighbour,Box)){ k+=1;}
    if (k==str.atom_pos.get_size()){cout<<"Warning: no match of neighbors found!"; exit;}
    e2=(str.atom_type[k]+1)%numEle;
    occ_j[j][e2]=1;
  }             
  for (int j=0;j<N_neighbor;j++){  
    for (int elem1=0;elem1<numEle;elem1++){        
    for (int elem2=0;elem2<numEle;elem2++){
      loc_avg_i_j[elem1][elem2] +=1.0/N_neighbor * occ_i[elem1]*occ_j[j][elem2];}}
  }               
}                 
for (int elem1=0;elem1<numEle;elem1++){        
for (int elem2=0;elem2<numEle;elem2++){
  glo_avg_i_j[elem1][elem2]= numEle*loc_avg_i_j[elem1][elem2]/str.atom_pos.get_size();}} //numEle*numEle*loc_avg_i_j[elem]/str.atom_pos.get_size()-1;
return glo_avg_i_j; //(glo_avg_i_j-0.25)/0.25;
}

int wlce::matchVector(Vector3d<double> atom1, Vector3d<double> atom2, double Box[3]){
 int boole[3];
 for (int i=0;i<3;i++){
   if (abs(atom1[i] -atom2[i])<1e-4 or abs(abs(atom1[i] -atom2[i])-Box[i])<1e-4) 
     {boole[i]=1;} 
   else {boole[i]=0;}}
 return boole[0]*boole[1]*boole[2];
}

void wlce::moveInBox(Vector3d<double> &Atom, double Box[3]){
for (int i=0;i<3;i++){
 while (Atom[i]<0) {Atom[i] += Box[i];}
 while (Atom[i]>Box[i]) {Atom[i] -= Box[i];}
}}

float wlce::cal_etot(Structure ideal_str,LinkedList<MultiCluster> clusterlist,SpaceGroup spacegroup, Array<Arrayint> labellookup, Array<Real> eci){

Real pred=0.;
  LinkedList<Array<MultiCluster> > eq_clusterlist;
  LinkedListIterator<MultiCluster> icluster(clusterlist);
  for ( ; icluster; icluster++) {  
    Array<MultiCluster> *pmulticlus=new Array<MultiCluster>;
    find_equivalent_clusters(pmulticlus, *icluster, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
    eq_clusterlist << pmulticlus;                                                                                                                   
  } 

char *corrfunc_label="trigo";
int multincl=0;
CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
pcorrfunc->init_from_site_type_list(labellookup);
LinkedListIterator<Array<MultiCluster> > icluster_eq(eq_clusterlist);
int ieci=0;
for ( ; icluster_eq; icluster_eq++, ieci++) {
  Real rho=calc_correlation(ideal_str, *icluster_eq, spacegroup.cell, *pcorrfunc);
//  if (strlen(ecifile)>0) {
    pred+=eci(ieci)*(multincl ? 1 : icluster_eq->get_size())*rho;
//  }        
//  else {   
//    if (!doconc) {cout << rho << delim;}
  }        
//if (strlen(ecifile)>0) {
  return float(pred);
//}
}

int wlce::genRandom(int range_from, int range_to){
std::random_device                  rand_dev;
std::mt19937                        generator(rand_dev());
std::uniform_int_distribution<int>  distr(range_from, range_to);
return distr(generator);
}

void wlce::writeLog(double gamma,vector<double> dos, vector<int> hist, vector<int> count, matrix sro, char *logfile){
ofstream logf(logfile);  
logf <<gamma<<'\n';                                                
logf <<"DOS:" <<'\n';                                                
for (int i = 0; i < sro.size(); i++){                                          
   logf<< dos[i] <<'\n';}                                                     
logf <<"HISTOGRAM:" <<'\n';                                          
for (int i = 0; i < sro.size(); i++){                                          
   logf<< hist[i] <<'\n';} 
logf <<"Short Range Parameter (sro):" <<'\n'; 
for (int i = 0; i < sro.size(); i++){
   for (int j=0; j < sro.begin()->size(); j++){                                          
     logf<< sro[i][j]/count[i]<<",";}
     logf<<count[i]<<'\n';}
}

void wlce::readLog(double gamma,vector<double> dos, vector<int> hist, vector<int> count, matrix sro, char *logfile){
const char * split = ",";
char * p;    
char ch[101];
int binSize=sro.size();
int numEle=sro.begin()->size();
dos.resize(binSize);
hist.resize(binSize);
count.resize(binSize);
fstream logf(logfile,ios::in);   
int row=0;
  if (!logf) {cout<<"Sorry bad file"; exit(0);}
  while (!logf.eof())    
   {                    
     row=row+1;         
     logf.getline(ch,100);
     //cout <<row<<" "<<ch<<endl;
     if (row==1) { p=strtok(ch,split);  gamma=atof(p); }
     if (row>2 and row<binSize+3) { p=strtok(ch,split);  dos[row-3]=atof(p); }
     if (row>binSize+3 and row<2*binSize+4) { p=strtok(ch,split);  hist[row-binSize-4]=atof(p); }
     if (row>2*binSize+4 and row<3*binSize+5) { 
        p=strtok(ch,split); int i=0;           
        while (i<numEle){                      
          sro[row-2*binSize-5][i]=atof(p);i++;} 
        count[row-2*binSize-5]=atoi(p);}
   }
}

// is flat?
int wlce::isFlat(vector<int> hist,int maskedN, int N_step, int maskedThreshold){
int binSize=hist.size();
vector<float> rest;
rest.resize(binSize);
float avg,maskedMax;
int maskedSize,sum;
maskedSize=0;
sum=0;
if (maskedN>500){
  for(int i = 0; i < binSize; i++){
    if (hist[i]>maskedThreshold){       
      sum += hist[i];
      maskedSize +=1;}}
  avg = float(sum) / maskedSize;
//  array<float,maskedSize> rest;
  for(int i = 0; i < binSize; i++){rest[i]=abs(hist[i]/avg-1);}
  maskedMax=abs(float(maskedThreshold)/avg-1);
//  cout<< maxArray(rest,maskedMax)<<" "<<maskedMax<<'\n';
  if (maxArray(rest,maskedMax)<0.5) {return 1;}
  else {return 0;} 
}
else{
  for(int i = 0; i < binSize; i++){sum += hist[i];}      
  avg = float(sum) / binSize;  
  for(int i = 0; i < binSize; i++){rest[i]=abs(hist[i]/avg-1);}
  if ((maxArray(rest,100)<0.5) and (avg>1)) {return 1;}
  else {return 0;} 
}
}

float wlce::maxArray(vector<float>rest, float maskedValue){
  int binSize=rest.size();
  float mx=0.0;
  for(int i=0;i<binSize;i++){
    if(mx<rest[i] and rest[i]<maskedValue){mx = rest[i];}}
  return mx;
}

// dos
int wlce::bin(double E, double lb, double ub, int binSize){
//return (E-lb)%((ub-lb)/binNum); //round((E-lb)/(ub-lb)*binNum);
double param=(E-lb)/(ub-lb)*binSize;
double fractpart, intpart;
fractpart = modf (param , &intpart);
return intpart;
}
//write structure
/*
void swap_write_struct(int atom_a, int atom_b, Structure &str, Array<AutoString> &atom_label, rMatrix3d &axes, char *file){
    Vector3d<double> mid;
    ofstream wfile(file); 
    mid=str.atom_pos(atom_a);
    str.atom_pos(atom_a)=str.atom_pos(atom_b);
    str.atom_pos(atom_b)=mid;
    write_structure(str,atom_label, axes, wfile);
}
*/
Structure wlce::swapped_str(int atom_a,int atom_b, Structure &str){
    Structure new_str=str;
    Vector3d<double> mid;
    mid=new_str.atom_pos(atom_a);
    new_str.atom_pos(atom_a)=new_str.atom_pos(atom_b);
    new_str.atom_pos(atom_b)=mid;
    return new_str;
}

void wlce::resizeVec( vector<vector<double> > &vec , const unsigned short ROWS , const unsigned short COLUMNS )
{                 
    vec.resize( ROWS );
    for( vector<vector<double> >::iterator it = vec.begin(); it != vec.end(); ++it)
    {             
        it->resize( COLUMNS );
    }             
} 
int wlce::countSameElement(Arrayint labels,int elementLabel){
  int count=0;
  for (int i=0;i<labels.get_size();i++){
    if (labels(i) == elementLabel){count +=1;}
  }
  return count;
} 
vector<int> wlce::random_pick_two_numbers(vector<int> components){
  int atom_a, atom_b;
  int component_a, component_b;
  int id_lb_a, id_lb_b, id_ub_a, id_ub_b;
  id_lb_a=0; id_lb_b=0; 
  int numEle = components.size();
  vector<int> two_int;
  component_a=genRandom(0,numEle-1);    
  component_b=genRandom(0,numEle-1);    
  while (component_a==component_b){
    component_b=genRandom(0,numEle-1);
  } 
  if (component_a==0) {id_ub_a=components[component_a]-1;}
  else {  
    for (int j=0; j<component_a; j++){    
      id_lb_a += components[j];}  
    id_ub_a =id_lb_a +components[component_a]-1;    
  }  
  if (component_b==0) {id_ub_b=components[component_b]-1;}
  else {                         
    for (int j=0; j<component_b; j++){    
      id_lb_b += components[j];}  
    id_ub_b =id_lb_b+components[component_b]-1;
  }  
  atom_a= genRandom(id_lb_a,id_ub_a);     
  atom_b= genRandom(id_lb_b,id_ub_b);     
  two_int.push_back(atom_a);
  two_int.push_back(atom_b);
  return two_int;
}

