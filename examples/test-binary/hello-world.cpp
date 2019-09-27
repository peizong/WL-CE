// my first program in C++
#include <iostream>
#include <vector>

void resizeVec( std::vector<std::vector<double> > &vec , const unsigned short ROWS , const unsigned short COLUMNS );

int main()
{
using namespace std;
std::cout << "Hello World!";
typedef std::vector< std::vector<double> > matrix;
matrix name(2, std::vector<double>(2));
resizeVec(name,3,3);
for (int i=0;i<3;i++){
  for (int j=0;j<3;j++){
    name[i][j]=i*3+j;
    cout << name[i][j] <<"\n";
  }
}
}

void resizeVec( std::vector<std::vector<double> > &vec , const unsigned short ROWS , const unsigned short COLUMNS )
{
    vec.resize( ROWS );
    for( std::vector<std::vector<double> >::iterator it = vec.begin(); it != vec.end(); ++it)
    {
        it->resize( COLUMNS );
    }
}
