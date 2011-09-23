/*crude monte carlo annealing*/
#include <set>
#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unistd.h>

using namespace std;

typedef vector<double> Coord;
typedef vector<Coord> Coords;
typedef vector<int> Label;

/*
 * routing2++.cppから、correspondだけを使う。単純に、格子点に近くない分子をさがす。


 *-dオプションは、1つの格子点に2もしくは0分子がassignされているケースを数える。
 *-vオプションは-dをVMRK形式で出力する。
 */


void
usage( string myname, string error, int code )
{
  cerr << "Error: " << error << endl;
  cerr << "Usage: " << myname << " reference.ar3a" << endl;
  cerr << "reference.ar3a Reference structure to be matched." << endl;
  exit(code);
}

void pbc( Coord& c, const Coord& box )
{
  for(int i=0;i<3;i++){
    c[i] -= floor(c[i] / box[i]) * box[i];
  }
}



double
ddistance( const Coord& a, const Coord& b, const Coord& box )
{
  double dd = 0;
  for(int j=0;j<3;j++){
    double x = a[j] - b[j];
    x -= rint( x / box[j] ) * box[j];
    dd += x*x;
  }
  return dd;
}



void
trial( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box, double beta )
{
  int n = fixed.size();
  int i = lrand48() % n;
  int j = lrand48() % n;
  if ( i != j ){
    double olddd = ddistance( fixed[i], atoms[label[i]], box ) + ddistance( fixed[j], atoms[label[j]], box );
    double newdd = ddistance( fixed[i], atoms[label[j]], box ) + ddistance( fixed[j], atoms[label[i]], box );
    if ( olddd > newdd  || exp(-beta*(newdd-olddd)) > drand48() ){
      int tmp = label[i];
      label[i] = label[j];
      label[j] = tmp;
      //cout << i << "x" << j << "|" << newdd-olddd << endl;
    }
  }
}



//液体分子位置に最も近い結晶格子のリストを返す。
Label
correspondence( const Coords& fixed, const Coords& atoms, Coord& box )
{
  Label correspond = Label(atoms.size());
  for(int i=0;i<atoms.size(); i++){
    int jmin = -1;
    double min = 1e10;
    for(int j=0;j<fixed.size(); j++){
      double dd = ddistance( atoms[i], fixed[j], box );
      if ( dd < min ){
	jmin = j;
	min = dd;
      }
    }
    correspond[i] = jmin;
    //cerr << "Corr " << i << "-->" << jmin << endl;
  }
  return correspond;
}



//全体の並進ずれを最小化する。
void
translation( const Coords& fixed, Coords& atoms, const Label& corr,Coord& box )
{
  for(int k=0;k<3;k++){
    double d = 0;
    for(int i=0;i<atoms.size(); i++){
      int j = corr[i];
      double x = atoms[i][k] - fixed[j][k];
      x -= rint( x / box[k] ) * box[k];
      d += x;
    }
    d /= atoms.size();
    for(int i=0;i<atoms.size(); i++){
      atoms[i][k] -= d;
    }
  }
}



//loop pathをみつけて、ひとつずらして解消する。
void
rondo( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box, const Label& correspond )
{
  Label done = Label(atoms.size());
  for(int i=0;i<atoms.size(); i++){
    Label members = Label(atoms.size());
    if ( done[i] == 0 ){
      int lattice = i;
      int size = 0;
      while(1){
	members[size] = lattice;
	int liquid = label[lattice];
	done[lattice] = 1;
	size ++;
	lattice = correspond[liquid];
	if ( lattice == i ){
	  break;
	}
	if ( done[lattice] ){
	  break;
	}
      }
      if ( lattice == i && size > 1 ){
	cerr << "RONDO " << endl;

	for(int m=0;m < size ;m++){
	  cerr << label[members[m]] << "-";
	}
	cerr << endl;

	int tmp = label[members[size-1]];
	for(int m=size-2; m>=0 ;m--){
	  label[members[m+1]] = label[members[m]];
	}
	label[members[0]] = tmp;

	for(int m=0;m < size ;m++){
	  cerr << label[members[m]] << "+";
	}
	cerr << endl;
      }
    }
  }	
}



//チェーンを見付け、末端同士がもし近ければチェーンを解消して末端を直結する。rondoの亜種。
void
linedance( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box, const Label& correspond )
{
  Label accept(atoms.size());
  Label donate(atoms.size());
  for(int i=0;i<atoms.size(); i++){
    int liquid = label[i];
    int lattice = correspond[liquid];
    if ( lattice != i ){
      accept[lattice] ++;
      donate[i] ++;
    }
  }
  /*
  for(int i=0;i<atoms.size(); i++){
    if ( accept[i] || donate[i] ){
      cerr << accept[i] << " " << donate[i] << " ";
      cerr << i << endl;
    }
  }
  */
  //acceptしておらず、donateしていれば、それが末端である。
  //環になっていれば末端はない。
  for(int i=0;i<atoms.size(); i++){
    //cerr << i<<":"<<atoms[i][0] << endl;
    if ( accept[i] == 0 && donate[i] != 0 ){
      Coord d(3);
      Label members(atoms.size());
      int size = 0;
      int head = i;
      double pathlen = 0;
      while(1){
	members[size] = head;
	size++;
	int next=correspond[label[head]];
	if ( next == head ){
	  break;
	}
        double sum = 0.0;
	for(int k=0;k<3;k++){
	  double x = atoms[next][k] - atoms[head][k];
	  x -= rint(x / box[k] ) * box[k];
          sum += x*x;
	  d[k] += x;
	}
        pathlen += sum;
        head = next;
      }
      bool warp=false;
      for(int k=0;k<3;k++){
	if ( rint(d[k] / box[k]) != 0.0 ){
	  warp = true;
	}
      }
      double sum = 0.0;
      for(int k=0;k<3;k++){
        double x = atoms[head][k] - atoms[i][k];
        sum += x*x;
      }
      //もしpathの両端が周期境界をまたいでいるなら
      //あるいは鎖の末端距離が、鎖の長さより短い場合はMak
      //実はあとの処理はrondo()と同じ。
      if ( warp || sum < pathlen ){
	cerr << "LINEDANCE" << endl;

	for(int m=0;m < size ;m++){
	  cerr << label[members[m]] << "-";
	}
	cerr << endl;

	int tmp = label[members[size-1]];
	for(int m=size-2; m>=0 ;m--){
	  label[members[m+1]] = label[members[m]];
	}
	label[members[0]] = tmp;
	for(int m=0;m < size ;m++){
	  cerr << label[members[m]] << "+";
	}
	cerr << endl;
      }
      else{
	cerr << "NO LINEDANCE" << endl;

	for(int m=0;m < size ;m++){
	  cerr << label[members[m]] << "-";
	}
	cerr << endl;
      }
    }
  }	
}




double
totalscore( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box )
{
  int n = fixed.size();
  double dd = 0.0;
  for(int i=0;i<n;i++){
    dd += ddistance( fixed[i], atoms[label[i]], box );
  }
  return dd;
}






Coords load_AR3A( istream& fin )
{
  Coords coord;
  string buf;
  getline( fin, buf );
  int nmol = atoi( buf.c_str() );
  for(int i=0; i<nmol; i++){
    getline( fin, buf );
    Coord pos(3);
    sscanf( buf.c_str(), "%lf %lf %lf", &pos[0], &pos[1], &pos[2] );
    coord.push_back( pos );
  }
  return coord;
}



int main(int argc, char* argv[] )
{
  Coords fixed;
  Coord box(3);
  string buf;
  srand48(getpid());
  //Label label;
  int arg = 1;
  //int outerloop = 20;
  //mode==0: @VMRK
  //mode==1: count off-lattice molecules
  //mode==2: count doubly-occupied/vacant voronoi cells
  //mode==3: Same as mode==2 but in VMRK format.
  //mode==4: 融解構造側ではなく参照構造側での分子番号でVMRK出力する。
  int mode = 0;
  if ( argc <= arg ){
    usage(string(argv[0]), "Milssing argument.", 1);
  }
  if ( string(argv[arg]) == "-c" ){
    arg ++;
    mode = 1;
  }
  else if ( string(argv[arg]) == "-d" ){
    arg ++;
    mode = 2;
  }
  else if ( string(argv[arg]) == "-v" ){
    arg ++;
    mode = 3;
  }
  else if ( string(argv[arg]) == "+v" ){
    arg ++;
    mode = 4;
  }
  ifstream probefile( argv[arg] );
  while( getline( probefile, buf )){
    if ( buf.size() > 4 && (
			    buf.substr( 0,5 ) == "@NX4A" ||
			    buf.substr( 0,5 ) == "@AR3A" )){
      fixed = load_AR3A( probefile );
      break;
    }
  }
  while( getline( cin, buf ) ){
    if ( buf.size() > 4 && buf.substr( 0,5 ) == "@BOX3" ){
      getline( cin, buf );
      sscanf( buf.c_str(), "%lf %lf %lf", &box[0], &box[1], &box[2] );
    }
    else if ( buf.size() > 4 && (
				 buf.substr( 0,5 ) == "@NX4A" ||
				 buf.substr( 0,5 ) == "@AR3A" )){
      Coords atoms = load_AR3A( cin );
      int n = atoms.size();
      if ( mode == 0 ){
        cout << "@VMRK" << endl;
        cout << n << endl;
      }
      Label correspond = correspondence( fixed, atoms, box );
      Label location = Label(atoms.size());
      int count = 0;
      for(int i=0;i<n; i++){
	double dd = ddistance( atoms[i], fixed[correspond[i]], box );
	if ( mode == 1 )
	  count += (dd>0.7*0.7);
        if ( mode == 0 ){
          cout << (dd>0.7*0.7) << " " << dd << endl;
        }
	if ( mode == 2 || mode == 3 || mode == 4 ){
	  location[correspond[i]]++;
	}
      }
      if ( mode == 2 ){
	count = 0;
	for ( int i=0; i<n; i++ ){
	  if ( location[i] != 1 ){
	    count ++;
	  }
	}
      }
      if ( mode ==1 || mode == 2 ){
        cout << count << " " << n << endl;
      }
      else if ( mode == 3 ){
	cout << "@VMRK" << endl;
	cout << n << endl;
	for( int i=0; i<n; i++ ){
	  cout << location[correspond[i]] << endl;
	}
      }
      else if ( mode == 4 ){
	cout << "@VMRK" << endl;
	cout << n << endl;
	for( int i=0; i<n; i++ ){
	  cout << location[i] << endl;
	}
      }
    }
  }
}
