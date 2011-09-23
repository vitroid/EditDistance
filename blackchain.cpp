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
  yellowchainと同じように、-yオプションで黒鎖の末端に印を描く。
  末端の個数を数える必要がある。欠陥濃度の推定。
  routing4と同じように、VMRKはオプションに変更。
*/

void
usage( string myname, string error, int code )
{
  cerr << "Error: " << error << endl;
  cerr << "Usage: " << myname << " [-c] reference.ar3a" << endl;
  cerr << " -c       Count the chain length (default: output sum square displacement.)" << endl;
  exit(code);
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




/*黒鎖は分岐の心配はないはず。*/
string
blackchainlength( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box, const Label& correspond, int mode )
{
  string buf;
  stringstream sout;
  Label accept(atoms.size());
  Label donate(atoms.size());
  for(int i=0;i<atoms.size(); i++){
    int liquid = label[i];
    int lattice = correspond[liquid];
    if ( lattice != i ){
      accept[lattice] ++;
      donate[i] ++;
      //cerr << i << ":" << lattice << endl;
    }
  }
  //acceptしておらず、donateしていれば、それが末端である。
  //環になっていれば末端はない。
  double maxpathlen = 0;
  for(int i=0;i<atoms.size(); i++){
    //cerr << i<<":"<<atoms[i][0] << endl;
    if ( accept[i] == 0 && donate[i] != 0 ){
      Label members(atoms.size());
      Label done(atoms.size());
      int size = 0;
      int head = i;
      double pathlen = 0;
      while(1){
	members[size] = head;
        done[head] = 1;
	size++;
	int next=correspond[label[head]];
	if ( next == head ){
	  break;
	}
        //avoid loop
        if ( done[next] ){
          size = 0;
          break;
        }
        if ( mode == 2 || mode == 3 || mode == 4 ){
          pathlen += 1;
        }
        else{
          double sum = 0.0;
          for(int k=0;k<3;k++){
            double x = atoms[next][k] - atoms[head][k];
            x -= rint(x / box[k] ) * box[k];
            sum += x*x;
          }
          pathlen += sum;
        }
        head = next;
      }
      if ( mode == 3 && size > 1 ){
	sout << "@ 0" << endl;
	sout << "y 1" << endl;
	sout << "c " << fixed[members[0]][0] << " " << fixed[members[0]][1] << " " << fixed[members[0]][2] << endl;
	sout << "@ 2" << endl;
	sout << "y 2" << endl;
	sout << "c " << atoms[members[size-1]][0] << " " << atoms[members[size-1]][1] << " " << atoms[members[size-1]][2] << endl;
      }
      else if ( mode == 4 ){
        maxpathlen += pathlen;
      }
      else {
	if ( maxpathlen < pathlen ){
	  maxpathlen = pathlen;
	  //cerr << maxpathlen << endl;
	  stringstream tmp;
	  tmp << maxpathlen;
	  buf = tmp.str();
	  //char cbuf[100];
	  //sprintf(cbuf, "%f", maxpathlen);
	  //buf = string(cbuf);
	}
      }
    }
  }
  if ( mode == 4 ){
    stringstream tmp;
    tmp << maxpathlen;
    buf = tmp.str();
  }
  //cerr << buf << endl;
  if ( buf == "" ){
    buf = sout.str();
  }
  return buf;
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
  Label label;
  int arg = 1;
  int outerloop = 20;
  int mode = 1; // 1: sum square displacement;  2: chain count; 3: yaplot
  while( arg < argc && argv[arg][0] == '-' ){
    buf = string(argv[arg]);
    if ( buf == "-c" ){
      mode = 2;
      arg++;
    }
    if ( buf == "-n" ){
      //黒チェインの総数を数える。(長いものに限定しない)
      mode = 4;
      arg++;
    }
    else if ( buf == "-y" ){
      mode = 3;
      arg++;
    }
    else if ( buf == "-v" ){
      arg++;
      if ( argc <= arg ){
	usage(string(argv[0]), "Missing argument.", 1);
      }
      ifstream labelfile( argv[arg] );
      while( getline( labelfile, buf ) ){
	if ( buf.size() > 4 && buf.substr( 0,5 ) == "@VMRK" ){
	  getline( labelfile, buf );
	  int nmol = atoi( buf.c_str() );
	  label.resize(nmol);
	  for(int i=0; i<nmol; i++){
	    getline( labelfile, buf );
	    Coord pos(3);
	    sscanf( buf.c_str(), "%d", &label[i] );
	  }
	}
      }
      arg++;
    }
  }
  if ( argc <= arg ){
    usage(string(argv[0]), "Missing argument.", 1);
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
      if ( label.size() == 0 ){
	label.resize(n);
	for(int i=0;i<n; i++){
	  label[i] = i;
	}
      }
      //atomsをbox内におさめる。
      for(int i=0;i<n;i++){
	atoms[i][0] -= rint(atoms[i][0] / box[0]) * box[0];
	atoms[i][1] -= rint(atoms[i][1] / box[1]) * box[1];
	atoms[i][2] -= rint(atoms[i][2] / box[2]) * box[2];
	fixed[i][0] -= rint(fixed[i][0] / box[0]) * box[0];
	fixed[i][1] -= rint(fixed[i][1] / box[1]) * box[1];
	fixed[i][2] -= rint(fixed[i][2] / box[2]) * box[2];
      }
      //結晶のi番目の格子点に一番近い液体分子の番号を記す。
      Label correspond = correspondence( fixed, atoms, box );
      cout << blackchainlength( fixed, atoms, label, box, correspond, mode ) << endl;
    }
    else if ( buf.size() > 4 && buf.substr( 0,5 ) == "@VMRK" ){
      cerr << "WARNING VMRK from stdin is just ignored." << endl;
    }
  }
}
