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
#include <cstdlib>

using namespace std;

typedef vector<double> Coord;
typedef vector<Coord> Coords;
typedef vector<int> Label;

/*
 * routing2.cpp�˲ä������Ū�˴Ĥ�õ����Ԥ������Υ��饹�����르�ꥺ��ȸ��äƤ�褤��
 * �Ĥ����դ���Ф��餹���Ȥǵ�®�˲����塼�˥󥰤��롣
 */


void
usage( string myname, string error, int code )
{
  cerr << "Error: " << error << endl;
  cerr << "Usage: " << myname << " [-y][-l loop] reference.ar3a" << endl;
  cerr << " -y       Output in yaplot format instead of @AR3A format." << endl;
  cerr << " -l loop  Specify the number of annealing/tempering cycles.(default=20)" << endl;
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



//����ʬ�Ұ��֤˺Ǥ�ᤤ�뾽�ʻҤΥꥹ�Ȥ��֤���
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



//���Τ��¿ʤ����Ǿ������롣
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



//loop path��ߤĤ��ơ��ҤȤĤ��餷�Ʋ�ä��롣
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



//����������դ�����üƱ�Τ��⤷�ᤱ��Х���������ä�����ü��ľ�뤹�롣rondo�ΰ��
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
  //accept���Ƥ��餺��donate���Ƥ���С����줬��ü�Ǥ��롣
  //�ĤˤʤäƤ������ü�Ϥʤ���
  for(int i=0;i<atoms.size(); i++){
    //cerr << i<<":"<<atoms[i][0] << endl;
    if ( accept[i] == 0 && donate[i] != 0 ){
      Coord d(3);
      Label members(atoms.size());
      Label done(atoms.size());
      int size = 0;
      int first = i;
      while(1){
	members[size] = first;
        done[first] = 1;
	size++;
	int next=correspond[label[first]];
        if ( size > 20 ){
          cerr << size << "("<<first<<":"<<next<<")";
        }
	if ( next == first ){
	  break;
	}
        //avoid loop
        if ( done[next] ){
          size = 0;
          break;
        }
	for(int k=0;k<3;k++){
	  double x = atoms[next][k] - atoms[first][k];
	  x -= rint(x / box[k] ) * box[k];
	  d[k] += x;
	}
        first = next;
      }
      
      bool warp=false;
      for(int k=0;k<3;k++){
	if ( rint(d[k] / box[k]) != 0.0 ){
	  warp = true;
	}
      }
      //�⤷path��ξü������������ޤ����Ǥ���ʤ�
      //�¤Ϥ��Ȥν�����rondo()��Ʊ����
      if ( warp ){
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




bool yaplot = 0;


string
output( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box, double total)
{
  stringstream sout;
  int n=atoms.size();
  if ( yaplot ){
    sout << "# " << total << endl;
    sout << "r 0.1" << endl;
    for(int i=0; i<n; i++){
      Coord delta = fixed[i];
      for(int j=0;j<3;j++){
	delta[j] -= rint( delta[j] / box[j] ) * box[j];
      }
      //line start
      sout << "y 1" << endl;
      sout << "l";
      for(int j=0;j<3;j++){
	sout << " " << delta[j];
      }

      for(int j=0;j<3;j++){
	double d = atoms[label[i]][j] - fixed[i][j];
	d -= rint( d / box[j] ) * box[j];
        delta[j] += d;
      }
      for(int j=0;j<3;j++){
	sout << " " << delta[j];
      }
      sout << endl;
    }
    sout << endl;
  }
  else{
    sout << "@ETOT" << endl;
    sout << total << endl;
    sout << "@BOX3" << endl;
    sout << box[0] << " " << box[1] << " " << box[2] << endl;
    sout << "@AR3A" << endl << n << endl;
    for(int i=0; i<n; i++){
      for(int j=0;j<3;j++){
	sout << atoms[label[i]][j] << " ";
      }
      sout << endl;
    }
    sout << "@VMRK" << endl << n << endl;
    for(int i=0; i<n; i++){
      sout << label[i]<< endl;
    }
  }
  return sout.str();
}


double mine=100000;
string last;

double
full_optimize( const Coords& fixed, const Coords& atoms, Label& label, const Coord& box )
{
  double beta = 3.0;
  int loop = 50000;
  /*tempering*/
  while( beta > 0.2 ){
    double total = totalscore(fixed, atoms, label,box );
    if ( total < mine ){
      mine = total;
      last = output( fixed, atoms, label, box, total );
      cout << last << flush;
    }
    cerr << beta << " " << total << endl;
    for(int i=0;i<loop;i++){
      trial( fixed, atoms, label, box, beta );
    }
    beta -= 0.1;
  }
  /*annealing*/
  while( beta < 3.0 ){
    double total = totalscore(fixed, atoms, label,box );
    if ( total < mine ){
      mine = total;
      last = output( fixed, atoms, label, box, total );
      cout << last << flush;
    }
    cerr << beta << " " << total << endl;
    for(int i=0;i<loop;i++){
      trial( fixed, atoms, label, box, beta );
    }
    beta += 0.1;
  }
  return mine;
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
  Coords atoms;
  Coord box(3);
  string buf;
  srand48(getpid());
  Label label;
  int arg = 1;
  int outerloop = 20;
  while( arg < argc && argv[arg][0] == '-' ){
    buf = string(argv[arg]);
    if ( buf == "-y" ){
      yaplot = 1;
      arg++;
    }
    else if ( buf == "-l" ){
      arg++;
      if ( argc <= arg ){
	usage(string(argv[0]), "Missing argument.", 1);
      }
      outerloop = atoi(argv[arg]);
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
      atoms = load_AR3A( cin );
    }
    else if ( buf.size() > 4 && buf.substr( 0,5 ) == "@VMRK" ){
      getline( cin, buf );
      int nmol = atoi( buf.c_str() );
      label.resize(nmol);
      for(int i=0; i<nmol; i++){
        getline( cin, buf );
        Coord pos(3);
        sscanf( buf.c_str(), "%d", &label[i] );
      }
    }
  }
  int n = atoms.size();
  if ( label.size() == 0 ){
    label.resize(n);
    for(int i=0;i<n; i++){
      label[i] = i;
    }
  }
  //�뾽��i���ܤγʻ����˰��ֶᤤ����ʬ�Ҥ��ֹ�򵭤���
  Label correspond = correspondence( fixed, atoms, box );
  if ( outerloop > 0 ){
    for(int loop=0;loop<outerloop;loop++){
      full_optimize( fixed, atoms, label, box );
      rondo( fixed, atoms, label, box, correspond );
      linedance( fixed, atoms, label, box, correspond );
      translation( fixed, atoms, correspond, box );
    }
  }
  else{
    //loop������ͤ����ꤵ�줿���ϡ�-loop�󡢺Ǿ��ͤ���������ʤ����˽�λ���롣
    //-5�����꤬��������
    double last = 0;
    int count = -outerloop;
    while( count ){
      double min = full_optimize( fixed, atoms, label, box );
      rondo( fixed, atoms, label, box, correspond );
      linedance( fixed, atoms, label, box, correspond );
      translation( fixed, atoms, correspond, box );
      if ( min != last ){
        count = -outerloop;
        last = min;
      }
      else{
        count --;
      }
    }
  }
  cout << "#####" << endl;
  cout << last << flush;
}
