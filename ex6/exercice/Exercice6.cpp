#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "ConfigFile.tpp"


using namespace std;

void coutBigFatVec(vector<double> const& vec, string name);

// Resolution d'un systeme d'equations lineaires par elimination de Gauss-Jordan:
template <class T>
vector<T> solve(vector<T> const& diag,
                vector<T> const& lower,
                vector<T> const& upper,
                vector<T> const& rhs)
{
  vector<T> solution(diag.size());
  vector<T> new_diag(diag);
  vector<T> new_rhs(rhs);

  for(int i(1); i<diag.size(); ++i)
  {
    double pivot = lower[i-1]/new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  for(int i(diag.size()-2); i>=0; --i)
    solution[i] = (new_rhs[i] - upper[i]*solution[i+1]) / new_diag[i];

  return solution;
}


// Classe pour epsilon_r(r)
class Epsilonr {
public:
  Epsilonr(bool const& trivial_, double const& b_, double const& c_)
    : b(b_), R(c_), trivial(trivial_) {};

  inline double operator()(double const& r, bool const& left) {
  // Le booleen "left" indique s'il faut prendre la limite a gauche ou a droite en cas de discontinuite
    double eps(1e-12*b);
    if(trivial or r<=b-eps or (abs(r-b)<=eps and left))
      return 1.0;
    else
      return 8.0 - 6.0*(r-b)/(R-b);
  }

private:
  double b, R;
  bool trivial;
};

// Classe pour rho_lib(r)/epsilon_0
class Rho_lib {
public:
  Rho_lib(bool const& trivial_, double const& b_, double const& a0_)
    : b(b_), a0(a0_), trivial(trivial_) {};

  inline double operator()(double const& r) {
    // if(trivial or r>b)
    //   return 1.0;
    // else
    //   return a0*(1.0-pow(r/b,2));
    if(trivial){
      return 1.0;
    }
    else{
      if(r>b)
        return 0.0;
      else
        return a0*(1.0 - pow(r/b,2));
        // return a0*sin(M_PI*(1.0 - pow(r/b, 2)));
    }
  }

private:
  double b, a0;
  bool trivial;
};

int main(int argc, char* argv[])
{

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice6 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice6 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Fichier de sortie :
  string output = configFile.get<string>("output");

  // Domaine :
  const double b(configFile.get<double>("b"));
  const double R(configFile.get<double>("R"));

  // Conditions aux bords :
  const double V0(configFile.get<double>("V0"));

  // Instanciation des objets :
  Epsilonr epsilonr(configFile.get<bool>("trivial"), b, R);
  Rho_lib rho_lib(configFile.get<bool>("trivial"), b, configFile.get<double>("a0"));

  // Facultatif
  double p(configFile.get<double>("p"));

  // Discretisation du domaine :
  int N1 = configFile.get<int>("N1");
  int N2 = configFile.get<int>("N2");
  int ninters = N1 + N2;
  int npoints = ninters + 1;
  double h1 = b/N1;
  double h2 = (R-b)/N2;
  vector<double> r(npoints);
  for(int i(0); i<N1; ++i)
    r[i] = i*h1;
  for(int i(0); i<=N2; ++i)
    r[N1+i] = b + i*h2;
  vector<double> h(ninters);
  for(int i(0); i<ninters; ++i)
    h[i] = r[i+1] - r[i];

  vector<double> diag(npoints,0.);  // Diagonale
  vector<double> lower(ninters,0.); // Diagonale inferieure
  vector<double> upper(ninters,0.); // Diagonale superieure
  vector<double> rhs(npoints,0.);   // Membre de droite

  //Matrice A
  //diagonale
  diag[0]=1./(2.*h[0])*(p*(r[1]*epsilonr(r[1], true)+r[0]*epsilonr(r[0], true))+(1-p)*(r[0]+r[1])*epsilonr((r[0]+r[1])/2.,true));
  for (size_t i(1); i<diag.size()-1; ++i){
    // if (i==N1-1){
    //   diag[i] = 1./(2.*h[i])*(r[i+1]*epsilonr(r[i+1], false)+r[i]*epsilonr(r[i], true)) + 1./(2.*h[i-1])*(r[i-1]*epsilonr(r[i-1], true)+r[i]*epsilonr(r[i], true));
    // } else
    if (i==N1){
      diag[i] = 1./(2.*h[i])*(p*(r[i+1]*epsilonr(r[i+1], true)+r[i]*epsilonr(r[i], false))+(1-p)*(r[i]+r[i+1])*epsilonr((r[i]+r[i+1])/2.,true)) + 1./(2.*h[i-1])*(p*(r[i-1]*epsilonr(r[i-1], true)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i-1])*epsilonr((r[i]+r[i-1])/2.,true));
    } else if (i==N1+1){
      diag[i] = 1./(2.*h[i])*(p*(r[i+1]*epsilonr(r[i+1], true)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i+1])*epsilonr((r[i]+r[i+1])/2.,true)) + 1./(2.*h[i-1])*(p*(r[i-1]*epsilonr(r[i-1], false)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i-1])*epsilonr((r[i]+r[i-1])/2.,true));
    } else {
      diag[i] = 1./(2.*h[i])*(p*(r[i+1]*epsilonr(r[i+1], true)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i+1])*epsilonr((r[i]+r[i+1])/2.,true)) + 1./(2.*h[i-1])*(p*(r[i-1]*epsilonr(r[i-1], true)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i-1])*epsilonr((r[i]+r[i-1])/2.,true));
    }
  };
  //à rk+1=b:
  // diag[N1-1] = 1./(2.*h[N1-1])*(r[N1]*epsilonr(r[N1], false)+r[N1-1]*epsilonr(r[N1-1], true)) + 1./(2.*h[N1-2])*(r[N1-2]*epsilonr(r[N1-2], true)+r[N1-1]*epsilonr(r[N1-1], true));
  // diag[N1] = 1./(2.*h[N1])*(r[N1+1]*epsilonr(r[N1+1], false)+r[N1]*epsilonr(r[N1], true)) + 1./(2.*h[N1-2])*(r[N1-2]*epsilonr(r[N1-2], true)+r[N1-1]*epsilonr(r[N1-1], true));
  // for (size_t i(N1+1); i<diag.size()-1; ++i){
  //   diag[i] = 1./(2.*h[i])*(r[i+1]*epsilonr(r[i+1], false)+r[i]*epsilonr(r[i], false)) + 1./(2.*h[i-1])*(r[i-1]*epsilonr(r[i-1], false)+r[i]*epsilonr(r[i], false));
  // }

  //sous et sur-diagonale
  for (size_t i(0); i<lower.size(); ++i){
    // if (i==N1-1){
    //   lower[i] = -1./(2.*h[i])*(r[i+1]*epsilonr(r[i+1], false)+r[i]*epsilonr(r[i], true));
    // } else
    if (i==N1){
      lower[i] = -1./(2.*h[i])*(p*(r[i+1]*epsilonr(r[i+1], true)+r[i]*epsilonr(r[i], false))+(1-p)*(r[i]+r[i+1])*epsilonr((r[i]+r[i+1])/2.,true));
    } else {
      lower[i] = -1./(2.*h[i])*(p*(r[i+1]*epsilonr(r[i+1], true)+r[i]*epsilonr(r[i], true))+(1-p)*(r[i]+r[i+1])*epsilonr((r[i]+r[i+1])/2.,true));
    }
    upper[i] = lower[i];
  }

  //à rk+1=b:
  // lower[N1-1] = -1./(2.*h[N1-1])*(r[N1]*epsilonr(r[N1], false)+r[N1-1]*epsilonr(r[N1-1], true));
  // upper[N1-1] = lower[N1-1];
  // for (size_t i(N1); i<lower.size(), ++i){
  //   lower[i] = -1./(2.*h[i])*(r[i+1]*epsilonr(r[i+1], false)+r[i]*epsilonr(r[i], false));
  //   upper[i] = lower[i];
  // }

  // sans epsilon_0
  rhs[0] = .5*(p*r[0]*rho_lib(r[0])*h[0]+(1-p)*h[0]*.5*(r[0]+r[1])*rho_lib(.5*(r[0]+r[1])));
  for (size_t i(1); i<rhs.size()-1; ++i){
    rhs[i] = .5*(p*r[i]*rho_lib(r[i])*(h[i-1]+h[i])+(1-p)*(h[i]*.5*(r[i]+r[i+1])*rho_lib(.5*(r[i]+r[i+1]))+h[i-1]*.5*(r[i]+r[i-1])*rho_lib(.5*(r[i]+r[i-1]))));
  }

  //Condition au bord
  diag[diag.size()-1]=1;
  lower[lower.size()-1]=0;
  rhs[rhs.size()-1]=V0;

  // coutBigFatVec(diag,"diag");
  // coutBigFatVec(lower,"lower");
  // coutBigFatVec(rhs,"rhs");

  //couts
  // for (size_t i=0; i<diag.size(); ++i){
  //   cout << diag[i] << " ";
  // }
  // cout << endl;

  // Resolution:
  vector<double> phi(solve(diag,lower,upper,rhs));

  // Export des resultats:
  // 1. phi
  ofstream ofs((output+"_phi.out").c_str());
  ofs.precision(15);
  for(int i(0); i<npoints; ++i)
    ofs << r[i] << " " << phi[i] << endl;
  ofs.close();

  // 2. E_r et D_r
  vector<double> rmid(ninters);
  vector<double> Er(ninters);
  vector<double> Dr(ninters);
  for(int i(0); i<ninters; ++i)
  {
    rmid[i] = 0.5*r[i] + 0.5*r[i+1];
    Er[i] = -(phi[i+1]-phi[i])/h[i];
    Dr[i] = epsilonr(rmid[i],true)*Er[i];
  }
  coutBigFatVec(Er, "Er");
  coutBigFatVec(Dr, "Dr");
  ofs.open((output+"_Er_Dr.out").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters; ++i)
    ofs << rmid[i] << " " << Er[i] << " " << Dr[i] << endl;
  ofs.close();

  // 3. rho_lib, div(E_r) et div(D_r)
  vector<double> rmidmid(ninters-1);
  vector<double> div_Er(ninters-1);
  vector<double> div_Dr(ninters-1);
  for(int i(0); i<ninters-1; ++i)
  {
    rmidmid[i] = 0.5*rmid[i] + 0.5*rmid[i+1];
    div_Er[i] = (rmid[i+1]*Er[i+1]-rmid[i]*Er[i])/(rmid[i+1]-rmid[i])/rmidmid[i]; // DF en cylindrique (d'où le 1/r)
    div_Dr[i] = (rmid[i+1]*Dr[i+1]-rmid[i]*Dr[i])/(rmid[i+1]-rmid[i])/rmidmid[i];//epsilonr(rmidmid[i],true)*div_Er[i];
  }
  ofs.open((output+"_rholib_divEr_divDr.out").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters-1; ++i)
    ofs << rmidmid[i] << " " << rho_lib(rmidmid[i]) << " " << div_Er[i] << " " << div_Dr[i] << endl;
  ofs.close();

  return 0;
}

void coutBigFatVec(vector<double> const& vec, string name){
  cout << name << "=( ";
  for (size_t i = 0; i < vec.size() - 1; i++) {
    cout.precision(5);
    cout << vec[i] << ", ";
  }
  cout << vec[vec.size() - 1] << " )" << endl;
}
