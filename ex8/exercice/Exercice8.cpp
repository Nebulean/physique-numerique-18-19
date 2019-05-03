#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.tpp"

using namespace std;
typedef vector<complex<double> > vec_cmplx;

template <class T>
void coutBigFatVec(vector<T> const& vec, string name);


// Fonction resolvant le systeme d'equations A * solution = rhs
// ou A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution)
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // forward elimination
  for(int i(1); i<diag.size(); ++i)
  {
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution.resize(diag.size());

  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  // backward substitution
  for(int i = diag.size() - 2; i >= 0; --i)
  {
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) :
double V(double const& x, double const& omega, double const& delta)
{
  return .5*omega*omega*min((x-delta)*(x-delta),(x+delta)*(x+delta));
}


// Declaration des diagnostics de la particule d'apres sa fonction d'onde psi :
//  - prob calcule la probabilite de trouver la particule entre les points nL.dx et nR.dx,
//  - E calcule son energie,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carre moyenne,
//  - pmoy calcule sa quantite de mouvement moyenne,
//  - p2moy calcule sa quantite de mouvement au carre moyenne.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);
double delx(vec_cmplx const& psi, vector<double> const& x, double const& dx);
double delp(vec_cmplx const& psi, double const& dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.


int main(int argc,char **argv)
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres physiques :
  double hbar    = 1.;
  double m       = 1.;
  double tfin    = configFile.get<double>("tfin");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  double omega   = configFile.get<double>("omega");
  double delta   = configFile.get<double>("delta");
  double x0      = configFile.get<double>("x0");
  double k0      = 2. * M_PI * configFile.get<int>("n") / (xR-xL);
  double sigma0  = configFile.get<double>("sigma_norm") * (xR-xL);

  // Parametres numeriques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  for(int i(0); i<Npoints; ++i)
    psi[i] = exp(complex_i*k0*x[i])*exp(-(x[i]-x0)*(x[i]-x0)/(2*sigma0*sigma0));
  // Modifications des valeurs aux bords :
  psi[0] = complex<double> (0.,0.);
  psi[Npoints-1] = complex<double> (0.,0.);
  // Normalisation :
  psi = normalize(psi, dx);

  // Incertitudes de la position et qu. de mouvement
  // double delx(sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)));
  // double delp(sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)));
  // cout << p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx) << endl;
  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'equation (4.90)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'equation (4.90)

  complex<double> a, b;

  // Ces matrices sont stockees sous forme tridiagonale, d:diagonale, c et a: diagonales superieures et inferieures
  for (size_t i = 0; i < Npoints; i++) {
    dH[i] = pow(hbar, 2)/( m*pow(dx, 2) ) + V(x[i], omega, delta);
    dA[i] = 1.0 + complex_i*dt/(2.0*pow(hbar, 2)) * dH[i];
    dB[i] = 1.0 - complex_i*dt/(2.0*pow(hbar, 2)) * dH[i];
  }
  for (size_t i = 0; i < Ninters; i++) {
    aH[i] = -pow(hbar, 2)/(2.0*m*pow(dx, 2));
    cH[i] = -pow(hbar, 2)/(2.0*m*pow(dx, 2));
    aA[i] = complex_i*dt/(2.0*pow(hbar, 2)) * aH[i];
    cA[i] = complex_i*dt/(2.0*pow(hbar, 2)) * cH[i];
    aB[i] = -complex_i*dt/(2.0*pow(hbar, 2)) * aH[i];
    cB[i] = -complex_i*dt/(2.0*pow(hbar, 2)) * cH[i];
  }

  // Conditions aux limites: psi nulle aux deux bords
  // Bord gauche
  dA[0] = 1;
  cA[0] = 0;
  cB[0] = 0;
  // cout << "size"; //TODO : Pourquoi y'a pas de seg. fault si on change les bord droits avec +10 au lieu de -1 par exemple ?
  // Bord droit
  dA[Npoints-1] = 1;
  aA[Ninters-1] = 0;
  aB[Ninters-1] = 0;

  // cout << "size2" << dA.size() << endl;

  // coutBigFatVec(aA, "dA");
  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_potentiel((output + "_pot.out").c_str());
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega, delta) << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((output + "_psi2.out").c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((output + "_obs.out").c_str());
  fichier_observables.precision(15);

  // Boucle temporelle :
  double t;
  for(t=0.; t+dt/2.<tfin; t+=dt)
  {
    // Ecriture de |psi|^2 :
    for(int i(0); i<Npoints; ++i)
      fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";
    fichier_psi << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,Ninters*xL/(xL-xR),dx) << " "       // probabilite que la particule soit en x < 0
                        << prob(psi,Ninters*xL/(xL-xR),Ninters,dx) << " " // probabilite que la particule soit en x > 0
                        << E(psi,dH,aH,cH,dx) << " "                      // Energie
                        << xmoy(psi,x,dx) << " "                          // Position moyenne
                        << x2moy(psi,x,dx) << " "                         // Position^2 moyenne
                        << pmoy(psi,dx) << " "                            // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << " "                         // (Quantite de mouvement)^2 moyenne
                        << delx(psi, x, dx) << " "
                        << delp(psi, dx) << endl;

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints,0.);

    // Multiplication psi_tmp = B * psi :
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i)
    {
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }

    // Resolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);

  } // Fin de la boucle temporelle

  for(int i(0); i<Npoints; ++i)
    fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";

  fichier_observables << t << " "
                      << prob(psi,0,Ninters*xL/(xL-xR),dx) << " " // Proba pour la partie droite
                      << prob(psi,Ninters*xL/(xL-xR),Ninters,dx) << " " // proba pour la partie gauche.
                      << E(psi,dH,aH,cH,dx) << " "
                      << xmoy(psi,x,dx) << " "
                      << x2moy(psi,x,dx) << " "
                      << pmoy(psi,dx) << " "
                      << p2moy(psi,dx) << " "                         // (Quantite de mouvement)^2 moyenne
                      << delx(psi, x, dx) << " "
                      << delp(psi, dx) << endl;

  fichier_observables.close();
  fichier_psi.close();

}


double prob(vec_cmplx const& psi, int nL, int nR, double dx)
{
  //calcule la probabilite de trouver la particule entre les points nL.dx et nR.dx
  double pr(0.);
  //max nR = Npoints-2
  for (size_t i = nL; i < nR; i++) {
    pr+=abs(psi[i])*abs(psi[i])+abs(psi[i+1])*abs(psi[i+1]);
  }
  pr*=dx/2.;
  return pr;
}


double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx)
{
  // vec_cmplx psi_tmp(psi.size());
  complex<double> E(0.0, 0.0);
  size_t Npoints = psi.size();

  // TODO: calculer la moyenne de l'Hamiltonien

  // H(psi)
  // On utilise la matrice H calculée plus haut
  vec_cmplx psiH(Npoints, 0.0); // H*psi
  for (size_t i = 1; i < Npoints - 1; i++) {
    psiH[i] = lowerH[i-1]*psi[i-1] + diagH[i]*psi[i] + upperH[i]*psi[i+1];
  }
  // Bords
  psiH[0] = diagH[0]*psi[0] + upperH[0]*psi[1];
  psiH[Npoints - 1] = lowerH[Npoints - 2]*psi[Npoints - 2] + diagH[Npoints - 1]*psi[Npoints - 1];

  // coutBigFatVec(psi, "psi");

  // Integrale de psi* H(psi) dx
  for (size_t i = 0; i < Npoints - 1; i++) {
    E += conj(psi[i])*psiH[i] + conj(psi[i+1])*psiH[i+1];
  }
  E *= dx/2.0;
  // cout << "E = " << real(E) << " " << imag(E) << endl;
  return real(E);
}


double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  //calcule la moyenne de la position
  complex<double> res(0.0, 0.0);
  for (size_t i = 0; i < psi.size()-1; i++) {
    res+=psi[i]*x[i]*conj(psi[i])+psi[i+1]*x[i+1]*conj(psi[i+1]);
  }
  // cout << "res = " << real(res) << " " << imag(res) << endl;
  return real(res*dx/2.0);
}


double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  //calcule la moyenne du x^2
  complex<double> res(0.0, 0.0);
  for (size_t i = 0; i < psi.size()-1; i++) {
    res+=psi[i]*x[i]*x[i]*conj(psi[i])+psi[i+1]*x[i+1]*x[i+1]*conj(psi[i+1]);
  }
  // cout << "res = " << real(res) << " " << imag(res) << endl;
  return real(res*dx/2.0);
}


double pmoy(vec_cmplx const& psi, double const& dx)
{
  //calcule la moyenne de p
  int Npoints(psi.size());
  double hbar(1.);
  complex<double> res(0,0);
  res+=conj(psi[0])*2.*(psi[1]-psi[0])+conj(psi[1])*(psi[2]-psi[1]);
  res+=conj(psi[Npoints-2])*(psi[Npoints-2]-psi[Npoints-3])+conj(psi[Npoints-1])*2.*(psi[Npoints-1]-psi[Npoints-2]);

  for (size_t i = 1; i < psi.size()-2; i++) {
    res+=conj(psi[i])*(psi[i+1]-psi[i-1])+conj(psi[i+1])*(psi[i+2]-psi[i]);
  }

  // cout << "pm= " << res;

  double pm(imag(res));
  return pm*hbar/4.;
}


double p2moy(vec_cmplx const& psi, double const& dx)
{
  //calcule la moyenne du p^2
  double hbar(1.);
  complex<double> res(0,0);
  for (size_t i = 1; i < psi.size()-2; i++) {
    res+=conj(psi[i])*(psi[i+1]-2.*psi[i]+psi[i-1])+conj(psi[i+1])*(psi[i+2]-2.*psi[i+1]+psi[i]);
  }

  // cout << "p2m= " << res;

  double p2m(real(res));
  return -p2m*hbar*hbar/(2.*dx);
}

vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi,0,psi.size()-1,dx));
  for(unsigned int i(0); i<psi.size(); ++i)
    psi_norm[i] = psi[i]/norm;
  return psi_norm;
}

double delx(vec_cmplx const& psi, vector<double> const& x, double const& dx)
{
  return sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx));
}

double delp(vec_cmplx const& psi, double const& dx)
{
  return sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx));
}

template <class T>
void coutBigFatVec(vector<T> const& vec, string name){
  cout << name << "=( ";
  for (size_t i = 0; i < vec.size() - 1; i++) {
    cout.precision(5);
    cout << vec[i] << ", ";
  }
  cout << vec[vec.size() - 1] << " )" << endl;
}
