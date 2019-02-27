#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include "ConfigFile.tpp"

using namespace std;

double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2);

vector<vector<double> > deriv(vector<vector<double> > const& Told, vector<vector<double> > const& Tnew, double step, vector<vector<bool> > const& f);

double max(vector<vector<double> > const& vec);

int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice5 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice5 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Geometrie:
  double L  = configFile.get<double>("L");
  double xa = configFile.get<double>("xa");
  double xb = configFile.get<double>("xb");
  double xc = configFile.get<double>("xc");
  double xd = configFile.get<double>("xd");
  double ya = configFile.get<double>("ya");
  double yb = configFile.get<double>("yb");

  // Temperatures:
  double Tc = configFile.get<double>("Tc");
  double Tf = configFile.get<double>("Tf");
  double Tb = configFile.get<double>("Tb");

  // Physical coefficients
  double kappa = configFile.get<double>("kappa");
  double rho = configFile.get<double>("rho");
  double Ccoef = configFile.get<double>("Ccoef");
  double Dcoef = kappa/(rho*Ccoef);

  // Duree de la simulation:
  double tfin = configFile.get<double>("tfin");
  double eps = configFile.get<double>("eps"); // Condition d'arret si etat stationnaire

  // Discretisation:
  int N = configFile.get<int>("N"); // Nombre d'intervalles dans chaque dimension
  double dt = configFile.get<double>("dt");
  double h = L/N;
  double alpha = Dcoef * dt / h / h;

  // Fichiers de sortie:
  string output = configFile.get<string>("output");
  ofstream output_T((output+"_T.out").c_str()); // Temperature au temps final
  ofstream output_P((output+"_P.out").c_str()); // Puissance au cours du temps
  output_T.precision(15);
  output_P.precision(15);

  // Tableaux:
  vector<vector<bool> > flag(N+1,vector<bool>(N+1));
  vector<vector<double> > T(N+1,vector<double>(N+1));


  // TODO: Initialisation des tableaux
  //////////////////////////////////////
  for (size_t i = 0; i < flag.size(); i++) {
    for (size_t j = 0; j < flag.size(); j++) {
      // si on se trouve sur le bord
      if (i==0 || i==flag.size()-1 || j==0 || j==flag.size()-1) {
        flag[i][j] = true;
        T[i][j] = Tb;
      }
      // si on se trouve dans la source chaude
      else if (xa <= h*i && xb >= h*i && ya <= h*j && yb >= h*j) {
        flag[i][j] = true;
        T[i][j] = Tc;
      }
      // si on se trouve dans la source froide
      else if (xc <= h*i && xd >= h*i && ya <= h*j && yb >= h*j) {
        flag[i][j] = true;
        T[i][j] = Tf;
      }
      // sinon on se trouve ailleurs
      else {
        flag[i][j] = false;
        T[i][j] = Tb;
      }
    }
  }

  // Iterations:
  //////////////////////////////////////
  vector<vector<double>> Told(T);
  vector<vector<double> > Tstar(T);
  double m(eps+1);

  //Initial step
  // if (!flag[i][j]) {
  //   for (size_t i = 0; i < T.size(); i++) {
  //     for (size_t j = 0; j < count; j++) {
  //       Tstar[i][j] = Told[i][j] + Dcoef*dt/(h*h)*(Told[i+1][j] + Told[i-1][j] + Told[i][j+1] + Told[i][j-1] + 4*Told[i][j]);
  //       T[i][j] = Told[i][j] + alpha*(Tstar[i][j] - Told[i][j]);
  //     }
  //   }
  // }

  for(int iter=0; iter*dt<tfin && m>=eps; ++iter)
  {
    // cout << "itération: " << iter << endl;
    // TODO: Schema a 2 niveaux et calcul de max(|dT/dt|)
    // schema a 2 niveaux
    for (size_t i = 0; i < T.size(); i++) {
      for (size_t j = 0; j < T.size(); j++) {
        if (!flag[i][j]) {
          // Jacobi sans surrelaxation
          T[i][j] = Told[i][j] + Dcoef*dt/(h*h)*(Told[i+1][j] + Told[i-1][j] + Told[i][j+1] + Told[i][j-1] - 4*Told[i][j]);
          // GS sans surrelaxation
          // T[i][j] = T[i][j] + Dcoef*dt/(h*h)*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]);
          // Jacobi avec surrelaxation
          // Tstar[i][j] = Told[i][j] + Dcoef*dt/(h*h)*(Told[i+1][j] + Told[i-1][j] + Told[i][j+1] + Told[i][j-1] - 4*Told[i][j]);
          // T[i][j] = Told[i][j] + alpha*(Tstar[i][j] - Told[i][j]); // ATTENTION: Alpha est peut-être faux, c'est pas le même que dans le cours
        }
      }
    }

    m = max(deriv(Told,T,dt,flag));
    // cout << "m=" << m << endl;
    Told = T;
    Tstar = T;

    // Diagnostiques:
    output_P << iter*dt << " " << puissance(T, kappa, h, xa, xb, ya, yb)
                        << " " << puissance(T, kappa, h, xc, xd, ya, yb)
                        << " " << puissance(T, kappa, h, xa, xd, ya, yb) << endl;
  }
  output_P.close();


  // Ecriture de la temperature finale:
  for(int i(0);i<N+1;++i)
    for(int j(0);j<N+1;++j)
      output_T << i*h << " " << j*h << " " << T[i][j] << endl;
  output_T.close();
  return 0;
}

double max(vector<vector<double> > const& vec){
  double res(-300);
  for (size_t i = 0; i < vec.size(); i++) {
    for (size_t j = 0; j < vec.size(); j++) {
      if (vec[i][j] > res) {
        res = vec[i][j];
      }
    }
  }

  return res;
}

vector<vector<double> > deriv(vector<vector<double> > const& Told, vector<vector<double> > const& Tnew, double step, vector<vector<bool> > const& f){
  vector<vector<double>> res(Tnew.size(),vector<double>(Tnew.size()));
  for (size_t i = 0; i < res.size(); i++){
    for (size_t j = 0; j < res.size(); j++){
      if (!f[i][j])
        res[i][j] = 1./step * (Tnew[i][j]-Told[i][j]);
      else
        res[i][j] = 0;
    }
  }
  return res;
}

// TODO: Calculer la puissance calorifique emise/recue par le rectangle allant de (x1,y1) a (x2,y2)
double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2)
{
  return 0;
}
