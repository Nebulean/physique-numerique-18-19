#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include "ConfigFile.tpp"

using namespace std;

double puissance(vector<vector<double> > const& T, vector<vector<double> >& jx, vector<vector<double> >& jy, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2, double const& L);

vector<vector<double> > deriv(vector<vector<double> > const& Told, vector<vector<double> > const& Tnew, double step, vector<vector<bool> > const& f);

double max(vector<vector<double> > const& vec);

void flux(vector<vector<double> > const& T, vector<vector<bool> > const& flag, vector<vector<double> >& jx, vector<vector<double> >& jy, vector<vector<double> >& jcx, vector<vector<double> >& jcy, double kappa, double h);

double sideOfSurfaceInt(vector<vector<double> > const& vec, double const h, bool const horizontal, size_t const where, size_t const lowerIndex, size_t const upperIndex, bool const positive);

size_t getIndexT(double const h, double const n, double const max);

template <typename type>
void printflag(vector<vector<type> > const& vec);

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
  // double alpha = Dcoef * dt / h / h;

  // Fichiers de sortie:
  string output = configFile.get<string>("output");
  ofstream output_T((output+"_T.out").c_str()); // Temperature au temps final
  ofstream output_P((output+"_P.out").c_str()); // Puissance au cours du temps
  ofstream output_F((output+"_F.out").c_str()); // Flux de chaleur
  output_T.precision(15);
  output_P.precision(15);
  output_F.precision(15);

  // Tableaux:
  vector<vector<bool> > flag(N+1,vector<bool>(N+1));
  vector<vector<double> > T(N+1,vector<double>(N+1));

  // j for the edges
  vector<vector<double> > jx(N, vector<double>(N+1, 0.0));
  vector<vector<double> > jy(N+1, vector<double>(N, 0.0));

  vector<vector<double> > jcx(N, vector<double>(N, 0.0));
  vector<vector<double> > jcy(N, vector<double>(N, 0.0));

  // vector<vector<vector<double> > > jcx(N+1, vector<vector<double> >(N, vector<double>(2, 0.0)));
  // vector<vector<double> > j(N+1, vector<double>(N+1));

  // j on the middle of the sides of the square
  // vector<vector<vector<double> > > jxc(N+1, vector<vector<double> >(N, vector<double>(2, 0.0)));
  // vector<vector<vector<double> > > jyc(N+1, vector<vector<double> >(N, vector<double>(2, 0.0)));

  // j on the middle of square
  // vector<vector<vector<double> > > jc(N+1, vector<vector<double> >(N, vector<double>(2, 0.0)));



  // TODO: Initialisation des tableaux
  //////////////////////////////////////
  for (size_t i = 0; i < flag.size(); i++) {
    for (size_t j = 0; j < flag.size(); j++) {
      // si on se trouve sur le bord
      if (i==0 || i==flag.size()-1 || j==0 || j==flag[i].size()-1) {
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

  printflag(flag);
  // printflag(flag);
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

    // flux
    flux(T, flag, jx, jy, jcx, jcy, kappa, h);

    // Diagnostiques:
    // output_P << iter*dt << " " << puissance(T, jx, jy, kappa, h, xa, xb, ya, yb, L)
    //                     << " " << puissance(T, jx, jy, kappa, h, xc, xd, ya, yb, L)
    //                     << " " << puissance(T, jx, jy, kappa, h, xa, xd, ya, yb, L) << endl;
  }
  output_P.close();

  // printflag(jx);

  // flux:
  // flux(T, flag, jx, jy, jcx, jcy, kappa, h);
  for(int i(0);i<N;++i)
    for(int j(0);j<N;++j)
      output_F << i*h + h/2 << " " << j*h + h/2 << " " << jcx[i][j] << " " << jcy[i][j] << " " << jx[i][j] << " " << jy[i][j] << endl;
  output_F.close();


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
// double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2)
double puissance(vector<vector<double> > const& T, vector<vector<double> >& jx, vector<vector<double> >& jy, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2, double const& L)
{
  // le résultat
  double power(0.0);

  // Les indices que l'on doit trouver
  size_t yBottom(y1/h-1);
  size_t yTop(y2/h);
  size_t xLeft(x1/h-1);
  size_t xRight(x2/h);

  // bottom
  for (size_t i = xLeft; i <= xRight; i++) {
    power -= jy[i][yBottom]*h;
  }

  // Top
  for (size_t i = xLeft; i <= xRight; i++) {
    power += jy[i][yTop]*h;
  }

  // Left
  for (size_t j = yBottom; j <= yTop; j++) {
    power -= jx[xLeft][j]*h;
  }

  // right
  for (size_t j = yBottom; j <= yTop; j++) {
    power += jx[xRight][j]*h;
  }

  // size_t yBottom(getIndexT(h, y1, L));
  // size_t yTop(getIndexT(h, y2, L));
  // size_t xLeft(getIndexT(h, x1, L));
  // size_t xRight(getIndexT(h, x2, L));
  // cout << yBottom << " " << yTop << " " << xLeft << " " << xRight << endl;
  // cout << yBottom*h+h/2 << " " << yTop*h+h/2 << " " << xLeft*h+h/2 << " " << xRight*h+h/2 << endl;

  // On calcul l'intégrale de surface.
  // bottom surface
  // power += sideOfSurfaceInt(jy, h, true, yBottom, xLeft, xRight, false);
  // for (size_t i = xLeft; i <= xRight; i++) {
  //   power -= jy[i][yBottom]*h;
  // }
  // cout << "Bottom: " << sideOfSurfaceInt(jy, h, true, yBottom, xLeft, xRight, false) << endl;

  // top surface
  // power += sideOfSurfaceInt(jy, h, true, yTop, xLeft, xRight, true);
  // for (size_t i = xLeft; i <= xRight; i++) {
  //   power += jy[i][yTop]*h;
  // }
  // cout << "Top: " << sideOfSurfaceInt(jy, h, true, yTop, xLeft, xRight, true) << endl;

  // left surface
  // power += sideOfSurfaceInt(jx, h, false, xLeft, yBottom, yTop, false);
  // for (size_t j = yBottom; j <= yTop; j++) {
  //   power -= jx[xLeft][j]*h;
  // }
  // cout << "Left: " << sideOfSurfaceInt(jx, h, false, xLeft, yBottom, yTop, false) << endl;

  // right surface
  // for (size_t j = yBottom; j <= yTop; j++) {
  //   power += jx[xRight][j]*h;
  // }
  // power += sideOfSurfaceInt(jx, h, false, xRight, yBottom, yTop, true);
  // cout << "Right: " << sideOfSurfaceInt(jx, h, false, xRight, yBottom, yTop, true) << endl;

  return power;
  // return 0;
}

// calcul la somme de tous les points sur une ligne. vec est la ligne à sommer, h est le paramètre pour l'élément de surface et positive permet de dire si on doit faire une somme de positifs (true) ou de négatifs (false)
double sideOfSurfaceInt(vector<vector<double> > const& vec, double const h, bool const horizontal, size_t const where, size_t const lowerIndex, size_t const upperIndex, bool const positive){
  double res(0.0);

  // cout << "Adding";

  if (horizontal) {
    for (size_t i = lowerIndex; i <= upperIndex; i++) {
      res += vec[i][where]*h;
      // cout << " " << vec[i][where];
    }
  } else {
    for (size_t j = lowerIndex; j <= upperIndex; j++) {
      res += vec[where][j]*h;
      // cout << " " << vec[where][j];
    }
  }

  if(!positive){
    return -res;
  }

  return res;
}

// Donne la position de la maille de T qui se trouve le plus proche de "n".
size_t getIndexT(double const h, double const n, double const max){
  double pos(h/2.0);
  size_t index(0);

  while (pos < max) {
    // cout << "Is " << pos << " < " << n << " < " << pos+h << " ? " << endl;
    if (pos+h > n && n > pos) {
      // cout << "yes" << endl;
      return index;
    }
    // cout << "no" << endl;
    pos += h;
    ++index;
  }
  cout << "error in index getter" << endl;
  return -1;

  // ancien code qui fait pas du tout ce que je veux, j'ai mal compris les paramètres dimensionnels.
  // double pos(h/2.0);
  // // double pos(0.0);
  // size_t index(0);
  //
  // while(pos < max){
  //   if (n < pos && pos < n+h) {
  //       return index;
  //   }
  //   pos += h;
  //   index += 1;
  // }
  // return -1; // error
}

void flux(vector<vector<double> > const& T, vector<vector<bool> > const& flag, vector<vector<double> >& jx, vector<vector<double> >& jy, vector<vector<double> >& jcx, vector<vector<double> >& jcy, double kappa, double h)
{
  // compute the flux - horizontal
  for (size_t i = 0; i < jx.size(); i++) {
    for (size_t j = 0; j < jx[i].size(); j++) {
      // if (!flag[i][j]) {
        jx[i][j] = -kappa/h *(T[i+1][j] - T[i][j]);
        // cout << "i=" << i << "=" << i*h+h/2 << " j=" << j << "=" << j*h+h/2 << " val=" << jx[i][j] << endl;
      // }
    }
  }

  // compute the flux - vertical
  for (size_t i = 0; i < jy.size(); i++) {
    for (size_t j = 0; j < jy[i].size(); j++) {
      // if (!flag[i][j]) {
        jy[i][j] = -kappa*(T[i][j+1] - T[i][j])/h;
      // }
    }
  }

  // compute the flux on the middle of the square
  for (size_t i = 0; i < jcx.size(); i++) {
    for (size_t j = 0; j < jcy.size(); j++) {
      // if (!flag[i][j]) {
        jcx[i][j] = (jx[i][j] + jx[i][j+1])/2;
        jcy[i][j] = (jy[i][j] + jy[i+1][j])/2;
      // }
    }
  }
}

template <typename type>
void printflag(vector<vector<type> > const& vec){
  cout << "Printing vector" << endl << "-----" << endl;
  for (size_t j = vec[1].size(); j > 0; j--) {
    for (size_t i = 0; i < vec.size(); i++) {
      cout << vec[i][j-1] << " ";
    }
    cout << endl;
  }
  cout << "-----" << endl;
}
