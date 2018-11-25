#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template
#include <valarray>
#include <vector>
#include <cmath>

using namespace std;



class Engine
{
private:
  double t, tFin, dt;
  int sampling;
  int last;
  ofstream *outputFile;
  valarray<double> p;   // vector containing all the informations about position and speed.
                        // Convention: (r1x, r1y, ... , rNx, rNy, v1x, v1y, ..., vNx, vNy)
                        // If there is Apollo 13, it is always the LAST body of the list.
                        // Earth is ALWAYS the first body.
  valarray<double> m;   // vector containing all the masses
  bool atm;             // false if no atmosphere, true is atmosphere.

  // some constants
  const double Rt     = 6378.1 * 1000;        // Earth's radius
  const double G      = 6.67408*pow(10,-11);  // Gravitationnal constant
  const double rho0   = 1.2;                  // Air density at see level
  const double lambda = 7238.2;               // The caracteristical width
  const double S      = 0;                    // Section area               //TODO: IL FAUT LE DEFINIR
  const double Cx     = 0.3;                  // Drag coefficient

  void printOut(bool force){
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << endl; // tous les trucs à print.
      last = 1;
    }
    else{
      last++;
    }
  }

  double grav(double m1, double m2, double r){
    return -G * m1 * m2 / pow(r,2);
  }

  double a(size_t body){ // "body" is the mass where the force is applied. Body starts at 0 !!!
    //TODO : CE N'EST PAS UNE VERSION CORRECTE, C'EST JUSTE UNE PERMIERE ETAPE DE SA REALISATION
    double result(0);
    double dist(0);
    double grav_effect(0);
    double atm_effect(0);


    // the positions are defined, it simplifies the notations.
    // the corps with index 'body' stays the same.
    double x1(p[2*body]);
    double x2(0);
    double y1(p[2*body + 1]);
    double y2(0);
    double vx1(p[p.size()/2. + 2*body]);
    double vx2(0);
    double vy1(p[p.size()/2. + 2*body + 1]);
    double vy2(0);

    // ====================== GRAVITATIONNAL EFFECT ============================
    // This for loop will computate the gravitationnal effect on the body.
    for (size_t i = 0; i < m.size(); i++) {
      // check that the body is not affected by itself
      if (i == body)
        i++;
      if (i >= m.size())
        return grav_effect;

      // computation of the distance between both objects
      x2 = p[2*i]; y2 = p[2*i + 1];
      dist = sqrt(pow(x2-x1,2) + pow(y2-y1,2));

      // Computation of the effect of gravity
      grav_effect += grav(m[body], m[i], dist);
    }
    // TODO: On a calculé la norme, il faut encore le projeter.

    // ======================== ATMOSPHERE EFFECT ==============================
    // we add the effect of atmosphere on Apollo 13, which is always the LAST body in p.
    // If this effect does not apply, or there is no Apollo 13,
    if (atm && (body == m.size() - 1)) {
      // computation of the distance between earth and Apollo 13.
      // first, we compute the distance between both bodies
      x2 = p[0];
      y2 = p[1];
      dist = sqrt(pow(x2-x1,2) + pow(y2-y1,2));

      // then, the difference of speed
      vx2 = p[p.size()/2.];
      vy2 = p[p.size()/2. + 1];
      double diff_speed = sqrt(pow(vx2 - vx1,2) + pow(vy2 - vy1,2));

      // then, rho
      double rho( rho0 * exp( - (dist - Rt)/lambda) );

      // then, the effect of earth's atmosphere
      atm_effect = -0.5 * rho * S * Cx * diff_speed;
    }
    // TODO: On a calculé la norme, il faut encore le projeter.
  }



  void step(){
    // RK4
    // some initialisations
    vector<double> k1, k2, k3, k4;
    size_t mid = p.size()/2.;
    size_t full = p.size();

    // We start by computing the changes k.
    // k1 - positions.
    for (size_t i = mid; i < full; i++) {
      k1.push_back( dt*p[i] );
    }

    // k1 - speed.
    for (size_t i = mid; i < full; i++) {
      k1.push_back( dt*a(i) );
    }

    // k2 - positions.
    for (size_t i = mid; i < full; i++) {
      k2.push_back( dt*( p[i]+0.5*k1[i] ) ); // C'est bizarre, mais ça semble être ce qu'il faut corriger par rapport au rapport 2.
    }

    // k2 - speed.
    for (size_t i = mid; i < full; i++) {
      // k2.push_back( dt*( a() ) ); // Finir l'acceleration avant d'écrire.
    }

    // k3 - positions.
    for (size_t i = mid; i < full; i++) {
      k3.push_back( dt*( p[i]+0.5*k2[i] ));
    }

    // k3 - speed.
    for (size_t i = mid; i < full; i++) {
      // k3.push_back( dt*( a() ) ); // Finir l'acceleration avant d'écrire.
    }

    // k4 - positions.
    for (size_t i = mid; i < full; i++) {
      k4.push_back( dt*(p[i]+k3[i]) );
    }

    // k4 - speed.
    for (size_t i = mid; i < full; i++) {
      // k4.push_back( dt*( a() )); // Finir l'acceleration avant d'écrire.
    }

    // We apply to the current vector p.
    for (size_t i = 0; i < full; i++) {
      p[i] += 1./6.*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
    }
  }

public:
  Engine(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);
      //TODO: Il faut que l'on implémente nos propres variables.
    tFin        = configFile.get<double>("tFin");
    dt          = configFile.get<double>("dt");
    atm         = configFile.get<bool>("atm");
    // m           = configFile.get<valarray<double> >("m"); //TODO: Comment faire pour avoir une liste d'éléments ?
    // d        = configFile.get<double>("d");
    // Omega    = configFile.get<double>("Omega");
    // kappa    = configFile.get<double>("kappa");
    // m        = configFile.get<double>("m");
    // g        = configFile.get<double>("g");
    // L        = configFile.get<double>("L");
    // theta    = configFile.get<double>("theta0");
    // thetadot = configFile.get<double>("thetadot0");
    // sampling = configFile.get<int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  }

  ~Engine(){
    outputFile->close();
    delete outputFile;
  }

  void run()
  {
    t = 0.;
    last = 0;
    printOut(true);
    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  }
};



int main(int argc, char* argv[]){
  Engine engine(argc, argv);
  engine.run();
  return 0;
}
