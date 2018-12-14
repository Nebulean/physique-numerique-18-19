#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template
#include <valarray>
#include <vector>
#include <cmath>
#include <string>

using namespace std;



class Engine
{
private:
  double t, tFin, dt, epsilon; //epsilon la précision pour le dt adaptatif
  int sampling;
  int last;
  ofstream *outputFile;
  valarray<double> p;   // vector containing all the informations about position and speed.
                        // Convention: (x1, y1, vx1, vy1, ..., x3, y3, vx3, vy3)
  double m1, m2, m3;    // mass of bodies
  bool atm;             // false if no atmosphere, true is atmosphere.
  bool dtad;            // indicates whether run() should use the adaptative dt or not

  // some constants
  const double Rt     = 6378.1 * 1000;        // Earth's radius
  const double G      = 6.674e-11;  // Gravitational constant
  const double rho0   = 1.2;                  // Air density at see level
  const double lambda = 7238.2;               // The caracteristical width
  const double S      = 11.9459060516;                    // Section area               //TODO: IL FAUT LE DEFINIR
  const double Cx     = 0.3;                  // Drag coefficient

  // prints result in output file
  void printOut(bool force){
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << t << " " << dt; // tous les trucs à print. 1 2
      for (size_t i = 0; i < p.size(); i++) { // 3-14
        *outputFile << " " << p[i];
      }
      // we print the acceleration of the third body (apollo in general), use in 4.3.
      valarray<double> accel(2);
      accel = a(2, p);
      *outputFile << " " << accel[0] << " " << accel[1]; // 15 16

      // we print the power of drag.
      *outputFile << " " << powerOfDrag(p); // 17

      // prints the mechanical energy of Earth and Moon
      *outputFile << " " << emec(0) << " " << emec(1); // 18 19

      double qmvmt1(m1*norm(getVel(0,p))), qmvmt2(m2*norm(getVel(1,p))); // momentum of Earth and Moon
      *outputFile << " " << qmvmt1 << " " << qmvmt2; // 20 21

      *outputFile << " " << norm(getPos(0,p)-getPos(1,p)); // 22

      *outputFile << endl;
      last = 1;
    }
    else{
      last++;
    }
  }

  valarray<double> centerOfMass(size_t bod1, size_t bod2, valarray<double> const& vec){ //returns coordinates of the center of mass of bod1 and bod2
    valarray<double> r1(getPos(bod1, vec)), r2(getPos(bod2, vec));
    double M1(getMass(bod1)), M2(getMass(bod2));
    return (r1*M1+r2*M2)/(M1+M2);
  }

  double emec(size_t body){
    double pot = -G * getMass(body) * (getMass((body+1)%3) + getMass((body+2)%3))/norm(centerOfMass((body+1)%3,(body+2)%3,p) - getPos(body,p));
    double kin = 0.5 * getMass(body) * norm(getVel(body, p))*norm(getVel(body, p));
    return pot+kin;
  }

  double powerOfDrag(valarray<double> const& vec){
    if (atm) {
      double power(0.);
      valarray<double> F(0., 2);
      valarray<double> v(0., 2);

      F = drag(2, 0, vec);
      v = getVel(2, vec);
      power = F[0]*v[0] + F[1]*v[1];

      return power;
    }
    return 0.;
  }

  // returns the position of the body.
  valarray<double> getPos(size_t body, valarray<double> const& vec){
    valarray<double> r(2);
    // r[0] = vec[4*body];
    // r[1] = vec[4*body + 1];
    r = vec[slice(4*body,2,1)];
    return r;
  }

  // returns the velocity of the body.
  valarray<double> getVel(size_t body, valarray<double> const& vec){
    valarray<double> v(2);
    // v[0] = vec[4*body + 2];
    // v[1] = vec[4*body + 3];
    v = vec[slice(4*body+2,2,1)];
    return v;
  }

  double getMass(size_t body){
    switch (body) {
      case 0:
        return m1;
        break;
      case 1:
        return m2;
        break;
      case 2:
        return m3;
        break;
      default:
        return 0;
        break;
    }
  }

  // returns the euclidian norm of the array. (as if it were a vector)
  double norm(valarray<double> const& v){
    double res(0);
    for (auto const& el : v) {
      res += pow(el,2);
    }
    return sqrt(res);
  }

  // valarray<double> dist(size_t bod1, size_t bod2, valarray<double> const& vec){ // distance between bod1 and bod2
  //   return (getPos(bod1, vec) - getPos(bod2, vec));
  // }

  // returns the gravitational effect of mass 2 on mass 1.
  valarray<double> grav(size_t target, size_t actor, valarray<double> const& vec){
    valarray<double> g(2);

    valarray<double> r(2);
    r = getPos(target, vec) - getPos(actor, vec);
    // r = dist(target, actor, p);
    g = -G * getMass(target) * getMass(actor) / pow(norm(r), 3) * r;
    return g;
  }

  // returns the air density on earth.
  double rho(valarray<double> const& r){
    double n(norm(r));
    return rho0 * exp(-(n - Rt)/lambda);
  }

  // returns the drag.
  valarray<double> drag(size_t target, size_t actor, valarray<double> const& vec){
    if (actor == 0 && target == 2 && atm) { // Earth is always on first position and Apollo on third position. atm is a toggle.
      valarray<double> f(2);
      valarray<double> v(2);
      valarray<double> r(2);

      r = getPos(target, vec) - getPos(actor, vec);

      v = getVel(target, vec) - getVel(actor, vec);

      f = -0.5 * rho(r) * S * Cx * norm(v) * v;

      return f;
    }

    return {0,0};
  }

  valarray<double> a(size_t body, valarray<double> const& vec){
    valarray<double> res(0.,2);
    for(size_t i(0); i<3; ++i){
      if(i != body){
        res += grav(body, i, vec)/getMass(body);
        res += drag(body, i, vec)/getMass(body);
      }
    }
    return res;
  }

  valarray<double> step(valarray<double> const& v, double time_step){
    // RK4
    // some initialisations
    size_t full = v.size();
    valarray<double> k1(full), k2(full), k3(full), k4(full);

    // cout << "time step = " << time_step << endl;

    // We start by computing the changes k.
    // k1
    for (size_t body = 0; body < 3; body++) {
      // position
      k1[4*body] = time_step*v[4*body + 2];
      k1[4*body + 1] = time_step*v[4*body + 3];

      // speed
      valarray<double> vec(2);
      vec = a(body, v);
      k1[4*body + 2] = time_step*vec[0];
      k1[4*body + 3] = time_step*vec[1];
    }

    // k2
    for (size_t body = 0; body < 3; body++) {
      // position
      k2[4*body] = time_step*( v[4*body + 2] + 0.5*k1[4*body + 2] );
      k2[4*body + 1] = time_step*( v[4*body + 3] + 0.5*k1[4*body + 3] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + 0.5*k1);
      k2[4*body + 2] = time_step*vec[0];
      k2[4*body + 3] = time_step*vec[1];
    }


    // k3
    for (size_t body = 0; body < 3; body++) {
      // position
      k3[4*body] = time_step*( v[4*body + 2] + 0.5*k2[4*body + 2] );
      k3[4*body + 1] = time_step*( v[4*body + 3] + 0.5*k2[4*body + 3] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + 0.5*k2);
      k3[4*body + 2] = time_step*vec[0];
      k3[4*body + 3] = time_step*vec[1];
    }


    // k4
    for (size_t body = 0; body < 3; body++) {
      // position
      k4[4*body] = time_step*( v[4*body + 2] + k3[4*body + 2] );
      k4[4*body + 1] = time_step*( v[4*body + 3] + k3[4*body + 3] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + k3);
      k4[4*body + 2] = time_step*vec[0];
      k4[4*body + 3] = time_step*vec[1];
    }


    valarray<double> res(v);

    res += (k1 + 2.0*k2 + 2.0*k3 + k4)/6;

    return res;
  }

  // changes the dt with a better one, and returns the computed step.
  valarray<double> stepDtAdapt(){
    valarray<double> p1(step(p, dt));
    valarray<double> ptemp(step(p, dt/2.));
    valarray<double> p2(step(ptemp, dt/2.));

    double d = abs(p1-p2).max();


    if(d<=epsilon){
      t += dt;
      dt *= pow(epsilon/d, 1./5.); // power 1/(n+1) with n the order of convergence
      return p2;
    } else {
      dt *= 0.99 * pow(epsilon/d, 1./5.);
      return stepDtAdapt();
    }
  }

  void coutBigFatVec(valarray<double> const& vec, string name){
    cout << name << "=( ";
    for (size_t i = 0; i < vec.size() - 1; i++) {
      cout << vec[i] << ", ";
    }
    cout << vec[vec.size() - 1] << " )" << endl;
  }

public:
  Engine(int argc, char* argv[]) {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin        = configFile.get<double>("tFin");
    dt          = configFile.get<double>("dt");
    atm         = configFile.get<bool>("atm");
    epsilon     = configFile.get<double>("epsilon");
    dtad        = configFile.get<bool>("dtad");
    m1          = configFile.get<double>("m1");
    m2          = configFile.get<double>("m2");
    m3          = configFile.get<double>("m3");
    sampling    = configFile.get<int>("sampling");
    p           = {configFile.get<double>("x1"), configFile.get<double>("y1"), configFile.get<double>("vx1"), configFile.get<double>("vy1"), configFile.get<double>("x2"), configFile.get<double>("y2"), configFile.get<double>("vx2"), configFile.get<double>("vy2"), configFile.get<double>("x3"), configFile.get<double>("y3"), configFile.get<double>("vx3"), configFile.get<double>("vy3")};
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  }

  ~Engine(){
    outputFile->close();
    delete outputFile;
  }

  void run() {
    t = 0.;
    last = 0;
    printOut(true);

    while( t < tFin-0.5*dt )
    {

      if(dtad){
        p = stepDtAdapt();
        // cout << "On est sorti de stepDtAdapt" << endl;
        //t += dt; //done in dtadapt
      } else {
        p = step(p, dt);
        t += dt;
      }

      printOut(false);
      // temporaire, afin d'arrêter la simulation quand l'engin s'écrase sur terre.
      if (norm(getPos(0, p) - getPos(2, p)) < Rt) {
        t = tFin;
      }
    }
    printOut(true);
  }
};



int main(int argc, char* argv[]){
  Engine engine(argc, argv);
  engine.run();
  return 0;
}
