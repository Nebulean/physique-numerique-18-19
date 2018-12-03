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
  const double G      = 6.67408e-11;  // Gravitational constant
  const double rho0   = 1.2;                  // Air density at see level
  const double lambda = 7238.2;               // The caracteristical width
  const double S      = 0;                    // Section area               //TODO: IL FAUT LE DEFINIR
  const double Cx     = 0.3;                  // Drag coefficient

  void printOut(bool force){
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << t << " " << dt; // tous les trucs à print.
      for (size_t i = 0; i < p.size(); i++) {
        *outputFile << " " << p[i];
      }
      *outputFile << endl;
      last = 1;
    }
    else{
      last++;
    }
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

  // void setPos(size_t body, valarray<double> vectarg){
  //   // p[4*body] = vectarg[0];
  //   // p[4*body + 1] = vectarg[1];
  //   p[slice(4*body,)]
  // }
  //
  // void setVel(size_t body, valarray<double> vectarg){
  //   p[4*body + 2] = vectarg[0];
  //   p[4*body + 3] = vectarg[1];
  // }

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

  // returns the gravitational effect of mass 2 on mass 1.
  // valarray<double> grav(double mass1, double mass2, valarray<double> const& r){
  //   valarray<double> g(2);
  //   g = -G * mass1 * mass2 * (1./pow(norm(r),3)) * r;
  //   return g;
  // }
  valarray<double> grav(size_t target, size_t actor, valarray<double> const& vec){
    valarray<double> g(2);
    // cout << "vec = ( ";
    // for (size_t i = 0; i < vec.size(); i++) {
    //   cout << vec[i] << ", ";
    // }
    // cout << ")" << endl;
    // double x = getPos(target, vec)[0]-getPos(actor, vec)[0];
    // double y = getPos(target, vec)[1]-getPos(actor, vec)[1];
    //
    // if(x==0){
    //   g[0] = 0.;
    // } else {
    //   g[0] = -G * getMass(target) * getMass(actor) / pow(x,2);
    // }
    //
    // if(y==0){
    //   g[1] = 0.;
    // } else {
    //   g[1] = -G * getMass(target) * getMass(actor) / pow(y,2);
    // }

    valarray<double> r(2);
    r = getPos(target, vec) - getPos(actor, vec);
    g = -G * getMass(target) * getMass(actor) / pow(norm(r), 3) * r;

    // g = -G * getMass(target) * getMass(actor) / pow(norm(getPos(target, vec)-getPos(actor, vec)),3) * (getPos(target, vec)-getPos(actor, vec));

    // cout << "grav = " << g[0] << " " << g[1] << endl;
    return g;
  }


  // returns the air density on earth.
  double rho(valarray<double> const& r){
    double n(norm(r));
    return rho0 * exp(-(n - Rt)/lambda);
  }

  // returns the drag.
  valarray<double> drag(valarray<double> const& r, valarray<double> const& v){
    valarray<double> f(2);
    f = -0.5 * rho(r) * S * Cx * norm(v) * v;
    return f;
  }

  valarray<double> a(size_t body, valarray<double> const& vec){
    valarray<double> res(0.,2);
    for(size_t i(0); i<3; ++i){
      if(i != body){
        res+= grav(body, i, vec)/getMass(body);
        // cout << "grav of " << i << " on " << body << ": (" << grav(body, i, vec)[0] << ", " << grav(body, i, vec)[0] << ")" << endl;
      }
    }
    // cout << "a =" << res[0] << " " << res[1] << endl;
    return res;
  }

  valarray<double> step(valarray<double> const& v, double time_step){
    // RK4
    // some initialisations
    size_t full = v.size();
    valarray<double> k1(full), k2(full), k3(full), k4(full);

    // We start by computing the changes k.
    // k1
    for (size_t body = 0; body < 3; body++) {
      // position
      k1[4*body] = time_step*v[4*body + 2];
      k1[4*body + 1] = time_step*v[4*body + 3];
      // k1[slice(4*body,2,1)] = time_step*v[slice(4*body + 2,2,1)];

      // speed
      valarray<double> vec(2);
      vec = a(body, v);
      k1[4*body + 2] = time_step*vec[0];
      k1[4*body + 3] = time_step*vec[1];
      // k1[slice(4*body+2,2,1)] = time_step*a(body, v);
    }

    // cout << "k1 = ( ";
    // for (size_t i = 0; i < k1.size(); i++) {
    //   cout << k1[i] << ", ";
    // }
    // cout << ")" << endl;

    // k2
    for (size_t body = 0; body < 3; body++) {
      // position
      k2[4*body] = time_step*( v[4*body + 2] + 0.5*k1[4*body + 2] );
      k2[4*body + 1] = time_step*( v[4*body + 3] + 0.5*k1[4*body + 3] );
      // k2[slice(4*body,2,1)] = time_step*( v[slice(4*body + 2,2,1)] + 0.5*k1[slice(4*body + 2,2,1)] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + 0.5*k1);
      k2[4*body + 2] = time_step*vec[0];
      k3[4*body + 3] = time_step*vec[1];
      // k2[slice(4*body + 2,2,1)] = time_step*a(body, v + 0.5*k1);
    }

    // cout << "k2 = ( ";
    // for (size_t i = 0; i < k2.size(); i++) {
    //   cout << k2[i] << ", ";
    // }
    // cout << ")" << endl;

    // k3
    for (size_t body = 0; body < 3; body++) {
      // position
      k3[4*body] = time_step*( v[4*body + 2] + 0.5*k2[4*body + 2] );
      k3[4*body + 1] = time_step*( v[4*body + 3] + 0.5*k2[4*body + 3] );
      // k3[slice(4*body,2,1)] = time_step*( v[slice(4*body + 2,2,1)] + 0.5*k2[slice(4*body + 2,2,1)] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + 0.5*k2);
      k3[4*body + 2] = time_step*vec[0];
      k3[4*body + 3] = time_step*vec[1];
      // k3[slice(4*body + 2,2,1)] = time_step * a(body, v + 0.5*k2);
    }

    // cout << "k3 = ( ";
    // for (size_t i = 0; i < k3.size(); i++) {
    //   cout << k3[i] << ", ";
    // }
    // cout << ")" << endl;

    // k4
    for (size_t body = 0; body < 3; body++) {
      // position
      k4[4*body] = time_step*( v[4*body + 2] + k3[4*body + 2] );
      k4[4*body + 1] = time_step*( v[4*body + 3] + k3[4*body + 3] );
      // k4[slice(4*body,2,1)] = time_step*( v[slice(4*body + 2,2,1)] + k3[slice(4*body + 2,2,1)] );

      // speed
      valarray<double> vec(2);
      vec = a(body, v + k3);
      k4[4*body + 2] = time_step*vec[0];
      k4[4*body + 3] = time_step*vec[1];
      // k4[slice(4*body + 2,2,1)] = time_step * a(body, v + k3);
    }

    // cout << "k4 = ( ";
    // for (size_t i = 0; i < k4.size(); i++) {
    //   cout << k4[i] << ", ";
    // }
    // cout << ")" << endl;

    valarray<double> res(v);

    res += 1./6. * (k1 + 2.*k2 + 2.*k3 + k4);

    // cout << "res = (";
    // for (size_t i = 0; i < res.size(); i++) {
    //   cout << res[i] << ", ";
    // }
    // cout << ")" << endl;

    return res;
  }

  // changes the dt with a better one, and returns the computed step.
  valarray<double> stepDtAdapt(){
    valarray<double> p1(step(p, dt));
    valarray<double> ptemp(step(p, dt/2));
    valarray<double> p2(step(ptemp, dt/2));

    double d = abs(p1-p2).max();
    coutBigFatVec(p1, "p1");
    coutBigFatVec(p2, "p2");

    if(d<=epsilon){
      t += dt;
      dt *= pow(epsilon/d, 1./5.); // power 1/(n+1) with n the order of convergence
      return p2;
    } else {
      dt *= 0.99 * dt * pow(epsilon/d, 1./5.);
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
    // atm         = configFile.get<bool>("atm");
    epsilon     = configFile.get<double>("epsilon");
    dtad        = configFile.get<bool>("dtad");
    m1          = configFile.get<double>("m1");
    m2          = configFile.get<double>("m2");
    m3          = configFile.get<double>("m3");
    sampling    = configFile.get<int>("sampling");
    p           = {configFile.get<double>("x1"), configFile.get<double>("y1"), configFile.get<double>("vx1"), configFile.get<double>("vy1"), configFile.get<double>("x2"), configFile.get<double>("y2"), configFile.get<double>("vx2"), configFile.get<double>("vy2"), configFile.get<double>("x3"), configFile.get<double>("y3"), configFile.get<double>("vx3"), configFile.get<double>("vy3")};
    sampling    = configFile.get<double>("sampling");
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
        //t += dt; //done in dtadapt
      } else {
        p = step(p, dt);
        t += dt;
      }

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
