#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;



class Engine
{
private:
  double t, tFin, dt;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force){
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << endl; // tous les trucs à print.
      last = 1;
    }
    else{
      last++;
    }
  }

  void step(){
    // implémenter RK4.
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
    tFin     = configFile.get<double>("tFin");
    dt       = configFile.get<double>("dt");
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
