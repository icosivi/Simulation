void Compilatore(){
  gSystem->CompileMacro("Punto.cxx","kgf");
  gSystem->CompileMacro("Retta.cxx","kgf");
  gSystem->CompileMacro("Neutron.cxx","kgf");
  gSystem->CompileMacro("Generatore.cxx","kgf");
  gSystem->CompileMacro("Propagatore.cxx","kgf");
  gSystem->CompileMacro("Rivelatore.cxx","kgf");
  gSystem->CompileMacro("Simulatore.C","kgf");
}