/*
Se trabajo en el Codigo base realizado por el grupo de instructores de la materia 
de Tecnicas de Simulacion por computadoras.
Derechos a los respectivos creadores del codigo
 */

#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"
#include "display_tools.h"

int main(){
    vector<Matrix> localKs;
    vector<Vector> localbs;
    
    Matrix K;
    Vector b;
    Vector U;
    cout<< "TAREA LABo 3: Encontrar los valores para U de la ecuacion dada\n";
    mesh m;

    leerMallayCondiciones(m);
    crearSistemasLocales(m,localKs,localbs);

    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    
    applyNeumann(m,b);    
    applyDirichlet(m,K,b);

    zeroes(U,b.size());
    calculate(K,b,U);

    cout << "La respuesta es: \n";
    showVector(U);

    return 0;
}