/*
********************************************************************************************
********************************************************************************************
********************************************************************************************
 
    PROGRAMA VOLUNTARIO SCHROEDINGER

    Input: Headers con los modulos que contienen los modulos principales
    Output: Archivos de salida de la funcion de onda y sus observables a lo largo del tiempo

********************************************************************************************
********************************************************************************************
********************************************************************************************
*/

#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>

#include "complex.h"    //Header con las estructura de complejos

const double PI = 3.14159265358979323846;
const double e = 2.7182818284;

using namespace std;

//Defino las funciones utilizadas
void Inicio_F_Onda_Gaussiana (vector <fcomplex>& , double , double , double );
void Inicio_Potencial_Armonico (vector <double> & , vector <double> & , double , double , double );
void Normalizar_F_Onda (vector <fcomplex>& , double);
void Probabilidad_F_Onda (vector <fcomplex> & , double , double & , double &);
void Valor_Medio_Posicion (vector <fcomplex> & , double , double & , double & );
void Valor_Medio_Posicion_Gauss (vector <fcomplex> & , double , double & ); 
void Valor_Medio_Posicion_Cuadrado (vector <fcomplex> & , double , double & , double & );
void Valor_Medio_Momento (vector <fcomplex> & , double , double & , double & );
void Valor_Medio_Momento_Centrado (vector <fcomplex> & , double , double & , double & );
void Valor_Medio_Momento_Cuadrado (vector <fcomplex> & , double , double & , double & );
void Valor_Medio_Potencial (vector <fcomplex> & , vector <double> &, double ,  double & , double & );


/*
    BLOQUE PRINCIPAL:

    Bloque de la funcion principal donde usamos las funciones definidas para poder inicializar el programa,
    realizar los cálculos y posteriormente escribir los resultados en los ficheros de salida

*/
int main()
{
    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
                                            BLOQUE DECLARACIONES
    /////////////////////////////////////////////////////////////////////////////////////////////////////*/

    //Declaro unas variables auxiliares complejas
    fcomplex zero = {0.0, 0.0};
    fcomplex i = {0.0,1.0};
    fcomplex unidad = {1.0,0.0};

    //Declaro los ficheros de salida
    ofstream Funcion_Onda ("./sigma116_centro05/Gauss_Funcion_Onda.txt");
    ofstream Observable ("./sigma116_centro05/Gauss_Observables.txt");
    ofstream Condiciones_Iniciales ("./sigma116_centro05/Condiciones_Iniciales.txt");

    //Escribo los encabezados de los ficheros
    
    Funcion_Onda << "Tiempo" << "\t" << "Posicion" << "\t" << "Probabilidad" << "\t" << "Parte_Real" << "\t" << "Parte_Imaginaria" <<  "\n";
    Observable << "Tiempo" << "\t" << "Probabilidad" << "\t" << "Valor_Esperado_x" << "\t" << "Valor_Esperado_x_Gauss" << "\t" << "Valor_Esperado_p" << "\t" << "Incertidumbre_x" << "\t" << "Incertidumbre_p" << "\t" << "Valor_Esperado_Energia" << "\t" << "Heisenberg"<< "\n";   
    Condiciones_Iniciales << "Longitud_red" << "\t" << "Paso_Espacial" << "\t" << "Paso_Temporal" << "\t" << "Sigma" << "\t" << "Omega" << "\t" << "Centro_Gauss" << "\n";

    //Aseguro que el programa lea/escriba en los ficheros
    if (!Funcion_Onda || !Observable || !Condiciones_Iniciales)
    {
        cout << "Error al abrir los ficheros de entrada/salida !!!" << endl;
        return 0;
    }

    //Declaro los parámetros iniciales
    double L = 1;                   //Longitud Red
    int S = 1000;                   //Numero de pasos intermedios
    double centro = L/2;            //Guarda el centro del potencial (es simetrico respecto el centro)
    double dx = L/(S*1.0);          //Paso espacial
    double omega = 200;             //Frecuencia del potencial
    double dt = 0.0001;             //Paso temporal
    double x_0;                     //Valor centro gaussiana
    double sigma;                   //Anchura de la Gaussiana
    int n = 0;                      //Contador temporal

    //Variables iniciales gaussiana
    sigma = 1/16;
    x_0 = 0.5;

    Condiciones_Iniciales << L << "\t" << dx << "\t" << dt << "\t" << sigma << "\t" << omega << "\t" << x_0 << "\n";

    //Declaro los parametros necesarios para el algoritmo
    double stilde = dt/(dx*dx);
    vector <double> V (S, 0); 
    vector <double> Vtilde (S, 0); 

    //Declaro los observables
    double x_simpson, x_trapecio, x_simpson2, x_trapecio2;
    double p_simpson, p_trapecio, p_simpson2, p_trapecio2;
    double prob_trapecio, prob_simpson;
    double x_gauss;
    double V_trapecio, V_simpson;
    double E_trapecio, E_simpson;
    double Delta_x_trapecio, Delta_x_simpson, Delta_p_trapecio, Delta_p_simpson;
    double Heisenberg_trapecio, Heisenberg_simpson;


    //Declaro el vector que me guarda la funcion de onda
    vector <fcomplex> F ((int)S, zero); //Funcion de onda principal (de num_puntos pasos e inciada a 0)

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
                                            BLOQUE INICIALIZACION
    /////////////////////////////////////////////////////////////////////////////////////////////////////*/
    
    //Inicializo el potencial y la funcion de onda
    Inicio_Potencial_Armonico(V,Vtilde,omega,dx, centro);

    //Inicializo la funcion de onda con una autofuncion
    Inicio_F_Onda_Gaussiana(F,sigma, x_0, dx);

    //Normalizo la funcion de onda
    Normalizar_F_Onda(F, dx);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

     
    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
                                            BLOQUE EVOLUCION TEMPORAL 
    /////////////////////////////////////////////////////////////////////////////////////////////////////*/
    
    //Introducimos el algoritmo de evolucion temporal, extraido del código del trabajo obligatorio
    //////Escribimos los parametros iniciales////////
    //Recorro todos los valores de la funcion de onda para escribirlos en un fichero
    for (int j = 0; j< F.size(); j++)
    {
        double prob = 0;

        //Obtengo la densidad de probabilidad en un punto
        prob = Cmul(F[j],Conjg(F[j])).r;


        //Escribo en el fichero los observables de la funcion de onda
        Funcion_Onda << setprecision(10) << 0 << "\t" << setprecision(10) << (j*dx) << "\t" << prob << setprecision(10) << "\t" << setprecision(10) << F[j].r << "\t" << setprecision(10) << F[j].i <<  "\n";
    }

    //Calculo y Escritura de Observables//

    //Obtengo los valores de los observables (probabilidad y valor esperado de x, x^2, p, p^2)
    Probabilidad_F_Onda(F,dx,prob_trapecio,prob_simpson);
    Valor_Medio_Posicion(F,dx, x_trapecio, x_simpson);
    Valor_Medio_Posicion_Gauss(F,dx,x_gauss);
    Valor_Medio_Posicion_Cuadrado(F,dx,x_trapecio2,x_simpson2);
    Valor_Medio_Momento_Centrado(F,dx,p_trapecio,p_simpson);
    Valor_Medio_Momento_Cuadrado(F,dx,p_trapecio2,p_simpson2);
    Valor_Medio_Potencial(F,V,dx,V_trapecio,V_simpson);

    //Obtengo Valor Medio Energia
    E_trapecio = p_trapecio2 + V_trapecio;
    E_simpson = p_simpson2 + V_simpson;

    //Obtengo las incertidumbres de p y x mediante la definicion formal
    Delta_p_trapecio = sqrt(abs(p_trapecio2 -p_trapecio*p_trapecio));
    Delta_p_simpson = sqrt(abs(p_simpson2-p_simpson*p_simpson));
    Delta_x_trapecio = sqrt(abs(x_trapecio2-x_trapecio*x_trapecio));
    Delta_x_simpson = sqrt(abs(x_simpson2-x_simpson*x_simpson));

    Heisenberg_trapecio = Delta_p_trapecio*Delta_x_trapecio;
    Heisenberg_simpson = Delta_x_simpson*Delta_p_simpson;

    Observable << setprecision(10) << dt*0 << "\t" << setprecision(10) << prob_trapecio << "\t" << setprecision(15) << x_trapecio << "\t" << setprecision(15) << x_gauss << "\t" << setprecision(10) << p_trapecio << "\t" << setprecision(10) << Delta_x_trapecio << "\t" << setprecision(10) << Delta_p_trapecio << "\t" << setprecision(10) << E_trapecio << "\t" << setprecision(10) << Heisenberg_trapecio << "\n";     

    //Declaro las variables recursivas que me haran falta para la evolucion temporal
    vector <fcomplex> A_0 (S, zero); //Variable recursiva A_0 (de num_puntos pasos e inciada a 0)
    vector <fcomplex> alpha (S, zero); //Variable de iteracion alpha
    vector <fcomplex> gamma (S, zero); //Variable de iteracion gamma
    vector <fcomplex> G (S, zero); //Funcion de onda auxiliar (de num_puntos pasos e inciada a 0)
    vector <fcomplex> b (S, zero); //Variable recursiva b (de num_puntos pasos e inciada a 0)
    vector <fcomplex> beta (S, zero); //Variable recursiva beta (de num_puntos pasos e inciada a 0)

    //Aplico la condicion inicial de la funcion de onda es 0 en 0 y N (lo que significa que alpha tmb en esos puntos)
    alpha[0].r = alpha[S-1].r = alpha[0].i = alpha[S-1].i = 0;

    //Inicio los valores de A_0 (que son independientes del tiempo)
    for (int j = 0; j< S; j++)
    {
        double aux = 2/stilde;   //Variable auxiliar para multiplicar por el complejo
        A_0[j].r = -2-Vtilde[j];
        A_0[j].i = aux;
    }

    //Obtengo los valores de alpha y gamma (que son independientes del tiempo)
    for (int j = S-1 ; j>=0; j--)
    {
        fcomplex aux = Cadd(A_0[j],alpha[j]);  //Variable auxiliar del sumando del denominador (A+ = 1)
        gamma[j] = Cdiv({1,0},aux);         //Obtengo la gamma en j
        if ((j-1)>0)    //Filtro de j-1 es mayor que 0 para meter el valor de alpha, ya que alpha en 0 es 0
        {
            alpha[j-1] = Csub({0,0},gamma[j]);
        }
    }
    
    n=1;    //Empiezo el ciclo con el tiempo en 1 pq ya he hecho el incial de 0 y lo he escrito

    while (n<=5000)
    {
        beta[0].r = beta[S-1].r = beta[0].i = beta[S-1].i = 0;    //Condicion de contorno de beta
        
        //Recorro todo el espacio para inicializar b y beta (de N a 0)
        for (int j = S-1; j>=0; j--)
        {
            b[j].i = 4/stilde;
            b[j].r = 0;
            b[j] = Cmul(b[j],F[j]);
            
            if ((j-1)>0)    //Filtro de j-1 es mayor que 0 para meter el valor de beta, ya que beta en 0 es 0
            {
                beta[j-1] = Csub(b[j],beta[j]);
                beta[j-1] = Cmul(beta[j-1],gamma[j]);
            }
        }

        //Impongo la condicion de contorno de G[0] = 0
        G[0] = G[G.size()-1] = {0.0};

        //Recorro todo el espacio menos los extremos, y obtengo la funcion de onda auxiliar (de 0 a N)
        for (int j = 0; j<S; j++)
        {
            //Declaro la variable probabilidad
            double prob = 0;

            if (j>0)
            {
                fcomplex aux = Cadd(A_0[j],alpha[j]);   //Variable auxiliar que me guarda el denominador comun de la expresion
                fcomplex suma = Csub(b[j],beta[j]);     //Variable auxiliar que guarda el sumando derecho
                suma = Cdiv(suma,aux);
            
                fcomplex prod = {-1,0};     //Variable auziliar que me guarda el factor de G[j-1]
                prod = Cdiv(prod,aux);
                prod = Cmul(prod,G[j-1]);

                G[j] = Cadd(prod,suma);

                //Aprovecho el bucle para actualizar el valor de la funcion de onda
                F[j] = Csub(G[j],F[j]);
            }

            //Obtengo la densidad de probabilidad en un punto
            prob = Cmul(F[j],Conjg(F[j])).r;

            //Escribo en el fichero los observables de la funcion de onda
            Funcion_Onda << setprecision(10) << (n*dt)  << "\t" << setprecision(10) << (j*dx) << "\t" << setprecision(10) << prob << "\t" << setprecision(10) << F[j].r << "\t" << setprecision(10) << F[j].i <<  "\n";
        }

        //Calculo y Escritura de Observables//

        //Obtengo los valores de los observables (probabilidad y valor esperado de x, x^2, p, p^2)
        Probabilidad_F_Onda(F,dx,prob_trapecio,prob_simpson);
        Valor_Medio_Posicion(F,dx, x_trapecio, x_simpson);
        Valor_Medio_Posicion_Gauss(F,dx,x_gauss);
        Valor_Medio_Posicion_Cuadrado(F,dx,x_trapecio2,x_simpson2);
        Valor_Medio_Momento_Centrado(F,dx,p_trapecio,p_simpson);
        Valor_Medio_Momento_Cuadrado(F,dx,p_trapecio2,p_simpson2);
        Valor_Medio_Potencial(F,V,dx,V_trapecio,V_simpson);

        //Obtengo Valor Medio Energia
        E_trapecio = p_trapecio2 + V_trapecio;
        E_simpson = p_simpson2 + V_simpson;

        //Obtengo las incertidumbres de p y x mediante la definicion formal
        Delta_p_trapecio = sqrt(abs(p_trapecio2 -p_trapecio*p_trapecio));
        Delta_p_simpson = sqrt(abs(p_simpson2-p_simpson*p_simpson));
        Delta_x_trapecio = sqrt(abs(x_trapecio2-x_trapecio*x_trapecio));
        Delta_x_simpson = sqrt(abs(x_simpson2-x_simpson*x_simpson));

        Heisenberg_trapecio = Delta_p_trapecio*Delta_x_trapecio;
        Heisenberg_simpson = Delta_x_simpson*Delta_p_simpson;

        Observable << setprecision(10) << dt*n << "\t" << setprecision(10) << prob_trapecio << "\t" << setprecision(15) <<  x_trapecio << "\t" << setprecision(15) <<  x_gauss << "\t" << setprecision(10) << p_trapecio << "\t" << setprecision(10) << Delta_x_trapecio << "\t" << setprecision(10) << Delta_p_trapecio << "\t" << setprecision(10) << E_trapecio << "\t" << setprecision(10) << Heisenberg_trapecio << "\n";   
        
        n++;
    }

    Funcion_Onda.close();
    Observable.close();
    Condiciones_Iniciales.close();

    cout << endl;
    cout << "He terminado :)" << endl;
    cout << endl;

    return 0;
}


/*
        INICIALIZAR FUNCION DE ONDA AUTOFUNCION

        Objetivo: Insertar los valores de la funcion de onda segun la funcion dada
        Input: Vector funcion de onda (por referencia)
        Output: Vector funcion de onda

*/
void Inicio_F_Onda_Gaussiana (vector <fcomplex>& F, double sigma, double x_0, double dx)
{
    double aux;

    //Recorro todo el epacio actualizando a los valores iniciales de la funcion de onda
    for (int j=1; j<F.size()-1; j++)
    {
        aux = -((j*dx)-x_0)*((j*dx)-x_0)/(2*sigma*sigma);
        
        F[j].r = exp(aux);      //Actualizamos solo la parte real de la funcion de onda pq la imaginaria es nula en t=0
    }

    //Impongo las condiciones iniciales
    F[0] = F[F.size()-1] = {0,0};

    return;
}

/*
        INICIALIZAR POTENCIAL ARMONICO

        Objetivo: Insertar los valores del potencial segun la funcion dada (para un potencial barrera)
        Input: Vector potencial (por referencia) // Vector potencial Vtilde (potencial reescalado, por referencia) //
               frecuencia potencial // Paso espacial // centro potencial
        Output: Vector potencial // Vector potencial Vtilde

*/
void Inicio_Potencial_Armonico (vector <double> & V, vector <double> & Vtilde, double omega, double dx, double centro)
{
    //Declaro un fichero para escribir los datos del potencial
    ofstream Potencial_Armonico ("./sigma116_centro05/Gauss_Potencial_Armonico.txt");

    double aux = omega*omega/4;

    //Recorro todo el espacio y inicializo los valores del potencial en cada punto
    for (int j = 0; j<V.size(); j++)
    {
        V[j] = aux *((dx*j - centro)*(dx*j - centro));  //Valor potencial
        Vtilde[j] = dx*dx*V[j];                         //Reescalamiento potencial

        //Escribo el valor del potencial en un fichero
        Potencial_Armonico << setprecision(10) << (dx*j) << "\t" << setprecision(10) << V[j] << "\n";
    }

    Potencial_Armonico.close();
    return;
}

/*
        NORMALIZAR FUNCION DE ONDA

        Objetivo: Normalizar la funcion de onda dada.
        Input: Vector funcion de onda (por referencia)
        Output: Vector funcion de onda normalizado
*/
void Normalizar_F_Onda (vector <fcomplex>& F, double dx)
{
    //Declaro una variable auxiliar donde guardo la norma (al cuadrado)
    double norma = 0;

    //Recorro toda la funcion de onda y sumo la probabilidad (obtengo la norma al cuadrado)
    for (int j = 0; j< F.size(); j++)
    {
        norma += Cmul(F[j],Conjg(F[j])).r;
    }

    //Pongo la norma por pantalla
    cout << endl;
    cout << "La norma ANTES de normalizar es: " << norma << endl;
    cout << endl;

    norma = norma * dx;
    
    //Divido cada componente de la funcion de onda por la norma
    for (int j = 0; j< F.size(); j++)
    {
        F[j] = Cdiv(F[j],{sqrt(norma),0});
    }

    // BLOQUE DE COMPROBACIÓN PARA QUEDARNOS TRANQUILOS
    double norma_comprobacion = 0.0;
    for (int j = 0; j< F.size(); j++)
    {
        norma_comprobacion += (F[j].r * F[j].r) + (F[j].i * F[j].i);
    }
    
    cout << "Norma DESPUES de normalizar: " << (norma_comprobacion * dx) << "\n\n";

    //Cierro la funcion
    return;
}

    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
                                            BLOQUE OBSERVABLES 
    /////////////////////////////////////////////////////////////////////////////////////////////////////*/

/*
        FUNCION PROBABILIDAD FUNCION ONDA

        Objetivo: Obtener la densidad de probabilidad integrada en todo el espacio, para cualquier valor temporal usando integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // Variables probabilidad 
        Output: Valor del valor esperado de la posicion mediante Simpson y Trapecio

*/
void Probabilidad_F_Onda (vector <fcomplex> & F, double dx, double & prob_trapecio, double & prob_simpson)  //Mediante Suma de Riemann
{
    prob_trapecio = 0;
    prob_simpson = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada), y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)    //Empiezo en 1 y S-1 pq f(0) y f(S) = 0
    {
        double aux;

        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);  //Auxiliar que me da el valor de las funciones de onda con el complejo conjugado

        //Algoritmo Trapecio (integracion normal)
        prob_trapecio += aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz.
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            prob_simpson += 2*(dx*j)*aux;
        }
        else   //Impar
        {
            prob_simpson += 4*(dx*j)*aux;
        }
    }

    prob_trapecio = prob_trapecio*dx;
    prob_simpson = prob_simpson*dx/3;

    //Devuelvo el valor de la media
    return;
}

/*
        FUNCION VALOR ESPERADO POSICION

        Objetivo: Obtener el valor esperado de la posicion de la funcion de onda, para cualquier valor temporal usando integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // 
        Output: Valor del valor esperado de la posicion mediante Simpson y Trapecio

*/
void Valor_Medio_Posicion (vector <fcomplex> & F, double dx, double & x_trapecio, double & x_simpson)  //Mediante Suma de Riemann
{
    x_trapecio = 0;
    x_simpson = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada), y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)    //Empiezo en 1 y S-1 pq f(0) y f(S) = 0
    {
        double aux;

        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);  //Auxiliar que me da el valor de las funciones de onda con el complejo conjugado

        //Algoritmo Trapecio (integracion normal)
        x_trapecio += (dx*j)*aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz.
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            x_simpson += 2*(dx*j)*aux;
        }
        else   //Impar
        {
            x_simpson += 4*(dx*j)*aux;
        }
    }

    x_trapecio = x_trapecio*dx;
    x_simpson = x_simpson*dx/3;

    //Devuelvo el valor de la media
    return;
}

/*
        FUNCION VALOR ESPERADO POSICION GAUSS

        Objetivo: Obtener el valor esperado de la posicion de la funcion de onda, para cualquier valor temporal buscando la semianchura de la funcion de onda
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // Valor Posicion Medio
        Output: Valor del valor esperado de la posicion mediante la busqueda de la semianchura

*/
void Valor_Medio_Posicion_Gauss (vector <fcomplex> & F, double dx, double & x_med)  
{

    int cont_max, cont_med;     //contadores
    double pos_max, max, med;     //Posiciones y valores maximos y medios
    max = med = pos_max = cont_max = cont_med =  0;

    for (int j = 0; j<F.size(); j++)    //algorotmo de busqueda del maximo
    {
        double aux;
        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);

        if (aux>=max)    //Actualizo los datos del maximo
        {
            max = aux;
            cont_max = j;
            pos_max = cont_max*dx;
        }
    }

    med = max/2;

    for (int j = cont_max; j < F.size(); j++)
    {
        double aux;
        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);

        if (aux<=med)
        {
            med = aux;
            cont_med = j;
            x_med = cont_med*dx-pos_max;
            break;
        }
    }

    //Devuelvo el valor de la media
    return;
}


/*
        FUNCION VALOR ESPERADO POSICION CUADRADO

        Objetivo: Obtener el valor esperado de la posicion de la funcion de onda, para cualquier valor temporal usando integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // 
        Output: Valor del valor esperado de la posicion mediante Simpson y Trapecio

*/
void Valor_Medio_Posicion_Cuadrado (vector <fcomplex> & F, double dx, double & x_trapecio2, double & x_simpson2)  //Mediante Suma de Riemann
{
    x_trapecio2 = 0;
    x_simpson2 = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada), y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)    //Empiezo en 1 y S-1 pq f(0) y f(S) = 0
    {
        double aux, x_cuad;

        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);  //Auxiliar que me da el valor de las funciones de onda con el complejo conjugado
        x_cuad = (dx*j)*(dx*j);

        //Algoritmo Trapecio (integracion normal)
        x_trapecio2 += x_cuad*aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz.
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            x_simpson2 += 2*x_cuad*aux;
        }
        else   //Impar
        {
            x_simpson2 += 4*x_cuad*aux;
        }
    }

    x_trapecio2 = x_trapecio2*dx;
    x_simpson2 = x_simpson2*dx/3;

    //Devuelvo el valor de la media
    return;
}

/*
        FUNCION VALOR MEDIO MOMENTO 

        Objetivo: Obtener el valor esperado del momento de la funcion de onda, para cualquier valor temporal y medainte integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // 
        Output: Valor del valor esperado del momento

*/
void Valor_Medio_Momento (vector <fcomplex> & F, double dx, double & p_trapecio, double & p_simpson)  //Mediante Suma de Riemann
{
    fcomplex dF;    //Variable auxiliar que guarda el valor de la derivada del momento

    p_trapecio = p_simpson = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada) excepto los extremos que al ser la funcion nula no contribuira, y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)
    {
        double aux;

        //Obtengo la derivada
        dF = Csub(F[j+1], F[j]);
        dF.r = dF.r/dx;
        dF.i = dF.i/dx;
        dF = Cmul(dF,{0,-1});
        
        aux = ((F[j].r*dF.r)+(F[j].i*dF.i));
        
        //Trapecio
        p_trapecio += aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz. 
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            p_simpson += 2*aux;
        }
        else   //Impar
        {
            p_simpson += 4*aux;
        }
    }

    //Devuelvo el valor de la media
    p_trapecio = p_trapecio*dx;
    p_simpson = p_simpson*dx/3;

    return;
}

//Funcion identica a la del momento, pero tomando la definicion de derivada centrada
void Valor_Medio_Momento_Centrado (vector <fcomplex> & F, double dx, double & p_trapecio, double & p_simpson)  //Mediante Suma de Riemann
{
    fcomplex dF;    //Variable auxiliar que guarda el valor de la derivada del momento

    p_trapecio = p_simpson = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada) excepto los extremos que al ser la funcion nula no contribuira, y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)
    {
        double aux;

        dF = Csub(F[j+1], F[j-1]);
        dF.r = dF.r/(2*dx);
        dF.i = dF.i/(2*dx);
        dF = Cmul(dF,{0,-1});

        aux = ((F[j].r*dF.r)+(F[j].i*dF.i));
            
        //Trapecio
        p_trapecio += aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz. 
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            p_simpson += 2*aux;
        }
        else   //Impar
        {
            p_simpson += 4*aux;
        }
    }

    //Devuelvo el valor de la media
    p_trapecio = p_trapecio*dx;
    p_simpson = p_simpson*dx/3;

    return;
}

/*
        FUNCION VALOR ESPERADO MOMENTO CUADRADO 

        Objetivo: Obtener el valor esperado del momento cuadrado de la funcion de onda, para cualquier valor temporal y medainte integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // 
        Output: Valor del valor esperado del momento cuadrado

*/
void Valor_Medio_Momento_Cuadrado (vector <fcomplex> & F, double dx, double & p_trapecio2, double & p_simpson2)  //Mediante Suma de Riemann
{
    fcomplex dF;    //Variable auxiliar que guarda el valor de la derivada del momento

    p_trapecio2 = p_simpson2 = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada) excepto los extremos que al ser la funcion nula no contribuira, y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)
    {
        double aux;

        dF = Csub(F[j+1], Cmul({2,0},F[j]));
        dF = Cadd(dF, F[j-1]);
        dF.r = dF.r/(dx*dx);
        dF.i = dF.i/(dx*dx);
        dF = Cmul(dF,{-1,0});
            
        aux = ((F[j].r*dF.r)+(F[j].i*dF.i));

        //Trapecio
        p_trapecio2 += aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz. 
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            p_simpson2 += 2*aux;
        }
        else   //Impar
        {
            p_simpson2 += 4*aux;
        }
    }

    //Devuelvo el valor de la media
    p_trapecio2 = p_trapecio2*dx;
    p_simpson2 = p_simpson2*dx/3;

    return;
}

/*
        FUNCION VALOR ESPERADO POTENCIAL (funcion reciclada del valor esperado de la posicion, cambiando x por V[j])

        Objetivo: Obtener el valor esperado del potencial de la funcion de onda, para cualquier valor temporal y medainte integracion Trapecio y Simpson
        Input: Funcion de onda (por referencia para no crear una copia de la funcion de onda) // Paso espacial // Potencial 
        Output: Valor del valor esperado del potencial

*/
void Valor_Medio_Potencial (vector <fcomplex> & F, vector <double> & V, double dx,  double & V_trapecio, double & V_simpson)  //Mediante Suma de Riemann
{
    //Obtengo el valor medio del potencial (reciclo el algoritmo del valor medio del )
    V_trapecio = 0;
    V_simpson = 0;

    //Recorro todo el espacio de la funcion de onda (ya normalizada), y multiplico el valor por la posicion (la integral de la posicion discretizada)
    for (int j=1; j<F.size()-1; j++)    //Empiezo en 1 y S-1 pq f(0) y f(S) = 0
    {
        double aux;

        aux = (F[j].r*F[j].r)+(F[j].i*F[j].i);  //Auxiliar que me da el valor de las funciones de onda con el complejo conjugado

        //Algoritmo Trapecio (integracion normal)
        V_trapecio += V[j]*aux;    //Como el modulo es la raiz de la suma de la parte real^2 + parte imaginaria^2, al elevar al cuadrado quito la raiz.
        
        //Algoritmo Simpson Compuesto
        if (j%2 == 0)   //Par
        {
            V_simpson += 2*V[j]*aux;
        }
        else   //Impar
        {
            V_simpson += 4*V[j]*aux;
        }
    }

    V_trapecio = V_trapecio*dx;
    V_simpson = V_simpson*dx/3;

    //Devuelvo el valor de la media
    return;
}
