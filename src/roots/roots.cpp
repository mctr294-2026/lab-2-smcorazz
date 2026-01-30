#include <functional>

/* Tries to find a zero crossing in f() in the interval [a,b] with the bisection method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root){
    double fa = f(a);
    double fb = f(b);
    double c=0;
    bool apos = false , bpos = false;
    int iterations = trunc(1.5*log2(abs(b-a)/(1e-6)));
    if (abs(fa)<1e-6){
        *root = a;
        return true;
    }
    if ((abs(fb)<1e-6)){
        *root = b;
        return true;
    }
    
    for (int i = 0; i<=iterations;i++){
        fa = f(a),fb = f(b);
        c = (a+b)/2;
        if (abs(f(c))<1e-6){
            *root = c;
            return true;
        }
        if (f(a) * f(c) > 0) {
        a = c;
        } else {
        b = c;
        }
    }
    return false;
        
    }
    
/* Tries to find a zero crossing in f() in the interval [a,b] with the
 * false positive / regula falsi method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root){
    double fa = f(a);
    double fb = f(b);
    double c=0;
    bool apos = false , bpos = false;
    int iterations = 1e6;
    if (abs(fa)<1e-6){
        *root = a;
        return true;
    }
    if ((abs(fb)<1e-6)){
        *root = b;
        return true;
    }
    
    for (int i = 0; i<=iterations;i++){
        fa = f(a),fb = f(b);
        if (abs(fa) < 1e-6) { *root = a; return true; }
        if (abs(fb) < 1e-6) { *root = b; return true; }
        c = a-(fa*(b-a)/(fb-fa));
        if (abs(f(c))<1e-6){
            *root = c;
            return true;
        }
        if (f(a) * f(c) > 0) {
        a = c;
        fa = f(c);
        fb *= 0.5;
        } else {
        b = c;
        fb = f(c);
        fa *= 0.5;
        }
        
    }
    return false;
        
    }

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the netwon-raphson method, given a function that computes the
 * derivative g() and a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or the derivative is zero.
 */
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root){
    if (abs(f(a))<1e-6){
        *root = a;
        return true;
    }
    else if ((abs(f(b))<1e-6)){
        *root = b;
        return true;
    }
    double x_n =c;
    double x_np1;
    for(int i=0;i<1e6;i++){
        if (x_n < a || x_n > b) {
            break;
        }
        if (abs(f(x_n))<1e-6){
            *root = x_n;
            return true;
        }
        if( (abs(g(x_n)) < 1e-12)){
            break;
        }
        
        x_np1 = x_n -(f(x_n)/g(x_n));
        x_n= x_np1;
    }
    return false;
    
                    }

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the secant method, given a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or derivative is zero.
 */
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root){
    if (abs(f(a))<1e-6){
        *root = a;
        return true;
    }
    else if ((abs(f(b))<1e-6)){
        *root = b;
        return true;
    }
    double x_nm1 = c;
    double x_n = c + 1e-6;
    double x_np1;
    for(int i=0;i<1e6;i++){
        if (x_n < a || x_n > b) {
            break;
        }
        if (abs(f(x_n))<1e-6){
            *root = x_n;
            return true;
        }
        
        if (abs(f(x_n)-f(x_nm1))<1e-12){
            break;
        }
        x_np1 = x_n -(f(x_n)*(x_n-x_nm1)/(f(x_n)-f(x_nm1)));
        x_n= x_np1;
    }
    
    return false;
            }