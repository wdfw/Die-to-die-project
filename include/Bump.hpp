#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <utility>

using namespace std ; 

enum DieType {
    DUMMY, 
    SIGNAL, 
    VDD,
    VSS,
    OFFSET
} ;

class Bump ; 
class Net ; 

ostream& operator<<(ostream& os, const Bump& bump) ;
ostream& operator<<(ostream& os, const Net& bump) ;

string DieType2Str(const DieType& type) ;
DieType Str2DieType(const string& str) ;

string Bump2Str(const Bump& bump) ;
string Bump2Net(const Net& net) ;

class Bump {
public:
    string name;
    DieType type;
    int id;
    double x, y;

    Bump(const string& name="", DieType type=DUMMY, int id=-1, double x=0.0, double y=0.0) : name(name), type(type), id(id), x(x), y(y) {};
    bool operator==(const Bump& bump) const ;  
} ;

class Net : public vector<tuple<double,double,double,double>> {
public:
    string name ;
    Net(const string& name="", const vector<tuple<double,double,double,double>>& sequence=vector<tuple<double,double,double,double>>()) : name(name), vector<tuple<double,double,double,double>>(sequence) {} ;
} ;

class Teardrop : public Bump {
public:
    double sx, sy, ex, ey ; 
    Teardrop(const string& name="", DieType type=DUMMY, int id=-1, double sx=0.0, double sy=0.0, double ex=0.0, double ey=0.0) : 
       Bump(name, type, id), sx(sx), sy(sy), ex(ex), ey(ey) {}

} ;
// Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);