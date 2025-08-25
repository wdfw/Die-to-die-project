#include "Utils.hpp"

string GetTrailingNumber(const string& str){
    string res ;
    for(auto c=str.crbegin(); c!=str.crend(); ++c){
        if(isdigit(*c))  res.push_back(*c) ;
        else break ;
    }
    reverse(res.begin(), res.end()) ;
    return res ;
}
