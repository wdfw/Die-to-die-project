#include "Utils.hpp"

string Strip(const string& str){
    int l=0, r=str.size(), f=0 ; 
    for(int i=0; i<str.size(); i++){
        if(isspace(str[i])){
            if(!f){
                l = i+1 ;  
            }
        }else{
            if(!f) f = 1 ;
            else r = i+1 ;
        }
    }
    return str.substr(l, r-l) ;
}

string GetTrailingNumber(const string& str){
    string res ;
    for(auto c=str.crbegin(); c!=str.crend(); ++c){
        if(isdigit(*c))  res.push_back(*c) ;
        else break ;
    }
    reverse(res.begin(), res.end()) ;
    return res ;
}
