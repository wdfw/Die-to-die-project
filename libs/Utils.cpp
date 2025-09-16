#include "Utils.hpp"

//------------------------------------------ Timer Method Begin ------------------------------------------ 

Timer::Timer() {
    SetClock();
}

void Timer::SetClock() {
    _clock_time = chrono::steady_clock::now();
}

clock_t Timer::GetDurationSeconds() const {
    return chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - _clock_time).count();
}

clock_t Timer::GetDurationMilliseconds() const {
    return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - _clock_time).count();
}

//------------------------------------------ Timer Method Begin ------------------------------------------ 


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
