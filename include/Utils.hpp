#pragma once

#include <iostream>
#include <algorithm>
#include <chrono>

using namespace std ; 

class Timer {
private:
    chrono::steady_clock::time_point _clock_time ;
public:
    Timer() ;
    void SetClock() ;
    clock_t GetDurationSeconds() const ;
    clock_t GetDurationMilliseconds() const ;
};


string Strip(const string& str) ; // 截斷字串左右兩邊的空白符(\n, ' ', '\t')
string GetTrailingNumber(const string& str) ; // 取得字串尾端的數字