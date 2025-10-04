#pragma once

#include <iostream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>

using namespace std ; 

class Timer ;
string Strip(const string& str) ; // 截斷字串左右兩邊的空白符(\n, ' ', '\t')
string GetTrailingNumber(const string& str) ; // 取得字串尾端的數字

template<typename T>
void FindElementsInEachRow(const vector<T>& elements, vector<vector<T>>& rowElements, double epsilonY = 1e-6){
    vector<T> sortedElements = elements ; sort(sortedElements.begin(), sortedElements.end(), [](const T& a, const T& b) {return a.y < b.y;});
    
    rowElements.clear() ; 
    for (const auto& element : sortedElements) {
        if(!rowElements.size()) rowElements.push_back({element}) ;
        else if(fabs(rowElements.back().back().y - element.y) >= epsilonY) rowElements.push_back({element}) ;
        else rowElements.back().push_back(element) ;
    }
}

template<typename T>
void FindElementsInEachRow(const vector<shared_ptr<T>>& elements, vector<vector<shared_ptr<T>>>& rowElements, double epsilonY = 1e-6){
    vector<shared_ptr<T>> sortedElements = elements ; sort(sortedElements.begin(), sortedElements.end(), [](const shared_ptr<T>& a, const shared_ptr<T>& b) {return a->y < b->y;});
    rowElements.clear() ; 
    for (const auto& element : sortedElements) {
        if(!rowElements.size()) rowElements.push_back({element}) ;
        else if(fabs(rowElements.back().back()->y - element->y) >= epsilonY) rowElements.push_back({element}) ;
        else rowElements.back().push_back(element) ;
    }
}

template<typename T>
void Matrixization(const vector<T>& elements, vector<vector<T>>& matrix, double epsilonY = 1e-6){
    FindElementsInEachRow(elements, matrix) ; 
    for(auto& row : matrix) sort(row.begin(), row.end(), [](const T& a, const T& b) {return a.x < b.x;}) ;
}

template<typename T>
void Matrixization(const vector<shared_ptr<T>>& elements, vector<vector<shared_ptr<T>>>& matrix, double epsilonY = 1e-6){
    FindElementsInEachRow(elements, matrix) ; 
    for(auto& row : matrix) sort(row.begin(), row.end(), [](const shared_ptr<T>& a, const shared_ptr<T>& b) {return a->x < b->x;}) ;
}

template<typename T>
void FindLeftMostInEachRow(const vector<T>& elements, vector<T>& leftMostElements, double epsilonY = 1e-6){
    vector<vector<T>> rowElements ;
    leftMostElements.clear() ; 
    
    FindElementsInEachRow(elements, rowElements, epsilonY) ; 
    for(const auto& row : rowElements) 
        leftMostElements.push_back( *min_element(row.begin(), row.end(), [](const T& a, const T& b) {return a.x < b.x;}) );
}


class Timer {
private:
    chrono::steady_clock::time_point _clock_time ;
public:
    Timer() ;
    void SetClock() ;
    clock_t GetDurationSeconds() const ;
    clock_t GetDurationMilliseconds() const ;
};




