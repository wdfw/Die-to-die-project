#include <iostream>
#include <algorithm>

using namespace std ; 

string Strip(const string& str) ; // 截斷字串左右兩邊的空白符(\n, ' ', '\t')
string GetTrailingNumber(const string& str) ; // 取得字串尾端的數字