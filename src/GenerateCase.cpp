#include "Generator.hpp"

int main(int argc, char *argv[]) {

    int num_signals = stoi(argv[3]) ; //總共Signal bump的數量
    string input_file = argv[1] ; //輸入參考檔的位置
    string output_path = argv[2] ; //輸出bump file的位置
    double offset = stoi(argv[4]) ; //Die1 Die2 的垂直差
    GenerateBumpCaseFiles(num_signals, input_file, output_path, offset);
    return 0 ;
}

//修改方向
//1. 可以調整任意兩點間的寬度
//2. 不用參考檔也能生成

//若沒有result資料夾輸入
// mkdir result/RDL1

// ./bin/GenerateCase ./case/d2d_case_bump.location ./result/RDL1 
// ./bin/ShowGraph ./result

