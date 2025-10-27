#include "Generator.hpp"

void GenerateBumpCaseFiles(int num_signals, const string &case_path, const string &via_path, double &offset){
    // offset = 2.5;
    const double Triangle_Length = 7.0;
    const int max_signal_rows = 4; // 可以是 4, 6, 8, 10, ... => 奇數的話, bump 的排列會跟原本不同
    const double Ver_SPACING_Y = Triangle_Length * sin(M_PI/3.0);
    const double Hor_SPACING_X = Triangle_Length; //兩個bump間的水平長度
    const double Die_gap = 3.0 * Hor_SPACING_X;

    const double Die1_START_X = 15.0;
    const double Die1_START_Y = 15.0;

    vector<tuple<double,double,int>> die1Signal;
    vector<pair<double,double>> die1VSS, die1VDD;

    int columns = static_cast<int>(ceil((double)num_signals / max_signal_rows));
    int col_last = columns - 1;

    // Die-1 VSS
    for (int col = 0; col < columns; ++col) {
        if (col % 2 == 1) {
            die1VSS.emplace_back(Die1_START_X + col * Hor_SPACING_X,
                                 Die1_START_Y);
        }
        die1VSS.emplace_back(Die1_START_X + col * Hor_SPACING_X + Hor_SPACING_X/2.0,
                             Die1_START_Y + Ver_SPACING_Y);
    }

    // Die-1 Signal
    for (int i = 0; i < num_signals; ++i) {
        int col = i / max_signal_rows;
        int row = i % max_signal_rows;
        double offset_x = (row % 2) * (Hor_SPACING_X/2.0);
        double x = Die1_START_X + col * Hor_SPACING_X + offset_x;
        double y = Die1_START_Y + (2 + row) * Ver_SPACING_Y;
        die1Signal.emplace_back(x, y, i + 1);
    }

    // Die-1 VDD
    for (int col = 0; col < columns; ++col) {
        if (col % 2 == 1) {
            die1VDD.emplace_back(Die1_START_X + col * Hor_SPACING_X,
                                 Die1_START_Y + (2 + max_signal_rows) * Ver_SPACING_Y);
        }
        die1VDD.emplace_back(Die1_START_X + col * Hor_SPACING_X + Hor_SPACING_X/2.0,
                             Die1_START_Y + (2 + max_signal_rows + 1) * Ver_SPACING_Y);
    }

    // Die-1 special cases
    {
        // left side VSS (row1, col0)
        double x0 = Die1_START_X + Hor_SPACING_X/2.0 - Hor_SPACING_X;
        double y0 = Die1_START_Y + Ver_SPACING_Y;
        die1VSS.emplace_back(x0, y0);

        // left side VSS (row3, col0)
        double x1 = Die1_START_X + Hor_SPACING_X/2.0 - Hor_SPACING_X;
        double y1 = Die1_START_Y + (2 + 1) * Ver_SPACING_Y;
        die1VSS.emplace_back(x1, y1);

        // right side VSS (row2, col_last)
        double x2 = Die1_START_X + col_last * Hor_SPACING_X + Hor_SPACING_X;
        double y2 = Die1_START_Y + (2 + 0) * Ver_SPACING_Y;
        die1VSS.emplace_back(x2, y2);

        // right side VDD (row4, col_last)
        double x3 = Die1_START_X + col_last * Hor_SPACING_X + Hor_SPACING_X;
        double y3 = Die1_START_Y + (2 + 2) * Ver_SPACING_Y;
        die1VDD.emplace_back(x3, y3);
    }

    // Die-2 start coordinates
    double signal_row4_col_last_x = Die1_START_X + col_last * Hor_SPACING_X;
    double Die2_START_X = signal_row4_col_last_x + 2.0 * Hor_SPACING_X + Die_gap;
    double Die2_START_Y = Die1_START_Y;

    vector<tuple<double,double,int>> die2Signal;
    vector<pair<double,double>> die2VSS, die2VDD;

    // Die-2 VSS
    for (int col = 0; col < columns; ++col) {
        if (col % 2 == 1) {
            die2VSS.emplace_back(Die2_START_X + col * Hor_SPACING_X,
                                 Die2_START_Y);
        }
        die2VSS.emplace_back(Die2_START_X + col * Hor_SPACING_X + Hor_SPACING_X/2.0,
                             Die2_START_Y + Ver_SPACING_Y);
    }

    // Die-2 Signal
    for (int i = 0; i < num_signals; ++i) {
        int col = i / max_signal_rows;
        int row = i % max_signal_rows;
        double offset_x = (row % 2) * (Hor_SPACING_X/2.0);
        double x = Die2_START_X + col * Hor_SPACING_X + offset_x;
        double y = Die2_START_Y + (2 + row) * Ver_SPACING_Y;
        die2Signal.emplace_back(x, y, i + 1);
    }

    // Die-2 VDD
    for (int col = 0; col < columns; ++col) {
        if (col % 2 == 1) {
            die2VDD.emplace_back(Die2_START_X + col * Hor_SPACING_X,
                                 Die2_START_Y + (2 + max_signal_rows) * Ver_SPACING_Y);
        }
        die2VDD.emplace_back(Die2_START_X + col * Hor_SPACING_X + Hor_SPACING_X/2.0,
                             Die2_START_Y + (2 + max_signal_rows + 1) * Ver_SPACING_Y);
    }

    // write helper lambdas, 寫入檔案時, 數值都會被放大 10 倍, 在 QT 顯示時, 點才不會擠在一起
    auto printBoundary = [&](ostream &os, double x, double y){
        os << fixed << setprecision(1) << x*10 << " "
           << fixed << setprecision(1) << y*10 << "\n";
    };
    auto printXY = [&](ostream &os, double x, double y){
        os << fixed << setprecision(1) << x*10 << " "
           << defaultfloat << setprecision(numeric_limits<double>::max_digits10) << y*10;
    };

    // write case file
    {
        ofstream f(case_path);
        if (!f) throw runtime_error("Failed to open file: " + case_path);
        double min_x = die1VSS[0].first, min_y = die1VSS[0].second;
        for (auto &p : die1VSS) { min_x = min(min_x, p.first); min_y = min(min_y, p.second); }
        double max_x = die2VDD[0].first, max_y = die2VDD[0].second;
        for (auto &p : die2VDD) { max_x = max(max_x, p.first); max_y = max(max_y, p.second); }
        max_x += 2.0 * Hor_SPACING_X;
        printBoundary(f, min_x, min_y);
        f << fixed << setprecision(1) << max_x*10 << " "
          << defaultfloat << setprecision(numeric_limits<double>::max_digits10) << max_y*10 << "\n";
        for (auto &[x, y, id] : die1Signal) {
            f << "DIE1 SIG " << (id - 1) << " ";
            printXY(f, x, y);
            f << "\n";
        }
        for (size_t i = 0; i < die1VDD.size(); ++i) {
            f << "DIE1 VDD " << i << " ";
            printXY(f, die1VDD[i].first, die1VDD[i].second);
            f << "\n";
        }
        for (size_t i = 0; i < die1VSS.size(); ++i) {
            f << "DIE1 VSS " << i << " ";
            printXY(f, die1VSS[i].first, die1VSS[i].second);
            f << "\n";
        }
        for (auto &[x, y, id] : die2Signal) {
            f << "DIE2 SIG " << (id - 1) << " ";
            printXY(f, x, y + offset);
            f << "\n";
        }
        for (size_t i = 0; i < die2VDD.size(); ++i) {
            f << "DIE2 VDD " << i << " ";
            printXY(f, die2VDD[i].first, die2VDD[i].second + offset);
            f << "\n";
        }
        for (size_t i = 0; i < die2VSS.size(); ++i) {
            f << "DIE2 VSS " << i << " ";
            printXY(f, die2VSS[i].first, die2VSS[i].second + offset);
            f << "\n";
        }
    }

    // write via file (same content)
    {
        filesystem::create_directories(filesystem::path(via_path).parent_path()); // 自動創建資料夾
        ofstream f(via_path);
        if (!f) throw runtime_error("Failed to open file: " + via_path);
        double min_x = die1VSS[0].first, min_y = die1VSS[0].second;
        for (auto &p : die1VSS) { min_x = min(min_x, p.first); min_y = min(min_y, p.second); }
        double max_x = die2VDD[0].first, max_y = die2VDD[0].second;
        for (auto &p : die2VDD) { max_x = max(max_x, p.first); max_y = max(max_y, p.second); }
        max_x += 2.0 * Hor_SPACING_X;
        printBoundary(f, min_x, min_y);
        f << fixed << setprecision(1) << max_x*10 << " "
          << defaultfloat << setprecision(numeric_limits<double>::max_digits10) << max_y*10 << "\n";
        for (auto &[x, y, id] : die1Signal) {
            f << "DIE1 SIG " << (id - 1) << " ";
            printXY(f, x, y);
            f << "\n";
        }
        for (size_t i = 0; i < die1VDD.size(); ++i) {
            f << "DIE1 VDD " << i << " ";
            printXY(f, die1VDD[i].first, die1VDD[i].second);
            f << "\n";
        }
        for (size_t i = 0; i < die1VSS.size(); ++i) {
            f << "DIE1 VSS " << i << " ";
            printXY(f, die1VSS[i].first, die1VSS[i].second);
            f << "\n";
        }
        for (auto &[x, y, id] : die2Signal) {
            f << "DIE2 SIG " << (id - 1) << " ";
            printXY(f, x, y + offset);
            f << "\n";
        }
        for (size_t i = 0; i < die2VDD.size(); ++i) {
            f << "DIE2 VDD " << i << " ";
            printXY(f, die2VDD[i].first, die2VDD[i].second + offset);
            f << "\n";
        }
        for (size_t i = 0; i < die2VSS.size(); ++i) {
            f << "DIE2 VSS " << i << " ";
            printXY(f, die2VSS[i].first, die2VSS[i].second + offset);
            f << "\n";
        }
    }
}
