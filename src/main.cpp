#include <iostream>

#include "Drawer.hpp"
#include "Bump.hpp"
#include "DesignRule.hpp"
#include "Generator.hpp"
#include "Parser.hpp"

//Drawer.hppㄧ定要include在最上層

using namespace std ; 
namespace bg = boost::geometry ;

int main(int argc, char *argv[]) {

    // return 0 ;
    if (argc<=3) {
        cerr << "Usage: " << argv[0] << " <input_file> <design_rule_file> <output_path> [#_of_signals] [offset]" << "\n";
        return 1 ;
    }

    string viaPath = argv[1] ; 
    string offsetViaPath ;
    string designRulePath = argv[2] ;
    string outputPath = argv[3] ;
    string netListPath ;
    string teardropPath ;
    int numOfSignal = 4<argc ? stoi(argv[4]) : 0 ; 
    double offset = 5<argc ? stoi(argv[5]) : 0 ; 
    
    //Step 1. Generate Bump File
    // GenerateBumpCaseFiles(numOfSignal, viaPath, outputPath + "RDL1" + "/via_layer_1", offset) ;

    //Step 2. D2D Routing 
    vector<Bump> bumps, offsetBumps1, offsetBumps2 ; // 在 bump layer 上的點
    vector<double> coordinates ;  // min_x, min_y, max_x, max_y
    vector<Net> nets ; 
    vector<tuple<Bump, double, double, double, double>> teardrops ;
    DesignRule designRule;
    
    viaPath = "result/RDL2/via_layer_2" ;
    offsetViaPath = "result/RDL2/offset_via_layer_2" ; 
    designRulePath = "case/design.rule" ;
    netListPath = "result/RDL2/netlist_2" ;
    teardropPath = "result/RDL2/teardrop_2" ;

    ParseBump(viaPath, bumps, coordinates) ;
    ParseDesignRule(designRulePath, designRule) ;
    ParseOffsetBump(offsetViaPath, offsetBumps1, offsetBumps2) ;
    ParseNet(netListPath, nets) ;
    ParseTeardrop(teardropPath, teardrops) ;
    // // ParseBump(viaPath, bumps, coordinates) ;
    // // ParseDesignRule(designRulePath, designRule);

    QApplication app(argc, argv) ;
    QGraphicsScene scene ; scene.setSceneRect(0, 0, 5000, 2000) ;
    QtView view(&scene) ; view.setStyleSheet("background-color: black;") ;
    QHBoxLayout *buttonLayout = new QHBoxLayout;
    QVBoxLayout *layout = new QVBoxLayout;
    QWidget mainWidget;
    QPushButton *button = new QPushButton("Bump");
    QGraphicsScene* cachedScene = new QGraphicsScene(); cachedScene->setSceneRect(0, 0, 5000, 2000);

    Drawer drawer(&designRule, cachedScene) ;
    unordered_map<int, QGraphicsScene*> sceneCache;

    vector<Bump> allBumps = bumps ; allBumps.insert(allBumps.end(), offsetBumps1.begin(), offsetBumps1.end()) ;

    for(int i=0; i<teardrops.size(); i++) drawer.DrawTeardrop(teardrops[i]) ;
    for(int i=0; i<bumps.size(); i++) drawer.DrawBump(bumps[i]) ; 
    for(int i=0; i<offsetBumps1.size(); i++) drawer.DrawOffsetBump(offsetBumps1[i], offsetBumps2[i]) ;
    for(int i=0; i<nets.size(); i++) drawer.DrawNets(nets[i]) ;

    sceneCache[0] = cachedScene ;
    view.resetTransform();
    view.setScene(sceneCache[0]);  // 在一開始顯示 bump layer 的場景
    view.centerOn(0, 0);
    
    QObject::connect(button, &QPushButton::clicked, [&]() { view.resetTransform(); view.setScene(sceneCache[0]); }); buttonLayout->addWidget(button);
    
    // 生成 RDL 按鈕及場景
    for (int i = 1; i <= 1; ++i) {
        QPushButton *button = new QPushButton("RDL " + QString::number(i));
        // 根據按鈕的編號綁定相應的場景繪製函數
        QObject::connect(button, &QPushButton::clicked, [i, &view, &sceneCache]() { view.resetTransform(); view.setScene(sceneCache[i]); });
        buttonLayout->addWidget(button);
    }

    // 創建主窗口布局

    layout->addWidget(&view);
    layout->addLayout(buttonLayout);
    mainWidget.setLayout(layout);

    // 顯示主窗口
    mainWidget.resize(1000, 1000);
    mainWidget.show();

    return app.exec();
    
}
//./bin/D2D case/d2d_case_bump.location case/design.rule d2d_result 48
//這裡的via N代表第N與N-1層的via