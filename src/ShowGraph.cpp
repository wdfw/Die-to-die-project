#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>

#include "Drawer.hpp" //Drawer.hppㄧ定要include在最上層
#include "Bump.hpp"
#include "DesignRule.hpp"
#include "Generator.hpp"
#include "Parser.hpp"
#include "Utils.hpp"

using namespace std ; 


void Reload(const vector<string>& RDLDirectories, const string& designRulePath, vector<QGraphicsScene*>& sceneCache){
    Drawer drawer ; 
    string viaPath ;
    string offsetViaPath ;
    string netListPath ;
    string teardropPath ;

    vector<Bump> bumps, offsetBumps1, offsetBumps2 ; // 在 bump layer 上的點
    vector<double> coordinate ;  // min_x, min_y, max_x, max_y
    vector<Net> nets ; 
    vector<tuple<Bump, double, double, double, double>> teardrops ;
    vector<Bump> allBumps ;
    DesignRule designRule;

    // if(RDLDirectories.size() != sceneCache.size()) throw runtime_error("Different size of scnces and directories.") ; 

    for(int i=0; i<sceneCache.size(); i++) sceneCache[i]->clear() ; 

    ParseDesignRule(designRulePath, designRule) ; drawer.SetDesignRule(&designRule) ;

    drawer.SetDesignScence(sceneCache[0]) ;
    ParseBump(RDLDirectories[0] + "/via_layer_" + GetTrailingNumber(RDLDirectories[0]), bumps, coordinate) ;
    for(int i=0; i<bumps.size(); i++) drawer.DrawBump(bumps[i]) ; 

    for(int i=0; i<RDLDirectories.size(); i++){
        string dir = RDLDirectories[i] ;
        string layerNumber = GetTrailingNumber(dir) ; //取得RDLX的數字部分(e.g. X的部分)

        viaPath = dir + "/via_layer_" + layerNumber ;
        offsetViaPath = dir + "/offset_via_layer_" + layerNumber ; 
        netListPath = dir + "/netlist_" + layerNumber ; 
        teardropPath = dir + "/teardrop_" + layerNumber ; 
        
        drawer.SetDesignScence(sceneCache[i+1]) ;

        if(filesystem::exists(teardropPath)){
            ParseTeardrop(teardropPath, teardrops) ;
            for(int i=0; i<teardrops.size(); i++) drawer.DrawTeardrop(teardrops[i]) ;
        }

        if(filesystem::exists(viaPath)){
            ParseBump(viaPath, bumps, coordinate) ;
            for(int i=0; i<bumps.size(); i++) drawer.DrawBump(bumps[i]) ; 
            // drawer.DrawDieBoundary(coordinate) ; //檔案中die boundary僅參考用, 所以不畫了
        }

        if(filesystem::exists(offsetViaPath)){
            ParseOffsetBump(offsetViaPath, offsetBumps1, offsetBumps2) ;
            for(int i=0; i<offsetBumps1.size(); i++) drawer.DrawOffsetBump(offsetBumps1[i], offsetBumps2[i]) ; 
        }

        if(filesystem::exists(netListPath)){
            ParseNet(netListPath, nets) ;
            for(int i=0; i<nets.size(); i++) drawer.DrawNet(nets[i]) ;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc!=3) { //!=3
        cerr << "Usage: " << argv[0] << " <result_folder> <design_rule_file>" << "\n";
        return 1 ;
    }

    QApplication app(argc, argv) ;
    QGraphicsScene scene ; scene.setSceneRect(0, 0, 5000, 2000) ;
    QtView view(&scene) ; view.setStyleSheet("background-color: black;") ;
    QHBoxLayout *buttonLayout = new QHBoxLayout;
    QVBoxLayout *layout = new QVBoxLayout;
    QWidget mainWidget;
    QPushButton *button = new QPushButton("Bump");
    vector<QGraphicsScene*> sceneCache ;

    string inputDirectory = argv[1] ; 
    string designRulePath = argv[2] ;
    vector<string> RDLDirectories ;
    
    //走訪放置RDLX的資料夾
    for(auto& dir : filesystem::directory_iterator(inputDirectory)){
        if(string(dir.path()).find("RDL")==string::npos) throw runtime_error("Invalid folder " + string(dir.path())) ;
        RDLDirectories.push_back(dir.path()) ; 
    }
    if(!RDLDirectories.size()) throw runtime_error("Empty folder " + inputDirectory) ;
    sort(RDLDirectories.begin(), RDLDirectories.end()) ;

    sceneCache.resize(RDLDirectories.size()+1) ;
    for(int i=0; i<RDLDirectories.size()+1; i++){
        sceneCache[i] = new QGraphicsScene() ; 
        sceneCache[i]->setSceneRect(0, 0, 5000, 2000); 
    }

    view.resetTransform();
    view.setScene(sceneCache[0]);  // 在一開始顯示 bump layer 的場景
    view.centerOn(0, 0);

    // 生成 RDL 按鈕及場景
    QObject::connect(button, &QPushButton::clicked, [&]() { view.resetTransform(); view.setScene(sceneCache[0]); }); buttonLayout->addWidget(button);
    for (int i = 1; i < sceneCache.size(); ++i) {
        QPushButton *button = new QPushButton("RDL " + QString::number(i));
        QObject::connect(button, &QPushButton::clicked, [i, &view, &sceneCache]() { view.resetTransform(); view.setScene(sceneCache[i]); });
        buttonLayout->addWidget(button);
    }

    // 生成執行reload的按鈕
    QPushButton *reloadButton = new QPushButton("Reload");
    QObject::connect(reloadButton, &QPushButton::clicked, [&RDLDirectories, &designRulePath, &sceneCache](){ Reload(RDLDirectories, designRulePath, sceneCache); });
    buttonLayout->addWidget(reloadButton);

    // 創建主窗口布局
    layout->addWidget(&view);
    layout->addLayout(buttonLayout);
    mainWidget.setLayout(layout);

    // 顯示主窗口
    mainWidget.resize(1000, 1000);
    mainWidget.show();
    
    // 先自動載入一次, 之後可以點擊reload來重載
    Reload(RDLDirectories, designRulePath, sceneCache);

    return app.exec();
    
}
// ./bin/ShowResult result/ case/design.rule 
//這裡的via N代表第N與N-1層的via