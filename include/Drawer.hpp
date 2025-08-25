#pragma once

#include <map>

#include <QApplication>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QWheelEvent>
#include <QKeyEvent>
#include <QScrollBar>
#include <QMouseEvent>
#include <QPushButton>
#include <QVBoxLayout>
#include <QWidget>
#include <QGraphicsProxyWidget>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/strategies/buffer.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include <boost/geometry/algorithms/closest_points.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include "Bump.hpp"
#include "DesignRule.hpp"

using point = boost::geometry::model::d2::point_xy<double>;
using polygon = boost::geometry::model::polygon<point>;
using multi_point = boost::geometry::model::multi_point<point>;
using multi_polygon = boost::geometry::model::multi_polygon<polygon>;
using linestring = boost::geometry::model::linestring<point>;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Triangulation;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using box = bg::model::box<point>;
using value = std::pair<box, size_t>;
using namespace std;

// 畫 bump layer
void QtDrawBump(QGraphicsScene *scene, vector<Bump> &bumps, vector<double> &routing_area_coordinate, DesignRule designRule);

// // 畫 VSS bump 和 VDD bump
// void drawPG(QGraphicsScene *scene, vector<Bump> &bumps);

// // 畫 routing area
// void drawRoutingArea(QGraphicsScene *scene, vector<double> &routing_area_coordinate);

// // 畫 RDL layer
// void QtDrawRDL(QGraphicsScene *scene, vector<Bump> &vias, vector<double> &routing_area_coordinate, DesignRule designRule, const string &netlistFile, vector<pair<double, double>>& triangleEdgeSource, vector<pair<double, double>>& triangleEdgeTarget, vector<Bump>& offset_vias, RoutingGraph& RDL, vector<wireLength>& wireLengths);

// // 畫 via offset (橢圓形結構)
// void QtDrawOffsetVia(QGraphicsScene *scene, vector<Bump>& offset_vias, DesignRule designRule);

class Drawer {
private:
    QGraphicsScene *scene ;
    DesignRule *designRule ;

    void _DrawTeardrop(const point& center1, const point& center2) ;
public:
    bg::strategy::buffer::distance_symmetric<double> viaRadiutStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::distance_symmetric<double> bufferDistanceStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::distance_symmetric<double> trackDistanceStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::join_round joinStrategy = bg::strategy::buffer::join_round(0) ;
    bg::strategy::buffer::end_round endStrategy = bg::strategy::buffer::end_round(0) ;
    bg::strategy::buffer::point_circle circleStrategy = bg::strategy::buffer::point_circle(0) ;
    bg::strategy::buffer::side_straight sideStrategy = bg::strategy::buffer::side_straight() ;

    Drawer(DesignRule *designRule, QGraphicsScene *scene) ;
    void SetDesignRule(DesignRule *designRule) ; 
    void SetDesignScence(QGraphicsScene *scene) ; 

    void DrawBump(Bump &bump, Qt::GlobalColor bumpColor=Qt::yellow) ; // yellow for default indiactor (will use other colors to show bumps)
    void DrawOffsetBump(Bump &bump, Bump &matchedBump) ;
    void DrawNets(Net& net, Qt::GlobalColor bumpColor=Qt::yellow) ; // 目前只能塗特定顏色, 因為net list資訊刮
    void DrawTeardrop(tuple<Bump, double, double, double, double>& teardrop) ; //Bump與最近的線段接Teardrop
 
    

    // void DrawNet(Bump &bump) ; 

    // void DrawOneRDL(vector<Bump> &bumps, vector<double> &coordinate) ; 
    // void DrawAllRDL(vector<Bump> &bumps, vector<double> &coordinate) ; 
} ;

// 視窗相關的設定
class QtView : public QGraphicsView {
public:
    QtView(QGraphicsScene *scene) : QGraphicsView(scene) {
        setRenderHint(QPainter::Antialiasing);
        setDragMode(QGraphicsView::ScrollHandDrag);
        setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
        setStyleSheet("background-color: black;");
    }

protected:
    bool rightMousePressed = false;
    QPoint lastMousePos;

    // 滑鼠按下事件
    void mousePressEvent(QMouseEvent *event) override {
        if (event->button() == Qt::RightButton) {
            rightMousePressed = true;
            lastMousePos = event->pos();
            setCursor(Qt::ClosedHandCursor);
        } else {
            QGraphicsView::mousePressEvent(event);
        }
    }

    // 滑鼠移動事件
    void mouseMoveEvent(QMouseEvent *event) override {
        if (rightMousePressed) {
            QPoint delta = event->pos() - lastMousePos;
            lastMousePos = event->pos();
            horizontalScrollBar()->setValue(horizontalScrollBar()->value() - delta.x());
            verticalScrollBar()->setValue(verticalScrollBar()->value() - delta.y());
        } else {
            QGraphicsView::mouseMoveEvent(event);
        }
    }

    // 滑鼠釋放事件
    void mouseReleaseEvent(QMouseEvent *event) override {
        if (event->button() == Qt::RightButton) {
            rightMousePressed = false;
            setCursor(Qt::ArrowCursor);
        } else {
            QGraphicsView::mouseReleaseEvent(event);
        }
    }


    void wheelEvent(QWheelEvent *event) override {
        if (event->modifiers() & Qt::ControlModifier) {
            double factor = (event->angleDelta().y() > 0) ? 1.1 : 0.9;
            scale(factor, factor);
        } else {
            QGraphicsView::wheelEvent(event);
        }
    }

    void keyPressEvent(QKeyEvent *event) override {
        if (event->key() == Qt::Key_Z) {
            resetTransform();  // 重置縮放和平移
            // centerOn(0, 0);     // 回到初始座標位置 (0, 0)
        } else {
            QGraphicsView::keyPressEvent(event);
        }
    }

};

