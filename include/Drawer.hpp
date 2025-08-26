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

class Drawer {
private:
    bg::strategy::buffer::distance_symmetric<double> viaRadiutStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::distance_symmetric<double> bufferDistanceStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::distance_symmetric<double> trackDistanceStrategy = bg::strategy::buffer::distance_symmetric<double>(0.0) ;
    bg::strategy::buffer::join_round joinStrategy = bg::strategy::buffer::join_round(0) ;
    bg::strategy::buffer::end_round endStrategy = bg::strategy::buffer::end_round(0) ;
    bg::strategy::buffer::point_circle circleStrategy = bg::strategy::buffer::point_circle(0) ;
    bg::strategy::buffer::side_straight sideStrategy = bg::strategy::buffer::side_straight() ;
    double teaddropRadius = 0.0 ;

    const DesignRule *designRule = nullptr ;
    QGraphicsScene *scene = nullptr ;

    void _DrawTeardrop(const point& center1, const point& center2) ; // 後端負責畫Teardrop的函數, 給定兩個座標點即可畫出由center1->center2的Teardrop
    
public:
    Drawer(const DesignRule *designRule=nullptr, QGraphicsScene *scene=nullptr) ;
    void SetDesignRule(const DesignRule *designRule) ; // 設定Design rule, 並根據規則生成緩衝區策略
    void SetDesignScence(QGraphicsScene *scene) ; // 設定要畫的QT場景

    void DrawBump(const Bump &bump, Qt::GlobalColor bumpColor=Qt::yellow) ; // 畫單個Bump, bumpColor可以指定顏色, 當bumpColor為黃色時會根據Bump的種類畫預設的顏色(e.g. 不能畫黃色的Bump)
    void DrawOffsetBump(const Bump &bump, const Bump &matchedBump) ; // 畫一組Bump與其Offset via, 只能畫預設顏色, 顏色較深的為Offset via, 反之為Bump
    void DrawNet(const Net& net, Qt::GlobalColor netColor=Qt::yellow) ; // 畫單條Net, netColor可以指定顏色, 當netColor為黃色時會根據Net的種類畫預設的顏色(e.g. 不能畫黃色的Net)
    void DrawTeardrop(const tuple<Bump, double, double, double, double>& teardrop) ; // 畫一顆Teardrop
    void DrawDieBoundary(const vector<double>& coordinate) ; // 畫Die的邊界, 目前邊界僅為參考用

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

