#include "Drawer.hpp"


map<DieType, Qt::GlobalColor> BumpColor = {{DUMMY, Qt::gray}, {SIGNAL, Qt::red}, {VDD, Qt::blue}, {VSS, Qt::green}} ;
map<DieType, Qt::GlobalColor> OffsetBumpColor = {{DUMMY, Qt::darkGray}, {SIGNAL, Qt::darkRed}, {VDD, Qt::darkBlue}, {VSS, Qt::darkGreen}} ;
map<string, Qt::GlobalColor> NetColor = {{"power", Qt::blue}, {"ground", Qt::green}, {"Net", Qt::red}} ;

//-------------------------------------------------------------------------------------
Drawer::Drawer(DesignRule *designRule, QGraphicsScene *scene){
    // 設定緩衝區策略
    const double via_radius_distance = designRule->viaOpeningDiameter / 2;
    const double via_buffer_distance = (designRule->viaPadDiameter - designRule->viaOpeningDiameter )/ 2;
    const double track_buffer_distance = designRule->minimumLineWidth / 2;
    const int points_per_circle = 36;
    
    viaRadiutStrategy = boost::geometry::strategy::buffer::distance_symmetric<double>(via_radius_distance);
    bufferDistanceStrategy = boost::geometry::strategy::buffer::distance_symmetric<double>(via_buffer_distance);
    trackDistanceStrategy = boost::geometry::strategy::buffer::distance_symmetric<double>(track_buffer_distance);
    joinStrategy = boost::geometry::strategy::buffer::join_round(points_per_circle);
    endStrategy = boost::geometry::strategy::buffer::end_round(points_per_circle);
    circleStrategy = boost::geometry::strategy::buffer::point_circle(points_per_circle);
    sideStrategy = boost::geometry::strategy::buffer::side_straight();

    this->scene = scene ; 
    this->designRule = designRule ; 
}

void Drawer::DrawBump(Bump &bump, Qt::GlobalColor bumpColor){
    multi_polygon buffer, final_outer_buffer;
    point pt(bump.x, bump.y) ; 
    
    QBrush outerBrush(Qt::yellow); 
    QPen outerPen(Qt::yellow); outerPen.setWidth(1);
    
    bumpColor = (bumpColor==Qt::yellow) ? BumpColor[bump.type] : bumpColor ;
    QBrush innerdBrush(bumpColor);
    QPen innerPen(bumpColor); innerPen.setWidth(1);

    boost::geometry::buffer(pt, buffer, viaRadiutStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);
    boost::geometry::buffer(buffer, final_outer_buffer, bufferDistanceStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);

    for (const auto &poly : final_outer_buffer) {
        QPolygonF qpoly;
        for (const auto &pt : poly.outer()) qpoly << QPointF(pt.x(), pt.y());
        scene->addPolygon(qpoly, outerPen, outerBrush);
    }

    for (const auto &poly : buffer) { 
        QPolygonF qred;
        for (const auto &pt : poly.outer()) qred << QPointF(pt.x(), pt.y());
        scene->addPolygon(qred, innerPen, innerdBrush);
    }

    // 繪製 bump id 文字
    if (bump.type == SIGNAL){
        QString label = QString::number(bump.id);
        QGraphicsSimpleTextItem* textItem = new QGraphicsSimpleTextItem(label);
        textItem->setBrush(Qt::white);
        textItem->setFont(QFont("Arial", 6));  // 字小一點
        textItem->setPos(bump.x - 5, bump.y - 5);  // bump 位置
        scene->addItem(textItem);  // 直接加進 scene，不加進 group
    }
}

void Drawer::DrawOffsetBump(Bump &bump, Bump &matchedBump){


    point center1(bump.x, bump.y);
    point center2(matchedBump.x, matchedBump.y);  // 兩個圓心相距 50，彼此不相交
    polygon convex_hull;
    multi_polygon buffer1, buffer2, merged_area, final_outer_buffer;
    QBrush outerYellowBrush(Qt::yellow);
    QPen outerYellowPen(Qt::yellow);  outerYellowPen.setWidth(1);

    boost::geometry::buffer(center1, buffer1, viaRadiutStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);
    boost::geometry::buffer(center2, buffer2, viaRadiutStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);
    boost::geometry::union_(buffer1, buffer2, merged_area);
    boost::geometry::convex_hull(merged_area, convex_hull);
    boost::geometry::buffer(convex_hull, final_outer_buffer, bufferDistanceStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);



    for (const auto &poly : final_outer_buffer) {
        QPolygonF qpoly;
        for (const auto &pt : poly.outer()) {
            qpoly << QPointF(pt.x(), pt.y());
        }
        scene->addPolygon(qpoly, outerYellowPen, outerYellowBrush);
    }

    DrawBump(bump) ; 
    DrawBump(matchedBump, OffsetBumpColor[bump.type]) ; 

}

void Drawer::DrawNets(Net& net, Qt::GlobalColor netColor){
    
    if(netColor==Qt::yellow){
        for(auto& [type, color] : NetColor){
            if(net.name.find(type)!=string::npos) netColor = color ;
        }
    }

    QBrush brush(netColor);
    QPen pen(netColor); pen.setWidth(1);

    for(auto& seg : net){
        linestring line ;
        multi_polygon buffer ;
        bg::append(line, point(get<0>(seg) , get<1>(seg))) ;
        bg::append(line, point(get<2>(seg) , get<3>(seg))) ;
        boost::geometry::buffer(line, buffer, trackDistanceStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy) ;
        for (auto &poly : buffer) {
            QPolygonF qpoly;
            for (auto &pt : poly.outer()) {
                qpoly << QPointF(pt.x(), pt.y());
            }
            scene->addPolygon(qpoly, pen, brush);
        }
    }
}

void Drawer::DrawTeardrop(tuple<Bump, double, double, double, double>& teardrop) {
    point bumpPt(get<1>(teardrop), get<2>(teardrop)) ;
    point segPt(get<3>(teardrop), get<4>(teardrop)) ;
    cout << get<0>(teardrop) << " " << get<1>(teardrop) << " " << get<2>(teardrop) << " " << get<3>(teardrop) << " " << get<4>(teardrop) << "\n" ;
    _DrawTeardrop(bumpPt, segPt) ;
}

void Drawer::_DrawTeardrop(const point& center1, const point& center2) {
    double radius = designRule->viaOpeningDiameter / 2 + (designRule->viaPadDiameter - designRule->viaOpeningDiameter )/ 2;
    multi_polygon circle_poly, final_outer_buffer;

    boost::geometry::buffer(center1, circle_poly, bufferDistanceStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);
    boost::geometry::buffer(circle_poly, final_outer_buffer, bufferDistanceStrategy, sideStrategy, joinStrategy, endStrategy, circleStrategy);

    double x1 = bg::get<0>(center1);
    double y1 = bg::get<1>(center1);
    double x2 = bg::get<0>(center2);
    double y2 = bg::get<1>(center2);

    double dx = x2 - x1;
    double dy = y2 - y1;
    double d  = std::sqrt(dx*dx + dy*dy);

    // 圓心->外部點的角度
    double theta = std::atan2(dy, dx);
    // 切線與圓心連線的夾角
    double alpha = std::acos(radius / d);

    // 兩條切線與圓心連線的角度
    double theta1 = theta + alpha;
    double theta2 = theta - alpha;

    // 對應切點的座標
    double tx1_1 = x1 + radius * std::cos(theta1);
    double ty1_1 = y1 + radius * std::sin(theta1);

    double tx1_2 = x1 + radius * std::cos(theta2);
    double ty1_2 = y1 + radius * std::sin(theta2);

    // 建立 Boost.Geometry point
    point T1(tx1_1, ty1_1);
    point T2(tx1_2, ty1_2);

    // ----------------------------------------------------
    // 4. 生成「尾巴」(三角形) 並和「圓形」做 union
    // ----------------------------------------------------
    polygon tail_poly;
    tail_poly.outer().push_back(T1);
    tail_poly.outer().push_back(center2);
    tail_poly.outer().push_back(T2);
    tail_poly.outer().push_back(T1); // 關閉多邊形

    multi_polygon merged_area;
    bg::union_(final_outer_buffer, tail_poly, merged_area);

    // (1) 用 Boost.Geometry 計算「圓之外 (尾巴部分)」
    multi_polygon tail_portion; 
    bg::difference(merged_area, final_outer_buffer, tail_portion);

    QBrush tailBrush(Qt::yellow);      // 其餘部分(尾巴)用藍色
    QPen   tailPen(Qt::yellow);
    tailPen.setWidth(1);

    // (4) 繪製「尾巴區域」(tail_portion)
    for (auto const &poly : tail_portion)
    {
        QPolygonF qpoly;
        for (auto const &pt : poly.outer())
        {
            qpoly << QPointF(bg::get<0>(pt), bg::get<1>(pt));
        }
        scene->addPolygon(qpoly, tailPen, tailBrush);
    }
}