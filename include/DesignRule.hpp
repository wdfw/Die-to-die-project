#pragma once

struct DesignRule {
    double viaOpeningDiameter; // 藍綠紅 圓的直徑(不含黃色)
    double viaPadDiameter; // 藍綠紅 圓的直徑(含黃色)
    double minimumViaPadSpacing; // 兩顆圓的距離(含黃色)
    double minimumViaSpacing; // via offset(橢圓形結構)中, 兩顆有顏色的圓的距離(不含黃色)
    double minimumLineWidth; // 線寬
    double minimumLineSpacing;  // 線與線之間的距離
    double minimumLineViaSpacing; // 線與 via 之間的距離
    double minimumTeardropDist; // teardrop 外面的點, 到 teardrop 圓心的距離
};