////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include <DGtal/geometry/tools/MelkmanConvexHull.h>

#include "DGtal/io/boards/Board2D.h"
///////////////////////////////////////////////////////////////////////////////


using namespace DGtal;
using namespace Z2i;

int main()
{
    // randomely make a set of 10 points
    std::vector<Z2i::Point> points;
    srand(time(0));
    for (int i=0; i<10; i++) {
        points.push_back(Point(rand()%100, rand()%100));
    }

    // make a convex hull
    typedef InHalfPlaneBySimple3x3Matrix<Z2i::Point, DGtal::int64_t> Functor;
    Functor functor;
    MelkmanConvexHull<Z2i::Point, Functor>  cvx( functor );
    for (auto &p : points)
        cvx.add(p);

    // scan the CVX points and draw the edges
    Board2D aBoard;
    auto it = cvx.begin();
    Z2i::Point p = *it;
    aBoard.setPenColor(Color::Red);
    for(it = cvx.begin()+1 ; it != cvx.end(); ++it)
    {
        Z2i::Point q=*it;
        aBoard.drawArrow(p[0]-0.5, p[1]-0.5, q[0]-0.5, q[1]-0.5);  //there is a little +1/2 shift in the board exporter
        p = q;
    }
    aBoard.drawArrow(p[0]-0.5, p[1]-0.5, (*cvx.begin())[0]-0.5, (*cvx.begin())[1]-0.5);

    // make a pdf file
    aBoard.saveCairo("convexHull.pdf",Board2D::CairoPDF);
}