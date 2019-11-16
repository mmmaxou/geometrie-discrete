///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//shape and digitizer
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"

//tracking grid curve
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/GridCurve.h"

#include "DGtal/io/boards/Board2D.h"
///////////////////////////////////////////////////////////////////////////////

using namespace DGtal;
using namespace Z2i;
using namespace std;

typedef ImplicitBall<Z2i::Space> Disk;
typedef AccFlower2D<Z2i::Space> Flower;
typedef Ellipse2D<Z2i::Space> Ellipse;

const Z2i::Point CENTER(0, 0);

const Flower makeFlower()
{
  return Flower(Z2i::Point(0, 0), 15, 5, 2, 0.2);
}

const Disk makeDisk()
{
  return Disk(Z2i::Point(0, 0), 500);
}

const Ellipse makeEllipse()
{
  int a = 200;
  return Ellipse(Z2i::Point(0, 0), a * 1.5, a / 2, 0.0);
}

void findValues(const std::vector<Z2i::Point> &boundaryPoints, const DigitalSet &shape)
{
  cout << "Perimeter is " << boundaryPoints.size() << endl;
  cout << "Area is " << shape.size() << endl;
}

int main()
{
  // define an Euclidean shape (disk)
  auto shape = makeEllipse();

  // Gauss discretization
  double h = 1;
  GaussDigitizer<Z2i::Space, Ellipse> dig;
  dig.attach(shape);
  dig.init(shape.getLowerBound() + Z2i::Vector(-1, -1),
           shape.getUpperBound() + Z2i::Vector(1, 1), h);

  Domain domain(dig.getLowerBound(), dig.getUpperBound());

  // make a Kovalevsky-Khalimsky space
  Z2i::KSpace ks;
  ks.init(dig.getLowerBound(), dig.getUpperBound(), true);
  // set an adjacency (4-connectivity)
  SurfelAdjacency<2> sAdj(true);

  // search for one boundary element
  Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel(ks, dig, 1000);
  // boundary tracking
  std::vector<Z2i::Point> boundaryPoints;
  Surfaces<Z2i::KSpace>::track2DBoundaryPoints(boundaryPoints, ks, sAdj, dig, bel);

  Z2i::Curve c;
  c.initFromVector(boundaryPoints);

  // display the perimeter and the area of the shape
  Z2i::DigitalSet aSet(domain);
  Shapes<Z2i::Domain>::digitalShaper(aSet, dig);
  findValues(boundaryPoints, aSet);

  // make a convex hull
  typedef InHalfPlaneBySimple3x3Matrix<Z2i::Point, DGtal::int64_t> Functor;
  Functor functor;
  MelkmanConvexHull<Z2i::Point, Functor> cvx(functor);
  for (auto &p : boundaryPoints)
    cvx.add(p);

  // draw convex hull
  // Compute the convex hull area and perimeter
  auto it = cvx.begin();
  Z2i::Point p = *it;
  double convexHullArea = 0.;
  double convexHullPerimeter = 0.;
  Board2D aBoard;
  aBoard << c;
  aBoard.setPenColor(Color::Red);
  for (it = cvx.begin() + 1; it != cvx.end(); ++it)
  {
    Z2i::Point q = *it;

    convexHullArea += abs(p[0] * (q[1] - CENTER[1]) + q[0] * (CENTER[1] - p[1]) + CENTER[0] * (p[1] - q[1])) / 2.; //abs(Ax(By - Cy) + Bx(Cy - Ay) + Cx(Ay - By)) / 2
    convexHullPerimeter += sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]));

    aBoard.drawArrow(p[0] - 0.5, p[1] - 0.5, q[0] - 0.5, q[1] - 0.5); //there is a little +1/2 shift in the board exporter
    p = q;
  }
  cout << "Convex hull area: " << convexHullArea << endl;
  cout << "Convex hull perimeter: " << convexHullPerimeter << endl;

  aBoard.drawArrow(p[0] - 0.5, p[1] - 0.5, (*cvx.begin())[0] - 0.5, (*cvx.begin())[1] - 0.5);

  aBoard.saveCairo("boundaryCurve.pdf", Board2D::CairoPDF);
}

///////////////////////////////////////////////////////////////////////////////
