#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/PGMReader.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/io/Color.h>
#include <DGtal/geometry/curves/GreedySegmentation.h>
#include <math.h>
#define MAXIMUM_SEARCH 100000

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char>::Type ImageType;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSetType;
typedef Object<DT4_8, DigitalSet> ObjectType48;
typedef Object<DT8_4, DigitalSet> ObjectType84;

typedef FreemanChain<int> Border4;
typedef ArithmeticalDSSComputer<Border4::ConstIterator, int, 4> DSS4;
typedef GreedySegmentation<DSS4> Decomposition4;

template <class T>
Curve boundary(T &o, bool is4_8)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(o.domain().lowerBound() - Point(2, 2), o.domain().upperBound() + Point(2, 2), true);

    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(is4_8);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, o.pointSet(), MAXIMUM_SEARCH);

    // boundary points
    vector<Point> t_BoundaryPoints;
    Surfaces<KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, o.pointSet(), bel);

    // obtain a curve
    Curve boundaryCurve;
    boundaryCurve.initFromVector(t_BoundaryPoints);

    return boundaryCurve;
}

template <class T>
double segmentation(T &o, Domain domain)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(o.domain().lowerBound() - Point(2, 2), o.domain().upperBound() + Point(2, 2), true);

    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(true);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, o.pointSet(), 100000);

    // boundary tracking
    std::vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, o.pointSet(), bel);

    Curve c;
    c.initFromVector(t_BoundaryPoints);

    // Construct the Freeman chain
    Border4 t_Contour(t_BoundaryPoints);
    // Segmentation
    Decomposition4 t_Decomposition(t_Contour.begin(), t_Contour.end(), DSS4());
    // Draw each segment
    string styleName = "";
    double perimeter = 0;
    double partialArea = 0;

    auto itEnd = t_Decomposition.end();
    auto firstPoint = t_Decomposition.begin().begin().get();
    Point lastPoint;
    for (Decomposition4::SegmentComputerIterator it = t_Decomposition.begin(); it != itEnd; ++it)
    {
        auto p = it.get().begin().get();
        auto q = it.get().end().get();
        perimeter += sqrt(pow(q[0] - p[0], 2) + pow(q[1] - p[1], 2));
        partialArea += p[0] * q[1] - p[1] * q[0];
        lastPoint = q;
    }

    perimeter += sqrt(pow(firstPoint[0] - lastPoint[0], 2) + pow(firstPoint[1] - lastPoint[1], 2));
    partialArea += lastPoint[0] * firstPoint[1] - lastPoint[1] * firstPoint[0];
    double area = abs(partialArea) * 0.5;
    double circularity = (4 * M_PI * area) / (perimeter * perimeter);

    cout << perimeter;
    cout << ';';
    // cout << area;
    // cout << ';';
    cout << circularity;
    cout << ';';
}

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        cout << "Please give me the picture name as argument" << endl;
        return 0;
    }
    // read an image
    // const string filestart = "/home/mmmaxou/Imac/imac3/Discrete geometry/TP2/RiceGrains/Rice_";
    const string filestart = "../RiceGrains/Rice_";
    const string filename(argv[1]);
    const string fileend = "_seg_bin.pgm";
    ImageType image = PGMReader<ImageType>::importPGM(filestart + argv[1] + fileend);

    // digital set
    DigitalSet set2d(image.domain());
    SetFromImage<DigitalSet>::append<ImageType>(set2d, image, 1, 255);

    // vector of digital object
    // (4,8) adjacency
    vector<ObjectType48> objects48;
    back_insert_iterator<vector<ObjectType48>> inserter48(objects48);

    // connected components
    // (4,8) adjacency
    ObjectType48 bdiamond48(dt4_8, set2d);
    bdiamond48.writeComponents(inserter48);

    Board2D aBoard;

    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;

    Domain domain = image.domain();

    cout << "Perimètre 1;Circularité 1;Perimètre 2;Circularité 2;" << endl;

    // Find limits
    int xLimit = image.domain().upperBound()[0] * 2;
    int yLimit = image.domain().upperBound()[1] * 2;

    for (auto &o : objects48)
    {
        Curve c = boundary(o, true);
        bool isIn = true;
        for (auto &p : c)
        {
            PointVector<2, Integer> point = p.preCell().coordinates;
            if (point[0] <= 0 || point[0] >= xLimit || point[1] <= 0 || point[1] >= yLimit)
            {
                isIn = false;
                break;
            }
        }
        if (isIn)
        {
            double circularity = (4 * M_PI * o.size()) / (c.size() * c.size());
            cout << c.size();
            cout << ';';
            // cout << o.size();
            // cout << ';';
            cout << circularity;
            cout << ';';
            segmentation(o, domain);
            cout << endl;
        }
    }

    //aBoard.saveCairo("TestGridCurve.pdf", Board2D::CairoPDF);
    return 0;
}
