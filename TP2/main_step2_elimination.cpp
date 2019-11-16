#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/PGMReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/io/Color.h>
#define MAXIMUM_SEARCH 100000

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char>::Type ImageType;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSetType;
typedef Object<DT4_8, DigitalSet> ObjectType48;
typedef Object<DT8_4, DigitalSet> ObjectType84;

template <class T>
Curve boundary(T &object, bool is4_8)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(object.domain().lowerBound() - Point(2, 2), object.domain().upperBound() + Point(2, 2), true);

    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(is4_8);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, object.pointSet(), MAXIMUM_SEARCH);

    // boundary points
    vector<Point> t_BoundaryPoints;
    Surfaces<KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, object.pointSet(), bel);

    // obtain a curve
    Curve boundaryCurve;
    boundaryCurve.initFromVector(t_BoundaryPoints);

    return boundaryCurve;
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
    // (8,4) adjacency
    vector<ObjectType84> objects84;
    back_insert_iterator<vector<ObjectType84>> inserter84(objects84);

    // connected components
    // (4,8) adjacency
    ObjectType48 bdiamond48(dt4_8, set2d);
    bdiamond48.writeComponents(inserter48);
    // (8,4) adjacency
    ObjectType84 bdiamond84(dt8_4, set2d);
    bdiamond84.writeComponents(inserter84);

    Board2D aBoard;

    // Find limits
    int xLimit = image.domain().upperBound()[0] * 2;
    int yLimit = image.domain().upperBound()[1] * 2;

    int count4_8 = 0;
    int count8_4 = 0;

    // draw kept grains
    for (auto &o : objects48)
    {
        Z2i::Curve c = boundary(o, true);
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
            aBoard << c;
            count4_8++;
        }
    }

    for (auto &o : objects84)
    {
        Z2i::Curve c = boundary(o, true);
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
            count8_4++;
        }
    }

    cout << "Nombre de grains de riz 4_8: " << endl;
    cout << count4_8 << endl;
    cout << "Nombre de grains de riz 8_4: " << endl;
    cout << count8_4 << endl;

    aBoard.saveCairo(("pdf/boundaryCurve_" + string(argv[1]) + "_NoBorder.pdf").c_str(), Board2D::CairoPDF);
    return 0;
}