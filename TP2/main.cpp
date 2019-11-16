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
typedef Object<DT4_8, DigitalSetType> ObjectType;

template <class T>
Curve boundary(T &object)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(object.domain().lowerBound() - Point(2, 2), object.domain().upperBound() + Point(2, 2), true);

    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(true);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, object.pointSet(), MAXIMUM_SEARCH);

    // boundary points
    vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, object.pointSet(), bel);

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
    Z2i::DigitalSet set2d(image.domain());
    SetFromImage<Z2i::DigitalSet>::append<ImageType>(set2d, image, 1, 255);

    // vector of digital object
    vector<ObjectType> objects;
    back_insert_iterator<vector<ObjectType>> inserter(objects);

    // connected components
    // (4,8) adjacency
    ObjectType bdiamond(dt4_8, set2d);
    bdiamond.writeComponents(inserter);

    // graph it
    Board2D aBoard;
    for (auto &o : objects)
    {
        aBoard << boundary(o);
    }

    aBoard.saveCairo(("pdf/boundaryCurve" + string(argv[1]) + ".pdf").c_str(), Board2D::CairoPDF);
    return 0;
}
