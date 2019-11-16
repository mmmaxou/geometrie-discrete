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
void drawObjectDSSAndCurve(T &o, Domain domain, Board2D &aBoard)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(o.domain().lowerBound() - Point(2, 2), o.domain().upperBound() + Point(2, 2), true);

    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(true);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, o.pointSet(), MAXIMUM_SEARCH);

    // boundary tracking
    std::vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, o.pointSet(), bel);

    Curve c;
    c.initFromVector(t_BoundaryPoints);

    aBoard << c;

    // Construct the Freeman chain
    Border4 boundaryBorder(t_BoundaryPoints);

    // Segmentation
    Decomposition4 boundaryDecomposition(boundaryBorder.begin(), boundaryBorder.end(), DSS4());

    // Draw each segment
    string styleName = "";
    for (Decomposition4::SegmentComputerIterator it = boundaryDecomposition.begin(), itEnd = boundaryDecomposition.end(); it != itEnd; ++it)
    {
        aBoard << SetMode("ArithmeticalDSS", "BoundingBox");
        aBoard << CustomStyle("ArithmeticalDSS/BoundingBox", new CustomPenColor(Color::Green));
        aBoard << it->primitive();
    }
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

    // graph it
    Board2D aBoard;
    Domain domain = image.domain();

    for (auto &o : objects48)
    {
        drawObjectDSSAndCurve(o, domain, aBoard);
    }
    aBoard.saveCairo("pdf/TestGridCurve.pdf", Board2D::CairoPDF);
    return 0;
}
