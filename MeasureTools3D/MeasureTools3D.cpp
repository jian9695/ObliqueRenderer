#include "MeasureTools3D.h"
#include "osgUtil/Tessellator"
#include "osgDB/WriteFile"
#include "osg/LineWidth"
#include <osg/Point>
#include <osg/BlendFunc>
#include <osg/Texture2D>
#include <osg/PointSprite>
#include <osg/PolygonMode>
#include <osg/Geometry>
#include "Qt/qfileinfo.h"
#include "ogrsf_frmts.h"
class ShapeFile
{
public:

	ShapeFile()
	{
		poDS = NULL;
	}
	ShapeFile(std::string filename, int update = 0) {
		g_mFileName = filename;
		poDS = NULL;
		poDS = OGRSFDriverRegistrar::Open(filename.data(), update);
		poLayer = poDS->GetLayer(0);
	}

	void close()
	{
		if (poDS)
			OGRDataSource::DestroyDataSource(poDS);
		poDS = NULL;
	}
	void create(const char* filename, OGRSpatialReference* spatialRef = NULL, OGRFeatureDefn *poFDefn = NULL, OGRwkbGeometryType geotype = wkbPoint)
	{


		g_mFileName = filename;
		if (poDS)
			OGRDataSource::DestroyDataSource(poDS);
		const char *pszDriverName = "ESRI Shapefile";
		OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
			pszDriverName);

		if (QFileInfo(filename).exists())
		{
			poDriver->DeleteDataSource(filename);
			//poDS = OGRSFDriverRegistrar::Open(filename, TRUE);
			//poLayer = poDS->GetLayer(0);
		}

		poDS = poDriver->CreateDataSource(filename, NULL);
		poLayer = poDS->CreateLayer(QFileInfo(filename).baseName().toLocal8Bit().data(), spatialRef, geotype, NULL);

		/*OGRSpatialReference oSRS;
		oSRS.SetWellKnownGeogCS("WGS84");
		poLayer = poDS->CreateLayer(QFileInfo(filename).baseName().toLocal8Bit().data(), &oSRS, tp, NULL);
		OGRSpatialReference* spatialRef = poLayer->GetSpatialRef();*/


		//else
		//{
		if (poFDefn)
		{
			for (int iField = 0; iField < poFDefn->GetFieldCount(); iField++)
			{
				OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
				poLayer->CreateField(poFieldDefn);
			}
		}
		//}
	}

	~ShapeFile() {
		close();
	}
public:
	OGRLayer       *poLayer;
	OGRDataSource   *poDS;
private:

	std::string g_mFileName;
};



MeasureTools3D::MeasureTools3D()
{
	OGRRegisterAll();
	osg::ref_ptr<osg::Program> program = new osg::Program;
	
	char vertexShaderSource[] =
	"void main(void)\n"
	"{\n"
	"   gl_Position    = gl_ModelViewProjectionMatrix *  gl_Vertex;\n"
	"}\n";
	char fragmentShaderSource[] =
	"void main(void) \n"
	"{\n"
	"     gl_FragColor = vec4(1,0,0,1);\n"
	"}\n";

	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentShaderSource));
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexShaderSource));
	g_pPointGeode = new osg::Geode;
	osg::StateSet* ss = g_pPointGeode->getOrCreateStateSet();

	ss->setAttribute(program.get(),osg::StateAttribute::ON |osg::StateAttribute::OVERRIDE);

	g_pPointGeode->setCullingActive(false);
	osg::ref_ptr<osg::Point> point = new osg::Point(5.0);
	//point->setDistanceAttenuation(osg::Vec3(0.0, 0.0000, 0.05f));
	ss->setAttribute(point.get());

}



MeasureTools3D::~MeasureTools3D(void)
{

}
void MeasureTools3D::loadFromShapeFile(std::string filename)
{
	OGRRegisterAll();
	g_mNormals.clear();
	g_mPoints.clear();
	ShapeFile shp(filename);
	OGRFeature *poFeature;
	shp.poLayer->ResetReading();
	int pointid = 0;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{

		OGRPoint* ogrp = (OGRPoint*)poFeature->GetGeometryRef();
		double x = poFeature->GetFieldAsDouble("x");
		double y = poFeature->GetFieldAsDouble("y");
		double z = poFeature->GetFieldAsDouble("z");
		double nx = poFeature->GetFieldAsDouble("nx");
		double ny = poFeature->GetFieldAsDouble("ny");
		double nz = poFeature->GetFieldAsDouble("nz");

		osg::Vec3d pos(x, y, z);
		osg::Vec3d normal(nx, ny, nz);
		g_mNormals.push_back(normal);
		g_mPoints.push_back(pos);
		OGRFeature::DestroyFeature(poFeature);

	}
	updateGeometry();

}
bool MeasureTools3D::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{

	if(ea.getEventType() == osgGA::GUIEventAdapter::RELEASE && ea.getButton() == osgGA::GUIEventAdapter::RIGHT_MOUSE_BUTTON)
	{
		osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
		if( !view) 
			return false;
		osgUtil::LineSegmentIntersector::Intersections intersections;

		if(view->computeIntersections(ea,intersections))
		{
			osg::Vec3d pos = intersections.begin()->getWorldIntersectPoint();
			osg::Vec3d normal = intersections.begin()->getWorldIntersectNormal();
			normal.normalize();
			pos = pos + normal * 0.01;

			printf("%f,%f,%f\n", pos.x(), pos.y(), pos.z());
			//pos.z() = g_mBound.zMin();
			g_mNormals.push_back(normal);
			g_mPoints.push_back(pos);
			updateGeometry();

		}
	}
	else if(ea.getKey() == osgGA::GUIEventAdapter::KEY_Delete  &&  ea.getEventType() == osgGA::GUIEventAdapter::KEYUP)
	{
		if(g_mPoints.size()>0)
		{
			g_mPoints.erase(--g_mPoints.end());
			g_mNormals.erase(--g_mNormals.end());
			updateGeometry();
		}
	}
	else if(ea.getKey() == osgGA::GUIEventAdapter::KEY_B  &&  ea.getEventType() == osgGA::GUIEventAdapter::KEYUP)
	{
		ShapeFile shp;
		shp.create(m_filename.data());

		OGRFeatureDefn *poFDefn = shp.poLayer->GetLayerDefn();

		
		OGRFieldDefn defx("x", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defx);
		int x = shp.poLayer->GetLayerDefn()->GetFieldIndex("x");

		OGRFieldDefn defy("y", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defy);
		int y = shp.poLayer->GetLayerDefn()->GetFieldIndex("y");

		OGRFieldDefn defz("z", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defz);
		int z = shp.poLayer->GetLayerDefn()->GetFieldIndex("z");


		OGRFieldDefn defnx("nx", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defnx);
		int nx = shp.poLayer->GetLayerDefn()->GetFieldIndex("nx");

		OGRFieldDefn defny("ny", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defny);
		int ny = shp.poLayer->GetLayerDefn()->GetFieldIndex("ny");

		OGRFieldDefn defnz("nz", OGRFieldType::OFTReal);
		shp.poLayer->CreateField(&defnz);
		int nz = shp.poLayer->GetLayerDefn()->GetFieldIndex("nz");


		for (size_t i = 0; i < g_mPoints.size(); i++)
		{

			osg::Vec3d p = g_mPoints[i];
			osg::Vec3d n = g_mNormals[i];
			OGRPoint po;
			po.setX(p.x());
			po.setY(p.y());
			OGRFeature* poFeaPoint = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
			poFeaPoint->SetGeometry(&po);
			poFeaPoint->SetField(x, p.x());
			poFeaPoint->SetField(y, p.y());
			poFeaPoint->SetField(z, p.z());

			poFeaPoint->SetField(nx, n.x());
			poFeaPoint->SetField(ny, n.y());
			poFeaPoint->SetField(nz, n.z());

			shp.poLayer->CreateFeature(poFeaPoint);
			OGRFeature::DestroyFeature(poFeaPoint);
		}

	}

	return false;
}

void MeasureTools3D::updateGeometry()
{
   g_pPointGeode->removeDrawables(0,g_pPointGeode->getNumDrawables());

   osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
   osg::ref_ptr<osg::Geometry> polyGeom = new osg::Geometry();
   for(std::vector<osg::Vec3d>::iterator iter = g_mPoints.begin(); iter != g_mPoints.end(); ++iter)
   {
	   vertices->push_back(*iter);
   }
   polyGeom->setVertexArray(vertices.get());
   polyGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS,0,vertices->size()));
   g_pPointGeode->addDrawable(polyGeom.get());

}

