

#include <Windows.h>
#include "qdir.h"
#include "SVFDatabasePager.h"
#include "SVFComputeTools.h"
#include "osg/ComputeBoundsVisitor"
#include "gdal_priv.h"
#include <time.h>
#include <chrono>
void setupCamera(osg::Camera* camera, osg::Vec3d& center,osg::Vec3d& bbsize)
{
	osg::Vec3d eye = center;
	double viewDist = 1000;
	eye.z() = center.z() + 1000;
	double znear = viewDist - 500;
	double zfar = viewDist + 500;
	camera->setViewMatrixAsLookAt(eye, center, osg::Vec3d(0,1,0));
	camera->setProjectionMatrixAsOrtho(-bbsize.x()*0.5, bbsize.x()*0.5, -bbsize.y()*0.5, bbsize.y()*0.5, znear, zfar);
}
osg::Camera* createOrtho(int texSizeX, int texSizeY)
{

	osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;

	tex->setTextureSize(texSizeX, texSizeY);
	tex->setResizeNonPowerOfTwoHint(false);
	tex->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
	tex->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
	tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
	tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
	tex->setInternalFormat(GL_RGBA);
	tex->setSourceFormat(GL_RGBA);
	tex->setSourceType(GL_UNSIGNED_BYTE);

	osg::ref_ptr<osg::Image> img = new osg::Image;
	//img->allocateImage(texSize, texSize,1,GL_ALPHA,GL_FLOAT);
	img->allocateImage(texSizeX, texSizeY, 1, GL_RGBA, GL_UNSIGNED_BYTE);
	//img->allocateImage(texSize, texSize,1,GL_RGBA,GL_FLOAT); 
	tex->setImage(img.get());

	osg::Camera* orthoCamera = new osg::Camera;
	orthoCamera->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	orthoCamera->setClearColor(osg::Vec4(0, 0, 0, 0));

	orthoCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
	orthoCamera->setViewport(0, 0, texSizeX, texSizeY);
	orthoCamera->setRenderOrder(osg::Camera::PRE_RENDER);
	orthoCamera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
	orthoCamera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, tex.get());
	orthoCamera->attach(osg::Camera::COLOR_BUFFER0, img.get());
	return orthoCamera;

}
struct PageLODControl
{
	osg::Node* NodePtr;
	osg::BoundingBox BB;
	PageLODControl(osg::Node* node)
	{
		osg::ComputeBoundsVisitor cbv;
		node->accept(cbv);
		BB = cbv.getBoundingBox();
		NodePtr = node;
	}
	bool intersect(double& xmin, double& ymin, double& xmax, double& ymax)
	{
		if (xmin > BB.xMax() || xmax < BB.xMin() || ymin > BB.yMax() || ymax < BB.yMin())
		{
			NodePtr->setNodeMask(false);
			return false;
		}
		//printf("%s\n", NodePtr->getName().data());
		NodePtr->setNodeMask(true);
		return true;
	}
};

class PageLODManager
{
public:
	std::vector<PageLODControl> PageLODControls;
	PageLODManager(osg::Group* pagelodGroup)
	{
		for (size_t i = 0; i < pagelodGroup->getNumChildren(); i++)
		{
			PageLODControl control(pagelodGroup->getChild(i));
			PageLODControls.push_back(control);
		}
	}
	void intersect(double& xmin, double& ymin, double& xmax, double& ymax)
	{
		for (size_t i = 0; i < PageLODControls.size(); i++)
		{
			PageLODControls[i].intersect(xmin, ymin, xmax, ymax);
		}

	}

};
//--------------------------------------------------------------------------
osgViewer::Viewer* createViewer(int topleftx,int toplefty,int w,int h)
{
	osgViewer::Viewer* viewer = new osgViewer::Viewer;
	viewer->setUpViewInWindow(topleftx, toplefty, w, h);
	//viewer->setUpViewInWindow(0, 0, 1, 1);
	//osg::ref_ptr<osgGA::CameraManipulator> manip = new osgGA::TrackballManipulator;
	//viewer->setCameraManipulator(manip.get());
	viewer->setThreadingModel(osgViewer::ViewerBase::ThreadingModel::ThreadPerContext);
	//ObliqueCameraWrapper* orthoCamera = new ObliqueCameraWrapper(viewer.getCamera(),bb, 45);
	viewer->getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	viewer->getCamera()->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//viewer.getCamera
	//viewer->addEventHandler(new KeyboardEventHandler());
	// add the state manipulator
	//viewer->addEventHandler( new osgGA::StateSetManipulator(viewer->getCamera()->getOrCreateStateSet()) );

	//databasePager->pause();
	// add the thread model handler
	viewer->addEventHandler(new osgViewer::ThreadingHandler);

	// add the window size toggle handler
	viewer->addEventHandler(new osgViewer::WindowSizeHandler);

	// add the stats handler
	viewer->addEventHandler(new osgViewer::StatsHandler);

	// add the LOD Scale handler
	viewer->addEventHandler(new osgViewer::LODScaleHandler);
	return viewer;
}
osg::Node* createTextureRect(osg::ref_ptr<osg::Texture2D>& tex)
{
	tex = new osg::Texture2D;
	tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::CLAMP_TO_BORDER);
	tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::CLAMP_TO_BORDER);
	//tex->setImage(img);
	osg::Geode* geode = new osg::Geode;
	//createGeometry( );
	//SetupStateSet(cubemap);
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	osg::Vec3Array* vertices = new osg::Vec3Array();
	osg::Vec3Array* normals = new osg::Vec3Array();
	vertices->push_back(osg::Vec3(-1, -1, 1));
	vertices->push_back(osg::Vec3(-1, 1, 1));
	vertices->push_back(osg::Vec3(1, -1, 1));
	vertices->push_back(osg::Vec3(1, -1, 1));
	vertices->push_back(osg::Vec3(-1, 1, 1));
	vertices->push_back(osg::Vec3(1, 1, 1));
	normals->push_back(osg::Vec3(0, 0, 1));
	geom->setVertexArray(vertices);
	geom->setNormalArray(normals);
	geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
	geom->setCullingActive(false);
	geode->setCullingActive(false);
	geode->getOrCreateStateSet()->setMode(GL_CULL_FACE,
		osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, 6));
	geode->addDrawable(geom.get());
	geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	geode->getOrCreateStateSet()->setTextureAttributeAndModes(0, tex.get(), osg::StateAttribute::ON);

	osg::Program* program = new osg::Program;

	char vertexSource[] =
		"void main(void)\n"
		"{\n"
		"   gl_TexCoord[0] = vec4(gl_Vertex.x,gl_Vertex.y,0,1.0);\n"
		"   gl_Position   = vec4(gl_Vertex.x,gl_Vertex.y,0,1);\n"
		"}\n";
	char fragmentSource[] =
		"uniform sampler2D texture0;\n"
		"void main(void) \n"
		"{\n"
		"    vec2 uv = gl_TexCoord[0].xy; \n"
		"    uv = uv * 0.5 + 0.5;\n"
		"    vec4 color =   texture2D( texture0, uv );\n"
		"    gl_FragColor = vec4(color.rgb,0.5);\n"
		"}\n";
	program->setName("sky_dome_shader");
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexSource));
	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentSource));
	geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
	geode->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF);
	geode->getOrCreateStateSet()->setRenderingHint(osg::StateSet::OPAQUE_BIN);
	return geode;
}
class OrthoCameraCB : public osg::NodeCallback
{
public:
	osg::Vec3d _bbcenter;
	osg::Vec3d _bbsize;
	osg::Camera* _camera;
	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{

		if (nv->getVisitorType() == osg::NodeVisitor::CULL_VISITOR)
		{
			osgUtil::CullVisitor* cv = static_cast<osgUtil::CullVisitor*>(nv);
		
			setupCamera(_camera, _bbcenter, _bbsize);
		}
		traverse(node, nv);
	}
};

void testPerformance3D(osgViewer::Viewer* viewer,double cellsize, int w,int h)
{

	osg::Group* root = new osg::Group;

	//"F:/Oblique_Photogrammetry/xinji/Data/all.desc")
	osg::Node* city = osgDB::readNodeFile("Z:/Jianming/xinji/Data/all.desc");
	city->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	osg::ComputeBoundsVisitor cbv;
	city->accept(cbv);
	osg::BoundingBox bb = cbv.getBoundingBox();
	osg::Vec3d bbsize(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	osg::Vec3d bbcenter = bb.center();
	printf("%f,%f,%f", bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	PageLODManager manager((osg::Group*)city);
	//double cellsize = 50;
	double halfcellsize = cellsize * 0.5;
	//double xmin = bb.center().x() - halfcellsize;
	//double xmax = bb.center().x() + halfcellsize;
	//double ymin = bb.center().y() - halfcellsize;
	//double ymax = bb.center().y() + halfcellsize;
	//manager.intersect(xmin, ymin, xmax, ymax);
	root->addChild(city);
	viewer->setSceneData(root);
	osg::Vec3d bbcellsize(cellsize, cellsize, cellsize);
	//OrthoCameraCB* cb = new OrthoCameraCB;
	//cb->_camera = viewer->getCamera();
	//cb->_bbcenter = bbcenter;
	//cb->_bbsize = bbsize;
	osg::Camera* camera = viewer->getCamera();
	osg::ref_ptr<osg::Texture2D> texture;
	osg::ref_ptr<osg::Image> image;
	texture = new osg::Texture2D;
	texture->setTextureSize(1024, 768);
	texture->setResizeNonPowerOfTwoHint(false);
	texture->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR);
	texture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
	texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
	texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
	//rtTexture->setDataVariance(osg::Object::DYNAMIC);
	texture->setInternalFormat(GL_RGBA);
	texture->setSourceFormat(GL_RGBA);
	texture->setSourceType(GL_UNSIGNED_BYTE);
	image = new osg::Image;
	image->allocateImage(w, h, 1, GL_RGBA, GL_UNSIGNED_BYTE);
	texture->setImage(image.get());

	/* camera->attach(osg::Camera::COLOR_BUFFER0, texture.get());*/
	camera->attach(osg::Camera::COLOR_BUFFER0, image.get());
	VGEDatabasePager* databasePager = new VGEDatabasePager;

	viewer->getScene()->setDatabasePager(databasePager);
	databasePager->cancel();
	double celly = bb.yMax();

	int numrows = 0;
	int numcols = 0;
	while (celly - cellsize >= bb.yMin())
	{
		celly = celly - cellsize;
		numrows++;
	}
	double cellx = bb.xMin();

	while (cellx + cellsize <= bb.xMax())
	{
		cellx = cellx + cellsize;
		numcols++;
	}
	clock_t begin = clock();
	LARGE_INTEGER frequency = { 1 }, counter = { 0 };
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	int64_t mFrequency = static_cast<int64_t>(frequency.QuadPart);
	int64_t mInitialTicks = static_cast<int64_t>(counter.QuadPart);
	int64_t mInvFrequency = 1.0 / static_cast<double>(mFrequency);
	//3D performance 
	//2253.624
	typedef std::chrono::high_resolution_clock Clock;
	//std::ofstream stats("stats25d.csv");
	//stats << 
	double visTime = 0;
	double loadTime = 0;
	int ncell = 0;
	for (int nrow = 0; nrow < numrows; nrow++)
	{
		celly = bb.yMax() - nrow * cellsize - halfcellsize;
		for (int ncol = 0; ncol < numcols; ncol++)
		{
			if (ncell % 50 == 0)
			{
				cellx = bb.xMin() + ncol * cellsize + halfcellsize;
				osg::Vec3d cellcenter(cellx + halfcellsize, celly - halfcellsize, 0);
				setupCamera(viewer->getCamera(), cellcenter, bbcellsize);
				viewer->frame();
				while (databasePager->getFileRequestListSize() > 0)
				{
					setupCamera(viewer->getCamera(), cellcenter, bbcellsize);
					//auto begin1 = std::chrono::high_resolution_clock::now();
					clock_t begin2 = clock();
					databasePager->frame();
					//auto end1 = std::chrono::high_resolution_clock::now();
					clock_t end2 = clock();
					//int64_t msecs = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count();
					//double secs = static_cast<double>(msecs) / 1000.0;
					//printf("loading time: %f,%f\n", secs, double(end2 - begin2) / CLOCKS_PER_SEC);
					loadTime += double(end2 - begin2) / CLOCKS_PER_SEC;
					begin2 = clock();
					viewer->frame();
					end2 = clock();
					visTime += double(end2 - begin2) / CLOCKS_PER_SEC;

				}
			}
		
			//std::stringstream ss;
			//ss << "B:/VGE/DSMRenderer_2016_02_14/PerformanceTest/output/" << nrow << "_" << ncol << ".png";
			//osgDB::writeImageFile(*image, ss.str());
			printf("row=%d,col=%d\n", nrow, ncol);
			ncell++;
		}
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("3D performance=%f", elapsed_secs);
	printf("loadtime=%f,vistime=%f", loadTime,visTime);
	getchar();
	//stats.close();
	//%73 loading
	//%27 visualization
}
void testPerformance25D(osgViewer::Viewer* viewer, double cellsize, int w, int h)
{
	//osg::Group* root = new osg::Group;
	//"F:/Oblique_Photogrammetry/xinji/Data/all.desc")
	osg::Node* city = osgDB::readNodeFile("Z:/Jianming/xinji/Data/all.desc");
	city->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	osg::ComputeBoundsVisitor cbv;
	city->accept(cbv);
	osg::BoundingBox bb = cbv.getBoundingBox();
	osg::Vec3d bbsize(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	osg::Vec3d bbcenter = bb.center();
	printf("%f,%f,%f", bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	double halfcellsize = cellsize * 0.5;
	osg::Vec3d bbcellsize(cellsize, cellsize, cellsize);
	int pixelspercol = (int)(cellsize / 0.25);
	osg::ref_ptr<osg::Texture2D> tex;
	osg::Node* screenQuad = createTextureRect(tex);
	viewer->setSceneData(screenQuad);
	double celly = bb.yMax();
	int numrows = 0;
	int numcols = 0;
	while (celly - cellsize >= bb.yMin())
	{
		celly = celly - cellsize;
		numrows++;
	}
	double cellx = bb.xMin();

	while (cellx + cellsize <= bb.xMax())
	{
		cellx = cellx + cellsize;
		numcols++;
	}
	clock_t begin = clock();

	//25D performance 
	std::string filename25d = "Z:/Jianming/OAP3D_XINJI/2D/90_0.tif";
	GDALDataset* pDataset25d = (GDALDataset*)GDALOpen(filename25d.data(), GA_ReadOnly);
	GDALRasterBand* pBand25d = pDataset25d->GetRasterBand(1);
	double resol = 0.25;

	int panBandMap[3] = { 1,2,3 };



	for (int nrow = 0; nrow < numrows; nrow++)
	{
		int pstarty = pixelspercol * nrow;
		for (int ncol = 0; ncol < numcols; ncol++)
		{
			int pstartx = pixelspercol * ncol;

			osg::ref_ptr<osg::Image> maptile = new osg::Image;
			maptile->allocateImage(pixelspercol, pixelspercol, 1, GL_RGB, GL_BYTE);
			pDataset25d->RasterIO(GF_Read, pstartx, pstarty, pixelspercol, pixelspercol, maptile->data(),
				pixelspercol, pixelspercol, GDT_Byte, 3, panBandMap, 3, pixelspercol * 3, 1);
			//tex->setImage(maptile.get());
			//viewer->frame();
			//viewer->frame();
			//viewer->frame();
			//viewer->frame();

			//pBand25d->RasterIO(GF_Read,pstartx,pstarty,pixelspercol,pixelspercol,)
			//pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Float64, 0, 0);
			//pBand->SetNoDataValue(-9999);
			//std::stringstream ss;
			//ss << "B:/VGE/DSMRenderer_2016_02_14/PerformanceTest/output/" << nrow << "_" << ncol << ".png";
			//osgDB::writeImageFile(*maptile, ss.str());

			printf("row=%d,col=%d\n", nrow, ncol);
		}
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("25D performance=%f", elapsed_secs);

}

void copyDir()
{

	osg::Group* root = new osg::Group;

	//"F:/Oblique_Photogrammetry/xinji/Data/all.desc")
	osg::Node* city = osgDB::readNodeFile("F:/Oblique_Photogrammetry/gongan/1/all.desc");
	city->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	osg::ComputeBoundsVisitor cbv;
	city->accept(cbv);
	osg::BoundingBox bb = cbv.getBoundingBox();
	osg::Vec3d bbsize(bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	osg::Vec3d bbcenter = bb.center();
	printf("%f,%f,%f", bb.xMax() - bb.xMin(), bb.yMax() - bb.yMin(), bb.zMax() - bb.zMin());
	PageLODManager manager((osg::Group*)city);
	double halfsize = 100;
	double xmin = bbcenter.x() - halfsize;
	double xmax = bbcenter.x() + halfsize;
	double ymin = bbcenter.y() - halfsize;
	double ymax = bbcenter.y() + halfsize;
	std::string outroot = "F:/Oblique_Photogrammetry/OAP3D/";
	std::ofstream ofs(outroot + "index.desc");
	for (size_t i = 0; i < manager.PageLODControls.size(); i++)
	{
		PageLODControl& control = manager.PageLODControls[i];
		if (!control.intersect(xmin, ymin, xmax, ymax))
			continue;
		QDir qdir = QFileInfo(control.NodePtr->getName().data()).absoluteDir();
		std::string dirname = QFileInfo(control.NodePtr->getName().data()).baseName().toLocal8Bit().data();
		ofs << dirname +"\\" + QFileInfo(control.NodePtr->getName().data()).fileName().toLocal8Bit().data() << std::endl;
		std::string outdir = outroot + dirname + "/";
		std::string indir = qdir.absolutePath().toLocal8Bit().data();
		std::stringstream ss;
		std::string command = "xcopy " + indir + " " + outdir;
		QDir(outdir.data()).mkpath(".");
		for (size_t nchar = 0; nchar < command.size(); nchar++)
		{
			if (command[nchar] == '/')
				command[nchar] = '\\';
		}
		system(command.data());
	}
	ofs.close();
}
int main(int argc, char **argv)
{
	copyDir();
	GDALAllRegister();
    // parse arguments
    osg::ArgumentParser arguments(&argc,argv);
	int winwidth = 1024;
	int winheight = 768;
	double cellsize = 50;
	osgViewer::Viewer* viewer = createViewer(100,100, winwidth, winheight);
	//testPerformance3D(viewer,cellsize,winwidth,winheight);
	testPerformance25D(viewer, cellsize, winwidth, winheight);
	 getchar();
//2.5D performance 
	
	 while (true)
	 {
		 viewer->frame();
	 }
	 //viewer->getCamera()->addUpdateCallback(cb);
	 //root->addCullCallback(cb);
	 //CameraBuffer* orthoCamera =  CameraBuffer::create(1024, 1024,osg::Vec3d(0,0,1),osg::Vec3d(0,1,0),"");
	 //orthoCamera->addChild(city);
	 //setupCamera(orthoCamera, bbcenter, bbsize);
	 //root->addChild(orthoCamera);


	
	// add the screen capture handler
	//viewer->addEventHandler(new osgViewer::ScreenCaptureHandler);
	 //viewer->addEventHandler(new SceneEventHandler(SVFCameras, manip));


	

	//osgDB::writeImageFile(*orthoCamera->image, "image.png");

    return viewer->run();
}



