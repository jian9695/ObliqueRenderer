#pragma once
#include <Windows.h>
#include "osg/Camera"
#include "osg/Texture2D"
#include "osg/Geode"
#include <osg/MatrixTransform>
#include <osg/TextureCubeMap>
#include <osgGA/TrackballManipulator>
#include <osgGA/FlightManipulator>
#include <osgGA/DriveManipulator>
#include <osgGA/KeySwitchMatrixManipulator>
#include <osgGA/StateSetManipulator>
#include <osgGA/AnimationPathManipulator>
#include <osgGA/TerrainManipulator>
#include <osgGA/SphericalManipulator>
#include <osg/GLExtensions>
#include <osgViewer/Renderer>
#include <osgGA/TrackballManipulator>
#include <osgDB/WriteFile>
#include <osgViewer/ViewerEventHandlers>
#include <osg/ClampColor>
/*#include "DynaLightingRenderer.h"*/
#include <osgDB/ReadFile>
#include "osg/LineWidth"
#include <osgSim/LineOfSight>
#include "SVFDatabasePager.h"
class CameraBuffer : public osg::Camera
{
public:

	osg::ref_ptr<osg::Texture2D> texture;
	osg::ref_ptr<osg::Image> image;
	std::string name;
	osg::Vec3d pos;
	osg::Vec3d dir;
	osg::Vec3d up;
	CameraBuffer()
		:osg::Camera()
	{

	}
	void setupBuffer(int w, int h, osg::Vec3d _dir, osg::Vec3d _up, std::string _name)
	{
		texture = new osg::Texture2D;
		texture->setTextureSize(w, h);
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

		//camera = new osg::Camera;

		setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		setClearColor(osg::Vec4(0, 0, 0, 0));
		//camera->setClearDepth(0);
		//camera->setClearColor(osg::Vec4(0.53f, 0.85f, 1.0f, 0.9f));				// Background
		setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
		setViewport(0, 0, w, h);
		setRenderOrder(osg::Camera::PRE_RENDER);
		setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);

		attach(osg::Camera::COLOR_BUFFER0, texture.get());
		attach(osg::Camera::COLOR_BUFFER0, image.get());

		dir = _dir;
		up = _up;
		name = _name;
	}
	void update();
	static CameraBuffer* create(int w, int h, osg::Vec3d _dir, osg::Vec3d _up, std::string _name)
	{
		CameraBuffer* camera = new CameraBuffer;
		camera->setupBuffer(w, h, _dir, _up, _name);
		return camera;
	}
};

class CameraCB : public osg::NodeCallback
{
public:
	osg::Group* cameraBuffers;
	CameraCB(osg::Group* _cameraBuffers) { cameraBuffers = _cameraBuffers; }
	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{

		if (nv->getVisitorType() == osg::NodeVisitor::CULL_VISITOR)
		{
			osgUtil::CullVisitor* cv = static_cast<osgUtil::CullVisitor*>(nv);
			for (size_t i = 0; i < cameraBuffers->getNumChildren(); i++)
			{
				CameraBuffer* cameraBuffer = (CameraBuffer*)cameraBuffers->getChild(i);
			    cameraBuffer->update();
			}

		}
		traverse(node, nv);
	}
};
class GeomCB : public osg::NodeCallback
{
public:

	GeomCB() {}
	virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
	{

		if (nv->getVisitorType() == osg::NodeVisitor::CULL_VISITOR)
		{
			osgUtil::CullVisitor* cv = static_cast<osgUtil::CullVisitor*>(nv);


		}
		traverse(node, nv);
	}
};


class SceneEventHandler : public osgGA::GUIEventHandler
{
public:
	osg::Group* cameraBuffers;
	osg::ref_ptr<osgGA::CameraManipulator> manip;
	SceneEventHandler(osg::Group* _cameraBuffers, osg::ref_ptr<osgGA::CameraManipulator> _manip)
	{
		cameraBuffers = _cameraBuffers;
		manip = _manip;
	}
	void printfVec3d(osg::Vec3d v)
	{
		printf("%f,%f,%f\n", v.x(), v.y(), v.z());
	}
	void printfVec3(osg::Vec3 v)
	{
		printf("%f,%f,%f\n", v.x(), v.y(), v.z());
	}
	void computeMouseIntersection(osgUtil::LineSegmentIntersector* ray, osgViewer::Viewer* viewer)
	{


		//osg::Vec3d start = eye + osg::Vec3d(0, 0, 1) * 10000;
		//osg::Vec3d end = eye + osg::Vec3d(0, 0, -1) * 10000;
		//   dir.normalize();
		//eye = eye + dir * 0.1;
		//osgSim::LineOfSight los;
		//los.addLOS(start, end);

		//los.computeIntersections(g_pSceneRoot);
		//osg::Vec3d intersection = los.getIntersections(0)[0];
		//double terrainZ = intersection.z();
		//double eyeZ = terrainZ + 1.4 + 1;
		//eye.z() = eyeZ;
		//center = eye + dir * 100;
		//g_pCameraPosition = eye;

		osg::Vec3d orieye, oricenter, oriup;
		//osgViewer::Viewer* viewer = g_pViewer;
		viewer->getCamera()->getViewMatrixAsLookAt(orieye, oricenter, oriup);


		osg::Vec3d dir = orieye - oricenter;
		dir.normalize();

		osg::Vec3d curcenter = ray->getFirstIntersection().getWorldIntersectPoint();

		osg::Vec3d cureye = ray->getFirstIntersection().getWorldIntersectPoint() + dir * 50;

		/*	osgGA::CameraManipulator* manip = viewer->getCameraManipulator();*/
		viewer->setCameraManipulator(NULL, false);
		/*	manip->setHomePosition(cureye, curcenter, oriup);
		manip->home(0);*/

		//viewer->getCamera()->setViewMatrixAsLookAt(cureye, curcenter, oriup);
		viewer->frame();
		VGEDatabasePager* databasePager = (VGEDatabasePager*)viewer->getDatabasePager();
		databasePager->pause();
		databasePager->frame();
		while (databasePager->getFileRequestListSize() > 0)
		{
			viewer->frame();
			databasePager->frame();
			//printf("page,");
		}
		osgUtil::IntersectionVisitor visitor(ray);
		viewer->getCamera()->accept(visitor);
		printfVec3(ray->getFirstIntersection().getWorldIntersectPoint());
		osg::Vec3d observer = ray->getFirstIntersection().getWorldIntersectPoint();
		osg::Vec3d observerNormal = ray->getFirstIntersection().getWorldIntersectNormal();
		observer = observer + observerNormal * 0.5;
		osg::Vec3d observerup = osg::Vec3(0, 0, 1);
		/*	g_pManip->setHomePosition(eye, center, up);
		g_pManip->home(0);*/
		for (size_t i = 0; i < cameraBuffers->getNumChildren(); i++)
		{
			CameraBuffer* cameraBuffer = (CameraBuffer*)cameraBuffers->getChild(i);
			cameraBuffer->pos = observer;
		}
		cameraBuffers->setNodeMask(true);

		/*for (size_t i = 0; i < cameraBuffers.size(); i++)
		{
		CameraBuffer* cameraBuffer = cameraBuffers[i];*/
		//for (size_t j = 0; j < cameraBuffers.size(); j++)
		//{
		//	if (j != i)
		//	{
		//		cameraBuffer->camera->setNodeMask(false);
		//	}
		//	else
		//	{
		//		cameraBuffer->camera->setNodeMask(true);
		//	}
		//}
		viewer->frame();
		databasePager->frame();
		while (databasePager->getFileRequestListSize() > 0)
		{
			//cameraBuffer->update(g_pCameraPosition);
			viewer->frame();
			databasePager->frame();
			//printf("page,");
		}
		viewer->frame();
		//}



		cameraBuffers->setNodeMask(false);
		viewer->getCamera()->setViewMatrixAsLookAt(orieye, oricenter, oriup);
		viewer->setCameraManipulator(manip.get(), false);
		//manip->setHomePosition(orieye, oricenter, oriup);
		//manip->home(0);
		//viewer->frame();
		//while (databasePager->getFileRequestListSize() > 0)
		//{
		//	viewer->frame();
		//	databasePager->frame();
		//}
		databasePager->resume();
	}
	bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
	{
		osgViewer::Viewer* viewer = dynamic_cast<osgViewer::Viewer*>(&aa);
		if (!viewer) return false;



		if (ea.getEventType() == osgGA::GUIEventAdapter::KEYDOWN)
		{

		}

		if (ea.getEventType() == osgGA::GUIEventAdapter::KEYUP)
		{

		}


		//if (ea.getEventType() == osgGA::GUIEventAdapter::DOUBLECLICK)
		//{
		//	osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
		//	osg::ref_ptr<osgUtil::LineSegmentIntersector> ray = new osgUtil::LineSegmentIntersector(osgUtil::Intersector::PROJECTION, ea.getXnormalized(), ea.getYnormalized());
		//	osgUtil::IntersectionVisitor visitor(ray);
		//	view->getCamera()->accept(visitor);
		//	printfVec3(ray->getFirstIntersection().getWorldIntersectPoint());
		//}

		if ((ea.getModKeyMask() & ea.MODKEY_CTRL) != 0 && ea.getEventType() == osgGA::GUIEventAdapter::PUSH)
		{
			osgViewer::Viewer* view = dynamic_cast<osgViewer::Viewer*>(&aa);
			osg::ref_ptr<osgUtil::LineSegmentIntersector> ray = new osgUtil::LineSegmentIntersector(osgUtil::Intersector::PROJECTION, ea.getXnormalized(), ea.getYnormalized());
			osgUtil::IntersectionVisitor visitor(ray);
			view->getCamera()->accept(visitor);
			printfVec3(ray->getFirstIntersection().getWorldIntersectPoint());

			computeMouseIntersection(ray.get(), view);
			//return true;
		}

		return false;
	}
};

class SVFComputeTools
{
public:
	SVFComputeTools();
	~SVFComputeTools();
	osg::Group* createSVFCameras(osg::Node* city)
	{
		osg::Group* cameras = new osg::Group;
		int w = 512;
		int h = 512;
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(1, 0, 0), osg::Vec3(0, 0, 1), "POS_X"));
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(-1, 0, 0), osg::Vec3(0, 0, 1), "NEG_X"));
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(0, 1, 0), osg::Vec3(0, 0, 1), "POS_Y"));
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(0, -1, 0), osg::Vec3(0, 0, 1), "NEG_Y"));
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(0, 0, 1), osg::Vec3(0, -1, 0), "POS_Z"));
		cameras->addChild(CameraBuffer::create(w, h, osg::Vec3(0, 0, -1), osg::Vec3(0, 1, 0), "NEG_Z"));

		for (size_t i = 0; i < cameras->getNumChildren(); i++)
		{
			CameraBuffer* cameraBuffer = (CameraBuffer*)cameras->getChild(i);
			//cameras->addChild(cameraBuffer->camera.get());
			cameraBuffer->addChild(city);
		}

		osg::StateSet* ss = cameras->getOrCreateStateSet();
		CameraCB* cameraCB = new CameraCB(cameras);
		cameras->addCullCallback(cameraCB);

		return cameras;
	}

	osg::Node* cubemap2hemispherical(osg::Group* _cameraBuffers)
	{
		//osg::TextureCubeMap* cubemap = SkyDome::loadCubeMapTextures("E:/OpenSceneGraphSVF/OpenSceneGraphSVF/images_WEIHAI", ".png");
		enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

		osg::TextureCubeMap* cubeMap = new osg::TextureCubeMap;
		cubeMap->setInternalFormat(GL_RGBA);

		cubeMap->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);
		cubeMap->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
		cubeMap->setWrap(osg::Texture::WRAP_S, osg::Texture::CLAMP_TO_EDGE);
		cubeMap->setWrap(osg::Texture::WRAP_T, osg::Texture::CLAMP_TO_EDGE);

		cubeMap->setImage(osg::TextureCubeMap::NEGATIVE_X, ((CameraBuffer*)_cameraBuffers->getChild(NEG_X))->image.get());
		cubeMap->setImage(osg::TextureCubeMap::POSITIVE_X, ((CameraBuffer*)_cameraBuffers->getChild(POS_X))->image.get());
		cubeMap->setImage(osg::TextureCubeMap::NEGATIVE_Y, ((CameraBuffer*)_cameraBuffers->getChild(POS_Z))->image.get());
		cubeMap->setImage(osg::TextureCubeMap::POSITIVE_Y, ((CameraBuffer*)_cameraBuffers->getChild(NEG_Z))->image.get());
		cubeMap->setImage(osg::TextureCubeMap::NEGATIVE_Z, ((CameraBuffer*)_cameraBuffers->getChild(NEG_Y))->image.get());
		cubeMap->setImage(osg::TextureCubeMap::POSITIVE_Z, ((CameraBuffer*)_cameraBuffers->getChild(POS_Y))->image.get());
		osg::Geode* geode = new osg::Geode;
		//createGeometry( );
		//SetupStateSet(cubemap);
		osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
		osg::Vec3Array* vertices = new osg::Vec3Array();
		osg::Vec3Array* normals = new osg::Vec3Array();
		vertices->push_back(osg::Vec3(-1, -1, -1));
		vertices->push_back(osg::Vec3(-1, 1, -1));
		vertices->push_back(osg::Vec3(1, -1, -1));
		vertices->push_back(osg::Vec3(1, -1, -1));
		vertices->push_back(osg::Vec3(-1, 1, -1));
		vertices->push_back(osg::Vec3(1, 1, -1));
		normals->push_back(osg::Vec3(0, 0, -1));
		geom->setVertexArray(vertices);
		geom->setNormalArray(normals);
		geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
		geode->getOrCreateStateSet()->setMode(GL_CULL_FACE,
			osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
		//geom->setInitialBound(osg::BoundingBox(-1000000, -1000000, -1000000, 1000000, 1000000, 1000000));
		//geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
		geom->setCullingActive(false);
		geode->setCullingActive(false);
		geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, 6));
		geode->addDrawable(geom.get());
		geode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		geode->getOrCreateStateSet()->setTextureAttributeAndModes(0, cubeMap, osg::StateAttribute::ON);
		//geode->getOrCreateStateSet()->setTextureAttributeAndModes(1, g_pPanoTexture, osg::StateAttribute::ON);
		osg::Program* program = new osg::Program;

		char vertexSource[] =
			"void main(void)\n"
			"{\n"
			"   gl_TexCoord[0] = vec4(gl_Vertex.x,gl_Vertex.y,0,1.0);\n"
			"   //gl_Position   = vec4(gl_Vertex.x,gl_Vertex.y,0,1);\n"
			"   gl_Position   = vec4(gl_Vertex.x*0.5+0.5,gl_Vertex.y*0.5+0.5,0,1);\n"
			"}\n";
		char fragmentSource[] =
			"uniform vec4 color;\n"
			"uniform float alpha;\n"
			"uniform samplerCube uEnvironmentMap;\n"
			"uniform float rotateAngle;\n"
			"\n"
			"vec3 spherical2Cartisian(float lon, float lat)\n"
			"{\n"
			"float theta = lon * 0.0174533;\n"
			"float phi =   lat* 0.0174533;\n"
			"return vec3(cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi));\n"
			"}\n"

			"vec2 rotate(vec2 uv,float angle)\n"
			"{\n"
			"    angle = angle * 0.0174533;\n"
			"    float sin_factor = sin(angle);\n"
			"    float cos_factor = cos(angle);\n"
			"    uv = (uv - 0.5) * mat2(cos_factor, sin_factor, -sin_factor, cos_factor);\n"
			"    uv += 0.5;\n"
			"    return uv;\n"
			"}\n"
			"void main(void) \n"
			"{\n"
			"    float radius = length(gl_TexCoord[0].xy);\n"
			"    vec3 tex = vec3(-gl_TexCoord[0].x, gl_TexCoord[0].y, -(1-radius));\n"
			"    vec4 fisheye1 = textureCube( uEnvironmentMap, tex.xzy );\n"
			"    vec2 uv = gl_TexCoord[0].xy; \n"
			"    uv = uv * 0.5 + 0.5;\n"
			"    uv = rotate(uv,rotateAngle);\n"
			"    gl_FragColor = vec4(fisheye1.rgb,1);\n"
			"    if(radius > 1)\n"
			"    {\n"
			"        gl_FragColor = vec4(0,0,0,0);\n"
			"    }\n"
			"    else if(fisheye1.a < 0.5 )\n"
			"    {\n"
			"        gl_FragColor = vec4(1,1,1,0.5);\n"
			"    }\n"
			"}\n";

		geode->setCullCallback(new GeomCB);
		program->setName("sky_dome_shader");
		program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexSource));
		program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentSource));
		geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
		geode->getOrCreateStateSet()->addUniform(new osg::Uniform("uEnvironmentMap", 0));
		geode->getOrCreateStateSet()->addUniform(new osg::Uniform("rotateAngle", 0.0f));
		geode->getOrCreateStateSet()->addUniform(new osg::Uniform("alpha", 0.0f));

		//geode->getOrCreateStateSet()->setRenderBinDetails(5000, "RenderBin");
		geode->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF);
		//geode->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
		geode->getOrCreateStateSet()->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);
		//g_pPanoState = geode->getOrCreateStateSet();

		osg::Projection * HUDProjectionMatrix = new osg::Projection;

		HUDProjectionMatrix->setMatrix(osg::Matrix::ortho2D(0, 512, 0, 512));
		osg::MatrixTransform* HUDModelViewMatrix = new osg::MatrixTransform;
		HUDModelViewMatrix->setMatrix(osg::Matrix::translate(512, 512, 0));
		// above it in the scene graph:
		HUDModelViewMatrix->setReferenceFrame(osg::Transform::ABSOLUTE_RF);

		HUDProjectionMatrix->addChild(HUDModelViewMatrix);
		HUDModelViewMatrix->addChild(geode);

		return HUDProjectionMatrix;
		//return geode;
	}
	osg::Node* createTextureRect(std::string texfile)
	{


		osg::ref_ptr<osg::Texture2D> tex = new osg::Texture2D;
		tex->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::CLAMP_TO_BORDER);
		tex->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::CLAMP_TO_BORDER);
		tex->setImage(osgDB::readImageFile(texfile));

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
			"   //gl_Position   = vec4(gl_Vertex.x*0.5-0.5,gl_Vertex.y*0.5+0.5,0,1);\n"
			"   gl_Position   = vec4(gl_Vertex.x*0.5+0.5,gl_Vertex.y*0.5+0.5,0,1);\n"
			"}\n";
		char fragmentSource[] =
			"uniform sampler2D texture0;\n"
			"uniform float rotateAngle;\n"
			"vec2 rotate(vec2 uv,float angle)\n"
			"{\n"
			"    angle = angle * 0.0174533;\n"
			"    float sin_factor = sin(angle);\n"
			"    float cos_factor = cos(angle);\n"
			"    uv = (uv - 0.5) * mat2(cos_factor, sin_factor, -sin_factor, cos_factor);\n"
			"    uv += 0.5;\n"
			"    return uv;\n"
			"}\n"
			"void main(void) \n"
			"{\n"
			"    vec2 uv = gl_TexCoord[0].xy; \n"
			"    uv = uv * 0.5 + 0.5;\n"
			"    uv = rotate(uv,rotateAngle);\n"
			"    vec4 color =   texture2D( texture0, uv );\n"
			"    gl_FragColor = vec4(color.rgb,0.5);\n"
			"}\n";
		program->setName("sky_dome_shader");
		program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexSource));
		program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentSource));
		geode->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
		geode->getOrCreateStateSet()->addUniform(new osg::Uniform("texture0", 0));
		geode->getOrCreateStateSet()->addUniform(new osg::Uniform("rotateAngle", 0.0f));

		//geode->getOrCreateStateSet()->setRenderBinDetails(100000, "RenderBin");
		geode->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF);
		geode->getOrCreateStateSet()->setRenderingHint(osg::StateSet::OPAQUE_BIN);
		//g_pPanoState = geode->getOrCreateStateSet();
		return geode;
	}
};

