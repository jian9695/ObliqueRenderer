#include "ScreenOverlay.h"
#include <osgDB/ReadFile>
#include <osgGA/GUIEventHandler>
class BlendFactorEventHandler :
	public osgGA::GUIEventHandler
{
private:
	osg::Uniform* _blendFactorUniform;
	bool _KeyB_IsDown;
	float _blendFactor;
public:

	BlendFactorEventHandler(osg::Uniform* blendFactorUniform)
	{
		_blendFactorUniform = blendFactorUniform;
		_blendFactor = 1;
		_KeyB_IsDown = false;
	}

	bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
	{
		if (ea.getKey() == osgGA::GUIEventAdapter::KEY_B)
		{
			if (ea.getEventType() == osgGA::GUIEventAdapter::KEYDOWN)
			{
				_KeyB_IsDown = true;
			}
			else if (ea.getEventType() == osgGA::GUIEventAdapter::KEYUP)
			{
				_KeyB_IsDown = false;
			}

		}

		if (_KeyB_IsDown && ea.getEventType() == osgGA::GUIEventAdapter::KEYDOWN)
		{
			if (ea.getKey() == osgGA::GUIEventAdapter::KEY_Minus)
			{
				_blendFactor -= 0.02;
				if (_blendFactor < 0)
					_blendFactor = 0;
				_blendFactorUniform->set(_blendFactor);
			}
			else if (ea.getKey() == 61 || ea.getKey() == 64451)
			{
				_blendFactor += 0.02;
				if (_blendFactor > 1)
					_blendFactor = 1;
				_blendFactorUniform->set(_blendFactor);
			}

		}
		return false;
	}

};

ScreenOverlay::ScreenOverlay(osgViewer::Viewer* viewer)
{
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
	this->setCullingActive(false);
	this->getOrCreateStateSet()->setMode(GL_CULL_FACE,
		osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, 6));
	this->addDrawable(geom.get());
	this->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	//this->getOrCreateStateSet()->setTextureAttributeAndModes(0, tex.get(), osg::StateAttribute::ON);

	osg::Program* program = new osg::Program;

	char vertexSource[] =
		"void main(void)\n"
		"{\n"
		"   gl_TexCoord[0] = vec4(gl_Vertex.x,gl_Vertex.y,0,1.0);\n"
		"   gl_Position   = vec4(gl_Vertex.x*0.5+0.5,gl_Vertex.y*0.5+0.5,0,1);\n"
		"}\n";
	char fragmentSource[] =
		"uniform sampler2D texture0;\n"
		"uniform float blendFactor;\n"
		"void main(void) \n"
		"{\n"
		"    vec2 uv = gl_TexCoord[0].xy; \n"
		"    uv = uv * 0.5 + 0.5;\n"
		"    vec4 color =   texture2D(texture0, uv);\n"
		"    gl_FragColor = vec4(color.rgb, blendFactor);\n"
		"}\n";
	program->setName("overlay");
	program->addShader(new osg::Shader(osg::Shader::VERTEX, vertexSource));
	program->addShader(new osg::Shader(osg::Shader::FRAGMENT, fragmentSource));
	this->getOrCreateStateSet()->setAttributeAndModes(program, osg::StateAttribute::ON);
	this->getOrCreateStateSet()->addUniform(new osg::Uniform("texture0", 0));
	this->getOrCreateStateSet()->addUniform(new osg::Uniform("blendFactor", 1.0f));

	//this->getOrCreateStateSet()->setRenderBinDetails(100000, "RenderBin");
	this->getOrCreateStateSet()->setMode(GL_BLEND, osg::StateAttribute::OFF);
	//this->getOrCreateStateSet()->setRenderingHint(osg::StateSet::OPAQUE_BIN);

	viewer->addEventHandler(new BlendFactorEventHandler(this->getOrCreateStateSet()->getUniform("blendFactor")));
}

ScreenOverlay::~ScreenOverlay()
{

}

void ScreenOverlay::setTextureLayer(osg::Texture* tex)
{
	getOrCreateStateSet()->setTextureAttributeAndModes(0, tex, osg::StateAttribute::ON);
}