#include "GeometryVisitor.h"
#include "osgDB\ReadFile"
#include <map>

class FaceVisitor
{
public:
	FaceVisitor(osg::PrimitiveSet* pset)
	{
		pSet = pset;
	}
	std::vector<unsigned int> getFaceIndices()
	{

		std::vector<unsigned int> indices;
		unsigned int idx = pSet->getNumIndices();
		unsigned int numofprims;


		if(pSet->getMode() == osg::PrimitiveSet::TRIANGLE_FAN)
		{
			numofprims = idx-2;
			for (unsigned int i=0;i<numofprims;i++)
			{
				indices.push_back(pSet->index(0));
				indices.push_back(pSet->index(i + 1));
				indices.push_back(pSet->index(i + 2));
			}
		}
		else if (pSet->getMode() == osg::PrimitiveSet::TRIANGLES)
		{
			for (unsigned int i=0;i<idx;i++)
			{
				indices.push_back(pSet->index(i));
			}
		}
		else if (pSet->getMode() == osg::PrimitiveSet::TRIANGLE_STRIP)
		{
			numofprims = (idx - 2) / 2;
			for (unsigned int i=0;i<numofprims;i++)
			{
				indices.push_back(pSet->index(i * 2));
				indices.push_back(pSet->index(i * 2 + 1));
				indices.push_back(pSet->index(i * 2 + 2));

				indices.push_back(pSet->index(i * 2 + 2));		
				indices.push_back(pSet->index(i * 2 + 1));
				indices.push_back(pSet->index(i * 2 + 3));
			}

			if(numofprims + 2 % 2 != 0)
			{
				indices.push_back(pSet->index(idx-3));		
				indices.push_back(pSet->index(idx-2));
				indices.push_back(pSet->index(idx-1));
			}
		}
		return indices;
	}

private:
	osg::PrimitiveSet* pSet;
};
void clampVector(osg::Vec2& val)
{
	if(val.x()<0)
		val.x()=0;
	if(val.y()<0)
		val.y()=0;
	if(val.x()>1)
		val.x()=1;
	if(val.y()>1)
		val.y()=1;
}


void GeometryVisitor::apply( osg::Geode& node )
{
	
	osg::MatrixList matlist = node.getWorldMatrices();
	osg::Matrix matWorld = osg::Matrix::identity();

	/*计算goede的坐标变换 ？*/
	for (unsigned int i=0; i<matlist.size(); i++)
	{
		matWorld = matlist[i] * matWorld;
	}

	for (int i=0; i<node.getNumDrawables(); i++)
	{
		osg::Drawable* drawable = node.getDrawable(i);
		osg::Geometry* geom = dynamic_cast<osg::Geometry*>(drawable);
		if(!geom) continue;

		//统计面积
		osg::Vec3Array* vertices = (osg::Vec3Array *) geom->getVertexArray();
		for(int j=0; j<geom->getNumPrimitiveSets(); ++j)
		{
			osg::PrimitiveSet* primitiveSet = geom->getPrimitiveSet(j);
			FaceVisitor visitor(primitiveSet);
			std::vector<unsigned int> indices = visitor.getFaceIndices();

			for(unsigned int k=0;k< indices.size();k++) 
			{
	             osg::Vec3 pos=(*vertices)[indices[k]];
				 osg::Vec2 uv((pos.x()-g_mExtent.xMin())/(g_mExtent.xMax()-g_mExtent.xMin())
					         ,(pos.y()-g_mExtent.yMin())/(g_mExtent.yMax()-g_mExtent.yMin()));
				 //"gl_TexCoord[2] = vec4((pos.x-min.x)/(max.x-min.x),(pos.y-min.y)/(max.y-min.y),0,1);\n"

				clampVector(uv);
				osg::Vec4 color = g_pDSM->getColor(uv);
				if(color.r() > -3000)
				   pos.z() = color.r();
				(*vertices)[indices[k]] = pos;
				//osg::Vec3 p1 = (*vertices)[k1] * matWorld;
			}
		}

		
	}

	
}

