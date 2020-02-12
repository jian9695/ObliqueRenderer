#include "SVFComputeTools.h"



void CameraBuffer::update()
{
	osg::Matrix localOffset;
	localOffset.makeLookAt(pos, pos + dir * 100, up);
	osg::Matrix viewMatrix = localOffset;
	setReferenceFrame(osg::Camera::ABSOLUTE_RF);
	float nearDist = 0.1;
	setProjectionMatrixAsFrustum(-nearDist, nearDist, -nearDist, nearDist, nearDist, 10000.0);
	setViewMatrix(viewMatrix);
	setClearColor(osg::Vec4(0, 0, 0, 0));
}

SVFComputeTools::SVFComputeTools()
{
}

SVFComputeTools::~SVFComputeTools()
{
}
