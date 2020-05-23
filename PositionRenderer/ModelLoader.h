#pragma once
#include <Windows.h>
#include <string>
#include <osg/Node>
class ModelLoader
{
public:
	ModelLoader();
	~ModelLoader();
	static osg::Node* Load3DTiles(std::string indir);
};

