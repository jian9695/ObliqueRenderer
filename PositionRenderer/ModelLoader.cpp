#include "ModelLoader.h"
#include <osgDB/ReadFile>
#include <vector>

ModelLoader::ModelLoader()
{
}


ModelLoader::~ModelLoader()
{
}
std::vector<std::string> findSubdirs(std::string dir)
{
	WIN32_FIND_DATAA FindFileData;
	HANDLE hFind;
	std::string sPath;
	std::vector<std::string> MyVect;//"C:\\Documents and Settings\\yugesh\\Desktop\\*"
	sPath.assign(dir);
	hFind = FindFirstFileA(sPath.data(), &FindFileData);
	do
	{
		if (FindFileData.dwFileAttributes == 16)
		{
			if (std::string(FindFileData.cFileName) != "." && std::string(FindFileData.cFileName) != "..")
			{
				MyVect.push_back(FindFileData.cFileName);
			}

		}
	} while (FindNextFileA(hFind, &FindFileData));
	FindClose(hFind);
	//for (int i = 0; i<MyVect.size(); i++)
	//	cout << MyVect.at(i).data() << endl;
	return MyVect;
}
osg::Node* loadModels(std::string file)
{
	return osgDB::readNodeFile(file);
}
osg::Node* loadModels(std::vector<std::string> files)
{
	//if (files.size() == 1)
	//{
	//	osg::Node* nd = osgDB::readNodeFile(files[0]);
	//	return nd;
	//}

	osg::Group* scene = new osg::Group;
	for (int i = 0; i<files.size(); i++)
	{
		osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(files[i]);
		if (!node || !node.valid())
			continue;
		node->setName(files[i]);
		scene->addChild(node.get());
	}
	//scene->setInitialBound(bb);
	//scene->setCenter(bb.center());
	return scene;
}

osg::Node* ModelLoader::Load3DTiles(std::string indir)
{
	
	std::vector<std::string> files;
	std::vector<std::string> dirs = findSubdirs(indir + "*");
	for (size_t i = 0; i < dirs.size(); i++)
	{
		if (dirs[i] == "." || dirs[i] == "..")
			continue;
		files.push_back(indir + dirs[i] + "\\" + dirs[i] + ".osgb");
	}

	return loadModels(files);
}