#include "ModelLoader.h"
#include <osgDB/ReadFile>
#include <vector>
#include <sstream>
#include <qdir.h>
#include <qfile.h>
#include "osg/ComputeBoundsVisitor"
#include "osgUtil/SmoothingVisitor"
#include "osgDB/writeFile"
#include "osgDB/readFile"
#include "ShapeFile.h"
#include "GrassSolar.h"
#include "GDAL_DS.h"

ModelLoader::ModelLoader()
{
}


ModelLoader::~ModelLoader()
{
}

//std::vector<std::string> findSubdirs(std::string dir)
//{
//	WIN32_FIND_DATAA FindFileData;
//	HANDLE hFind;
//	std::string sPath;
//	std::vector<std::string> MyVect;//"C:\\Documents and Settings\\yugesh\\Desktop\\*"
//	sPath.assign(dir);
//	hFind = FindFirstFileA(sPath.data(), &FindFileData);
//	do
//	{
//		if (FindFileData.dwFileAttributes == 16)
//		{
//			if (std::string(FindFileData.cFileName) != "." && std::string(FindFileData.cFileName) != "..")
//			{
//				MyVect.push_back(FindFileData.cFileName);
//			}
//
//		}
//	} while (FindNextFileA(hFind, &FindFileData));
//	FindClose(hFind);
//	//for (int i = 0; i<MyVect.size(); i++)
//	//	cout << MyVect.at(i).data() << endl;
//	return MyVect;
//}

std::string getBaseName(std::string path)
{
		char last = path[path.length() - 1];
		if (last == '/' || last == '\\')
		{
				path = path.substr(0, path.length() - 1);
		}
		return QFileInfo(path.data()).baseName().toLocal8Bit().data();
}

std::vector<std::string> findSubdirs(std::string dir)
{
		std::vector<std::string> subdirs;
		QDir rootdir(dir.data());
		rootdir.setFilter(QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
		rootdir.setSorting(QDir::Name);
		std::vector<std::string> files;
		QFileInfoList list = rootdir.entryInfoList();
		for (int i = 0; i < list.size(); ++i) {
				QFileInfo fileInfo = list.at(i);
				std::string dir = (fileInfo.absoluteFilePath() + "/").toLocal8Bit().data();
				subdirs.push_back(dir);
		}
		return subdirs;
}

std::vector<std::string> findFiles(std::string indir, std::string match)
{
		std::vector<std::string> files;
		QDir input_dir(indir.data());
		input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
		input_dir.setSorting(QDir::Name);
		indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
		QFileInfoList list = input_dir.entryInfoList();
		for (int i = 0; i < list.size(); ++i) {
				QFileInfo fileInfo = list.at(i);
				std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
				if (/*!fileInfo.fileName().contains(match.data(),Qt::CaseInsensitive) || */match != "" && !fileInfo.fileName().endsWith(match.data(), Qt::CaseInsensitive))
						continue;
				files.push_back(fileInfo.absoluteFilePath().toLocal8Bit().data());
		}
		return files;
}

std::vector<std::string> findLeafTileFiles(const std::vector<std::string>& tileFiles, const std::string& masterTileName)
{
		std::vector<int> levels;
		int maxLevel = -1;
		for (size_t i = 0; i < tileFiles.size(); i++)
		{
				std::string tilename = QFileInfo(tileFiles[i].data()).baseName().toLocal8Bit().data();
				if (tilename.length() - masterTileName.length() < 3)
				{
						levels.push_back(-1);
						continue;
				}

				tilename = tilename.substr(masterTileName.length() + 2, tilename.size() - masterTileName.length() - 2);
				std::string levelStr = "";
				for (size_t j = 0; j < tilename.length(); j++)
				{
						if (!isdigit(tilename[j]))
								break;
						levelStr += tilename[j];
				}

				int level = -1;
				std::stringstream ss;
				ss << levelStr;
				ss >> level;
				levels.push_back(level);
				if (maxLevel < level)
						maxLevel = level;
		}

		std::vector<std::string> leafTiles;
		for (size_t i = 0; i < tileFiles.size(); i++)
		{
				if (levels[i] == maxLevel)
				{
						leafTiles.push_back(tileFiles[i]);
				}
		}
		return leafTiles;
}

std::vector<std::string> findMasterTiles(std::string indir)
{
		std::vector<std::string> files;
		std::vector<std::string> dirs = findSubdirs(indir);
		for (size_t i = 0; i < dirs.size(); i++)
		{
				if (dirs[i].find("Tile_") == std::string::npos)
						continue;
				std::string masterTileName = getBaseName(dirs[i]);
				std::string masterTilePath = dirs[i] + masterTileName + ".osgb";
				files.push_back(masterTilePath);
		}
		return files;
}

OGRPolygon* toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb)
{

		OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
		OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
		linearRing->addPoint(bb.MinX, bb.MinY);
		linearRing->addPoint(bb.MinX, bb.MaxY);
		linearRing->addPoint(bb.MaxX, bb.MaxY);
		linearRing->addPoint(bb.MaxX, bb.MinY);
		linearRing->addPoint(bb.MinX, bb.MinY);

		poPolygon->addRing(linearRing);//also crashed
		return poPolygon;
}

OGRPolygon* toOGRPolygon(OGRLayer* layer, const osg::BoundingBoxd& bb)
{
		OGREnvelope ogrBB;
		ogrBB.MaxX = bb.xMax();
		ogrBB.MaxY = bb.yMax();
		ogrBB.MinX = bb.xMin();
		ogrBB.MinY = bb.yMin();
		return toOGRPolygon(layer, ogrBB);
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
		std::vector<std::string> files = findMasterTiles(indir);
	return loadModels(files);
}

osg::Node* ModelLoader::Load3DTiles(std::string indir, osg::BoundingBox mask, bool intersects)
{
		double area = (mask.xMax() - mask.xMin()) * (mask.yMax() - mask.yMin());
		std::vector<std::string> files;
		std::vector<std::string> dirs = findSubdirs(indir);
		for (size_t i = 0; i < dirs.size(); i++)
		{
				if (dirs[i].find("Tile_") == std::string::npos)
						continue;
				std::string masterTileName = getBaseName(dirs[i]);
				std::string masterTilePath = dirs[i] + masterTileName + ".osgb";

				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(masterTilePath);
				osg::ComputeBoundsVisitor cbs;
				node->accept(cbs);
				osg::BoundingBox masterTileBB = cbs.getBoundingBox();
				bool intersectsBB = masterTileBB.intersects(mask);
				osg::BoundingBox intersection = masterTileBB.intersect(mask);
				if ((intersection.xMax() - intersection.xMin()) * (intersection.yMax() - intersection.yMin()) < area * 0.1)
				{
						intersectsBB = false;
				}
				if (intersectsBB != intersects)
						continue;
				files.push_back(masterTilePath);
		}
		return loadModels(files);
}

osg::Node* ModelLoader::Load3DTiles(std::string indir, std::vector<std::string> maskTiles, bool include)
{

		std::set<std::string> tileset;
		for (size_t i = 0; i < maskTiles.size(); i++)
		{
				tileset.insert(maskTiles[i]);
		}

		std::vector<std::string> files;
		std::vector<std::string> dirs = findSubdirs(indir);
		for (size_t i = 0; i < dirs.size(); i++)
		{
				if (dirs[i].find("Tile_") == std::string::npos)
						continue;
				std::string masterTileName = getBaseName(dirs[i]);
				std::string masterTilePath = dirs[i] + masterTileName + ".osgb";
				bool found = false;
				if (tileset.find(masterTileName) != tileset.end())
				{
						found = true;
				}
				if (found != include)
						continue;
				files.push_back(masterTilePath);
		}
		return loadModels(files);
}

void ModelLoader::CopyLeafTiles(std::string indir, std::string outdir, osg::BoundingBox bb)
{
		bb.zMin() = -10000;
		bb.zMax() = 100000;
		QDir(outdir.data()).mkpath(".");
		QString qoutdir = QDir(outdir.data()).absolutePath() + "/";
		std::vector<std::string> subdirs = findSubdirs(indir);
		for each (std::string subdir in subdirs)
		{
				if (subdir.find("Tile_") == std::string::npos)
						continue;

				std::string masterTileName = getBaseName(subdir);
				std::string masterTilePath = subdir + masterTileName + ".osgb";
				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(masterTilePath);
				osg::ComputeBoundsVisitor cbs;
				node->accept(cbs);
				osg::BoundingBox masterTileBB = cbs.getBoundingBox();
				if (!masterTileBB.intersects(bb))
						continue;

				std::vector<std::string> allTiles = findFiles(subdir, ".osgb");
				std::vector<std::string> leafTiles = findLeafTileFiles(allTiles, masterTileName);
				for each (std::string leafTileFile in leafTiles)
				{
						osg::ComputeBoundsVisitor tileCBS;
						node = osgDB::readNodeFile(leafTileFile);
						node->accept(tileCBS);
						if (!tileCBS.getBoundingBox().intersects(bb))
								continue;
						QFile::copy(leafTileFile.data(), qoutdir + QFileInfo(leafTileFile.data()).fileName());
						printf("%s\n", QFileInfo(leafTileFile.data()).fileName().toLocal8Bit().data());
				}
		}
}

void ModelLoader::CopyLeafTiles(std::string indir, std::vector<std::string> tilenames, std::string outfile)
{
		osg::ref_ptr<osg::Group> group = new osg::Group;
		for each (std::string tilename in tilenames)
		{
				std::string subdir = indir + tilename + "/";
				if (subdir.find("Tile_") == std::string::npos)
						continue;

				//std::string masterTileName = getBaseName(subdir);
				//std::string masterTilePath = subdir + masterTileName + ".osgb";
			//	osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(masterTilePath);
				std::vector<std::string> allTiles = findFiles(subdir, ".osgb");
				std::vector<std::string> leafTiles = findLeafTileFiles(allTiles, tilename);
				for each (std::string leafTileFile in leafTiles)
				{
						osg::ComputeBoundsVisitor tileCBS;
						osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(leafTileFile);
						osgDB::writeNodeFile(*node, "ljm.osg");
						node = osgDB::readNodeFile("ljm.osg");
						group->addChild(node.get());
						//QFile::copy(leafTileFile.data(), qoutdir + QFileInfo(leafTileFile.data()).fileName());
						printf("%s\n", QFileInfo(leafTileFile.data()).fileName().toLocal8Bit().data());
				}
		}
		osgUtil::SmoothingVisitor smooth;
		group->accept(smooth);
		osgDB::writeNodeFile(*group, outfile);
}

osg::BoundingBoxd ModelLoader::CalBound(std::string indir, const std::vector<std::string>& tiles)
{
		osg::BoundingBoxd bb;
		bb.init();
		for (size_t i = 0; i < tiles.size(); i++)
		{
				std::string masterTilePath = indir + tiles[i] + "/" + tiles[i] + ".osgb";
				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(masterTilePath);
				osg::ComputeBoundsVisitor cbs;
				node->accept(cbs);
				osg::BoundingBoxd tileBB = cbs.getBoundingBox();
				bb.expandBy(tileBB);
		}
		return bb;
}

void ModelLoader::ConvertAspectQGIS2Grass(std::string infile, std::string outfile)
{
		GDAL_DS<float>* dt = new GDAL_DS<float>();
		dt->open(infile);
		GDAL_DS<float>* newdt = new GDAL_DS<float>();

		memcpy(newdt->adfGeoTransform, dt->adfGeoTransform, sizeof(double) * 6);
		//newdt->adfGeoTransform[0] = bb.MinX;
		//newdt->bound = destBB;
		newdt->numbands = 1;
		newdt->ncols = dt->ncols;
		newdt->nrows = dt->nrows;
		newdt->projection = dt->projection;
		newdt->pszFormat = dt->pszFormat;
		newdt->slice = dt->slice;
		newdt->create(outfile);
		double nodata = dt->getNoData(1);
		float* data = dt->readData(1);

		for (size_t i = 0; i < newdt->slice; i++)
		{
				float aspect = data[i];
				if (aspect == nodata)
						continue;
				aspect = 360 - aspect;
				aspect = aspect + 90;
				if (aspect > 360)
						aspect -= 360;
				data[i] = aspect;
		}
		newdt->writeData(1, data, nodata);
		delete newdt;
		delete dt;

}

void getBlockStats(std::string filename, BBWrapper bb)
{
		int startCol, startRow, startIndex, endCol, endRow, endIndex;
		GDAL_DS<float>* dsm = new GDAL_DS<float>();
		dsm->open(filename, GA_Update);
		std::tie(startCol, startRow, startIndex) = dsm->getCellIndex(bb.xMin(), bb.yMax());
		std::tie(endCol, endRow, endIndex) = dsm->getCellIndex(bb.xMax(), bb.yMin());
		if (startCol < 0)
				startCol = 0;
		if (startRow < 0)
				startRow = 0;
		if (endCol > dsm->ncols - 1)
				endCol = dsm->ncols - 1;
		if (endRow > dsm->nrows - 1)
				endRow = dsm->nrows - 1;

		double blockXSize = endCol - startCol + 1;
		double blockYSize = endRow - startRow + 1;
		float* blockBuf = new float[blockXSize * blockYSize];
		dsm->m_dataset->GetRasterBand(1)->RasterIO(GF_Read, startCol, startRow, blockXSize, blockYSize, blockBuf, blockXSize, blockYSize, dsm->getType(), 0, 0);
		float min = 1000000;
		float max = -min;
		float* pValue = blockBuf;
		double nodata = dsm->getNoData(1);
		for (size_t i = 0; i < blockXSize * blockYSize; i++)
		{
				float val = *pValue;
				if (val == nodata)
				{
						pValue++;
						continue;
				}
		   
				if (min > val)
						min = val;
				if (max < val)
						max = val;

				pValue++;
		}
		printf("min=%f,max=%f\n", min, max);
		delete dsm;
}

void ModelLoader::Test()
{
		AABBox aabb(osg::Vec3d(-0.5, -0.5, 0.0), osg::Vec3d(0.5, 0.5, 0.0));
		Ray ray(osg::Vec3d(0.0, 0.0, 0.0), osg::Vec3d(1.0, 1.0, 0.0));
		osg::Vec3d dir(osg::Vec3d(0.5, 0.5, 0.0) - osg::Vec3d(1.0, 0.5, 0.0));
		dir.normalize();
		ray = Ray(osg::Vec3d(1.0, 0.5, 0.0), dir);
		double dist;
		bool intersects = aabb.intersect(ray, dist);
		intersects = aabb.intersect(Ray(osg::Vec3d(1.0, 0.5, 5.0), osg::Vec3d(5.0, 5.0, 5.0)));
		std::vector<std::string> tilenames;
		tilenames.push_back("Tile_-002_-016");
		tilenames.push_back("Tile_-002_-017");
		tilenames.push_back("Tile_-003_-016");
		tilenames.push_back("Tile_-003_-017");
		BBWrapper bb = CalBound("E:/Data/weihai/Data/", tilenames);
		
		double bbScale = 4;
		BBWrapper expandedBB(bb.center().x() - bb.xhalfsize() * bbScale, bb.center().y() - bb.yhalfsize() * bbScale, bb.zMin(),
				bb.center().x() + bb.xhalfsize() * bbScale, bb.center().y() + bb.yhalfsize() * bbScale, bb.zMax());

		//getBlockStats("E:/Code/Weihai_DSM_025.tif", expandedBB);
		//getBlockStats("E:/Code/Weihai_DSM_050.tif", expandedBB);
		//getBlockStats("E:/Code/Weihai_DSM_075.tif", expandedBB);
		//getBlockStats("E:/Code/Weihai_DSM_1m.tif", expandedBB);
		int startCol, startRow, startIndex, endCol, endRow, endIndex;
		GDAL_DS<float>* dsm = new GDAL_DS<float>();
		dsm->open("E:/Code/Weihai_DSM_025.tif", GA_Update);
		std::tie(startCol, startRow, startIndex) = dsm->getCellIndex(expandedBB.xMin(), expandedBB.yMax());
		std::tie(endCol, endRow, endIndex) = dsm->getCellIndex(expandedBB.xMax(), expandedBB.yMin());
		if (startCol < 0)
				startCol = 0;
		if (startRow < 0)
				startRow = 0;
		if (endCol > dsm->ncols - 1)
				endCol = dsm->ncols - 1;
		if (endRow > dsm->nrows - 1)
				endRow = dsm->nrows - 1;

		//dsm->downSampleMaximum("E:/Code/Weihai_DSM_075.tif", 3);
	 double blockXSize = endCol - startCol + 1;
		double blockYSize = endRow - startRow + 1;
		float* blockBuf = new float[blockXSize * blockYSize];
		dsm->m_dataset->GetRasterBand(1)->RasterIO(GF_Read, startCol, startRow, blockXSize, blockYSize, blockBuf, blockXSize, blockYSize, dsm->getType(), 0, 0);
		double sum = 0;
		double count = 0;
		float* pValue = blockBuf;
		double nodata = dsm->getNoData(1);
		for (size_t i = 0; i < blockXSize * blockYSize; i++)
		{
				float val = *pValue;
				if (val == nodata)
				{
						pValue++;
						continue;
				}
				sum += val;
				count += 1;
				pValue++;
		}

		printf("%f\n", sum / count);

		//dsm->close();
		//std::vector<std::string> allTiles = findFiles("E:/Code/WeihaiLeafTiles/", ".osgb");
		//osg::ref_ptr<osg::Group> group = new osg::Group;
		//for each (std::string leafTileFile in allTiles)
		//{
		//		osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(leafTileFile);
		//		osgDB::writeNodeFile(*node, "ljm.osg");
		//		node = osgDB::readNodeFile("ljm.osg");
		//		group->addChild(node.get());
		//}
		//osgUtil::SmoothingVisitor smooth;
		//group->accept(smooth);
		//osgDB::writeNodeFile(*group, "WeihaiLeafTiles.osgb");


		//GDAL_DS<unsigned char>* dom = new GDAL_DS<unsigned char>();
		//dom->open("E:/Code/Weihai_DOM_025.tif", GA_Update);
		////dom->setGeoTransform(dsm->adfGeoTransform);
		//dsm->close();
		//dom->close();
}

void ModelLoader::TileBoundary2Shapefile(std::string indir, std::string outfile)
{
		ShapeFile tileMap;
		tileMap.create(outfile);
		OGRFeatureDefn *poFDefn = tileMap.poLayer->GetLayerDefn();
		int tileNameIndex = tileMap.getOrCreateField("Name", OGRFieldType::OFTString);
		int tileMinXIndex = tileMap.getOrCreateField("MinX", OGRFieldType::OFTReal);
		int tileMaxXIndex = tileMap.getOrCreateField("MaxX", OGRFieldType::OFTReal);
		int tileMinYIndex = tileMap.getOrCreateField("MinY", OGRFieldType::OFTReal);
		int tileMaxYIndex = tileMap.getOrCreateField("MaxY", OGRFieldType::OFTReal);
		int tileFilePathIndex = tileMap.getOrCreateField("Path", OGRFieldType::OFTString);

		std::vector<std::string> files;
		std::vector<std::string> dirs = findSubdirs(indir);
		for (size_t i = 0; i < dirs.size(); i++)
		{
				if (dirs[i].find("Tile_") == std::string::npos)
						continue;
				std::string masterTileName = getBaseName(dirs[i]);
				std::string masterTilePath = dirs[i] + masterTileName + ".osgb";
				osg::ref_ptr<osg::Node> node = osgDB::readNodeFile(masterTilePath);
				osg::ComputeBoundsVisitor cbs;
				node->accept(cbs);
				osg::BoundingBoxd masterTileBB = cbs.getBoundingBox();
				if (!masterTileBB.valid())
						continue;
				if(masterTileBB.xMin() < -100000000 || masterTileBB.xMax() > 100000000 || masterTileBB.yMin() < -100000000 || masterTileBB.yMax() > 100000000)
						continue;
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(tileMap.poLayer->GetLayerDefn());
				poFeaPolygon->SetField(tileNameIndex, masterTileName.data());
				poFeaPolygon->SetField(tileMinXIndex, masterTileBB.xMin());
				poFeaPolygon->SetField(tileMaxXIndex, masterTileBB.xMax());
				poFeaPolygon->SetField(tileMinYIndex, masterTileBB.yMin());
				poFeaPolygon->SetField(tileMaxYIndex, masterTileBB.yMax());
				poFeaPolygon->SetField(tileFilePathIndex, masterTilePath.data());
				OGRPolygon *poPolygon = toOGRPolygon(tileMap.poLayer, masterTileBB);
				poFeaPolygon->SetGeometry(poPolygon);
				tileMap.poLayer->CreateFeature(poFeaPolygon);
				OGRFeature::DestroyFeature(poFeaPolygon);
		}
		tileMap.close();
}

//int main(int argc, char **argv)
//{
//		AABBox box(Vec3f(-1), Vec3f(1));
//		gen.seed(0);
//		for (uint32_t i = 0; i < 16; ++i) {
//				Vec3f randDir(2 * dis(gen) - 1, 2 * dis(gen) - 1, 2 * dis(gen) - 1);
//				randDir.normalize();
//				Ray ray(Vec3f(0), randDir);
//				float t;
//				if (box.intersect(ray, t)) {
//						Vec3f Phit = ray.orig + ray.dir * t;
//						std::cerr << ray.orig << " " << Phit << std::endl;
//				}
//		}
//		return 0;
//}

ShadowCaster::~ShadowCaster()
{
		for each (GDAL_DS<float>* ds in m_rasters)
		{
				delete ds;
		}
		m_rasters.clear();

		delete[] m_slopeData;
		delete[] m_aspectData;
		m_slopeData = nullptr;
		m_aspectData = nullptr;
}

bool ShadowCaster::Intersects(const Ray& r, const BBWrapper& aabb, double& t)
{
		double tmin, tmax, tymin, tymax, tzmin, tzmax;

		tmin = (aabb.bounds[r.sign[0]].x() - r.orig.x()) * r.invdir.x();
		tmax = (aabb.bounds[1 - r.sign[0]].x() - r.orig.x()) * r.invdir.x();
		tymin = (aabb.bounds[r.sign[1]].y() - r.orig.y()) * r.invdir.y();
		tymax = (aabb.bounds[1 - r.sign[1]].y() - r.orig.y()) * r.invdir.y();

		if ((tmin > tymax) || (tymin > tmax))
				return false;

		if (tymin > tmin)
				tmin = tymin;
		if (tymax < tmax)
				tmax = tymax;

		tzmin = (aabb.bounds[r.sign[2]].z() - r.orig.z()) * r.invdir.z();
		tzmax = (aabb.bounds[1 - r.sign[2]].z() - r.orig.z()) * r.invdir.z();

		if ((tmin > tzmax) || (tzmin > tmax))
				return false;

		if (tzmin > tmin)
				tmin = tzmin;
		if (tzmax < tmax)
				tmax = tzmax;

		t = tmin;

		if (t < 0) {
				t = tmax;
				if (t < 0) return false;
		}

		return true;
}

void ShadowCaster::Initialize(std::vector<std::string> rasterFiles, std::vector<osg::BoundingBoxd> bounds, std::string slopeFile, std::string aspectFile)
{
		if (m_slopeData)
				delete[] m_slopeData;
		if (m_aspectData)
				delete[] m_aspectData;
		for each (GDAL_DS<float>* ds in m_rasters)
		{
				delete ds;
		}
		m_rasters.clear();
		m_bounds.clear();
		if (rasterFiles.size() != bounds.size())
				return;

		GDAL_DS<float>* slopeDS = new GDAL_DS<float>();
		slopeDS->open(slopeFile);
		m_slopeData = slopeDS->readData(1);
		delete slopeDS;
		GDAL_DS<float>* aspectDS = new GDAL_DS<float>();
		aspectDS->open(aspectFile);
	 m_aspectData = aspectDS->readData(1);
		delete aspectDS;

		for (size_t i = 0; i < rasterFiles.size(); i++)
		{
				GDAL_DS<float>* ds = new 	GDAL_DS<float>();
				ds->open(rasterFiles[i]);
				osg::BoundingBoxd bb = bounds[i];
				ds->readCache(bb);
				m_rasters.push_back(ds);
				m_bounds.push_back(bb);
		}
}

void ShadowCaster::CastShadow(osg::BoundingBoxd bb, std::string outfile, osg::Vec3d sundir)
{

}

struct Color3
{
		unsigned char m_r, m_g, m_b;
		Color3()
		{
				m_r = 0; m_g = 0; m_b = 0;
		}

		Color3(unsigned char r, unsigned char g, unsigned char b)
		{
				m_r = r; m_g = g; m_b = b;
		}
};

void ShadowCaster::ComputeViewshed(double px, double py)
{
		int viewshedSize = 32;
		osg::ref_ptr<osg::Image> colorImage = new osg::Image;
		colorImage->allocateImage(viewshedSize, viewshedSize, 1, GL_RGB, GL_UNSIGNED_BYTE);
		Color3* colorData = (Color3*)colorImage->data();
		double resol = 1.0 / viewshedSize;
		double totalarea = 0;
		double skyarea = 0;
		int npixel = 0;
		for (unsigned int row = 0; row < viewshedSize; row++)
		{
				double y = 0.5 - row * resol - 0.5 * resol;
				for (unsigned int col = 0; col < viewshedSize; col++)
				{
						double x = -0.5 + col * resol + 0.5 * resol;
						double dist2ori = sqrt(x * x + y * y);
						if (dist2ori >= 1.0)
						{
								*colorData = Color3(0, 0, 255);
								colorData++;
								continue;
						}
						double zenithD = sqrt(x * x + y * y) * 90.0;//in degrees
						if (zenithD <= 0.000000001)
								zenithD = 0.000000001;
						double zenithR = zenithD * 3.1415926 / 180.0;
						double x2 = 0.0;
						double y2 = 1.0;
						double	cosa = (x * x2 + y * y2) / sqrt((x * x + y * y) * (x2 * x2 + y2 * y2));
						double lon = acos(cosa) * 180.0 / 3.1415926;
						lon = 360.0 - lon;
						lon = 1.0 - (lon / 360.0);
						double lat = dist2ori;
				}
		}
}

std::tuple<int, int, int, int, double, BBWrapper> ShadowCaster::GetCellIndex(double x, double y)
{
		for (long i = 0; i < m_rasters.size(); i++)
		{
				GDAL_DS<float>* ds = (GDAL_DS<float>*)m_rasters[i];
				int col, row, index, dsIndex;
				std::tie(col, row, index) = ds->getCellIndexCached(x, y);
				if (index >= 0)
				{
						double z = ds->m_cache[index];
						if (z == ds->m_nodata)
								continue;
						double xmin = x - ds->adfGeoTransform[1];
						double xmax = x + ds->adfGeoTransform[1];
						double ymin = y - abs(ds->adfGeoTransform[5]);
						double ymax = y + abs(ds->adfGeoTransform[5]);
						double zmin = 0;
						double zmax = max(z,0);
						BBWrapper bb = BBWrapper(xmin, ymin, zmin, xmax, ymax, zmax);
						return std::make_tuple(col, row, index, i, ds->adfGeoTransform[1], bb);
				}
		}
		return std::make_tuple(-1, -1, -1, -1, -1, BBWrapper());
}

bool ShadowCaster::Intersects(const Ray& ray, double& t, double& slope, double& aspect)
{
		slope = 0;
		aspect = 0;
		//Ray r = Ray(r.orig,r.dir);
		Ray groundRay = Ray(osg::Vec2d(ray.orig.x(), ray.orig.y()), osg::Vec2d(ray.dir.x(), ray.dir.y()));
		osg::Vec3d pos = groundRay.orig;
		GDAL_DS<float>* ds = (GDAL_DS<float>*)m_rasters[0];
		double cellSize = ds->adfGeoTransform[1];
		pos = groundRay.orig + groundRay.dir * (cellSize * 0.5);
		int step = 0;
		while (true)
		{
				int col, row, index, dsIndex;
				BBWrapper bb;
				t = 0;
				std::tie(col, row, index, dsIndex, cellSize, bb) = GetCellIndex(pos.x(), pos.y());
				if (col == -1)
						return false;
				double t2;
				AABBox aabb(bb._min, bb._max);
				aabb.intersect(ray, t2);
				if (Intersects(ray, bb, t))
				{
						osg::Vec3d hit = ray.orig + ray.dir * t;
						if (dsIndex == 0) 
						{
								slope = m_slopeData[index];
								aspect = m_aspectData[index];
						}
						return true;
				}
				step++;
				pos = pos + groundRay.dir * (cellSize * 0.5);
		}
		return false;
}

void ShadowCaster::ComputeHorizonAngleMaps(std::string infile, std::string outfile, double startAngle, double endAngle, double stepAngle)
{
		GDAL_DS<float>* ds = new GDAL_DS<float>();
		ds->open(infile);
		float* data = ds->readData(1);
		std::vector<double> angles;
		double curAngle = startAngle;
		while (curAngle <= endAngle)
		{
				if (curAngle >= endAngle)
				{
						if (curAngle >= 360)
								curAngle -= 360;
						if (curAngle == startAngle && angles.size() > 0)
								break;
						angles.push_back(curAngle);
						break;
				}
				angles.push_back(curAngle);
				curAngle += stepAngle;
		}
		double nodata = ds->getNoData(1);

		GDAL_DS<unsigned short>* dsOutput = new GDAL_DS<unsigned short>();
		dsOutput->setDSInfo((GDAL_DSInfo*)ds);
		dsOutput->numbands = angles.size();
		dsOutput->create(outfile);

		for (size_t i = 0; i < angles.size(); i++)
		{
				curAngle = angles[i];
				osg::Vec3d sundir = GrassSolar::solarAngle2Vector(0, curAngle);
				unsigned short* output = ComputeHorizonAngleMap(ds, data, nodata, osg::Vec2d(sundir.x(), sundir.y()));
				dsOutput->writeData(i + 1, output, 100);
				dsOutput->m_dataset->FlushCache();
		}
		delete dsOutput;
		delete ds;
		delete[] data;
}

unsigned short* ShadowCaster::ComputeHorizonAngleMap(void* pDS, float*& data, const double& nodata, const osg::Vec2d& dir)
{
		GDAL_DS<float>* ds = (GDAL_DS<float>*)pDS;
		unsigned short* output = new unsigned short[ds->slice];
		unsigned short* poutput = output;
		for (int row = 0; row < ds->nrows; row++)
		{
				double y = ds->bound.MaxY - ds->adfGeoTransform[1] * row - ds->adfGeoTransform[1] * 0.5;
				for (int col = 0; col < ds->ncols; col++)
				{
						double x = ds->bound.MinX + ds->adfGeoTransform[1] * col + ds->adfGeoTransform[1] * 0.5;
						int index = col + row * ds->ncols;
						float elev = data[index];
						unsigned short encodedHorizon = 100;
						if (elev > -500 && elev < 10000)
						{
								Ray groundRay(osg::Vec2d(x, y), dir);
								double horizonAngle = ds->calHorizonAngle(groundRay, data);
								if (horizonAngle < 0)
										horizonAngle = 0;
								double intPart;
								double fractPart = modf(horizonAngle, &intPart);
								encodedHorizon = (unsigned short)(intPart * 100) + (unsigned short)(fractPart * 100);
						}
						*poutput = encodedHorizon;
						poutput++;
				}
				printf("%d/%d\n", row, ds->nrows);
		}
		return output;
}

//void ShadowCaster::ComputeHorizonAngleMapsFast(std::string infile, std::string outfile, double startAngle, double endAngle, double stepAngle)
//{
//		GDAL_DS<float>* ds = new GDAL_DS<float>();
//		ds->open(infile);
//		float* data = ds->readData(1);
//		QuadTree* tree = new QuadTree;
//		tree->setInfo(*ds);
//		tree->setData(data);
//		tree->buildChild();
//		tree->m_child->buildChild();
//		std::vector<QuadTree*> trees;
//		trees.push_back(tree->m_child->m_child);
//		trees.push_back(tree->m_child);
//		trees.push_back(tree);
//
//
//
//		std::vector<double> angles;
//		double curAngle = startAngle;
//		while (curAngle <= endAngle)
//		{
//				if (curAngle >= endAngle)
//				{
//						if (curAngle >= 360)
//								curAngle -= 360;
//						if (curAngle == startAngle && angles.size() > 0)
//								break;
//						angles.push_back(curAngle);
//						break;
//				}
//				angles.push_back(curAngle);
//				curAngle += stepAngle;
//		}
//		double nodata = ds->getNoData(1);
//
//		GDAL_DS<unsigned short>* dsOutput = new GDAL_DS<unsigned short>();
//		dsOutput->setDSInfo((GDAL_DSInfo*)ds);
//		dsOutput->numbands = angles.size();
//		dsOutput->create(outfile);
//
//		for (size_t i = 0; i < angles.size(); i++)
//		{
//				curAngle = angles[i];
//				osg::Vec3d sundir = GrassSolar::solarAngle2Vector(0, curAngle);
//				unsigned short* output = ComputeHorizonAngleMap(ds, data, nodata, osg::Vec2d(sundir.x(), sundir.y()));
//				dsOutput->writeData(i + 1, output, 100);
//				dsOutput->m_dataset->FlushCache();
//		}
//
//		delete dsOutput;
//		delete ds;
//		delete[] data;
//}
//
//void QuadTree::buildChild()
//{
//		m_child = new QuadTree;;
//		double resolX = this->adfGeoTransform[1] * 2;
//		double resolY = this->adfGeoTransform[5] * 2;
//		m_child->ncols = (int)(ceil((double)ncols / 2.0));
//		m_child->nrows = (int)(ceil((double)nrows / 2.0));
//		m_child->bound.MinX = bound.MinX;
//		m_child->bound.MaxY = bound.MaxY;
//		m_child->bound.MaxX = bound.MinX + m_child->ncols * resolX;
//		m_child->bound.MinY = bound.MaxY + m_child->nrows * resolY;
//		memcpy(m_child->adfGeoTransform, adfGeoTransform, sizeof(double) * 6);
//		m_child->adfGeoTransform[1] = resolX;
//		m_child->adfGeoTransform[5] = resolY;
//		m_child->numbands = 1;
//		m_child->projection = this->projection;
//		m_child->pszFormat = this->pszFormat;
//		m_child->slice = m_child->ncols * m_child->nrows;
//		m_child->nodata = this->nodata;
//		if (m_child->nodata == 0)
//				m_child->nodata = -9999;
//		m_child->m_data.resize(m_child->slice);
//		float block[4];
//		int index = 0;
//		for (int row = 0; row < m_child->nrows; row++)
//		{
//				for (int col = 0; col  < m_child->ncols; col++)
//				{
//						block[0] = getValue(col * 2,     row * 2);
//						block[1] = getValue(col * 2 + 1, row * 2);
//						block[2] = getValue(col * 2,     row * 2 + 1);
//						block[3] = getValue(col * 2 + 1, row * 2 + 1);
//						m_child->m_data[index] = sample(block);
//						index++;
//				}
//		}
//}
