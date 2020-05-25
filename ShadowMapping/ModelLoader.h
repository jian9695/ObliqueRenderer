#pragma once
#include <Windows.h>
#include <string>
#include <osg/Node>
#include <osg/BoundingBox>
#include <cstdlib>
#include <random>
//#include "GDAL_DS.h"

//class QuadTree : GDAL_DSInfo
//{
//public:
//		enum Sampler
//		{
//				Min,
//				Max,
//				Mean
//		};
//
//		QuadTree() 
//		{ 
//				m_sampler = Sampler::Min;
//				m_child = nullptr; 
//		}
//
//		~QuadTree() 
//		{ 
//				if (m_child)
//						delete m_child;
//		}
//		void setInfo(const GDAL_DSInfo& info)
//		{
//				this->ncols = info.ncols;
//				this->nrows = info.nrows;
//				memcpy(this->adfGeoTransform, info.adfGeoTransform, sizeof(double) * 6);
//				this->numbands = info.numbands;
//				this->projection = info.projection;
//				this->pszFormat = info.pszFormat;
//				this->slice = info.slice;
//				this->nodata = info.nodata;
//				if (this->nodata == 0)
//						this->nodata = -9999;
//		}
//		void setData(const float* data)
//		{
//				this->m_data.resize(slice);
//				memcpy(&m_data[0], data, slice * sizeof(float));
//		}
//		void buildChild();
//		QuadTree* m_child;
//		std::vector<float> m_data;
//		Sampler m_sampler;
//		float sample(float block[4])
//		{
//				float min = block[0];
//				for (int i = 0; i < 3; i++)
//				{
//						const float& val = block[i];
//						if (min == nodata)
//								min = val;
//						else if(val != nodata && min > val)
//								min = val;
//				}
//				return min;
//		}
//
//		float getValue(const double& x, const double& y)
//		{
//				int col = getCol(x);
//				int row = getRow(y);
//				return m_data[col + row * ncols];
//		}
//
//		float getValue(const int& col, const int& row)
//		{
//				return m_data[col + row * ncols];
//		}
//
//		int getCol(const double& x)
//		{
//				if (x < bound.MinX || x > bound.MaxX)
//						return -1;
//				int col = (int)((x - bound.MinX) / adfGeoTransform[1]);
//				if (col < 0)
//						col = 0;
//				else if (col > ncols - 1)
//						col = ncols - 1;
//				return col;
//		}
//
//		int getRow(const double& y)
//		{
//				if (y < bound.MinY || y > bound.MaxY)
//						return -1;
//				int row = (int)((y - bound.MaxY) / adfGeoTransform[5]);
//				if (row < 0)
//						row = 0;
//				else if (row > nrows - 1)
//						row = nrows - 1;
//				return row;
//		}
//
//		//double calHorizonAngle(const Ray& ray, int& row, int& col)
//		//{
//		//		double horizonAngle = 0;
//		//		if (m_data.size() == 0)
//		//				return horizonAngle;
//		//		osg::Vec3d pos = ray.orig;
//		//		int col = getCol(pos.x());
//		//		int row = getRow(pos.y());
//		//		if (col < 0 || row < 0)
//		//				return horizonAngle;
//		//		int index = col + row * ncols;
//		//		double baseHeight = (double)m_data[index];
//		//		double cellSize = adfGeoTransform[1];
//		//		double delta = cellSize * 4;
//		//		pos = ray.orig + ray.dir * delta;
//		//		double run = delta;
//		//		int lastCol = col;
//		//		int lastRow = row;
//		//		while (true)
//		//		{
//		//				col = getCol(pos.x());
//		//				row = getRow(pos.y());
//		//				if (col < 0 || row < 0)
//		//						return horizonAngle;
//		//				if (col != lastCol || row != lastRow)
//		//				{
//		//						index = col + row * ncols;
//		//						const float elev = m_data[index];
//		//						double rise = elev - baseHeight;
//		//						double angle = atan2(rise, run) * 180 / M_PI;
//		//						horizonAngle = max(angle, horizonAngle);
//		//				}
//		//				lastCol = col;
//		//				lastRow = row;
//		//				pos = pos + ray.dir * delta;
//		//				run += delta;
//
//		//		}
//		//		return horizonAngle;
//		//}
//
//		//std::vector<std::pair<unsigned short, unsigned short>> computeHorizonAngleMapIndices(const osg::Vec2d& dir)
//		//{
//		//		std::vector<std::pair<unsigned short, unsigned short>> output;
//		//		for (int row = 0; row < nrows; row++)
//		//		{
//		//				double y = bound.MaxY - adfGeoTransform[1] * row - adfGeoTransform[1] * 0.5;
//		//				for (int col = 0; col < ncols; col++)
//		//				{
//		//						double x = bound.MinX + adfGeoTransform[1] * col + adfGeoTransform[1] * 0.5;
//		//						int index = col + row * ncols;
//		//						float elev = m_data[index];
//		//						unsigned short encodedHorizon = 100;
//		//						int rowIdx = row;
//		//						int colIdx = col;
//		//						if (elev > -500 && elev < 10000)
//		//						{
//		//								Ray groundRay(osg::Vec2d(x, y), dir);
//		//								double horizonAngle = calHorizonAngle(groundRay, rowIdx, colIdx);
//		//						}
//		//						output.push_back(std::pair<unsigned short, unsigned short>(rowIdx, colIdx));
//		//				}
//		//		}
//		//		return output;
//		//}
//
//		//std::vector<std::pair<unsigned short, unsigned short>> computeHorizonAngleMapIndicesFromChild(QuadTree* child, std::vector<std::pair<unsigned short, unsigned short>> childMaxIndices)
//		//{
//		//		std::vector<std::pair<unsigned short, unsigned short>> output;
//		//		return output;
//		//		unsigned short blockRowIndices[4];
//		//		unsigned short blockColIndices[4];
//		//		float block[4];
//
//		//		for (int row = 0; row < child->nrows; row++)
//		//		{
//		//				double y = bound.MaxY - adfGeoTransform[1] * row - adfGeoTransform[1] * 0.5;
//		//				for (int col = 0; col < child->ncols; col++)
//		//				{
//		//						double x = bound.MinX + adfGeoTransform[1] * col + adfGeoTransform[1] * 0.5;
//		//						int index = col + row * ncols;
//		//						float elev = m_data[index];
//		//						unsigned short encodedHorizon = 100;
//		//						int rowIdx = row;
//		//						int colIdx = col;
//		//						blockRowIndices[0] = col * 2;
//		//						blockRowIndices[1] = col * 2 + 1;
//		//						blockRowIndices[2] = col * 2;
//		//						blockRowIndices[3] = col * 2 + 1;
//
//		//						blockColIndices[0] = row * 2;
//		//						blockColIndices[1] = row * 2;
//		//						blockColIndices[2] = row * 2 + 1;
//		//						blockColIndices[3] = row * 2 + 1;
//
//		//						block[0] = getValue(col * 2, row * 2);
//		//						block[1] = getValue(col * 2 + 1, row * 2);
//		//						block[2] = getValue(col * 2, row * 2 + 1);
//		//						block[3] = getValue(col * 2 + 1, row * 2 + 1);
//
//		//						double maxHorizon = 0;
//		//						int cellIndex = 0;
//		//						for (int cell = 0; cell < 4; cell++)
//		//						{
//		//								if (blockRowIndices[cell] > nrows - 1 || blockColIndices[cell] > ncols - 1)
//		//										continue;
//		//								double baseHeight = (double)m_data[index];
//		//								//double delta = cellSize * 4;
//		//							//	pos = ray.orig + ray.dir * delta;
//		//								//double run = delta;
//		//						}
//		//				}
//		//		}
//		//}
//};

class BBWrapper : public osg::BoundingBoxd
{
public:
		BBWrapper() :
				osg::BoundingBoxd()
		{

		}

		BBWrapper(const osg::BoundingBoxd& bb) :
				osg::BoundingBoxd(bb)
		{
				bounds[0] = _min;
				bounds[1] = _max;
		}

		BBWrapper(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) :
				osg::BoundingBoxd(xmin, ymin, zmin, xmax, ymax, zmax)
		{
				bounds[0] = _min;
				bounds[1] = _max;
		}

		BBWrapper(osg::Vec3d min, osg::Vec3d max) :
				osg::BoundingBoxd(min, max)
		{
				bounds[0] = _min;
				bounds[1] = _max;
		}

		double xsize() { return xMax() - xMin(); }
		double ysize() { return yMax() - yMin(); }
		double zsize() { return zMax() - zMin(); }
		double xhalfsize() { return xsize() * 0.5; }
		double yhalfsize() { return ysize() * 0.5; }
		double zhalfsize() { return zsize() * 0.5; }
		osg::Vec3d size() { return _max - _min; }
		osg::Vec3d halfsize() { return size() * 0.5; }
		osg::Vec3d bounds[2];
};
//https://www.scratchapixel.com/code.php?id=10&origin=/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes&src=1
class Ray
{
public:


		Ray(const osg::Vec3d &orig, const osg::Vec3d &dir)
		{
				this->orig = orig; 
				this->dir = dir;
				this->dir.normalize();
				invdir = osg::Vec3d(1.0 / this->dir.x(), 1.0 / this->dir.y(), 1.0 / this->dir.z());
				sign[0] = (invdir.x() < 0);
				sign[1] = (invdir.y() < 0);
				sign[2] = (invdir.z() < 0);
		}


		Ray(const osg::Vec2d &orig, const osg::Vec2d &dir)
				:Ray(osg::Vec3d(orig.x(), orig.y(), 0.0), osg::Vec3d(dir.x(), dir.y(), 0.0))
		{
		}

		//bool intersect(const AABBox &bb, double &t) const
		//{
		//		double tmin, tmax, tymin, tymax, tzmin, tzmax;

		//		tmin = (bb.bounds[sign[0]].x() - orig.x()) * invdir.x();
		//		tmax = (bb.bounds[1 - sign[0]].x() - orig.x()) * invdir.x();
		//		tymin = (bb.bounds[sign[1]].y() - orig.y()) * invdir.y();
		//		tymax = (bb.bounds[1 - sign[1]].y() - orig.y()) * invdir.y();

		//		if ((tmin > tymax) || (tymin > tmax))
		//				return false;

		//		if (tymin > tmin)
		//				tmin = tymin;
		//		if (tymax < tmax)
		//				tmax = tymax;

		//		tzmin = (bb.bounds[sign[2]].z() - orig.z()) * invdir.z();
		//		tzmax = (bb.bounds[1 - sign[2]].z() - orig.z()) * invdir.z();

		//		if ((tmin > tzmax) || (tzmin > tmax))
		//				return false;

		//		if (tzmin > tmin)
		//				tmin = tzmin;
		//		if (tzmax < tmax)
		//				tmax = tzmax;

		//		t = tmin;

		//		if (t < 0) {
		//				t = tmax;
		//				if (t < 0) return false;
		//		}

		//		return true;
		//}

		osg::Vec3d orig, dir; // ray orig and dir 
		osg::Vec3d invdir;
		int sign[3];
};

class AABBox
{
public:
		AABBox(const osg::Vec3d &b0, const osg::Vec3d &b1) { bounds[0] = b0, bounds[1] = b1; }

		AABBox(const osg::Vec2d &b0, const osg::Vec2d &b1) { bounds[0] = osg::Vec3d(b0.x(), b0.y(), 0.0), bounds[1] = osg::Vec3d(b1.x(), b1.y(), 0.0); }

		bool intersect(Ray r) const
		{
				double t;
				return	intersect(r, t);
		}

		bool intersect(const Ray &r, double &t) const
		{
				double tmin, tmax, tymin, tymax, tzmin, tzmax;

				tmin = (bounds[r.sign[0]].x() - r.orig.x()) * r.invdir.x();
				tmax = (bounds[1 - r.sign[0]].x() - r.orig.x()) * r.invdir.x();
				tymin = (bounds[r.sign[1]].y() - r.orig.y()) * r.invdir.y();
				tymax = (bounds[1 - r.sign[1]].y() - r.orig.y()) * r.invdir.y();

				if ((tmin > tymax) || (tymin > tmax))
						return false;

				if (tymin > tmin)
						tmin = tymin;
				if (tymax < tmax)
						tmax = tymax;

				tzmin = (bounds[r.sign[2]].z() - r.orig.z()) * r.invdir.z();
				tzmax = (bounds[1 - r.sign[2]].z() - r.orig.z()) * r.invdir.z();

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
		osg::Vec3d bounds[2];
};

class ModelLoader
{
public:
	ModelLoader();
	~ModelLoader();
	static osg::Node* Load3DTiles(std::string indir);
	static osg::Node* Load3DTiles(std::string indir, osg::BoundingBox mask, bool intersects);
	static osg::Node* Load3DTiles(std::string indir, std::vector<std::string> maskTiles, bool include);
	static void CopyLeafTiles(std::string indir, std::string outdir, osg::BoundingBox bb);
	static void CopyLeafTiles(std::string indir, std::vector<std::string> tilenames, std::string outfile);
	static void Test();
	static void TileBoundary2Shapefile(std::string indir, std::string outfile);
	static osg::BoundingBoxd CalBound(std::string indir, const std::vector<std::string>& tiles);
	static void ConvertAspectQGIS2Grass(std::string infile, std::string outfile);
};


class ShadowCaster
{
public:
		ShadowCaster() 
		{
				m_slopeData = nullptr;
				m_aspectData = nullptr;
		};
		~ShadowCaster();
		void Initialize(std::vector<std::string> rasterFiles, std::vector<osg::BoundingBoxd> bounds, std::string slopeFile, std::string aspectFile);
		void CastShadow(osg::BoundingBoxd bb, std::string outfile, osg::Vec3d sundir);
		void ComputeViewshed(double px, double py);
		static bool Intersects(const Ray& r, const BBWrapper& aabb, double& t);
		std::tuple<int, int, int, int, double, BBWrapper> GetCellIndex(double x, double y);
		bool Intersects(const Ray& ray, double& t, double& slope, double& aspect);
		static unsigned short* ComputeHorizonAngleMap(void* pDS, float*& data, const double& nodata, const osg::Vec2d& dir);
		static void ComputeHorizonAngleMaps(std::string infile, std::string outfile, double startAngle, double endAngle, double stepAngle);
		//static void ComputeHorizonAngleMapsFast(std::string infile, std::string outfile, double startAngle, double endAngle, double stepAngle);
		std::vector<void*> m_rasters;
		std::vector<osg::BoundingBoxd> m_bounds;
		float* m_slopeData;
		float* m_aspectData;
private:
};