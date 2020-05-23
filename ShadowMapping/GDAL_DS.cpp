//#pragma once
//#include "GDAL_DS.h"
//
//template <class T>
//GDAL_DS<T>::GDAL_DS()
//{
//		dataset = NULL;
//		pszFormat = "GTiff";
//}
//
//template <class T>
//GDAL_DS<T>::~GDAL_DS() {
//		close();
//	}
//
//template <class T>
//	bool GDAL_DS<T>::isInside(double x, double y)
//	{
//		if (x >= bound.MinX && x <= bound.MaxX && y >= bound.MinY && y <= bound.MaxY)
//			return true;
//		return false;
//	}
//
//	template <class T>
//	bool GDAL_DS<T>::open(std::string fname, GDALAccess mode = GA_ReadOnly)
//	{
//		close();
//		filename = fname;
//		dataset = (GDALDataset*)GDALOpen(filename.data(), mode);
//		if (!dataset)
//			return false;
//		dataset->GetGeoTransform(adfGeoTransform);
//		ncols = dataset->GetRasterXSize();
//		nrows = dataset->GetRasterYSize();
//		slice = nrows*ncols;
//		numbands = dataset->GetRasterCount();
//		projection = dataset->GetProjectionRef();
//		setGeoTransform(adfGeoTransform);
//	}
//
//	template <class T>
//	void GDAL_DS<T>::create(std::string fname)
//	{
//		close();
//		filename = fname;
//		char **papszOptions = NULL;
//		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat.data());
//		poDriver->Delete(filename.data());
//		dataset = poDriver->Create(filename.data(), ncols, nrows, numbands, getType(), papszOptions);
//		dataset->SetGeoTransform(adfGeoTransform);
//		dataset->SetProjection(projection.data());
//	}
//
//	template <class T>
//	void GDAL_DS<T>::setGeoTransform(double* transform)
//	{
//		memcpy(adfGeoTransform, transform, sizeof(double) * 6);
//		bound.MinX = adfGeoTransform[0];
//		bound.MaxY = adfGeoTransform[3];
//		bound.MaxX = adfGeoTransform[0] + adfGeoTransform[1] * ncols;
//		bound.MinY = adfGeoTransform[3] + adfGeoTransform[5] * nrows;
//		if (dataset)
//				dataset->SetGeoTransform(transform);
//	}
//
//	template <class T>
//	T* GDAL_DS<T>::readData(int bandidx)
//	{
//		if (!dataset)
//			return NULL;
//		T* data = new T[slice];
//		GDALDataType gdt = getType();
//		dataset->GetRasterBand(bandidx)->RasterIO(GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, gdt, 0, 0);
//		return data;
//	}
//
//	template <class T>
//	void GDAL_DS<T>::readData(int bandidx, T* data)
//	{
//		if (!dataset)
//			return;
//		GDALDataType gdt = getType();
//		dataset->GetRasterBand(bandidx)->RasterIO(GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, gdt, 0, 0);
//	}
//
//	template <class T>
//	void GDAL_DS<T>::writeData(int bandidx, T* data, double nodata)
//	{
//		GDALRasterBand *pBand = dataset->GetRasterBand(bandidx);
//		pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, data, ncols, nrows, getType(), 0, 0);
//		pBand->SetNoDataValue(nodata);
//	}
//
//	template <class T>
//	double GDAL_DS<T>::sum(int bandindex)
//	{
//		T* data = readData(bandindex);
//		T nodata = (T)getNoData(bandindex);
//		T* pdata = data;
//		double sum = 0;
//		for (size_t i = 0; i < slice; i++)
//		{
//			if (*pdata != nodata)
//			{
//				sum += *pdata;
//			}
//			pdata++;
//		}
//		delete[] data;
//		return sum;
//	}
//
//	template <class T>
//	void GDAL_DS<T>::multiply(T factor, int bandindex)
//	{
//		T* data = readData(bandindex);
//		T nodata = (T)getNoData(bandindex);
//		T* pdata = data;
//		double sum = 0;
//		for (size_t i = 0; i < slice; i++)
//		{
//			if (*pdata != nodata)
//			{
//				*pdata = *pdata * factor;
//			}
//			pdata++;
//		}
//		writeData(bandindex, data, nodata);
//		delete[] data;
//	}
//
//	template <class T>
//	double GDAL_DS<T>::getNoData(int bandidx)
//	{
//		return dataset->GetRasterBand(bandidx)->GetNoDataValue();
//	}
//
//	template <class T>
//	void GDAL_DS<T>::setNoData(int bandidx,double nodata)
//	{
//		 dataset->GetRasterBand(bandidx)->SetNoDataValue(nodata);
//	}
//
//	template <class T>
//	void GDAL_DS<T>::close()
//	{
//		if (dataset)
//			GDALClose(dataset);
//		dataset = NULL;
//	}
//
//	template <class T>
//	GDALDataType GDAL_DS<T>::getType()
//	{
//		if (typeid(T) == typeid(unsigned char) || typeid(T) == typeid(char))
//			return GDT_Byte;
//		if (typeid(T) == typeid(unsigned short))
//			return GDT_UInt16;
//		if (typeid(T) == typeid(short))
//			return GDT_Int16;
//		if (typeid(T) == typeid(unsigned int))
//			return GDT_UInt32;
//		if (typeid(T) == typeid(int))
//			return GDT_Int32;
//		if (typeid(T) == typeid(float))
//			return GDT_Float32;
//		if (typeid(T) == typeid(double))
//			return GDT_Float64;
//		return GDT_Unknown;
//	}
//
//	template <class T>
//	void GDAL_DS<T>::setDSInfo(GDAL_DSInfo* info)
//	{
//		nrows = info->nrows;
//		ncols = info->ncols;
//		slice = info->slice;
//		numbands = info->numbands;
//		filename = info->filename;
//		pszFormat = info->pszFormat;
//		projection = info->projection;
//		setGeoTransform(info->adfGeoTransform);
//	}
//
//	template <class T>
//	void GDAL_DS<T>::crop(OGREnvelope destBB, std::string outfile)
//	{
//		int startCol = (int)((destBB.MinX - bound.MinX) / adfGeoTransform[1]);
//		int endCol = (int)((destBB.MaxX - bound.MinX) / adfGeoTransform[1]);
//		int numCols = endCol - startCol + 1;
//
//		int startRow = (int)((destBB.MaxY - bound.MaxY) / adfGeoTransform[5]);
//		int endRow = (int)((destBB.MinY - bound.MaxY) / adfGeoTransform[5]);
//
//		int numRows = endRow - startRow + 1;
//		destBB.MinX = bound.MinX + startCol * adfGeoTransform[1];
//		destBB.MaxY = bound.MaxY + startRow * adfGeoTransform[5];
//		GDAL_DS* newdt = new GDAL_DS<T>();
//
//		memcpy(newdt->adfGeoTransform, adfGeoTransform, sizeof(double) * 6);
//		newdt->adfGeoTransform[0] = destBB.MinX;
//		newdt->adfGeoTransform[3] = destBB.MaxY;
//		newdt->bound = destBB;
//		newdt->numbands = 1;
//		newdt->ncols = numCols;
//		newdt->nrows = numRows;
//		newdt->projection = this->projection;
//		newdt->pszFormat = this->pszFormat;
//		newdt->slice = numCols * numRows;
//		T* srcdata = this->readData(1);
//		T nodata = (T)this->readData(1);
//		T* destdata = new T[numRows*numCols];
//		int idx = 0;
//
//		for (int irow = 0; irow < numRows; irow++)
//		{
//			for (int icol = 0; icol < numCols; icol++)
//			{
//				destdata[idx++] = nodata;
//			}
//		}
//
//		for (int irow = 0; irow < numRows; irow++)
//		{
//			int srcRow = startRow + irow;
//			if (srcRow > nrows - 1 || srcRow < 0)
//				continue;
//			for (int icol = 0; icol < numCols; icol++)
//			{
//				int srcCol = startCol + icol;
//				if (srcCol > ncols - 1 || srcCol < 0)
//					continue;
//				destdata[icol + irow * numCols] = srcdata[srcCol + srcRow * this->ncols];
//			}
//		}
//		newdt->create(outfile);
//		newdt->writeData(1, destdata, nodata);
//
//		delete[] srcdata;
//		delete[] destdata;
//		delete newdt;
//	}
//
//	//<col, row, index>
//	template <class T>
//	std::tuple<int, int, int> GDAL_DS<T>::getCellIndex(double x, double y)
//	{
//			int col = (int)((x - bound.MinX) / adfGeoTransform[1]);
//			int row = (int)((y - bound.MaxY) / adfGeoTransform[5]);
//			int index = row * ncols + col;
//			return std::make_tuple(col, row, index);
//	}
//
//	//template <class T>
//	//void GDAL_DS<T>::downSampleMaximum(std::string outfile, int factor)
//	//{
//	//		double resolX = adfGeoTransform[1] / factor;
//	//		double resolY = adfGeoTransform[5] / factor;
//
//	//		GDAL_DS* newdt = new GDAL_DS<T>();
//	//		newdt->ncols = (int)(ceil((double)ncols / (double)factor));
//	//		newdt->nrows = (int)(ceil((double)nrows / (double)factor));
//	//		newdt->bound.MinX = bound.MinX;
//	//		newdt->bound.MaxY = bound.MinY;
//	//		newdt->bound.MaxX = bound.MinX + newdt->ncols * resolX;
//	//		newdt->bound.MinY = bound.MaxY + newdt->nrows * resolY;
//	//		memcpy(newdt->adfGeoTransform, adfGeoTransform, sizeof(double) * 6);
//	//		newdt->adfGeoTransform[1] = resolX;
//	//		newdt->adfGeoTransform[5] = resolY;
//	//		newdt->setGeoTransform(adfGeoTransform);
//	//		newdt->numbands = 1;
//	//		newdt->projection = this->projection;
//	//		newdt->pszFormat = this->pszFormat;
//	//		newdt->slice = newdt->ncols * newdt->nrows;
//
//	//		int bufXSize = factor * 1000;
//	//		int bufYSize = factor * 1000;
//	//		int bufSize = bufXSize * bufYSize;
//	//		T* srcBlockBuf = new T[bufSize];
//	//		T* destBlockBuf = new T[1000 * 1000];
//	//		T* destdata = new T[newdt->ncols * newdt->nrows];
//	//		int yoffset = 0;
//	//		double nodata = this->getNoData(1);
//	//		while (true)
//	//		{
//	//				int xoffset = 0;
//	//				int blockYSize = MAX(bufYSize, nrows - yoffset);
//	//				while (true)
//	//				{
//	//						int blockXSize = MAX(bufXSize, ncols - xoffset);
//	//						dataset->GetRasterBand(1)->RasterIO(GF_Read, xoffset, yoffset, blockXSize, blockYSize, srcBlockBuf, blockXSize, blockYSize, dsm->getType(), 0, 0);
//	//						float* pDestBlockBuf = destBlockBuf;
//	//						for (int x = 0; x < 1000; x++)
//	//						{
//	//								for (int y = 0; y < 1000; y++)
//	//								{
//	//										float windowMax = -1000000;
//	//										for (int windowx = 0; windowx < factor; windowx++)
//	//										{
//	//												for (int windowy = 0; windowy < factor; windowy++)
//	//												{
//	//														float val = srcBlockBuf[(y * factor + windowy) * bufXSize + (x * factor + windowx)];
//	//														if (windowMax < val)
//	//																windowMax = val;
//	//												}
//	//										}
//	//										if (windowMax <= -1000000)
//	//												windowMax = nodata;
//	//										*pDestBlockBuf = windowMax;
//	//										pDestBlockBuf++;
//	//								}
//	//						}
//
//	//						pDestBlockBuf = destBlockBuf;
//	//						for (int y = yoffset / blockYSize; y < yoffset / blockYSize + 1000; y++)
//	//						{
//	//								if (y >= newdt->nrows)
//	//										break;
//	//								for (int x = xoffset / blockXSize; x < xoffset / blockXSize + 1000; x++)
//	//								{
//	//										if (x >= newdt->ncols)
//	//										{
//	//												pDestBlockBuf++;
//	//												break;
//	//										}
//	//										destdata[y * newdt->ncols + x] = *pDestBlockBuf++;
//	//								}
//	//						}
//
//	//						xoffset += bufXSize;
//	//						if (xoffset >= ncols)
//	//								break;
//	//				}
//	//				yoffset += bufYSize;
//	//				if (yoffset >= nrows)
//	//						break;
//	//		}
//
//	//		newdt->create(outfile);
//	//		newdt->writeData(1, destdata, nodata);
//
//	//		delete[] srcBlockBuf;
//	//		delete[] destBlockBuf;
//	//		delete[] destdata;
//	//		delete newdt;
//	//}