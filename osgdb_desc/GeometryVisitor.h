/* -*-c++-*- OpenSceneGraph - Copyright (C) 1998-2006 Robert Osfield
 *
 * This library is open source and may be redistributed and/or modified under
 * the terms of the OpenSceneGraph Public License (OSGPL) version 0.0 or
 * (at your option) any later version.  The full license is in LICENSE file
 * included with this distribution, and on the openscenegraph.org website.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * OpenSceneGraph Public License for more details.
*/

#ifndef _GeometryVisitor
#define _GeometryVisitor

#include <osg/NodeVisitor>
#include <osg/Matrix>
#include <osg/Geometry>
#include <osg/Transform>
#include <osg/Texture2D>

#include <osgUtil/Export>

#include <set>
#include <osg/Geode>
#include <osg/PagedLOD>

/** Helper base class for implementing Optimizer techniques.*/
class GeometryVisitor : public osg::NodeVisitor
{
    public:
		GeometryVisitor():osg::NodeVisitor(TRAVERSE_ALL_CHILDREN){setNodeMaskOverride(0xffffffff);}
		virtual void apply(osg::Geode& node);
		osg::Image* DSM() const { return g_pDSM; }
		void DSM(osg::Image* val) { g_pDSM = val; }
	    void Extent(osg::BoundingBox val) { g_mExtent = val; }
		osg::BoundingBox Extent() const { return g_mExtent; }
private:
	osg::Image* g_pDSM;
	osg::BoundingBox g_mExtent;
};

#endif
