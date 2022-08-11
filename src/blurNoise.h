#pragma once

#include <maya/MTypeId.h> 

#include <maya/MMatrix.h>
#include <maya/MDoubleArray.h>

#include <maya/MPxDeformerNode.h> 
#include <maya/MItGeometry.h>
#include <maya/MFnNurbsSurface.h>

#include "OpenSimplex2F.h"

#define DEFORMER_NAME "blurNoise"

class blurNoise : public MPxDeformerNode
{
public:
    static void*   creator();
    static MStatus initialize();

    virtual MStatus deform(
		MDataBlock& block,
		MItGeometry& vertIter,
		const MMatrix& mat,
		unsigned int multiIndex
	);

    void getSurfaceCVParams(const MFnNurbsSurface &fnSurf, MDoubleArray &uParams, MDoubleArray &vParams) const;

    static const MTypeId id;
    OpenSimplexEnv *ose;
    OpenSimplexGradients *osg;

	blurNoise();
	~blurNoise() override;
};
