#pragma once

#include <maya/MTypeId.h> 

#include <maya/MMatrix.h>
#include <maya/MDoubleArray.h>

#include <maya/MPxDeformerNode.h> 
#include <maya/MItGeometry.h>
#include <maya/MFnNurbsSurface.h>

#include <FastNoise/FastNoise.h>

#define DEFORMER_NAME "blurNoise"

class blurNoise : public MPxDeformerNode
{
public:
	blurNoise();
	~blurNoise() override;

    static void*   creator();
    static MStatus initialize();

    // Rich crontrol over the amplitude 
    static MObject aAmp; // float
    static MObject aAmpOffset; // float
    static MObject aAmpClampHigh; // float
    static MObject aAmpClampLow; // float
    static MObject aAmpUseClampHigh; // bool
    static MObject aAmpUseClampLow; // bool

    // The Time offset for the noise
    // I'm not providing a seed because its easier
    // to just let the artist shift the time value
    // to get different noise
    static MObject aTime; // time

    // Transform channels for noise control
    // Make sure these connect directly to a
    // transform node so the artist can just
    // hook that up
    static MObject aTranslate;
        static MObject aTranslateX;
        static MObject aTranslateY;
        static MObject aTranslateZ;
    static MObject aRotate;
        static MObject aRotateX;
        static MObject aRotateY;
        static MObject aRotateZ;
    static MObject aScale;
        static MObject aScaleX;
        static MObject aScaleY;
        static MObject aScaleZ;
    static MObject aShear;
        static MObject aShearXY;
        static MObject aShearXZ;
        static MObject aShearYZ;
    static MObject aRotateOrder;

    // The number of octaves
    // Then the difference between octaves scale and range is N**(octNum)
    // Where N is one of the bases, and octNum is the current octave index
    static MObject aOctaveCount; // int
    static MObject aOctaveScaleBase; // float > 1
    static MObject aOctaveRangeBase; // float > 1

    static const MTypeId id;
private:
    virtual MStatus deform(
		MDataBlock& block,
		MItGeometry& vertIter,
		const MMatrix& mat,
		unsigned int multiIndex
	);

    void getSurfaceCVParams(
        const MFnNurbsSurface &fnSurf,
        MDoubleArray &uParams,
        MDoubleArray &vParams
    ) const;

    FastNoise::SmartNode<FastNoise::Simplex> fnSimplex;
    FastNoise::SmartNode<FastNoise::FractalFBm> fnFractal;
};
