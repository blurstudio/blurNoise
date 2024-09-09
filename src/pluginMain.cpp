#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include "blurNoise.h"
#include "version.h"

MStatus initializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj, "Blur Studio", VERSION_STRING, "Any");
	result = plugin.registerNode(
        DEFORMER_NAME,
        blurNoise::id,
        blurNoise::creator,
        blurNoise::initialize,
        MPxNode::kDeformerNode
    );
    MGlobal::executeCommand("makePaintable -attrType \"multiFloat\" -sm \"deformer\" \"" DEFORMER_NAME "\" \"weights\";");
	return result;
}

MStatus uninitializePlugin(MObject obj) {
	MStatus result;
	MFnPlugin plugin(obj);
	result = plugin.deregisterNode(blurNoise::id);
	return result;
}

