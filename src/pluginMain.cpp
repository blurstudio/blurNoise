#include <maya/MFnPlugin.h>
#include <maya/MTypeId.h> 
#include "blurNoise.h"

// standard initialization procedures
//

MStatus initializePlugin( MObject obj )
{
    MStatus result;

    MFnPlugin plugin( obj, PLUGIN_COMPANY, "3.0", "Any");
    result = plugin.registerNode(
        "blurNoise" ,
        blurNoise::id ,
        &blurNoise::creator ,
        &blurNoise::initialize ,
        MPxNode::kDeformerNode
        );

    return result;
}

MStatus uninitializePlugin( MObject obj )
{
    MStatus result;

    MFnPlugin plugin( obj );
    result = plugin.deregisterNode( blurNoise::id );

    return result;
}
