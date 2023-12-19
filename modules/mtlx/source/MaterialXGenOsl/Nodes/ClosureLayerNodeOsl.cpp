//
// Copyright Contributors to the MaterialX Project
// SPDX-License-Identifier: Apache-2.0
//

#include <MaterialXGenOsl/Nodes/ClosureLayerNodeOsl.h>

#include <MaterialXGenShader/ShaderGenerator.h>
#include <MaterialXGenShader/GenContext.h>
#include <MaterialXGenShader/TypeDesc.h>

MATERIALX_NAMESPACE_BEGIN

ShaderNodeImplPtr ClosureLayerNodeOsl::create()
{
    return std::make_shared<ClosureLayerNodeOsl>();
}

void ClosureLayerNodeOsl::emitFunctionCall(const ShaderNode& _node, GenContext& context, ShaderStage& stage) const
{
    const ShaderGenerator& shadergen = context.getShaderGenerator();

    ShaderNode& node = const_cast<ShaderNode&>(_node);

    ShaderInput* topInput = node.getInput(TOP);
    ShaderInput* baseInput = node.getInput(BASE);
    ShaderOutput* output = node.getOutput();

    // Make sure the layer is fully connected.
    if (!(topInput->getConnection() && baseInput->getConnection()))
    {
        // Just declare the output variable with default value.
        emitOutputVariables(node, context, stage);
        return;
    }

    ShaderNode* top = topInput->getConnection()->getNode();
    ShaderNode* base = baseInput->getConnection()->getNode();

#ifdef MATERIALX_OSL_LEGACY_CLOSURES

    ClosureContext* cct = context.getClosureContext();

    // Evaluate top and base nodes and combine their result
    // according to throughput.
    //
    // TODO: In the BSDF over BSDF case should we emit code
    //       to check the top throughput amount before calling
    //       the base BSDF?

    // Make sure the connections are sibling nodes and not the graph interface.
    if (top->getParent() == node.getParent())
    {
        // If this layer node has closure parameters set,
        // we pass this on to the top component only.
        ScopedSetClosureParams setParams(&node, top, cct);
        shadergen.emitFunctionCall(*top, context, stage);
    }
    if (base->getParent() == node.getParent())
    {
        shadergen.emitFunctionCall(*base, context, stage);
    }

    // Get the result variables.
    const string& topResult = topInput->getConnection()->getVariable();
    const string& baseResult = baseInput->getConnection()->getVariable();

    // Calculate the layering result.
    emitOutputVariables(node, context, stage);
    if (base->getOutput()->getType() == Type::VDF)
    {
        // Combining a surface closure with a volumetric closure is simply done with the add operator in OSL.
        shadergen.emitLine(output->getVariable() + ".response = " + topResult + ".response + " + baseResult, stage);
        // Just pass the throughput along.
        shadergen.emitLine(output->getVariable() + ".throughput = " + topResult + ".throughput", stage);
    }
    else
    {
        shadergen.emitLine(output->getVariable() + ".response = " + topResult + ".response + " + baseResult + ".response * " + topResult + ".throughput", stage);
        shadergen.emitLine(output->getVariable() + ".throughput = " + topResult + ".throughput * " + baseResult + ".throughput", stage);
    }

#else

    // Emit the function call for top and base layer.
    // Make sure the connections are sibling nodes and not the graph interface.
    if (top->getParent() == node.getParent())
    {
        shadergen.emitFunctionCall(*top, context, stage);
    }
    if (base->getParent() == node.getParent())
    {
        shadergen.emitFunctionCall(*base, context, stage);
    }

    // Get the result variables.
    const string& topResult = topInput->getConnection()->getVariable();
    const string& baseResult = baseInput->getConnection()->getVariable();

    // Emit the layer closure call.
    shadergen.emitLineBegin(stage);
    shadergen.emitOutput(output, true, false, context, stage);
    shadergen.emitString(" = layer(" + topResult + ", " + baseResult + ")", stage);
    shadergen.emitLineEnd(stage);

#endif // MATERIALX_OSL_LEGACY_CLOSURES
}

MATERIALX_NAMESPACE_END
