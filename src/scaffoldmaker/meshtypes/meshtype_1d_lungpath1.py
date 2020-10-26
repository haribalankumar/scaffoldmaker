"""
Generates a 1-D lung path mesh.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds, remapEftNodeValueLabelWithNodes
from scaffoldmaker.utils import interpolation as interp

class MeshType_1d_lungpath1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '1D lung Path 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Coordinate dimensions' : 3,
            'Length' : 1.0,
            'Number of nodes per curve' : 2,
            'Species' : 'Mouse'
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'Length',
            'Number of nodes per curve',
            'Species'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Coordinate dimensions'] < 1) :
            options['Coordinate dimensions'] = 1
        elif (options['Coordinate dimensions'] > 3) :
            options['Coordinate dimensions'] = 3
        if (options['Length'] < 0) :
            options['Length'] = 1
        if (options['Number of nodes per curve'] < 2) :
            options['Number of nodes per curve'] = 2
        if (options['Number of nodes per curve'] > 2) :
            options['Number of nodes per curve'] = 2

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        coordinateDimensions = options['Coordinate dimensions']
        length = options['Length']
        elementsCountpercurve = options['Number of nodes per curve']
        speciestype = options['Species']
        print('creating 1d path for species: ', speciestype)
        # centralPath = options['Control curves']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm, components_count=coordinateDimensions)
        cache = fm.createFieldcache()

        #################
        # Create nodes
        #################
        if speciestype == 'Mouse':
            #Left lung (diaph medial, diaph lateral, accessory edge, posterior lateral, anterior edge, posterior edge)
            elementsCountLeft = (elementsCountpercurve+1)*6 - 9 #6 curves for left lobe

            #Right lung  using Intersection ∩ operator
            # posteriorEdge ∩ RML
            # Apical (aka RUL)
            # AnteriorEdge ∩ RLL
            # (basemedial ∩ RML)
            # (basemedial ∩ RLL)
            # (baselateral ∩ RML)
            # (baselateral ∩ RLL)
            # (ObliqueFissureLateral ∩ RML)
            # (ObliqueFissureMedial ∩ RML)
            # (ObliqueFissureMedial ∩ RUL)
            # (ObliqueFissureLateral ∩ RUL)
            # (HorizontalFissureLateral)
            # (HorizontalFissureMedial)
            elementsCountRight = (elementsCountpercurve-1)*13 - 4 #13 curves for right lung
            elementsCount = elementsCountLeft + elementsCountRight

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        # nodeIdentifier = 1
        # zero = [ 0.0, 0.0, 0.0 ]
        # for n in range(elementsCount + 1):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        #     nodeIdentifier = nodeIdentifier + 1


        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ length/elementsCount, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0, 0.0 ]
        d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCount):
            x[0] = length*n/elementsCount
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

#print('dofs',Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3)
#OUTPUT: 1,2,3,5

        # nodeIdentifier = 1
        # for n in range(elementsCount + 1):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3)
        #     nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################
        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
        # setEftScaleFactorIds(eft1, [1], [])
        elementtemplateX = mesh.createElementtemplate()
        elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result2 = elementtemplateX.defineField(coordinates, -1, eft1)

        elementIdentifier = 1

        if speciestype == 'Mouse':
            # eft1 = eftfactory.createEftNoCrossDerivatives()
            ##############################
            ##LEFT LUNG
            ##############################
            if (elementsCountpercurve == 2):
                # anterior and posterior
                #------------------------
                for e in range(2):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    d2Map = (0, 1, 0)
                    for i in range(2):
                        remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                                           (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                                           d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + 1, e + 2])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                for e in range(2):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    for i in range(2):
                        if (e > 0):
                            d2Map = (0, 1, 0)
                            remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3), d2Map))
                        else:
                            d2Map = (0, (-1)**i, 0)
                            remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3), d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + 4, e + 3])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                # Diaphragm medial
                #-----------------------
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
                result = elementtemplate.defineField(coordinates, -1, eft)

                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [1, 6])
                elementIdentifier = elementIdentifier + 1

                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [6, 5])
                elementIdentifier = elementIdentifier + 1

                ##Diaphragm lateral
                #--------------------------
                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                d2Map = (0, 0, -1)
                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, [1, 7])
                element.setScaleFactors(eft1, [-1.0])
                elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, [7, 5])
                element.setScaleFactors(eft1, [-1.0])
                elementIdentifier = elementIdentifier + 1

                # accessory edge
                #------------------
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
                elementtemplate.defineField(coordinates, -1, eft)

                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [2, 8])
                elementIdentifier = elementIdentifier + 1

                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, [8, 4])
                elementIdentifier = elementIdentifier + 1

                ###lateral smooth edge
                #-----------------------
                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                d2Map = (-1, 0, 0)
                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                       derivativeSignsToExpressionTerms(
                                           (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                                           d2Map))
                elementtemplateX = mesh.createElementtemplate()
                elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_LINE)
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, [2, 9])
                element.setScaleFactors(eft1, [-1.0])
                elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                d2Map = (0, 0, 1)
                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                       derivativeSignsToExpressionTerms(
                                           (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                                           d2Map))
                elementtemplateX = mesh.createElementtemplate()
                elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_LINE)
                elementtemplateX.defineField(coordinates, -1, eft1)
                elementtemplate = elementtemplateX
                element = mesh.createElement(elementIdentifier, elementtemplate)
                result2 = element.setNodesByIdentifier(eft1, [9, 4])
                element.setScaleFactors(eft1, [-1.0])
                elementIdentifier = elementIdentifier + 1

                leftLungendNode = 9
                ##############################
                ##RIGHT LUNGS
                ##############################
                # anterior and posterior
                #=========================####
                for e in range(elementsCountpercurve-1):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    d2Map = (0, 1, 0)
                    for i in range(2):
                        remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + leftLungendNode + 1, e + leftLungendNode+2])
                    elementIdentifier = elementIdentifier + 1

                prevedge_endnode = elementsCountpercurve + leftLungendNode
                for e in range(elementsCountpercurve-1):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    d2Map = (0, 1, 0)
                    for i in range(2):
                        remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + prevedge_endnode, e + prevedge_endnode+1])
                    elementIdentifier = elementIdentifier + 1

                prevedge_endnode = elementsCountpercurve +1 + leftLungendNode
                for e in range(elementsCountpercurve-1):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    for i in range(2):
                        d2Map = (0, (-1)*i, 0)
                        remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + prevedge_endnode+1, e + prevedge_endnode])
                    elementIdentifier = elementIdentifier + 1

                prevedge_endnode = elementsCountpercurve + 2 + leftLungendNode
                for e in range(elementsCountpercurve-1):
                    eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                    setEftScaleFactorIds(eft1, [1], [])
                    d2Map = (0, 1, 0)
                    for i in range(2):
                        remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft1, [e + prevedge_endnode+1, e + prevedge_endnode])
                    elementIdentifier = elementIdentifier + 1


                # Diaphragm medial
                #--------------------
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
                result = elementtemplate.defineField(coordinates, -1, eft)

                for e in range(elementsCountpercurve-1):
                    nstart = leftLungendNode+1 if(e==0) else (e+leftLungendNode+2*elementsCountpercurve+1)
                    nend = (e+leftLungendNode+3*elementsCountpercurve)
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, [nstart, nend])
                    elementIdentifier = elementIdentifier + 1

                for e in range(elementsCountpercurve-1):
                    nstart = (e+leftLungendNode+3*elementsCountpercurve)
                    nend = (leftLungendNode+2*elementsCountpercurve+1) if(e==elementsCountpercurve-2) else (e+leftLungendNode+4*elementsCountpercurve)
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, [nstart, nend])
                    elementIdentifier = elementIdentifier + 1

                prevedge_endnode = leftLungendNode+5*elementsCountpercurve-4
                ##Diaphragm lateral
                #-------------------
                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = leftLungendNode+1 if(e==0) else (prevedge_endnode+e)
                    nend = prevedge_endnode+e+1
                    d2Map = (-1, 0, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                                d2Map))
                    d2Map = (1, 0, 0)
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                                d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = (e+leftLungendNode+5*elementsCountpercurve-3)
                    nend = (leftLungendNode+3*elementsCountpercurve-1) if(e==elementsCountpercurve-2) else (e+leftLungendNode+5*elementsCountpercurve)
                    d2Map = (1, 0, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    d2Map = (0, 0, -1)
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                # Oblique lateral edge
                #----------------------
                for e in range(elementsCountpercurve-1):
                    nstart = (e + leftLungendNode + 5 * elementsCountpercurve - 4)
                    nend = (leftLungendNode + 6 * elementsCountpercurve - 4)
                    d2Map = (0, 1, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = (e+leftLungendNode+6*elementsCountpercurve-4)
                    nend = (leftLungendNode+elementsCountpercurve) if(e==elementsCountpercurve-2) else (e+leftLungendNode+5*elementsCountpercurve)
                    d2Map = (0, 1, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    d2Map = (-1, 0, 0)
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                # Oblique laterl edge
                #----------------------
                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = (e + leftLungendNode + 6 * elementsCountpercurve - 5)
                    nend = (e + leftLungendNode + 7 * elementsCountpercurve - 5)
                    d2Map = (0, 1, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = (e + leftLungendNode + 7 * elementsCountpercurve - 5)
                    nend = (leftLungendNode+elementsCountpercurve) \
                        if(e==elementsCountpercurve-2) \
                        else (e+leftLungendNode+5*elementsCountpercurve)
                    d2Map = (0, 1, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    d2Map = (0, 0, -1)
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

                ### Medial Horizontal fissure edge
                #---------------------------------
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
                elementtemplate.defineField(coordinates, -1, eft)
                for e in range(elementsCountpercurve - 1):
                    nstart = (e+leftLungendNode+6*elementsCountpercurve-4)
                    nend = (e + leftLungendNode + 2 * elementsCountpercurve)
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, [nstart,nend])
                    elementIdentifier = elementIdentifier + 1

                eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
                setEftScaleFactorIds(eft1, [1], [])
                for e in range(elementsCountpercurve-1):
                    nstart = (e + leftLungendNode + 7 * elementsCountpercurve - 5)
                    nend = (leftLungendNode+2*elementsCountpercurve) \
                        if(e==elementsCountpercurve-2) else (e+leftLungendNode+5*elementsCountpercurve)
                    d2Map = (1, 0, 0)
                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                        (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                        d2Map))
                    d2Map = (0, 0, -1)
                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
                            (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
                            d2Map))
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    elementtemplate = elementtemplateX
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, [nstart,nend])
                    element.setScaleFactors(eft1, [-1.0])
                    elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []

def extractxyzPathParametersFromRegion(region):
    '''
    Returns parameters of all nodes in region in identifier order.
    Assumes nodes in region have field coordinates (1 to 3 components).
    Currently limited to nodes with exactly value, d_ds1, d_ds2, d2_ds3,
    same as path 1 scaffold.
    :return: cx, cd1, cd2, cd3 (all padded with zeroes to 3 components)
    '''
    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName('coordinates').castFiniteElement()
    componentsCount = coordinates.getNumberOfComponents()
    assert componentsCount in [1, 2, 3], 'extractxyzPathParametersFromRegion.  Invalid coordinates number of components'
    cache = fm.createFieldcache()
    cx = []
    cd1 = []
    cd2 = []
    cd3 = []

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodeIter = nodes.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        cache.setNode(node)
        result, x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, componentsCount)
        result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, componentsCount)
        result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, componentsCount)
        result, d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, componentsCount)
        for c in range(componentsCount, 3):
            x.append(0.0)
            d1.append(0.0)
            d2.append(0.0)
            d3.append(0.0)
        cx.append(x)
        cd1.append(d1)
        cd2.append(d2)
        cd3.append(d3)
        node = nodeIter.next()
    return cx, cd1, cd2, cd3


def derivativeSignsToExpressionTerms(valueLabels, signs):
    '''
    Return remap expression terms for summing derivative[i]*sign[i]
    :param valueLabels: List of node value labels to possibly include.
    :param signs: List of 1 (no scaling), -1 (scale by scale factor 1) or 0 (no term).
    '''
    expressionTerms = []
    for i in range(len(valueLabels)):
        if signs[i] is 1:
            expressionTerms.append((valueLabels[i], []))
        elif signs[i] is -1:
            expressionTerms.append((valueLabels[i], [1]))
    return expressionTerms






    # elementtemplateStandard = mesh.createElementtemplate()
    # elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    #
    # eftStandard = eftfactory.createEftBasic()
    # elementtemplateStandard.defineField(coordinates, -1, eftStandard)
    #
    # elementtemplate = mesh.createElementtemplate()
    # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    #
    # elementtemplateX = mesh.createElementtemplate()
    # elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    #
    # # result = elementtemplate.defineField(coordinates, -1, eft1)
    #
    #
    # mapDerivatives = True
    # mapEndDerivatives = True
    # mapStartDerivatives = False
     #
    #
    # #######################
    # # create elements
    # ########################
    # elementtemplate = elementtemplateStandard
    # eft1 = eftStandard
    #
    # # if mapDerivatives:
    # #     eft1 = eftfactory.createEftNoCrossDerivatives()
    # #     setEftScaleFactorIds(eft1, [1], [])
    # #
    # #     elementtemplateX.defineField(coordinates, -1, eft1)
    # #     elementtemplate = elementtemplateX
    # # else:
    # #     # eft1 = eft
    # #     elementtemplate = elementtemplateStandard
    #
    #     bni11 = firstNodeIdentifier + 3*elementsCountAround*(elementsCountAlong+1)*(elementsCountThroughWall+1)-1 #72
    #     bni12 = bni11 + (elementsCountAround*elementsCountAlong) #80
    #     bni21 = firstNodeIdentifier + elementsCountAround*(elementsCountAlong+1)*(elementsCountThroughWall+1)-1  #24
    #     bni22 = bni21 + (elementsCountAround*elementsCountAlong) #32
    #     onecircle = elementsCountAround
    #     onecirclem1 = elementsCountAround -1
    #     # nodeIdentifiers = [73, 80, 28, 25, 76, 82, 32, 29]
    #     nodeIdentifiers = [bni11+1, bni12, bni21+onecircle, bni21+1, bni11+onecircle, bni12+2, bni22, bni22-onecirclem1]
    #     eft1 = eftfactory.createEftNoCrossDerivatives()
    #     setEftScaleFactorIds(eft1, [1], [])
    #     if mapDerivatives:
    #         d2Map = (1,1,0)
    #         lns = [2,6]
    #         for n3 in range(2):
    #             ln = [lns[n3]]
    #             remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
    #                                derivativeSignsToExpressionTerms(
    #                                    (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3), d2Map))
    #         elementtemplateX.defineField(coordinates, -1, eft1)
    #         elementtemplate = elementtemplateX
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #     result3 = element.setScaleFactors(eft1, [-1.0])
    #     elementIdentifier = elementIdentifier + 1
    #     if nSegment == 0:
    #         for meshGroup in meshGroups3:
    #             meshGroup.addElement(element)
    #

#
# if (elementsCountpercurve == 4):
#     # anterior and posterior
#     for e in range(4):
#         eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#         setEftScaleFactorIds(eft1, [1], [])
#         d2Map = (0, 1, 0)
#         for i in range(2):
#             remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#                 (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#                 d2Map))
#         elementtemplateX.defineField(coordinates, -1, eft1)
#         elementtemplate = elementtemplateX
#         element = mesh.createElement(elementIdentifier, elementtemplate)
#         element.setNodesByIdentifier(eft1, [e + 1, e + 2])
#         element.setScaleFactors(eft1, [-1.0])
#         elementIdentifier = elementIdentifier + 1
#     for e in range(4):
#         eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#         setEftScaleFactorIds(eft1, [1], [])
#         for i in range(2):
#             if (e > 0):
#                 d2Map = (0, 1, 0)
#                 remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#                     (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#                     d2Map))
#             else:
#                 d2Map = (0, (-1) * i, 0)
#                 remapEftNodeValueLabel(eft1, [i + 1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#                     (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#                     d2Map))
#         elementtemplateX.defineField(coordinates, -1, eft1)
#         elementtemplate = elementtemplateX
#         element = mesh.createElement(elementIdentifier, elementtemplate)
#         element.setNodesByIdentifier(eft1, [e + 6, e + 5])
#         element.setScaleFactors(eft1, [-1.0])
#         elementIdentifier = elementIdentifier + 1
#
#     elementtemplate = mesh.createElementtemplate()
#     elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     elementtemplate.defineField(coordinates, -1, eft)
#
#     # accessory edge
#     # =-============
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [3, 10])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [10, 11])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [11, 12])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [12, 7])
#     elementIdentifier = elementIdentifier + 1
#
#     ###lateral smooth curve
#     # result2 = elementtemplateX.defineField(coordinates, -1, eft1)
#     eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#     setEftScaleFactorIds(eft1, [1], [])
#     d2Map = (0, 0, -1)
#     remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#         (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#         d2Map))
#     elementtemplateX = mesh.createElementtemplate()
#     elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     elementtemplateX.defineField(coordinates, -1, eft1)
#     elementtemplate = elementtemplateX
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     result2 = element.setNodesByIdentifier(eft1, [3, 13])
#     element.setScaleFactors(eft1, [-1.0])
#     elementIdentifier = elementIdentifier + 1
#
#     # element = mesh.createElement(elementIdentifier, elementtemplate)
#     # element.setNodesByIdentifier(eft, [ 3, 13])
#     # elementIdentifier = elementIdentifier + 1
#
#     elementtemplate = mesh.createElementtemplate()
#     elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     result = elementtemplate.defineField(coordinates, -1, eft)
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [13, 14])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [14, 15])
#     elementIdentifier = elementIdentifier + 1
#
#     eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#     setEftScaleFactorIds(eft1, [1], [])
#     d2Map = (0, 0, 1)
#     remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#         (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#         d2Map))
#     elementtemplateX = mesh.createElementtemplate()
#     elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     elementtemplateX.defineField(coordinates, -1, eft1)
#     elementtemplate = elementtemplateX
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     result2 = element.setNodesByIdentifier(eft1, [15, 7])
#     element.setScaleFactors(eft1, [-1.0])
#     elementIdentifier = elementIdentifier + 1
#
#     # Diaphragm medial
#     elementtemplate = mesh.createElementtemplate()
#     elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     result = elementtemplate.defineField(coordinates, -1, eft)
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [1, 16])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [16, 17])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [17, 18])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [18, 9])
#     elementIdentifier = elementIdentifier + 1
#
#     ##Diaphragm lateral
#     eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#     setEftScaleFactorIds(eft1, [1], [])
#     d2Map = (0, 0, -1)
#     remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#         (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#         d2Map))
#     elementtemplateX.defineField(coordinates, -1, eft1)
#     elementtemplate = elementtemplateX
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     result2 = element.setNodesByIdentifier(eft1, [1, 19])
#     element.setScaleFactors(eft1, [-1.0])
#     elementIdentifier = elementIdentifier + 1
#
#     elementtemplate = mesh.createElementtemplate()
#     elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
#     result = elementtemplate.defineField(coordinates, -1, eft)
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [19, 20])
#     elementIdentifier = elementIdentifier + 1
#
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     element.setNodesByIdentifier(eft, [20, 21])
#     elementIdentifier = elementIdentifier + 1
#
#     eft1 = mesh.createElementfieldtemplate(cubicHermiteBasis)
#     setEftScaleFactorIds(eft1, [1], [])
#     d2Map = (0, 0, 1)
#     remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, derivativeSignsToExpressionTerms(
#         (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
#         d2Map))
#     elementtemplateX.defineField(coordinates, -1, eft1)
#     elementtemplate = elementtemplateX
#     element = mesh.createElement(elementIdentifier, elementtemplate)
#     result2 = element.setNodesByIdentifier(eft1, [21, 9])
#     element.setScaleFactors(eft1, [-1.0])
#     elementIdentifier = elementIdentifier + 1
