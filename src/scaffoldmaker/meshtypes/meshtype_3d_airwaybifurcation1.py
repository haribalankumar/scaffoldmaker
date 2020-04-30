"""
Generates a 3-D bifurcation mesh with variable numbers of elements around and up.
"""

#######from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
##from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear

from scaffoldmaker.utils import tubebifurcationmesh
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_airwaybifurcation1(Scaffold_base):
    '''
    2-D Airway Bifurcation scaffold.
    '''

    @staticmethod
    def getName():
        return '2D Bifurcation template 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 2,
            'Number of elements around' : 4,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements around',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements along'] < 2):
            options['Number of elements along'] = 2
        if (options['Number of elements around'] > 4):
            options['Number of elements around'] = 4
            if (options['Number of elements around'] < 4):
                options['Number of elements around'] = 4
        elif (options['Number of elements around'] % 2) == 1:
            options['Number of elements around'] += 1

    @staticmethod
    def generateBaseMesh(region, options):
        '''
        Generate the base bicubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        '''

        elementsCountAlong = options['Number of elements along']
        elementsCountAround = options['Number of elements around']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

#        mesh = fm.findMeshByDimension(3)
#        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
#        eft = eftfactory.createEftBasic()

#        elementtemplate = mesh.createElementtemplate()
#        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
#        elementtemplate.defineField(coordinates, -1, eft)

        mesh = fm.findMeshByDimension(2)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        elementIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementAlong = (math.pi/4)/elementsCountAlong


        ########################################################################
        #### centerline details
        ## typically comes from a centerline definition. currently hardcoding it
        parentsegmentLength = 1.1 ##When using centerline this is trach length - junctionheight
        daughter1segmentLength = 0.8 ##When using centerline this is trach length - junctionheight
        daughter2segmentLength = 0.8 ##When using centerline this is trach length - junctionheight

        parentsegmentLength = parentsegmentLength * 0.85
        daughter1segmentLength = daughter1segmentLength * 0.85
        daughter2segmentLength = daughter2segmentLength * 0.85

        elementsCountAlongSegment = 2

        wallThickness = 0.1
        startPhase = 360.0
        ########################################################################
        #
        # x0 = [0,0,0]
        # x1 = [0,0,parentsegmentLength]
        # x2 = [-0.5,0,trachealLength-0.2]
        # x3 = [0.5,0.2,trachealLength-0.2]
        #
        # d10 = [0,0,-1]
        # d11 = [0,0,-1]

        startRadiusparent = 0.4
        endRadiusparent = 0.4
        startRadiusDaugh1 = 0.35
        endRadiusDaugh1 = 0.35
        startRadiusDaugh2 = 0.35
        endRadiusDaugh2 = 0.35

        startRadiusparentDerivative = 0
        endRadiusparentDerivative = 0
        startRadiusDaugh1Derivative = 0.1
        endRadiusDaugh1Derivative = 0.1
        startRadiusDaugh2Derivative = 0.1
        endRadiusDaugh2Derivative = 0.1

        #######################################################################
        segmentCount = 1  # Hardcoded for starters

        # Central path
        cx = [ [ 0.0, 0.0, 0.0 ], [ 0.0, 0.0, parentsegmentLength ] ]
        cd1 = [ [ 0.0, 0.0, parentsegmentLength ], [ 0.0, 0.0, parentsegmentLength ] ]
        cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
        cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]

        print('calling parent central path sampling')
        # Sample central path - PARENT
        sxparent, sd1parent, separent, sxiparent, ssfparent = \
            interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2parent = interp.interpolateSampleCubicHermite(cd2, cd12, separent, sxiparent, ssfparent)[0]
        print('called parent central path sampling', sxparent[0], sxparent[1])

        # Sample central path - DAUGHTER1
        cx = [ [ 0.0, 0.0, 0.0 ], [ daughter1segmentLength, 0.0, 0.0 ] ]
        cd1 = [ [ daughter1segmentLength, 0.0, 0.0 ], [ daughter1segmentLength, 0.0, 0.0 ] ]
        cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
        cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
        sxDaugh1, sd1Daugh1, seDaugh1, sxiDaugh1, ssfDaugh1 = \
            interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2Daugh1 = interp.interpolateSampleCubicHermite(cd2, cd12, seDaugh1, sxiDaugh1, ssfDaugh1)[0]

        # Sample central path - DAUGHTER2
        cx = [ [ 0.0, 0.0, 0.0 ], [ -daughter2segmentLength, 0.0, 0.0 ] ]
        cd1 = [ [ -daughter2segmentLength, 0.0, 0.0 ], [ -daughter2segmentLength, 0.0, 0.0 ] ]
        cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
        cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
        sxDaugh2, sd1Daugh2, seDaugh2, sxiDaugh2, ssfDaugh2 = \
            interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2Daugh2 = interp.interpolateSampleCubicHermite(cd2, cd12, seDaugh2, sxiDaugh2, ssfDaugh2)[0]

        # Find parameter variation along elementsCountAlongSegment
        radiusparentAlongSegment = []
        dRadiusparentAlongSegment = []

        radiusDaugh1AlongSegment = []
        dRadiusDaugh1AlongSegment = []

        radiusDaugh2AlongSegment = []
        dRadiusDaugh2AlongSegment = []

        for n2 in range(elementsCountAlongSegment + 1):
            xi = 1/elementsCountAlongSegment * n2
            print('n2 for xi value = ', n2, xi)
            radius = interp.interpolateCubicHermite([startRadiusparent], [startRadiusparentDerivative],
                                                    [endRadiusparent], [endRadiusparentDerivative], xi)[0]
            radiusparentAlongSegment.append(radius)
            dRadius = interp.interpolateCubicHermiteDerivative([startRadiusparent], [startRadiusparentDerivative],
                                                               [endRadiusparent], [endRadiusparentDerivative], xi)[0]
            dRadiusparentAlongSegment.append(dRadius)

            radius = interp.interpolateCubicHermite([startRadiusDaugh1], [startRadiusDaugh1Derivative],
                                                    [endRadiusDaugh1], [endRadiusDaugh1Derivative], xi)[0]
            radiusDaugh1AlongSegment.append(radius)
            dRadius = interp.interpolateCubicHermiteDerivative([startRadiusparent], [startRadiusparentDerivative],
                                                               [endRadiusparent], [endRadiusparentDerivative], xi)[0]
            dRadiusDaugh1AlongSegment.append(dRadius)

            radius = interp.interpolateCubicHermite([startRadiusDaugh2], [startRadiusDaugh2Derivative],
                                                    [endRadiusDaugh2], [endRadiusDaugh2Derivative], xi)[0]
            radiusDaugh2AlongSegment.append(radius)

            dRadius = interp.interpolateCubicHermiteDerivative([startRadiusparent], [startRadiusparentDerivative],
                                                               [endRadiusparent], [endRadiusparentDerivative], xi)[0]
            dRadiusDaugh2AlongSegment.append(dRadius)

        airwaysegmentTubeMeshInnerPoints = AirwaySegmentTubeMeshInnerPoints(
            region, elementsCountAround, elementsCountAlongSegment,
            parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
            wallThickness, radiusparentAlongSegment, radiusDaugh1AlongSegment,
            radiusDaugh2AlongSegment, dRadiusDaugh2AlongSegment,
            dRadiusDaugh1AlongSegment, dRadiusDaugh2AlongSegment, startPhase)

        # Create inner points
        nSegment = 0

        xParentInner, xDaugh1Inner, xDaugh2Inner, \
        d1ParentInner, d1Daugh1Inner, d1Daugh2Inner, \
        d2ParentInner, d2Daugh1Inner, d2Daugh2Inner, \
        transitElementList, segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,\
        ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ = \
            airwaysegmentTubeMeshInnerPoints.getAirwaySegmentTubeMeshInnerPoints(nSegment)

        contractedWallThicknessList = airwaysegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Warp segment points
        #####################
        # xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList, \
        # d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,\
        # d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,\
        # d3ParentWarpedUnitList, d3Daugh1WarpedList, d3Daugh2WarpedList \
        #     = tubebifurcationmesh.warpAirwaySegmentPoints(
        #     xParentInner, xDaugh1Inner, xDaugh2Inner,
        #     d1ParentInner, d1Daugh1Inner, d1Daugh2Inner,
        #     d2ParentInner, d2Daugh1Inner, d2Daugh2Inner,
        #     segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
        #     parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
        #     sxparent, sxDaugh1, sxDaugh2,
        #     sd1parent, sd1Daugh1, sd1Daugh2,
        #     sd2parent, sd2Daugh1, sd2Daugh2,
        #     elementsCountAround, elementsCountAlongSegment, nSegment,
        #     ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ)


#        # Create coordinates and derivatives
#        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList,
#            d2WarpedList, d3WarpedUnitList, contractedWallThicknessList,
#            elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall, transitElementList)

        # Create nodes and elements
        # nextNodeIdentifier, nextElementIdentifier = \
        #     tubebifurcationmesh.createSurfaceNodesAndElements\
        #         (region,
        #          xParentWarpedList, d1ParentWarpedList, d2ParentWarpedList,
        #          xDaugh1WarpedList, d1Daugh1WarpedList, d2Daugh1WarpedList,
        #          xDaugh2WarpedList, d1Daugh2WarpedList, d2Daugh2WarpedList,
        #          elementsCountAround, elementsCountAlongSegment,
        #         nodeIdentifier, elementIdentifier, useCrossDerivatives)

        print('calling create surf node = ', nodeIdentifier)
        print('NodeIdentifier = ', nodeIdentifier)
        nextNodeIdentifier, nextElementIdentifier = \
            tubebifurcationmesh.createSurfaceNodesAndElements\
                (region,
                 xParentInner, d1ParentInner, d2ParentInner,
                 xDaugh1Inner, d1Daugh1Inner, d2Daugh1Inner,
                 xDaugh2Inner, d1Daugh2Inner, d2Daugh2Inner,
                 elementsCountAround, elementsCountAlongSegment,
                nodeIdentifier, elementIdentifier, useCrossDerivatives)

        fm.endChange()
#        return annotationGroups

class AirwaySegmentTubeMeshInnerPoints:
    """
    Generates inner profile of a cylindrical segment for use by tubemesh.
    """
    def __init__(self, region, elementsCountAround, elementsCountAlongSegment,
                 parentsegmentLength, daugh1segmentLength, daugh2segmentLength,
                 wallThickness,
                 innerRadiusParentSegmentList, innerRadiusDaugh1SegmentList,
                 innerRadiusDaugh2SegmentList, dInnerRadiusParentSegmentList,
                 dInnerRadiusDaugh1SegmentList, dInnerRadiusDaugh2SegmentList, startPhase):

        self._region = region
        self._elementsCountAround = elementsCountAround
        self._elementsCountAlongSegment = elementsCountAlongSegment
        self._parentsegmentLength = parentsegmentLength
        self._daugh1segmentLength = daugh1segmentLength
        self._daugh2segmentLength = daugh2segmentLength
        self._wallThickness = wallThickness
        self._innerRadiusParentSegmentList = innerRadiusParentSegmentList
        self._innerRadiusDaugh1SegmentList = innerRadiusDaugh1SegmentList
        self._innerRadiusDaugh2SegmentList = innerRadiusDaugh2SegmentList
        self._dInnerRadiusParentSegmentList = dInnerRadiusParentSegmentList
        self._dInnerRadiusDaugh1SegmentList = dInnerRadiusDaugh1SegmentList
        self._dInnerRadiusDaugh2SegmentList = dInnerRadiusDaugh2SegmentList

        self._contractedWallThicknessList = []
        self._startPhase = startPhase


    def getAirwaySegmentTubeMeshInnerPoints(self, nSegment):

        # Unpack radius and rate of change of inner radius
        startRadiusParent = self._innerRadiusParentSegmentList[nSegment]
        startRadiusParentDerivative = self._dInnerRadiusParentSegmentList[nSegment]
        endRadiusParent = self._innerRadiusParentSegmentList[nSegment+1]
        endRadiusParentDerivative = self._dInnerRadiusParentSegmentList[nSegment+1]

        startRadiusDaugh1 = self._innerRadiusDaugh1SegmentList[nSegment]
        startRadiusDaugh1Derivative = self._dInnerRadiusDaugh1SegmentList[nSegment]
        endRadiusDaugh1 = self._innerRadiusDaugh1SegmentList[nSegment+1]
        endRadiusDaugh1Derivative = self._dInnerRadiusDaugh1SegmentList[nSegment+1]

        startRadiusDaugh2 = self._innerRadiusDaugh2SegmentList[nSegment]
        startRadiusDaugh2Derivative = self._dInnerRadiusDaugh2SegmentList[nSegment]
        endRadiusDaugh2 = self._innerRadiusDaugh2SegmentList[nSegment+1]
        endRadiusDaugh2Derivative = self._dInnerRadiusDaugh2SegmentList[nSegment+1]

        xParentInner, xDaugh1Inner, xDaugh2Inner, \
        d1ParentInner, d1Daugh1Inner, d1Daugh2Inner, \
        d2ParentInner, d2Daugh1Inner, d2Daugh2Inner, \
        transitElementList, contractedWallThicknessSegment, \
        segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2, \
        ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ\
            = getAirwaySegmentInnerPoints(self._region,
                                          self._elementsCountAround, self._elementsCountAlongSegment,
                                          self._parentsegmentLength, self._daugh1segmentLength,
                                          self._daugh2segmentLength, self._wallThickness,
                                          startRadiusParent, startRadiusDaugh1, startRadiusDaugh2,
                                          startRadiusParentDerivative, startRadiusDaugh1Derivative,
                                          startRadiusDaugh2Derivative,
                                          endRadiusParent, endRadiusDaugh1, endRadiusDaugh2,
                                          endRadiusParentDerivative,endRadiusDaugh1Derivative,
                                          endRadiusDaugh2Derivative,
                                          self._startPhase)

        startIdx = 0 if nSegment == 0 else 1

        # xi = xiSegment[startIdx:self._elementsCountAlongSegment + 1]
        # self._xiList += xi

        contractedWallThickness = contractedWallThicknessSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._contractedWallThicknessList += contractedWallThickness

        return xParentInner, xDaugh1Inner, xDaugh2Inner, \
        d1ParentInner, d1Daugh1Inner, d1Daugh2Inner, \
        d2ParentInner, d2Daugh1Inner, d2Daugh2Inner, transitElementList, \
               segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2, \
               ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ

    def getContractedWallThicknessList(self):
        return self._contractedWallThicknessList


def getAirwaySegmentInnerPoints(region, elementsCountAround, elementsCountAlongSegment,
                                parentsegmentLength, daugh1segmentLength, daugh2segmentLength,
                                wallThickness, startRadiusParent, startRadiusDaugh1, startRadiusDaugh2,
                                startRadiusParentDerivative, startRadiusDaugh1Derivative, startRadiusDaugh2Derivative,
                                endRadiusParent, endRadiusDaugh1, endRadiusDaugh2,
                                endRadiusParentDerivative, endRadiusDaugh1Derivative, endRadiusDaugh2Derivative,
                                     startPhase):
    """
    Generates a 3-D cylindrical segment mesh with variable numbers of elements
    around, along the central path, and through wall.
    :param elementsCountAround: Number of elements around.
    :param elementsCountAlongSegment: Number of elements along cylindrical segment.
    :param segmentLength: Length of a cylindrical segment.
    :param wallThickness: Thickness of wall.
    :param startRadius: Inner radius at proximal end.
    :param startRadiusDerivative: Rate of change of inner radius at proximal end.
    :param endRadius: Inner radius at distal end.
    :param endRadiusDerivative: Rate of change of inner radius at distal end.
    :param startPhase: Phase at start.
    :return coordinates, derivatives on inner surface of a cylindrical segment.
    :return transitElementList: stores true if element around is an element that
    transits between a big and small element.
    :return faceMidPointsZ: z-coordinate of midpoints for each element group
    along segment.
    """

    transitElementList = [0] * elementsCountAround

    contractedWallThicknessList = []

    # create nodes - PARENT
    #######################
    segmentAxisParent = [0.0, 0.0, 1.0]

    xParentFinal = []
    d1ParentFinal = []
    d2ParentFinal = []
    sRadiusParentAlongSegment = []


    for n2 in range(elementsCountAlongSegment + 1):
        phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
        xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
        radius = interp.interpolateCubicHermite([startRadiusParent], [startRadiusParentDerivative],
                                                [endRadiusParent], [endRadiusParentDerivative], xi)[0]
        sRadiusParentAlongSegment.append(radius)
        z = parentsegmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * parentsegmentLength

        xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
                                           elementsCountAround, startRadians=0.0)
        xParentFinal = xParentFinal + xLoop
        d1ParentFinal = d1ParentFinal + d1Loop

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2 * elementsCountAround + n1
            nx.append(xParentFinal[n])
            nd2.append(segmentAxisParent)
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        radius = sRadiusParentAlongSegment[n2]
        for n1 in range(elementsCountAround):
            d2ParentFinal.append(smoothd2Raw[n1][n2])

    # Calculate z mid-point for each element set along the segment
    faceParentMidPointsZ = []
    lengthToFirstPhase = startPhase / 360.0 * parentsegmentLength
    for n2 in range(elementsCountAlongSegment + 1):
        faceParentMidPointsZ += [lengthToFirstPhase +
                                 n2 * parentsegmentLength / elementsCountAlongSegment]

    # create nodes - DAUGHTER1
    ############################
    segmentAxisDaugh1 = [1.0, 0.0, 0.0]

    xDaugh1Final = []
    d1Daugh1Final = []
    d2Daugh1Final = []
    sRadiusDaugh1AlongSegment = []

    for n2 in range(elementsCountAlongSegment + 1):
        phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
        xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
        radius = interp.interpolateCubicHermite([startRadiusDaugh1], [startRadiusDaugh1Derivative],
                                                [endRadiusDaugh1], [endRadiusDaugh1Derivative], xi)[0]
        sRadiusDaugh1AlongSegment.append(radius)
        z = daugh1segmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * daugh1segmentLength

        xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
                                           elementsCountAround, startRadians=0.0)
        xDaugh1Final = xDaugh1Final + xLoop
        d1Daugh1Final = d1Daugh1Final + d1Loop

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2 * elementsCountAround + n1
            nx.append(xDaugh1Final[n])
            nd2.append(segmentAxisDaugh1)
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        radius = sRadiusDaugh1AlongSegment[n2]
        for n1 in range(elementsCountAround):
            d2Daugh1Final.append(smoothd2Raw[n1][n2])

    # Calculate z mid-point for each element set along the segment
    faceDaughter1MidPointsZ = []
    lengthToFirstPhase = startPhase / 360.0 * daugh1segmentLength
    for n2 in range(elementsCountAlongSegment + 1):
        faceDaughter1MidPointsZ += [lengthToFirstPhase +
                                 n2 * daugh1segmentLength / elementsCountAlongSegment]

    # create nodes - DAUGHTER2
    ############################
    segmentAxisDaugh2 = [-1.0, 0.0, 0.0]

    xDaugh2Final = []
    d1Daugh2Final = []
    d2Daugh2Final = []
    sRadiusDaugh2AlongSegment = []

    for n2 in range(elementsCountAlongSegment + 1):
        phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
        xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
        radius = interp.interpolateCubicHermite([startRadiusDaugh2], [startRadiusDaugh2Derivative],
                                                [endRadiusDaugh2], [endRadiusDaugh2Derivative], xi)[0]
        sRadiusDaugh2AlongSegment.append(radius)
        z = daugh2segmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * daugh2segmentLength

        xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
                                           elementsCountAround, startRadians=0.0)
        xDaugh2Final = xDaugh2Final + xLoop
        d1Daugh2Final = d1Daugh2Final + d1Loop

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2 * elementsCountAround + n1
            nx.append(xDaugh2Final[n])
            nd2.append(segmentAxisDaugh2)
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        radius = sRadiusDaugh2AlongSegment[n2]
        for n1 in range(elementsCountAround):
            d2Daugh2Final.append(smoothd2Raw[n1][n2])

    # Calculate z mid-point for each element set along the segment
    faceDaughter2MidPointsZ = []
    lengthToFirstPhase = startPhase / 360.0 * daugh1segmentLength
    for n2 in range(elementsCountAlongSegment + 1):
        faceDaughter2MidPointsZ += [lengthToFirstPhase +
                                 n2 * daugh2segmentLength / elementsCountAlongSegment]

    # WALL THICKNESS
    contractedWallThickness = wallThickness
    contractedWallThicknessList.append(contractedWallThickness)

    return xParentFinal, xDaugh1Final, xDaugh2Final, \
           d1ParentFinal, d1Daugh1Final, d1Daugh2Final, \
           d2ParentFinal, d2Daugh1Final, d2Daugh2Final, \
           transitElementList, contractedWallThicknessList,\
           segmentAxisParent, segmentAxisDaugh1, segmentAxisDaugh2, \
           faceParentMidPointsZ, faceDaughter1MidPointsZ, faceDaughter2MidPointsZ
