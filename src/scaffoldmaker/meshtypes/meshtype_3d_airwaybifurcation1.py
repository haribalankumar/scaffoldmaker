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

##from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_airwaybifurcation1(Scaffold_base):
    '''
    3-D Airway Bifurcation scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Bifurcation template 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 2,
            'Number of elements around' : 4,
            'Number of elements through wall': 1,
            'Daughter1 angle': 25,
            'Daughter2 angle': 55,
            'Daughter interradius factor': 1.0,
            'Use cross derivatives' : False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along segment': 1,
            'Refine number of elements through wall': 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements around',
            'Number of elements through wall',
            'Use cross derivatives',
            'Use linear through wall',
            'Daughter1 angle',
            'Daughter2 angle',
            'Daughter interradius factor',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall'
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
        if (options['Daughter1 angle'] < 20):
            options['Daughter1 angle'] = 20
        if (options['Daughter1 angle'] > 60):
            options['Daughter1 angle'] = 60

        if (options['Daughter2 angle'] < 20):
            options['Daughter2 angle'] = 20
        if (options['Daughter2 angle'] > 60):
            options['Daughter2 angle'] = 60

        if (options['Daughter interradius factor'] < 0.5):
            options['Daughter interradius factor'] = 0.5
        if (options['Daughter interradius factor'] > 2):
            options['Daughter interradius factor'] = 2

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
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])

        daughter1angle = options['Daughter1 angle']
        daughter2angle = options['Daughter2 angle']

        daughter_xrad = options['Daughter interradius factor']

        nodeIdentifier = 1
        elementIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementAlong = (math.pi/4)/elementsCountAlong

        ########################################################################
        #### centerline details
        ## typically comes from a centerline definition. currently hardcoding it
        parentsegmentLength = 1.0 ##When using centerline this is trach length - junctionheight
        daughter1segmentLength = 0.8 ##When using centerline this is trach length - junctionheight
        daughter2segmentLength = 0.6 ##When using centerline this is trach length - junctionheight

        elementsCountAlongSegment = 2

        wallThickness = 0.02
        startPhase = 360.0
        ########################################################################

        startRadiusparent = 0.3
        endRadiusparent = 0.3
        startRadiusDaugh1 = 0.26
        endRadiusDaugh1 = 0.26
        startRadiusDaugh2 = startRadiusDaugh1 * daughter_xrad
        endRadiusDaugh2 = startRadiusDaugh2

        tracheaoriginx = 0
        tracheaoriginy = 0
        tracheaoriginz = 0

        startRadiusparentDerivative = 0
        endRadiusparentDerivative = 0

        startRadiusDaugh1Derivative = 0.0
        endRadiusDaugh1Derivative = 0.0

        startRadiusDaugh2Derivative = 0.0
        endRadiusDaugh2Derivative = 0.0

        #######################################################################
        segmentCount = 1  # Hardcoded for starters


        xlensegmentparent = 0.9


        #Split ratio - decide where branching starts in daughter branches
        #################################################################
        cosangled1 = math.cos(math.pi/180.0 * (daughter1angle))
        sinangled1 = math.sin(math.pi/180.0 * (daughter1angle))
        cosangled2 = math.cos(math.pi/180.0 * (daughter2angle))
        sinangled2 = math.sin(math.pi/180.0 * (daughter2angle))


        # xlensegment is the proportion of branch lengths1 at which the branchlength starts
        # This offset value will be calculated from branch angle, daughter radius and daughter lengths.
        # The offset ensures smooth transition from parent to daughter
        # Having some bugs in this algorithm. Need to revisit

        # segmentratioDaughter1 = 2.0*startRadiusDaugh1*cosangled1/(daughter1segmentLength*sinangled1)
        # segmentratioDaughter2 = 2.0*startRadiusDaugh2*cosangled2/(daughter2segmentLength*sinangled2)
        # print('segment ratios = ', segmentratioDaughter1, segmentratioDaughter2)

        cval = -math.pow(1.3*endRadiusparent,2) + math.pow(cosangled1*startRadiusDaugh1,2)
        aval = math.pow(daughter1segmentLength,2)*math.pow(sinangled1,2)
        bval = 2 * daughter1segmentLength * endRadiusDaugh1 * cosangled1 * sinangled1
        segmentratioDaughter11 = -0.5*(bval/aval)+math.sqrt(math.pow(bval,2)-4.0*aval*cval)/(2.0*aval)

        cval = -math.pow(1.3*endRadiusparent,2) + math.pow(cosangled2*endRadiusDaugh2,2)
        aval = math.pow(daughter2segmentLength,2)*math.pow(sinangled2,2)
        bval = 2 * daughter2segmentLength * endRadiusDaugh2 * cosangled2 * sinangled2
        segmentratioDaughter21 = -0.5*(bval/aval)+math.sqrt(math.pow(bval,2)-4.0*aval*cval)/(2.0*aval)
        print('segment ratios = ', segmentratioDaughter11, segmentratioDaughter21)

        # #FIND THE BIGGEST and pick that
        # xlensegmentd1 = segmentratioDaughter1 if(segmentratioDaughter1>segmentratioDaughter11) else segmentratioDaughter11
        # xlensegmentd2 = segmentratioDaughter2 if(segmentratioDaughter2>segmentratioDaughter21) else segmentratioDaughter21
        xlensegmentd1 = segmentratioDaughter11
        xlensegmentd2 = segmentratioDaughter21

        print('xlen seg ratios = ', xlensegmentd1, xlensegmentd2)


        ####################################################################################
        # Central path - SAMPLE
        cxparent = [[tracheaoriginx, tracheaoriginy, tracheaoriginz],
                    [tracheaoriginx, tracheaoriginy, tracheaoriginz + xlensegmentparent * parentsegmentLength]]
        cd1parent = [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]
        cd2parent = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
        cd12parent = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        ##DAUGHETER 1 - SAMPLE SEGMENT
        cxdaughter1 = [[tracheaoriginx + xlensegmentd1 * daughter1segmentLength * sinangled1, tracheaoriginy,
               tracheaoriginz + parentsegmentLength + xlensegmentd1 * daughter1segmentLength * cosangled1],
              [tracheaoriginx + daughter1segmentLength * sinangled1, tracheaoriginy,
               tracheaoriginz + parentsegmentLength + daughter1segmentLength * cosangled1]]
        cd1daughter1 = [[sinangled1, 0.0, cosangled1], [sinangled1, 0.0, cosangled1]]
        cd2daughter1 = [[-cosangled1, 0.0, sinangled1], [-cosangled1, 0.0, sinangled1]]
        cd12daughter1 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        ##DAUGHETER 2 - SAMPLE SEGMENT
        cxdaughter2 = [ [ tracheaoriginx-xlensegmentd2*daughter2segmentLength*sinangled2, tracheaoriginy,
                 tracheaoriginz+parentsegmentLength+xlensegmentd2*daughter2segmentLength*cosangled2 ],
               [ tracheaoriginx-daughter2segmentLength*sinangled2, tracheaoriginy,
                 tracheaoriginz+parentsegmentLength+daughter2segmentLength*cosangled2] ]
        cd1daughter2 = [ [ -sinangled2, 0.0, cosangled2 ], [ -sinangled2, 0.0, cosangled2 ] ]
        cd2daughter2 = [ [ cosangled2, 0.0, sinangled2 ], [ cosangled2, 0.0, sinangled2 ] ]
        cd12daughter2 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]

        # Sample central path - PARENT
        sxparent, sd1parent, separent, sxiparent, ssfparent = \
            interp.sampleCubicHermiteCurves(cxparent, cd1parent, elementsCountAlongSegment*segmentCount)
        sd2parent = interp.interpolateSampleCubicHermite(cd2parent, cd12parent, separent, sxiparent, ssfparent)[0]

        # Sample central path - DAUGHTER1
        sxDaugh1, sd1Daugh1, seDaugh1, sxiDaugh1, ssfDaugh1 = \
            interp.sampleCubicHermiteCurves(cxdaughter1, cd1daughter1, elementsCountAlongSegment*segmentCount)
        sd2Daugh1 = interp.interpolateSampleCubicHermite(cd2daughter1, cd12daughter1, seDaugh1, sxiDaugh1, ssfDaugh1)[0]

        # Sample central path - DAUGHTER2
        sxDaugh2, sd1Daugh2, seDaugh2, sxiDaugh2, ssfDaugh2 = \
            interp.sampleCubicHermiteCurves(cxdaughter2, cd1daughter2, elementsCountAlongSegment*segmentCount)
        sd2Daugh2 = interp.interpolateSampleCubicHermite(cd2daughter2, cd12daughter2, seDaugh2, sxiDaugh2, ssfDaugh2)[0]

        # Find parameter variation along elementsCountAlongSegment
        radiusparentAlongSegment = []
        dRadiusparentAlongSegment = []

        radiusDaugh1AlongSegment = []
        dRadiusDaugh1AlongSegment = []

        radiusDaugh2AlongSegment = []
        dRadiusDaugh2AlongSegment = []

        for n2 in range(elementsCountAlongSegment + 1):
            xi = 1/elementsCountAlongSegment * n2
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
            xlensegmentparent*parentsegmentLength, xlensegmentd1*daughter1segmentLength, xlensegmentd2*daughter2segmentLength,
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

        # Warp segment points
        #####################
        xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList, \
        d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,\
        d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,\
        d3ParentWarpedUnitList, d3Daugh1WarpedUnitList, d3Daugh2WarpedUnitList \
            = tubebifurcationmesh.warpAirwaySegmentPoints(
            xParentInner, xDaugh1Inner, xDaugh2Inner,
            d1ParentInner, d1Daugh1Inner, d1Daugh2Inner,
            d2ParentInner, d2Daugh1Inner, d2Daugh2Inner,
            segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
            xlensegmentparent*parentsegmentLength, xlensegmentd1*daughter1segmentLength, xlensegmentd2*daughter2segmentLength,
            sxparent, sxDaugh1, sxDaugh2,
            sd1parent, sd1Daugh1, sd1Daugh2,
            sd2parent, sd2Daugh1, sd2Daugh2,
            elementsCountAround, elementsCountAlongSegment, nSegment,
            ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ)

        contractedWallThicknessList = airwaysegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Form junction  points
        ###############################
        xjunctionOuter, xjunctionInner, d1junctionOuter, d2junctionOuter, d3junctionOuter,\
            d1junctionInner, d2junctionInner, d3junctionInner\
            = tubebifurcationmesh.createjunctionAirwaySegmentPoints(
            xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
            d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
            d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
            segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
            parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
            sxparent, sxDaugh1, sxDaugh2,
            sd1parent, sd1Daugh1, sd1Daugh2,
            sd2parent, sd2Daugh1, sd2Daugh2,
            elementsCountAround, elementsCountAlongSegment, nSegment)

        # Create coordinates and derivatives - PARENT AND DAUGHERS1,2
        xParentList, d1ParentList, d2ParentList, d3ParentList, \
        xDaughter1List, d1Daughter1List, d2Daughter1List, d3Daughter1List, \
        xDaughter2List, d1Daughter2List, d2Daughter2List, d3Daughter2List, \
        curvatureList = \
            tubebifurcationmesh.getAirwaySegmentCoordinatesFromInner(
                xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
                d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
                d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
                d3ParentWarpedUnitList, d3Daugh1WarpedUnitList, d3Daugh2WarpedUnitList,
                contractedWallThicknessList,
                elementsCountAround, elementsCountAlongSegment,
                elementsCountThroughWall, transitElementList)

        # Create coordinates and derivatives - JUNCTION
        xJunctionOuterList, d1JunctionOuterList, d2JunctionOuterList, d3JunctionOuterList, \
        xJunctionInnerList, d1JunctionInnerList, d2JunctionInnerList, d3JunctionInnerList, \
        curvatureList = \
            tubebifurcationmesh.getAirwayJunctionCoordinatesFromInner(
                xjunctionOuter,d1junctionOuter,d2junctionOuter,d3junctionOuter,
                xjunctionInner,d1junctionInner,d2junctionInner,d3junctionInner,
                contractedWallThicknessList,
                elementsCountAround, elementsCountAlongSegment,
                elementsCountThroughWall, transitElementList)

        ##Create nodes and elements
        ##############################
        nextNodeIdentifier, nextElementIdentifier = \
            tubebifurcationmesh.createAirwaySegmentNodesAndElements\
                (region,
                 xParentList, d1ParentList, d2ParentList, d3ParentList,
                 xDaughter1List, d1Daughter1List, d2Daughter1List, d3Daughter1List,
                 xDaughter2List, d1Daughter2List, d2Daughter2List, d3Daughter2List,
                 xJunctionOuterList, xJunctionInnerList, d1JunctionOuterList, d1JunctionInnerList,
                 d2JunctionOuterList, d2JunctionInnerList,
                 d3JunctionOuterList, d3JunctionInnerList,
                 elementsCountAround, elementsCountAlongSegment,
                 nodeIdentifier, elementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives)

        # nextNodeIdentifier, nextElementIdentifier = \
        #     tubebifurcationmesh.createAirwaySegmentSurfaceNodesAndElements\
        #         (region,
        #          xParentWarpedList, d1ParentWarpedList, d2ParentWarpedList,
        #          xDaugh1WarpedList, d1Daugh1WarpedList, d2Daugh1WarpedList,
        #          xDaugh2WarpedList, d1Daugh2WarpedList, d2Daugh2WarpedList,
        #          xjunctionOuter, xjunctionInner, d1junction,d2junction,
        #          elementsCountAround, elementsCountAlongSegment,
        #         nodeIdentifier, elementIdentifier, useCrossDerivatives)

        # fm.endChange()
       # return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along segment']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong)

        return meshrefinement.getAnnotationGroups()


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

        return xParentInner, xDaugh1Inner, xDaugh2Inner,\
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

        if n2 == elementsCountAlongSegment:
            z = parentsegmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * parentsegmentLength

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

        if n2 == elementsCountAlongSegment:
            z = daugh1segmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * daugh1segmentLength

        #xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
        #                                   elementsCountAround, startRadians=0.0)
        xLoop, d1Loop = createCirclePoints([z, 0.0, 0.0], [0.0, radius, 0.0], [0.0, 0.0, radius],
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

        if (n2 == elementsCountAlongSegment):
            z = daugh2segmentLength / elementsCountAlongSegment * (1.0*n2 - 0.25) + startPhase / 360.0 * daugh2segmentLength

        # xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
        #                                    elementsCountAround, startRadians=0.0)
        xLoop, d1Loop = createCirclePoints([-z, 0.0, 0.0], [0.0, radius, 0.0], [0.0, 0.0, radius],
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

    # WALL THICKNESS - variable thickness not coded yet
    for n2 in range(elementsCountAlongSegment + 1):
        contractedWallThickness = wallThickness
        contractedWallThicknessList.append(contractedWallThickness)

    return xParentFinal, xDaugh1Final, xDaugh2Final, \
           d1ParentFinal, d1Daugh1Final, d1Daugh2Final, \
           d2ParentFinal, d2Daugh1Final, d2Daugh2Final, \
           transitElementList, contractedWallThicknessList,\
           segmentAxisParent, segmentAxisDaugh1, segmentAxisDaugh2, \
           faceParentMidPointsZ, faceDaughter1MidPointsZ, faceDaughter2MidPointsZ



# def getAirwaySegmentInnerPoints(region, elementsCountAround, elementsCountAlongSegment,
#                                 parentsegmentLength, daugh1segmentLength, daugh2segmentLength,
#                                 wallThickness, startRadiusParent, startRadiusDaugh1, startRadiusDaugh2,
#                                 startRadiusParentDerivative, startRadiusDaugh1Derivative, startRadiusDaugh2Derivative,
#                                 endRadiusParent, endRadiusDaugh1, endRadiusDaugh2,
#                                 endRadiusParentDerivative, endRadiusDaugh1Derivative, endRadiusDaugh2Derivative,
#                                      startPhase):
#     """
#     Generates a 3-D cylindrical segment mesh with variable numbers of elements
#     around, along the central path, and through wall.
#     :param elementsCountAround: Number of elements around.
#     :param elementsCountAlongSegment: Number of elements along cylindrical segment.
#     :param segmentLength: Length of a cylindrical segment.
#     :param wallThickness: Thickness of wall.
#     :param startRadius: Inner radius at proximal end.
#     :param startRadiusDerivative: Rate of change of inner radius at proximal end.
#     :param endRadius: Inner radius at distal end.
#     :param endRadiusDerivative: Rate of change of inner radius at distal end.
#     :param startPhase: Phase at start.
#     :return coordinates, derivatives on inner surface of a cylindrical segment.
#     :return transitElementList: stores true if element around is an element that
#     transits between a big and small element.
#     :return faceMidPointsZ: z-coordinate of midpoints for each element group
#     along segment.
#     """
#
#     transitElementList = [0] * elementsCountAround
#
#     contractedWallThicknessList = []
#
#     # create nodes - PARENT
#     #######################
#     segmentAxisParent = [0.0, 0.0, 1.0]
#
#     xParentFinal = []
#     d1ParentFinal = []
#     d2ParentFinal = []
#     sRadiusParentAlongSegment = []
#
#
#     for n2 in range(elementsCountAlongSegment + 1):
#         phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
#         xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
#         radius = interp.interpolateCubicHermite([startRadiusParent], [startRadiusParentDerivative],
#                                                 [endRadiusParent], [endRadiusParentDerivative], xi)[0]
#         sRadiusParentAlongSegment.append(radius)
#         z = parentsegmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * parentsegmentLength
#
#         if n2 == elementsCountAlongSegment:
#             z = parentsegmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * parentsegmentLength
#
#         xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
#                                            elementsCountAround, startRadians=0.0)
#         xParentFinal = xParentFinal + xLoop
#         d1ParentFinal = d1ParentFinal + d1Loop
#
#     # Smooth d2 for segment
#     smoothd2Raw = []
#     for n1 in range(elementsCountAround):
#         nx = []
#         nd2 = []
#         for n2 in range(elementsCountAlongSegment + 1):
#             n = n2 * elementsCountAround + n1
#             nx.append(xParentFinal[n])
#             nd2.append(segmentAxisParent)
#         smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
#         smoothd2Raw.append(smoothd2)
#
#     # Re-arrange smoothd2
#     for n2 in range(elementsCountAlongSegment + 1):
#         radius = sRadiusParentAlongSegment[n2]
#         for n1 in range(elementsCountAround):
#             d2ParentFinal.append(smoothd2Raw[n1][n2])
#
#     # Calculate z mid-point for each element set along the segment
#     faceParentMidPointsZ = []
#     lengthToFirstPhase = startPhase / 360.0 * parentsegmentLength
#     for n2 in range(elementsCountAlongSegment + 1):
#         faceParentMidPointsZ += [lengthToFirstPhase +
#                                  n2 * parentsegmentLength / elementsCountAlongSegment]
#
#     # create nodes - DAUGHTER1
#     ############################
#     segmentAxisDaugh1 = [1.0, 0.0, 0.0]
#
#     xDaugh1Final = []
#     d1Daugh1Final = []
#     d2Daugh1Final = []
#     sRadiusDaugh1AlongSegment = []
#
#     for n2 in range(elementsCountAlongSegment + 1):
#         phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
#         xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
#         radius = interp.interpolateCubicHermite([startRadiusDaugh1], [startRadiusDaugh1Derivative],
#                                                 [endRadiusDaugh1], [endRadiusDaugh1Derivative], xi)[0]
#         sRadiusDaugh1AlongSegment.append(radius)
#         z = daugh1segmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * daugh1segmentLength
#
#         if n2 == elementsCountAlongSegment:
#             z = daugh1segmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * daugh1segmentLength
#
#         #xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
#         #                                   elementsCountAround, startRadians=0.0)
#         xLoop, d1Loop = createCirclePoints([z, 0.0, 0.0], [0.0, radius, 0.0], [0.0, 0.0, radius],
#                                           elementsCountAround, startRadians=0.0)
#         xDaugh1Final = xDaugh1Final + xLoop
#         d1Daugh1Final = d1Daugh1Final + d1Loop
#
#     # Smooth d2 for segment
#     smoothd2Raw = []
#     for n1 in range(elementsCountAround):
#         nx = []
#         nd2 = []
#         for n2 in range(elementsCountAlongSegment + 1):
#             n = n2 * elementsCountAround + n1
#             nx.append(xDaugh1Final[n])
#             nd2.append(segmentAxisDaugh1)
#         smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
#         smoothd2Raw.append(smoothd2)
#
#     # Re-arrange smoothd2
#     for n2 in range(elementsCountAlongSegment + 1):
#         radius = sRadiusDaugh1AlongSegment[n2]
#         for n1 in range(elementsCountAround):
#             d2Daugh1Final.append(smoothd2Raw[n1][n2])
#
#     # Calculate z mid-point for each element set along the segment
#     faceDaughter1MidPointsZ = []
#     lengthToFirstPhase = startPhase / 360.0 * daugh1segmentLength
#     for n2 in range(elementsCountAlongSegment + 1):
#         faceDaughter1MidPointsZ += [lengthToFirstPhase +
#                                  n2 * daugh1segmentLength / elementsCountAlongSegment]
#
#     # create nodes - DAUGHTER2
#     ############################
#     segmentAxisDaugh2 = [-1.0, 0.0, 0.0]
#
#     xDaugh2Final = []
#     d1Daugh2Final = []
#     d2Daugh2Final = []
#     sRadiusDaugh2AlongSegment = []
#
#     for n2 in range(elementsCountAlongSegment + 1):
#         phase = startPhase + n2 * 360.0 / elementsCountAlongSegment
#         xi = (phase if phase <= 360.0 else phase - 360.0) / 360.0
#         radius = interp.interpolateCubicHermite([startRadiusDaugh2], [startRadiusDaugh2Derivative],
#                                                 [endRadiusDaugh2], [endRadiusDaugh2Derivative], xi)[0]
#         sRadiusDaugh2AlongSegment.append(radius)
#         z = daugh2segmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * daugh2segmentLength
#
#         if (n2 == elementsCountAlongSegment):
#             z = daugh2segmentLength / elementsCountAlongSegment * (1.0*n2 - 0.25) + startPhase / 360.0 * daugh2segmentLength
#
#         # xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
#         #                                    elementsCountAround, startRadians=0.0)
#         xLoop, d1Loop = createCirclePoints([-z, 0.0, 0.0], [0.0, radius, 0.0], [0.0, 0.0, radius],
#                                            elementsCountAround, startRadians=0.0)
#
#         xDaugh2Final = xDaugh2Final + xLoop
#         d1Daugh2Final = d1Daugh2Final + d1Loop
#
#     # Smooth d2 for segment
#     smoothd2Raw = []
#     for n1 in range(elementsCountAround):
#         nx = []
#         nd2 = []
#         for n2 in range(elementsCountAlongSegment + 1):
#             n = n2 * elementsCountAround + n1
#             nx.append(xDaugh2Final[n])
#             nd2.append(segmentAxisDaugh2)
#         smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
#         smoothd2Raw.append(smoothd2)
#
#     # Re-arrange smoothd2
#     for n2 in range(elementsCountAlongSegment + 1):
#         radius = sRadiusDaugh2AlongSegment[n2]
#         for n1 in range(elementsCountAround):
#             d2Daugh2Final.append(smoothd2Raw[n1][n2])
#
#     # Calculate z mid-point for each element set along the segment
#     faceDaughter2MidPointsZ = []
#     lengthToFirstPhase = startPhase / 360.0 * daugh2segmentLength
#     for n2 in range(elementsCountAlongSegment + 1):
#         faceDaughter2MidPointsZ += [lengthToFirstPhase +
#                                  n2 * daugh2segmentLength / elementsCountAlongSegment]
#
#
#     # WALL THICKNESS - variable thickness not coded yet
#     for n2 in range(elementsCountAlongSegment + 1):
#         contractedWallThickness = wallThickness
#         contractedWallThicknessList.append(contractedWallThickness)
#
#
#     return xParentFinal, xDaugh1Final, xDaugh2Final, \
#            d1ParentFinal, d1Daugh1Final, d1Daugh2Final, \
#            d2ParentFinal, d2Daugh1Final, d2Daugh2Final, \
#            transitElementList, contractedWallThicknessList,\
#            segmentAxisParent, segmentAxisDaugh1, segmentAxisDaugh2, \
#            faceParentMidPointsZ, faceDaughter1MidPointsZ, faceDaughter2MidPointsZ
