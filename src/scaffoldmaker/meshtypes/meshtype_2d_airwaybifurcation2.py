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

class MeshType_2d_airwaybifurcation2(Scaffold_base):
    '''
    2-D Airway Bifurcation scaffold.
    '''

    @staticmethod
    def getName():
        return '2D Bifurcation template 2'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 2,
            'Number of elements around' : 4,
            'Daughter interradius factor': 1.0,
            'Use cross derivatives' : False,
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
            'Use cross derivatives',
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

        nodeIdentifier = 1
        elementIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementAlong = (math.pi/4)/elementsCountAlong

        ########################################################################
        #### centerline details

        elementsCountAlongSegment = 2

        wallThickness = 0.02
        startPhase = 360.0


        startRadiusparentlist = [0.7, 0.52, 0.72]
        endRadiusparentlist = startRadiusparentlist
        startRadiusDaugh1list = [0.6, 0.35, 0.6]
        endRadiusDaugh1list = startRadiusDaugh1list
        startRadiusDaugh2list = [0.74, 0.64, 0.34]
        endRadiusDaugh2list = startRadiusDaugh2list

        tracheaoriginx = 0
        tracheaoriginy = 0
        tracheaoriginz = 0

        startRadiusparentDerivative = 0
        endRadiusparentDerivative = 0

        startRadiusDaugh1Derivative = 0.0
        endRadiusDaugh1Derivative = 0.0

        startRadiusDaugh2Derivative = 0.0
        endRadiusDaugh2Derivative = 0.0

        elementsCountAlongSegment = 2

        wallThickness = 0.02
        startPhase = 360.0

        ##DAUGHETER 1 - SAMPLE SEGMENT
        daughter1anglelist = [17, 4, 46]
        daughter2anglelist = [12, 39, 10]
        dirvecparentlist = [[0,0,1],[0.3864,0.0000,0.9223],[-0.4847,0.0000,0.87468]]
        # dirvecdaughter1list = [[-0.01513,0.28833,0.95741],[0.4100,1.7100,0.93000],[-0.1900,-0.34000,2.52000]]
        # dirvecdaughter2list= [[-0.07799,-0.20846,0.97492],[-0.150,-1.940,1.55000],[-0.290,0.3700,2.180]]

        dirvecdaughter1list = [[0.3864,0.0000,0.922],[0.4100,1.7100,0.93000],[-0.1900,-0.34000,2.52000]]
        dirvecdaughter2list = [[-0.485,0.0000,0.875],[-0.150,-1.940,1.55000],[-0.290,0.3700,2.180]]

        radiusvecparentlist =  [[0.0, 1.0, 0.0],[-0.08206,5.1855,1.55829],[-0.27244,-3.39860,0.73252]]
        radiusvecdaughter1list = [[-0.91,0.00000,0.3]]
        radiusvecdaughter2list = [[0.93,0.0000,0.51055]]

        parentx0list = [[0,0,0],[0,0,2.380],[-0.108,-0.0330,2.380]]
        daughter1x0list = [[0,0,4.000],[-0.190,1.5300,7.57000],[-0.3800,-0.760,5.780]]
        daughter2x0list = [[0,0,4.000],[-0.190,1.5300,7.57000],[-0.3800,-0.760,5.780]]

        daughter1x1list = [[1.18830,0.00000,6.99000],[3.2400,0.220,8.500],[-0.570,-1.100,8.300]]
        daughter2x1list = [[-1.11170,0.00000,5.89000],[-2.70000,-0.530,7.33000],[1.90000,-0.48000,9.75000]]

        segmentlengthparentlist = [4, 3.23, 2.14]
        segmentlengthdaughter1list = [3.23, 1.0, 1.275]
        segmentlengthdaughter2list = [2.14, 1.11, 1.24]
        #####################################################################

        #SMOOTHING PARAMETERS FOR EVERY SEGMENT
        xlensegmentparentlist = [0.7, 0.35, 0.45]
        xlensegmentd1list = [0.3, 0.3, 0.3]
        xlensegmentd2list = [0.35, 0.34, 0.34]

        ####################################################################################
        # Central path - SAMPLE
        segmentCount = 1  # Hardcoded for starters

        for nSegment in range(segmentCount):

            startRadiusparent = startRadiusparentlist[nSegment]
            endRadiusparent = startRadiusparent
            startRadiusDaugh1 = startRadiusDaugh1list[nSegment]
            endRadiusDaugh1 = startRadiusDaugh1
            startRadiusDaugh2 = startRadiusDaugh2list[nSegment]
            endRadiusDaugh2 = startRadiusDaugh2list[nSegment]

            cx_p1 = parentx0list[nSegment]
            cx0_d1 = daughter1x0list[nSegment]
            cx0_d2 = daughter2x0list[nSegment]

            cx1_d1 = daughter1x1list[nSegment]
            cx1_d2 = daughter2x1list[nSegment]

            cdv_p1 = dirvecparentlist[nSegment]
            cdv_d1 = dirvecdaughter1list[nSegment]
            cdv_d2 = dirvecdaughter2list[nSegment]
            cradv_d1 = radiusvecdaughter1list[nSegment]
            cradv_d2 = radiusvecdaughter2list[nSegment]

            parentsegmentLength = segmentlengthparentlist[nSegment]
            daughter1segmentLength = segmentlengthdaughter1list[nSegment]
            daughter2segmentLength = segmentlengthdaughter2list[nSegment]

            xlensegmentparent = xlensegmentparentlist[nSegment]
            xlensegmentd1 = xlensegmentd1list[nSegment]
            xlensegmentd2 = xlensegmentd2list[nSegment]

            if nSegment == 0:
                cxparent = [[cx_p1[0],cx_p1[1],cx_p1[2]],
                            [cx_p1[0]+cdv_p1[0]*xlensegmentparent*parentsegmentLength,
                             cx_p1[1]+cdv_p1[1]*xlensegmentparent*parentsegmentLength,
                             cx_p1[2]+cdv_p1[2]*xlensegmentparent*parentsegmentLength]]

            print('parent xcoord =',cxparent)

            cd1parent = [[cdv_p1[0],cdv_p1[1],cdv_p1[2]],[cdv_p1[0],cdv_p1[1],cdv_p1[2]]]
            print('parent dcs =',cd1parent)
            cd2parent = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
            cd12parent = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

            #Daughter 1
            cxdaughter1 = [[cx0_d1[0] + cdv_d1[0] * xlensegmentd1 * daughter1segmentLength,
                            cx0_d1[1] + cdv_d1[1] * xlensegmentd1 * daughter1segmentLength,
                            cx0_d1[2] + cdv_d1[2] * xlensegmentd1 * daughter1segmentLength],
                           [cx1_d1[0],cx1_d1[1],cx1_d1[2]]]
            cd1daughter1 = [[cdv_d1[0],cdv_d1[1],cdv_d1[2]],[cdv_d1[0],cdv_d1[1],cdv_d1[2]]]
            cd2daughter1 = [[cradv_d1[0],cradv_d1[1],cradv_d1[2]], [cradv_d1[0],cradv_d1[1],cradv_d1[2]]]
            cd12daughter1 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

            # Daughter 2
            cxdaughter2 = [[cx0_d2[0] + cdv_d2[0] * xlensegmentd2 * daughter2segmentLength,
                            cx0_d2[1] + cdv_d2[1] * xlensegmentd2 * daughter2segmentLength,
                            cx0_d2[2] + cdv_d2[2] * xlensegmentd2 * daughter2segmentLength],
                           [cx1_d2[0],cx1_d2[1],cx1_d2[2]]]
            print ('cx of daugh2=',cxdaughter2)
            cd1daughter2 = [[cdv_d2[0],cdv_d2[1],cdv_d2[2]],[cdv_d2[0],cdv_d2[1],cdv_d2[2]]]
            print ('cd1 of daughter2=',cd1daughter2)
            cd2daughter2 = [[cradv_d2[0],cradv_d2[1],cradv_d2[2]], [cradv_d2[0],cradv_d2[1],cradv_d2[2]]]
            print ('cd2 of daughter2=',cd2daughter2)

            cd12daughter2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

            daughter1angle = daughter1anglelist[nSegment]
            daughter2angle = daughter2anglelist[nSegment]

            ###########################################################################################

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
                xlensegmentparent*parentsegmentLength,
                (1-xlensegmentd1)*daughter1segmentLength,
                (1-xlensegmentd2)*daughter2segmentLength,
                wallThickness, radiusparentAlongSegment, radiusDaugh1AlongSegment,
                radiusDaugh2AlongSegment, dRadiusDaugh2AlongSegment,
                dRadiusDaugh1AlongSegment, dRadiusDaugh2AlongSegment, startPhase)


            xParentInner, xDaugh1Inner, xDaugh2Inner, \
            d1ParentInner, d1Daugh1Inner, d1Daugh2Inner, \
            d2ParentInner, d2Daugh1Inner, d2Daugh2Inner, \
            transitElementList, segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,\
            ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ = \
                airwaysegmentTubeMeshInnerPoints.getAirwaySegmentTubeMeshInnerPoints(nSegment)

            contractedWallThicknessList = airwaysegmentTubeMeshInnerPoints.getContractedWallThicknessList()

            # Warp segment points
            #####################
            xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList, \
            d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,\
            d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,\
            d3ParentWarpedUnitList, d3Daugh1WarpedList, d3Daugh2WarpedList \
                = tubebifurcationmesh.warpAirwaySegmentPoints(
                xParentInner, xDaugh1Inner, xDaugh2Inner,
                d1ParentInner, d1Daugh1Inner, d1Daugh2Inner,
                d2ParentInner, d2Daugh1Inner, d2Daugh2Inner,
                segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
                xlensegmentparent*parentsegmentLength, (1-xlensegmentd1)*daughter1segmentLength, (1-xlensegmentd2)*daughter2segmentLength,
                sxparent, sxDaugh1, sxDaugh2,
                sd1parent, sd1Daugh1, sd1Daugh2,
                sd2parent, sd2Daugh1, sd2Daugh2,
                elementsCountAround, elementsCountAlongSegment, nSegment,
                ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ)

            # Form junction  points
            ###############################
            xjunctionOuter, xjunctionInner, \
            d1junctionOuter, d1junctionInner, \
            d2junctionOuter, d2junctionInner \
                = tubebifurcationmesh.createjunctionAirwaySurfaceSegmentPoints(
                xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
                d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
                d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
                segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
                parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
                sxparent, sxDaugh1, sxDaugh2,
                sd1parent, sd1Daugh1, sd1Daugh2,
                sd2parent, sd2Daugh1, sd2Daugh2,
                elementsCountAround, elementsCountAlongSegment, nSegment)

            nextNodeIdentifier, nextElementIdentifier = \
                tubebifurcationmesh.createAirwaySegmentSurfaceNodesAndElements\
                    (region,
                    xParentWarpedList, d1ParentWarpedList, d2ParentWarpedList,
                    xDaugh1WarpedList, d1Daugh1WarpedList, d2Daugh1WarpedList,
                    xDaugh2WarpedList, d1Daugh2WarpedList, d2Daugh2WarpedList,
                    xjunctionOuter, xjunctionInner,
                    d1junctionOuter, d1junctionInner,
                    d2junctionOuter, d2junctionInner,
                    elementsCountAround, elementsCountAlongSegment,
                    nodeIdentifier, elementIdentifier, useCrossDerivatives)


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
    # NOTE: THe circle of points axis goes opposite
    # this has been done to allow appending multiple units past this segment
    #########################################################################
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
        xLoop, d1Loop = createCirclePoints([-z, 0.0, 0.0], [0.0, -radius, 0.0], [0.0, 0.0, radius],
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
