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
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite

##from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_airwaybifurcation3(Scaffold_base):
    '''
    3-D Airway Bifurcation scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Bifurcation with lung'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 2,
            'Number of elements around' : 4,
            'Number of elements through wall': 1,
            'Use junction elements': False,
            'Number of generations downstream': 1,
            'Include left lung': False,
            'Number of lung elements up': 4,
            'Number of lung elements laterally': 4,
            'Lung Height': 10,
            'Lung Width': 10,
            'Lung Upcurve coefficient1': 1.0,
            'Lung Upcurve coefficient2': 1.0,
            'Lung Medialsurface coefficient': 0.2,
            'Lung Lateralsurface coefficient': 0.2,
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
            'Use junction elements',
            'Number of generations downstream',
            'Include left lung',
            'Number of lung elements up',
            'Number of lung elements laterally',
            'Lung Height',
            'Lung Width',
            'Lung Upcurve coefficient1',
            'Lung Upcurve coefficient2',
            'Lung Medialsurface coefficient',
            'Lung Lateralsurface coefficient',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements along'] < 2):
            options['Number of elements along'] = 2

        if (options['Number of elements around'] > 8):
            options['Number of elements around'] = 8
        elif (options['Number of elements around'] % 2) == 1:
            options['Number of elements around'] += 1

        if (options['Number of elements around'] < 4):
            options['Number of elements around'] = 4

        if (options['Number of generations downstream'] < 1):
            options['Number of generations downstream'] = 1
        if (options['Number of generations downstream'] > 2):
            options['Number of generations downstream'] = 2

        #LUNG PARAMS
        if options['Number of lung elements up'] < 3:
            options['Number of lung elements up'] = 3
        if options['Number of lung elements laterally'] < 3:
            options['Number of lung elements laterally'] = 3

        if options['Lung Width'] < 2:
            options['Lung Width'] = 2
        if options['Lung Height'] < 2:
            options['Lung Height'] = 2

        if options['Lung Upcurve coefficient1'] > 4:
            options['Lung Upcurve coefficient1'] = 4
        if options['Lung Upcurve coefficient2'] > 4:
                options['Lung Upcurve coefficient2'] = 4

        if options['Lung Medialsurface coefficient'] > 0.95:
            options['Lung Medialsurface coefficient'] = 0.95

        if options['Lung Lateralsurface coefficient'] > 0.95:
            options['Lung Lateralsurface coefficient'] = 0.95

        if options['Lung Medialsurface coefficient'] < 0.05:
            options['Lung Medialsurface coefficient'] = 0.05

        if options['Lung Lateralsurface coefficient'] < 0.05:
            options['Lung Lateralsurface coefficient'] = 0.05


    @staticmethod
    def generateBaseMesh(region, options):
        '''
        Generate the base Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        '''

        elementsCountAlong = options['Number of elements along']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        numberofGenerations = options['Number of generations downstream']

        useCrossDerivatives = options['Use cross derivatives']
        useJunctionElements = options['Use junction elements']
        includeLeftLung = options['Include left lung']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])

        ######LUNG PARAMETERS
        #################################
        elementsCountlateral = options['Number of lung elements laterally']
        elementsCountUp = options['Number of lung elements up']
        lungheight = options['Lung Height']
        lungwidth = options['Lung Width']
        posteriorcurveradiuscoeff = options['Lung Upcurve coefficient1']
        anteriorcurveradiuscoeff = options['Lung Upcurve coefficient2']
        medialsurfcoeff = options['Lung Medialsurface coefficient']
        lateralsurfcoeff = options['Lung Lateralsurface coefficient']


        nodeIdentifier = 1
        elementIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementAlong = (math.pi/4)/elementsCountAlong

        xjunctionOuter = []
        xjunctionInner = []
        d1junctionOuter = []
        d2junctionOuter = []
        d3junctionOuter = []
        d1junctionInner = []
        d2junctionInner = []
        d3junctionInner = []
        xJunctionOuterList = []
        d1JunctionOuterList = []
        d2JunctionOuterList = []
        d3JunctionOuterList = []
        xJunctionInnerList = []
        d1JunctionInnerList = []
        d2JunctionInnerList = []
        d3JunctionInnerList = []

        ########################################################################
        #### centerline details
        ## typically comes from a centerline definition. currently hardcoding it


        startRadiusparentlist = [0.7, 0.6, 0.72]
        endRadiusparentlist = startRadiusparentlist
        startRadiusDaugh1list = [0.6, 0.35, 0.6]
        endRadiusDaugh1list = startRadiusDaugh1list
        startRadiusDaugh2list = [0.72, 0.64, 0.35]
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
        startPhase = 0.0

        ##DAUGHETER 1 - SAMPLE SEGMENT
        daughter1anglelist = [17, 4, 46]
        daughter2anglelist = [12, 39, 10]
        dirvecparentlist = [[0,0,1],[0.3864,0.0000,0.922],[-0.485,0.0000,0.875]]
        dirvecdaughter1list = [[0.3864,0.0000,0.922],[0.87929,0.00000,0.47628],[-0.07125,0.00000,0.99746]]
        dirvecdaughter2list = [[-0.485,0.0000,0.875],[0.05673,0.00000,0.99839],[-0.9445,0.00000,0.32852]]

        # radiusvecparentlist =  [[0.0, 1.0, 0.0],[-0.91,0.00000,0.3],[0.93,0.0000,0.51055]]
        # radiusvecdaughter1list = [[-0.91,0.00000,0.3],[-0.47628,0.00000,0.87929]]
        # radiusvecdaughter2list = [[0.93,0.0000,0.51055],[0.99839,0.00000,-0.05673]]

        radiusvecparentlist =  [[0.0, 1.0, 0.0],[0.0, 1.0, 0.0],[0.0,1.0000,0.0]]
        radiusvecdaughter1list = [[0.0,1.00,0.0],[0.0,1.00,0.0],[0.0,1.00,0.0]]
        radiusvecdaughter2list = [[0.0,1.00,0.0],[0.0,1.00,0.0],[0.0,1.00,0.0]]

        parentx0list = [[0,0,0],[0,0,4.000],[0,0,4.000]]
        daughter1x0list = [[0,0,4.000],[2.50530,0.00000,9.98000],[-2.09470,0.00000,7.78000]]
        daughter2x0list = [[0,0,4.000],[2.50530,0.00000,9.98000],[-2.09470,0.00000,7.78000]]

        parentx1list = [[0,0,4.000],[2.50530,0.00000,9.98000],[-2.09470,0.00000,7.78000]]
        daughter1x1list = [[1.25265,0.00000,6.99000],[4.90530,0.00000,11.28000],[-2.2947,0.000,10.580]]
        daughter2x1list = [[-1.04735,0.00000,5.89000],[2.70530,0.00000,13.50000],[-4.3947,0.00, 8.580]]

        segmentlengthparentlist = [4, 3.24, 2.16]
        segmentlengthdaughter1list = [3.24, 2.7, 2.8]
        segmentlengthdaughter2list = [2.16, 3.5, 2.43]
        #####################################################################

        # SMOOTHING PARAMETERS FOR EVERY SEGMENT
        xlensegmentparentlist = [1.0, 1., 1.]
        xlensegmentd1list = [0.0, 0.0, 0.0]
        xlensegmentd2list = [0.0, 0.0, 0.0]

        if useJunctionElements:
            #SMOOTHING PARAMETERS FOR EVERY SEGMENT
            xlensegmentparentlist = [0.8, 0.6, 0.6]
            xlensegmentd1list = [0.6, 0.5, 0.4]
            xlensegmentd2list = [0.7, 0.4, 0.5]

        segmentCount = pow(2,numberofGenerations-1)

        ####################################################################################
        #-----------------------------------------------
        # hardcoded LUNG VALUES

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        TrachealRotationAngleRadians = 0
        thetaup =  math.pi/3.0
        thetalateral = math.pi/6.0
        posteriorcurveradius = lungheight/(2*math.sin(thetaup))
        anteriorcurveradius = posteriorcurveradius

        LMBcentre = 0.0
        RMBcentre = LMBcentre + lungwidth*0.1

        zlateralcentre = 5.0
        ywidthcentre = 35.0
        halflungwidth = lungwidth*0.5
        lungdepth=halflungwidth

        ycentreleftmedial = LMBcentre - halflungwidth
        ycentreleftlateral = LMBcentre + halflungwidth

        ymidLateralMedial = (ycentreleftmedial+ycentreleftlateral)*0.5

        lateralsurfradius =  lateralsurfcoeff*lungwidth
        internalsurfradius =  lateralsurfradius*2.25

        lRadiansUp = []


        ####################################################################################
        # Central path
        cxminus_d1 = []
        cdvminus_d1 = []
        xlensegmentparentminus = 0
        parentsegmentLengthminus = 0

        for nSegment in range(segmentCount):
            childcount = pow(2,nSegment)
            for nChild in range(childcount):

                if nSegment > 0:
                    nodeIdentifier = nextnodeIdentifier
                    elementIdentifier = nextelementIdentifier

                nChild_number = nSegment + nChild
                startRadiusparent = startRadiusparentlist[nChild_number]
                endRadiusparent = startRadiusparent
                startRadiusDaugh1 = startRadiusDaugh1list[nChild_number]
                endRadiusDaugh1 = startRadiusDaugh1
                startRadiusDaugh2 = startRadiusDaugh2list[nChild_number]
                endRadiusDaugh2 = startRadiusDaugh2list[nChild_number]

                cx0_p1 = parentx0list[nChild_number]
                cx1_p1 = parentx1list[nChild_number]

                cx0_d1 = daughter1x0list[nChild_number]
                cx0_d2 = daughter2x0list[nChild_number]
                cx1_d1 = daughter1x1list[nChild_number]
                cx1_d2 = daughter2x1list[nChild_number]

                cdv_p1 = dirvecparentlist[nChild_number]
                cdv_d1 = dirvecdaughter1list[nChild_number]
                cdv_d2 = dirvecdaughter2list[nChild_number]

                if nSegment > 0:
                    cxminus_d1 = daughter1x0list[nSegment-1]
                    cdvminus_d1 = dirvecdaughter1list[nSegment-1]

                cradv_p = radiusvecparentlist[nChild_number]
                cradv_d1 = radiusvecdaughter1list[nChild_number]
                cradv_d2 = radiusvecdaughter2list[nChild_number]

                parentsegmentLength = segmentlengthparentlist[nChild_number]
                daughter1segmentLength = segmentlengthdaughter1list[nChild_number]
                daughter2segmentLength = segmentlengthdaughter2list[nChild_number]

                xlensegmentparent = xlensegmentparentlist[nChild_number]
                xlensegmentd1 = xlensegmentd1list[nChild_number]
                xlensegmentd2 = xlensegmentd2list[nChild_number]

                if nSegment > 0:
                    xlensegmentparentminus = xlensegmentd1list[nSegment-1]
                    parentsegmentLengthminus = segmentlengthdaughter1list[nSegment-1]


                if nSegment == 0:
                    cxparent = [[cx0_p1[0],cx0_p1[1],cx0_p1[2]],
                            [cx0_p1[0]+cdv_p1[0]*xlensegmentparent*parentsegmentLength,
                             cx0_p1[1]+cdv_p1[1]*xlensegmentparent*parentsegmentLength,
                             cx0_p1[2]+cdv_p1[2]*xlensegmentparent*parentsegmentLength]]
                else:
                    cxparent = [[cx1_p1[0] - cdv_p1[0] * parentsegmentLength,
                                 cx1_p1[1] - cdv_p1[1]  * parentsegmentLength,
                                 cx1_p1[2] - cdv_p1[2]  * parentsegmentLength],
                                [cx1_p1[0],cx1_p1[1],cx1_p1[2]]]
                print('parent xcoord =',cxparent)
                cd1parent = [[cdv_p1[0],cdv_p1[1],cdv_p1[2]],[cdv_p1[0],cdv_p1[1],cdv_p1[2]]]
                cd2parent = [[cradv_p[0],cradv_p[1],cradv_p[2]],[cradv_p[0],cradv_p[1],cradv_p[2]]]
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
                           [cx0_d2[0] + cdv_d2[0] * daughter2segmentLength,
                            cx0_d2[1] + cdv_d2[1] * daughter2segmentLength,
                            cx0_d2[2] + cdv_d2[2] * daughter2segmentLength]]
                           # [cx1_d2[0],cx1_d2[1],cx1_d2[2]]]
                cd1daughter2 = [[cdv_d2[0],cdv_d2[1],cdv_d2[2]],[cdv_d2[0],cdv_d2[1],cdv_d2[2]]]
                cd2daughter2 = [[cradv_d2[0],cradv_d2[1],cradv_d2[2]], [cradv_d2[0],cradv_d2[1],cradv_d2[2]]]

                cd12daughter2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

                # daughter1angle = daughter1anglelist[nSegment]
                # daughter2angle = daughter2anglelist[nSegment]
        #######################################################################

                # Sample central path - PARENT
                # sxparent, sd1parent, separent, sxiparent, ssfparent = \
                #     interp.sampleCubicHermiteCurves(cxparent, cd1parent, elementsCountAlongSegment*segmentCount)
                sxparent, sd1parent, separent, sxiparent, ssfparent = \
                    interp.sampleCubicHermiteCurves(cxparent, cd1parent, elementsCountAlongSegment)
                sd2parent = interp.interpolateSampleCubicHermite(cd2parent, cd12parent, separent, sxiparent, ssfparent)[0]
                print('after sampling hermite sxparent=',sxparent)
                print('after sampling hermite sd1parent=',sd1parent)
                print('after sampling hermite sd2parent=',sd2parent)

                # Sample central path - DAUGHTER1
                # sxDaugh1, sd1Daugh1, seDaugh1, sxiDaugh1, ssfDaugh1 = \
                #     interp.sampleCubicHermiteCurves(cxdaughter1, cd1daughter1, elementsCountAlongSegment*segmentCount)
                sxDaugh1, sd1Daugh1, seDaugh1, sxiDaugh1, ssfDaugh1 = \
                    interp.sampleCubicHermiteCurves(cxdaughter1, cd1daughter1, elementsCountAlongSegment)
                sd2Daugh1 = interp.interpolateSampleCubicHermite(cd2daughter1, cd12daughter1, seDaugh1, sxiDaugh1, ssfDaugh1)[0]

                # Sample central path - DAUGHTER2
                # sxDaugh2, sd1Daugh2, seDaugh2, sxiDaugh2, ssfDaugh2 = \
                #     interp.sampleCubicHermiteCurves(cxdaughter2, cd1daughter2, elementsCountAlongSegment*segmentCount)
                sxDaugh2, sd1Daugh2, seDaugh2, sxiDaugh2, ssfDaugh2 = \
                    interp.sampleCubicHermiteCurves(cxdaughter2, cd1daughter2, elementsCountAlongSegment)
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
                    radiusDaugh2AlongSegment, dRadiusparentAlongSegment,
                    dRadiusDaugh1AlongSegment, dRadiusDaugh2AlongSegment, startPhase)

                #Create inner points
                # nSegment = 0

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
                    xlensegmentparent*parentsegmentLength,
                    (1-xlensegmentd1)*daughter1segmentLength,
                    (1-xlensegmentd2)*daughter2segmentLength,
                    sxparent, sxDaugh1, sxDaugh2,
                    sd1parent, sd1Daugh1, sd1Daugh2,
                    sd2parent, sd2Daugh1, sd2Daugh2,
                    elementsCountAround, elementsCountAlongSegment, nSegment,
                    ParentfaceMidPointsZ, Daughter1faceMidPointsZ, Daughter2faceMidPointsZ)

                contractedWallThicknessList = airwaysegmentTubeMeshInnerPoints.getContractedWallThicknessList()

                # Form junction  points
                ###############################
                if useJunctionElements:
                    xjunctionOuter, xjunctionInner, d1junctionOuter, d2junctionOuter, d3junctionOuter,\
                    d1junctionInner, d2junctionInner, d3junctionInner\
                        = tubebifurcationmesh.createjunctionAirwaySegmentPoints(
                        xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
                        d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
                        d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
                        segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
                        sxparent, sxDaugh1, sxDaugh2, sd1parent, sd1Daugh1, sd1Daugh2, sd2parent, sd2Daugh1, sd2Daugh2,
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

                if useJunctionElements:
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
            # if nSegment == 1:
            #     nodeIdentifier = 3*elementsCountAround*(1+elementsCountAlongSegment)*(elementsCountThroughWall+1)+10
            #     elementIdentifier = 3*elementsCountAround*(elementsCountAlongSegment)*elementsCountThroughWall+10

                print('checking node identified BEFORE create node and elems =', nodeIdentifier, elementIdentifier)
                nextnodeIdentifier, nextelementIdentifier = \
                    tubebifurcationmesh.createAirwaySegmentNodesAndElements(
                        region,
                        xParentList, d1ParentList, d2ParentList, d3ParentList,
                        xDaughter1List, d1Daughter1List, d2Daughter1List, d3Daughter1List,
                        xDaughter2List, d1Daughter2List, d2Daughter2List, d3Daughter2List,
                        xJunctionOuterList, xJunctionInnerList, d1JunctionOuterList, d1JunctionInnerList,
                        d2JunctionOuterList, d2JunctionInnerList,
                        d3JunctionOuterList, d3JunctionInnerList,
                        elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall,
                        nodeIdentifier, elementIdentifier,
                        useJunctionElements, useCubicHermiteThroughWall, useCrossDerivatives)

                print('checking node identifie after create =', nodeIdentifier, elementIdentifier)


        ###############################################################################################
        nodeIdentifier = nextnodeIdentifier
        # # Create node location and derivative 2 on left posterior edge upwards - CURVE
        lRadiansUp = []
        lRadiansLateral = []
        x = []
        nx = []
        d1 = []
        nd1 = []
        d2 = []
        nd2 = []

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        else:
            nodetemplate = nodetemplateApex

        cache = fm.createFieldcache()


        for n2 in range(3):
            radiansUp = -thetaup + (n2) / (2.0) * (2 * thetaup)
            lRadiansUp.append(radiansUp)

            radiansLateral = -thetalateral + (n2) / (2.0) * (2 * thetalateral)
            lRadiansLateral.append(radiansLateral)

        if includeLeftLung:
            # RIGHT LUNG
            lungwidth = lungwidth * 0.8

            nx = []
            nd1 = []
            for n2 in range(3):
                radiansUp = lRadiansUp[n2]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                x1 = [RMBcentre, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                    zlateralcentre + posteriorcurveradius * sinRadiansUp]
                d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
                nx.append(x1)
                nd1.append(d1)

            sPosteriorx, sPosteriorderiv2, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1,
                                                                                     elementsCountOut=elementsCountUp)

            nx = []
            nd1 = []
            for n2 in range(3):
                radiansUp = lRadiansUp[n2]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                x1 = [RMBcentre + 2, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                    zlateralcentre + posteriorcurveradius * sinRadiansUp]
                d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
                nx.append(x1)
                nd1.append(d1)

            sPosteriorOffsetx, sPosteriorOffsetderiv2, _, _, _ = \
                interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

            zero = [0, 0, 0]

            # -----------------------------------------
            ## Anterior edge
            ## -------------------------------------
            nx = []
            nd1 = []
            for n2 in range(3):
                radiansUp = lRadiansUp[n2]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                x1 = [RMBcentre, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + anteriorcurveradius * sinRadiansUp]
                d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
                nx.append(x1)
                nd1.append(d1)

            sAnteriorx, sAnteriorderiv2, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1,
                                                                                   elementsCountOut=elementsCountUp)

            # -----------------------------------------
            ## Anterior edge offset
            ## ---------------------------------------
            nx = []
            nd1 = []
            for n2 in range(3):
                radiansUp = lRadiansUp[n2]
                cosRadiansUp = math.cos(radiansUp)
                sinRadiansUp = math.sin(radiansUp)
                x1 = [RMBcentre + 0.5, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
                    zlateralcentre + anteriorcurveradius * sinRadiansUp]
                d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
                nx.append(x1)
                nd1.append(d1)

            sAnteriorOffsetx, sAnteriorderivOffset2, _, _, _ = \
                interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

            # Apply tracheal rotation angle
            # ---------------------------------
            xcen = sPosteriorx[0][0]
            ycen = sPosteriorx[0][1]
            zcen = sPosteriorx[0][2]
            for n2 in range(0, elementsCountUp):
                newx = (sPosteriorx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sPosteriorx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sPosteriorx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sPosteriorx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                sPosteriorx[n2][0] = newx + xcen
                sPosteriorx[n2][1] = newy + ycen

                newx = (sAnteriorx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sAnteriorx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sAnteriorx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sAnteriorx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                sAnteriorx[n2][0] = newx + xcen
                sAnteriorx[n2][1] = newy + ycen

                newx = (sPosteriorOffsetx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sPosteriorOffsetx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sPosteriorOffsetx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sPosteriorOffsetx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                sPosteriorOffsetx[n2][0] = newx + xcen
                sPosteriorOffsetx[n2][1] = newy + ycen

                newx = (sAnteriorOffsetx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sAnteriorOffsetx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sAnteriorOffsetx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sAnteriorOffsetx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                sAnteriorOffsetx[n2][0] = newx + xcen
                sAnteriorOffsetx[n2][1] = newy + ycen

            # # ADD nodes on posterior edge upwards
            # # ------------------------------------
            latUpNodeId = []
            for n2 in range(0, elementsCountUp):
                layerNodeId = []
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sPosteriorx[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sPosteriorderiv2[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
                latUpNodeId.append(layerNodeId)

            # # ADD nodes on anterior edge upwards
            # latUpNodeId = []
            for n2 in range(0, elementsCountUp):
                layerNodeId = []
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sAnteriorx[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sAnteriorderiv2[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
                latUpNodeId.append(layerNodeId)

            # # ADD nodes on posterior edge OFFSET upwards
            # # ------------------------------------
            latUpNodeId = []
            for n2 in range(0, elementsCountUp):
                layerNodeId = []
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sPosteriorOffsetx[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sPosteriorOffsetderiv2[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
                latUpNodeId.append(layerNodeId)

            # # ADD nodes on anterior edge upwards
            # latUpNodeId = []
            for n2 in range(0, elementsCountUp):
                layerNodeId = []
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sAnteriorOffsetx[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sAnteriorderivOffset2[n2])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                layerNodeId.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1
                latUpNodeId.append(layerNodeId)

            # ----------------------------------
            # # Create nodes on MEDIAL side
            # ----------------------------------
            for n2 in range(0, elementsCountUp):

                x1 = sPosteriorx[n2]
                x3 = sAnteriorx[n2]
                halflung = 0.5 * (x1[1] + x3[1])
                zlung = 0.5 * (x1[2] + x3[2])
                zratio = (0.25 - 0.4 * abs(zlung / lungheight - 0.5)) ** 0.5

                radiansLateral = lRadiansLateral[0]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)

                d21 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                       0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

                radiansLateral = lRadiansLateral[1]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)
                x2 = [x1[0] + 0.2 * lungdepth + 0.5 * lungdepth * medialsurfcoeff * (cosRadiansLateral - 1.0),
                      halflung,
                      x1[2]]

                d22 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                       -0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral,
                       0]

                radiansLateral = lRadiansLateral[2]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)
                d23 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                       0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

                nx = [x1, x2, x3]
                nd2 = [d21, d22, d23]
                sMedialx, sMedialderiv2, _, _, _  = \
                    interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)

                x1 = sPosteriorOffsetx[n2]
                x3 = sAnteriorOffsetx[n2]
                halflung = 0.5 * (x1[1] + x3[1])

                radiansLateral = lRadiansLateral[0]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)
                d21 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                       0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                       0]

                radiansLateral = lRadiansLateral[1]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)
                x2 = [x1[0] + zratio * lungdepth + 0.5 * lungdepth * lateralsurfcoeff * (cosRadiansLateral - 1.0),
                      halflung,
                      x1[2]]
                d22 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                       -0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                       0]

                radiansLateral = lRadiansLateral[2]
                cosRadiansLateral = math.cos(radiansLateral)
                sinRadiansLateral = math.sin(radiansLateral)
                d23 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                       0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                       0]

                nx = [x1, x2, x3]
                nd2 = [d21, d22, d23]
                sLateralx, sLateralderiv2,  _, _, _  = \
                    interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)
                #

                # Apply tracheal rotation angle
                # ---------------------------------
                for n3 in range(1, elementsCountlateral):
                    newx = (sMedialx[n3][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                           (sMedialx[n3][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                    newy = (sMedialx[n3][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                           (sMedialx[n3][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                    sMedialx[n3][0] = newx + xcen
                    sMedialx[n3][1] = newy + ycen

                    newx = (sLateralx[n3][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                           (sLateralx[n3][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                    newy = (sLateralx[n3][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                           (sLateralx[n3][1] - ycen) * math.cos(-TrachealRotationAngleRadians)
                    sLateralx[n3][0] = newx + xcen
                    sLateralx[n3][1] = newy + ycen

                for n3 in range(1, elementsCountlateral):
                    layerNodeId = []
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sMedialx[n3])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sMedialderiv2[n3])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    layerNodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1
                    latUpNodeId.append(layerNodeId)

                for n3 in range(1, elementsCountlateral):
                    layerNodeId = []
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sLateralx[n3])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sLateralderiv2[n3])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    layerNodeId.append(nodeIdentifier)
                    nodeIdentifier = nodeIdentifier + 1
                    latUpNodeId.append(layerNodeId)


            ###########################################################################
            # RIGHT LUNG
            ###########################################################################
            mesh = fm.findMeshByDimension(3)
            tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
            eft = tricubichermite.createEftBasic()

            tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

            # Regular elements
            elementtemplate = mesh.createElementtemplate()
            elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            elementtemplate.defineField(coordinates, -1, eft)

            elementIdentifier = nextelementIdentifier
            e0 = nodeIdentifier - 1

            delQ = 4 * elementsCountUp
            delX = 2 * elementsCountlateral - 2
            delA = 4 * elementsCountUp + (elementsCountlateral - 1)
            delB = 4 * elementsCountUp + (3 * elementsCountlateral - 3)
            delM = 2 * elementsCountUp
            delP = 4 * elementsCountUp + (elementsCountlateral - 2)

            for e2 in range(1, elementsCountUp):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                nodeIdentifiers = [e0 + e2,
                                   e0 + 1 + delQ + (delX) * (e2 - 1),
                                   e0 + e2 + 1,
                                   e0 + 1 + delQ + (delX) * (e2 - 1) + delX,
                                   e0 + e2 + delM,
                                   e0 + 1 + delA + (delX) * (e2 - 1),
                                   e0 + e2 + delM + 1,
                                   e0 + 1 + delB + (delX) * (e2 - 1)
                                   ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

            for e3 in range(1, elementsCountlateral - 1):
                for e2 in range(1, elementsCountUp):
                    e22 = e0 + 1 + delQ + (delX) * (e2 - 1) + (e3 - 1)
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    nodeIdentifiers = [e22,
                                       e22 + 1,
                                       e22 + delX,
                                       e22 + delX + 1,
                                       e22 + (elementsCountlateral - 1),
                                       e22 + (elementsCountlateral),
                                       e22 + (3 * elementsCountlateral - 3),
                                       e22 + (3 * elementsCountlateral - 3) + 1
                                       ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            for e2 in range(1, elementsCountUp):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                e22 = 1 + delQ + (delX) * (e2 - 1)
                nodeIdentifiers = [e0 + 1 + delP + (delX) * (e2 - 1),
                                   e0 + (e2) + elementsCountUp,
                                   e0 + 1 + delP + (delX) * (e2),
                                   e0 + (e2) + elementsCountUp + 1,
                                   e0 + 1 + delP + (delX) * (e2 - 1) + (elementsCountlateral - 1),
                                   e0 + e2 + (3 * elementsCountUp),
                                   e0 + 1 + delP + (delX) * (e2) + (elementsCountlateral - 1),
                                   e0 + e2 + (3 * elementsCountUp) + 1
                                   ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1



        fm.endChange()
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

        nSegmenttemporary = 0
        # Unpack radius and rate of change of inner radius
        startRadiusParent = self._innerRadiusParentSegmentList[nSegmenttemporary]
        startRadiusParentDerivative = self._dInnerRadiusParentSegmentList[nSegmenttemporary]
        endRadiusParent = self._innerRadiusParentSegmentList[nSegmenttemporary+1]
        endRadiusParentDerivative = self._dInnerRadiusParentSegmentList[nSegmenttemporary+1]

        startRadiusDaugh1 = self._innerRadiusDaugh1SegmentList[nSegmenttemporary]
        startRadiusDaugh1Derivative = self._dInnerRadiusDaugh1SegmentList[nSegmenttemporary]
        endRadiusDaugh1 = self._innerRadiusDaugh1SegmentList[nSegmenttemporary+1]
        endRadiusDaugh1Derivative = self._dInnerRadiusDaugh1SegmentList[nSegmenttemporary+1]

        startRadiusDaugh2 = self._innerRadiusDaugh2SegmentList[nSegmenttemporary]
        startRadiusDaugh2Derivative = self._dInnerRadiusDaugh2SegmentList[nSegmenttemporary]
        endRadiusDaugh2 = self._innerRadiusDaugh2SegmentList[nSegmenttemporary+1]
        endRadiusDaugh2Derivative = self._dInnerRadiusDaugh2SegmentList[nSegmenttemporary+1]

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
                                          endRadiusDaugh2Derivative, nSegment,
                                          self._startPhase)

        # startIdx = 0 if nSegment == 0 else 1

        startIdx = 0
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
                                nSegment, startPhase):
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

        ##HARI COMMENTED _ IMPORTANT NOTE
        # if n2 == elementsCountAlongSegment:
        #     z = parentsegmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * parentsegmentLength

        xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [0.0, radius, 0.0], [-radius, 0.0, 0.0],
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

        ##HARI COMMENTED _ IMPORTANT NOTE
        # if n2 == elementsCountAlongSegment:
        #     z = daugh1segmentLength / elementsCountAlongSegment * (1.0*n2-0.25) + startPhase / 360.0 * daugh1segmentLength


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
    # lengthToFirstPhase = startPhase / 360.0 * daugh1segmentLength
    # for n2 in range(elementsCountAlongSegment + 1):
    #     faceDaughter1MidPointsZ += [lengthToFirstPhase +
    #                              n2 * daugh1segmentLength / elementsCountAlongSegment]

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

        ##HARI COMMENTED _ IMPORTANT NOTE
        # if (n2 == elementsCountAlongSegment):
        #     z = daugh2segmentLength / elementsCountAlongSegment * (1.0*n2 - 0.25) + startPhase / 360.0 * daugh2segmentLength

        # xLoop, d1Loop = createCirclePoints([-z, 0.0, 0.0], [0.0, -radius, 0.0], [0.0, 0.0, radius],
        #                                    elementsCountAround, startRadians=0.0)
        xLoop, d1Loop = createCirclePoints([-z, 0.0, 0.0], [0.0, radius, 0.0], [0.0, 0.0, -radius],
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
    lengthToFirstPhase = startPhase / 360.0 * daugh2segmentLength
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



