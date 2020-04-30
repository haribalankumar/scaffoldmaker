'''
Utility function for generating tubular mesh from a central line
using a segment profile.
'''
from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

##from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector

def warpAirwaySegmentPoints(x1ListParent, x1ListDaugh1, x1ListDaugh2,
                            d1ListParent, d1ListDaugh1, d1ListDaugh2,
                            d2ListParent, d2ListDaugh1, d2ListDaugh2,
                            segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
                            parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
                            sxParent, sxDaugh1, sxDaugh2,
                            sd1Parent, sd1Daugh1, sd1Daugh2,
                            sd2Parent, sd2Daugh1, sd2Daugh2,
                            elementsCountAround, elementsCountAlongSegment,
                            nSegment,
                            ParentfaceMidPointZ, Daughter1faceMidPointZ, Daughter2faceMidPointZ):
    """
    Warps points in Airway segment to account for bending and twisting
    along central path defined by nodes sx and derivatives sd1 and sd2.
    :param x1List: coordinates of segment points.
    :param d1List: derivatives around axis of segment.
    :param d2List: derivatives along axis of segment.
    :param segmentAxis: axis perpendicular to segment plane.
    :param sx: coordinates of points on central path.
    :param sd1: derivatives of points along central path.
    :param sd2: derivatives representing cross axes.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :param faceMidPointZ: z-coordinate of midpoint for each element
    groups along the segment.
    :return coordinates and derivatives of warped points.
    """
    x1ParentWarpedList = []
    x1Daugh1WarpedList = []
    x1Daugh2WarpedList = []

    d1ParentWarpedList = []
    d1Daugh1WarpedList = []
    d1Daugh2WarpedList = []

    d2ParentWarpedList = []
    d2Daugh1WarpedList = []
    d2Daugh2WarpedList = []

    #smoothd2WarpedList = []
    d3ParentWarpedUnitList = []
    d3Daugh1WarpedUnitList = []
    d3Daugh2WarpedUnitList = []

    #parent
    #######
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        xMid = [0.0, 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Parent[n2])
        cp = vector.crossproduct3(segmentAxisParent, unitTangent)
        dp = vector.dotproduct(segmentAxisParent, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisParent, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                axisRot = [1.0, 0, 0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sxParent[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Parent[n2]
                    pt = [midRot[j] + sd2Parent[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Parent[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1ParentWarpedList.append(xTranslate)
            d1ParentWarpedList.append(d1Rot2)
            d2ParentWarpedList.append(d2Rot2)


    #Daughter1
    ##########
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        print('rotating daught1 to ', sd1Daugh1[n2])

        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

        xMid = [Daughter1faceMidPointZ[nAlongSegment], 0.0, 0.0]
        #xMid = [Daughter1faceMidPointZ[nAlongSegment], 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Daugh1[n2])
        cp = vector.crossproduct3(segmentAxisDaughter1, unitTangent)
        dp = vector.dotproduct(segmentAxisDaughter1, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisDaughter1, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            print('daughter1 warping')
        else: # path tangent parallel to segment axis (z-axis)
            print('daughter1 warping parallel to zxia')
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                #axisRot = [1.0, 0, 0]
                axisRot = [0.0, 0, 1.0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sxDaugh1[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Daugh1[n2]
                    pt = [midRot[j] + sd2Daugh1[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Daugh1[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1Daugh1WarpedList.append(xTranslate)
            d1Daugh1WarpedList.append(d1Rot2)
            d2Daugh1WarpedList.append(d2Rot2)


    #Daughter2
    ############
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

        xMid = [-Daughter2faceMidPointZ[nAlongSegment], 0.0, 0.0]
        #xMid = [-Daughter2faceMidPointZ[nAlongSegment], 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Daugh2[n2])
        cp = vector.crossproduct3(segmentAxisDaughter2, unitTangent)
        dp = vector.dotproduct(segmentAxisDaughter2, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisDaughter2, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                axisRot = [1.0, 0, 0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sxDaugh2[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Daugh2[n2]
                    pt = [midRot[j] + sd2Daugh2[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Daugh2[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1Daugh2WarpedList.append(xTranslate)
            d1Daugh2WarpedList.append(d1Rot2)
            d2Daugh2WarpedList.append(d2Rot2)

    # # Smooth d2 for segment
    # ########################
    # smoothd2Raw = []
    # for n1 in range(elementsCountAround):
    #     nx = []
    #     nd2 = []
    #     for n2 in range(elementsCountAlongSegment + 1):
    #         n = n2*elementsCountAround + n1
    #         nx.append(xWarpedList[n])
    #         nd2.append(d2WarpedList[n])
    #     smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    #     smoothd2Raw.append(smoothd2)
    #
    # # Re-arrange smoothd2
    # for n2 in range(elementsCountAlongSegment + 1):
    #     for n1 in range(elementsCountAround):
    #         smoothd2WarpedList.append(smoothd2Raw[n1][n2])

    # Calculate unit d3 Parent
    ###########################
    for n in range(len(x1ParentWarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1ParentWarpedList[n]),
                                                       vector.normalise(d2ParentWarpedList[n])))
        d3ParentWarpedUnitList.append(d3Unit)

    # Calculate unit d3 - Daughter1
    ##################################
    for n in range(len(x1Daugh1WarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1Daugh1WarpedList[n]),
                                                       vector.normalise(d2Daugh1WarpedList[n])))
        d3Daugh1WarpedUnitList.append(d3Unit)

    # Calculate unit d3 Daughter2
    #################################
    for n in range(len(x1Daugh2WarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1Daugh2WarpedList[n]),
                                                       vector.normalise(d2Daugh2WarpedList[n])))
        d3Daugh2WarpedUnitList.append(d3Unit)

    return x1ParentWarpedList, x1Daugh1WarpedList, x1Daugh2WarpedList, \
           d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList, \
           d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList, \
           d3ParentWarpedUnitList, d3Daugh1WarpedUnitList, d3Daugh2WarpedUnitList


def getAirwaySegmentCoordinatesFromInner(xInner, d1Inner, d2Inner, d3Inner,
    wallThicknessList, elementsCountAround,
    elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Generates coordinates from inner to outer surface using coordinates
    and derivatives of inner surface.
    :param xInner: Coordinates on inner surface
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3Inner: Derivatives on inner surface through wall
    :param wallThicknessList: Wall thickness for each element along tube
    :param elementsCountAround: Number of elements around tube
    :param elementsCountAlong: Number of elements along tube
    :param elementsCountThroughWall: Number of elements through tube wall
    :param transitElementList: stores true if element around is a transition
    element that is between a big and a small element.
    return nodes and derivatives for mesh, and curvature along inner surface.
    """

    xOuter = []
    curvatureAroundInner = []
    curvatureAlong = []
    curvatureList = []
    xList = []
    d1List = []
    d2List = []
    d3List = []

    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Inner[n]
            # Calculate outer coordinates
            x = [xInner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)
            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[n + elementsCountAround], d2Inner[n + elementsCountAround], vector.normalise(d3Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xInner[n - elementsCountAround], d2Inner[n - elementsCountAround], xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xInner[n - elementsCountAround], d2Inner[n - elementsCountAround], xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[n + elementsCountAround], d2Inner[n + elementsCountAround], vector.normalise(d3Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3Inner[n]
                innerx = xInner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xList.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Inner[n]]
                d1List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xInner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Inner[n]]
                d2List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3List.append(d3)

    return xList, d1List, d2List, d3List, curvatureList



def createSurfaceNodesAndElements(region,
                                  xParent, d1Parent, d2Parent,
                                  xDaugh1, d1Daugh1, d2Daugh1,
                                  xDaugh2, d1Daugh2, d2Daugh2,
                                  elementsCountAround, elementsCountAlong,
                                  firstNodeIdentifier, firstElementIdentifier,
                                  useCrossDerivatives):
    #    annotationGroups, annotationArray,
    """
    :param xList: coordinates of centerline points.
    :param d1List: derivatives along axis of segment.
    :param radius1List: derivatives along axis of segment.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [ 0.0, 0.0, 0.0 ]

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
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)


    #    mesh = fm.findMeshByDimension(3)
#    eftfactory = eftfactory_bicubichermite(mesh, useCrossDerivatives)
#    eft = eftfactory.createEftBasic()

    mesh = fm.findMeshByDimension(2)
    bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    if not useCrossDerivatives:
        for n in range(4):
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    cache = fm.createFieldcache()

    # radiansPerElementAround = 2.0 * math.pi / elementsCountAround
    # x = [0.0, 0.0, 0.0]
    # dx_ds1 = [0.0, 0.0, 0.0]
    # dx_ds2 = [0.0, 0.0, 1.0 / elementsCountAlong]
    # zero = [0.0, 0.0, 0.0]

    # #parent
    # radius = radius1list[0]
    # for n2 in range(elementsCountAlong + 1):
    #     x[2] = n2 / elementsCountAlong
    #     for n1 in range(elementsCountAround):
    #         radiansAround = n1 * radiansPerElementAround
    #         cosRadiansAround = math.cos(radiansAround)
    #         sinRadiansAround = math.sin(radiansAround)
    #         x[0] = radius * cosRadiansAround
    #         x[1] = radius * sinRadiansAround
    #         dx_ds1[0] = radiansPerElementAround * radius * -sinRadiansAround
    #         dx_ds1[1] = radiansPerElementAround * radius * cosRadiansAround
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
    #         if useCrossDerivatives:
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
    #         nodeIdentifier = nodeIdentifier + 1

    # # Central path
    # cx = [ [ 0.0, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
    # cd1 = [ [ segmentLength, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
    # cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
    # cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
    #
    # # Sample central path
    # sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    # sd2 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)[0]

    # ###daughter 1
    # radius = radius1list[1]
    # for n2 in range(elementsCountAlong + 1):
    #     x[2] = n2 / elementsCountAlong
    #     for n1 in range(elementsCountAround):
    #         radiansAround = n1 * radiansPerElementAround
    #         cosRadiansAround = math.cos(radiansAround)
    #         sinRadiansAround = math.sin(radiansAround)
    #         x[0] = radius * cosRadiansAround
    #         x[1] = radius * sinRadiansAround
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #         if useCrossDerivatives:
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
    #         nodeIdentifier = nodeIdentifier + 1

    print('coming to write nodes in create surface node')

    # Create nodes - PARENT
    # Coordinates field
    for n in range(len(xParent)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xParent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Parent[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes - DauGH1
    # Coordinates field
    for n in range(len(xDaugh1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh1[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes - DAUGH2
    # Coordinates field
    for n in range(len(xDaugh2)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh2[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1


    # # create elements
    # ##################
    # for e2 in range(elementsCountAlong):
    #     for e1 in range(elementsCountAround):
    #         element = mesh.createElement(elementIdentifier, elementtemplate)
    #         bni1 = e2*elementsCountAround + e1 + 1
    #         bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
    #         nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
    #         result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #         elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return nodeIdentifier,elementIdentifier

