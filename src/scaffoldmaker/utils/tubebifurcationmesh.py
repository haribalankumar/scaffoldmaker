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

from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite

from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector

def createjunctionAirwaySegmentPoints(
        xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
        d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
        d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
        segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
        parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
        sxparent, sxDaugh1, sxDaugh2,
        sd1parent, sd1Daugh1, sd1Daugh2,
        sd2parent, sd2Daugh1, sd2Daugh2,
        elementsCountAround, elementsCountAlongSegment, nSegment):

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
    groups along the segment.
    :return coordinates and derivatives of warped points.
    """
    xjunctionOuterList = []
    xjunctionInnerList = []
    d2junctionList = []
    d1junctionList = []

    n2 = elementsCountAlongSegment * nSegment + elementsCountAlongSegment

    #JUNCTION OUTER
    ##################
    xParentAlongSegment = xParentWarpedList[9]
    xDaugh1AlongSegment = xDaugh1WarpedList[2]

    xjunction1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j])/2.0 for j in range(3)]
    xjunctionOuterList.append(xjunction1)

    xParentAlongSegment = xParentWarpedList[11]
    xDaugh2AlongSegment = xDaugh2WarpedList[2]
    xjunction1 = [(xParentAlongSegment[j]+xDaugh2AlongSegment[j])/2.0 for j in range(3)]
    xjunctionOuterList.append(xjunction1)

    xDaugh1AlongSegment = xDaugh1WarpedList[0]
    xDaugh2AlongSegment = xDaugh2WarpedList[0]
    xjunction1 = [(xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j]+sxparent[n2][j])/3.0 for j in range(3)]
    xjunctionOuterList.append(xjunction1)

    #JUNCTION INNER
    ##################

    xParentAlongSegment = xParentWarpedList[8]
    xDaugh1AlongSegment = xDaugh1WarpedList[1]
    xDaugh2AlongSegment = xDaugh2WarpedList[1]

    xjunction1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunction1)

    xParentAlongSegment = xParentWarpedList[10]
    xDaugh1AlongSegment = xDaugh1WarpedList[3]
    xDaugh2AlongSegment = xDaugh2WarpedList[3]
    xjunction1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunction1)


    return xjunctionOuterList, xjunctionInnerList, d1junctionList, d2junctionList



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


def getAirwaySegmentCoordinatesFromInner(
        xParentInner, xDaugh1Inner, xDaugh2Inner,
        d1ParentInner, d1Daugh1Inner, d1Daugh2Inner,
        d2ParentInner, d2Daugh1Inner, d2Daugh2Inner,
        d3ParentInner, d3Daugh1Inner, d3Daugh2Inner,
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

    xParentList = []
    d1ParentList = []
    d2ParentList = []
    d3ParentList = []

    xDaughter1List = []
    d1Daughter1List = []
    d2Daughter1List = []
    d3Daughter1List = []

    xDaughter2List = []
    d1Daughter2List = []
    d2Daughter2List = []
    d3Daughter2List = []

    #PARENT
    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3ParentInner[n]

            # Calculate outer coordinates
            x = [xParentInner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)

            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xParentInner[prevIdx], d1ParentInner[prevIdx], xParentInner[n], d1ParentInner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xParentInner[n], d1ParentInner[n], xParentInner[nextIdx], d1ParentInner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xParentInner[n], d2ParentInner[n], xParentInner[n + elementsCountAround], d2ParentInner[n + elementsCountAround], vector.normalise(d3ParentInner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xParentInner[n - elementsCountAround], d2ParentInner[n - elementsCountAround], xParentInner[n], d2ParentInner[n], vector.normalise(d3ParentInner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xParentInner[n - elementsCountAround], d2ParentInner[n - elementsCountAround], xParentInner[n], d2ParentInner[n], vector.normalise(d3ParentInner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xParentInner[n], d2ParentInner[n], xParentInner[n + elementsCountAround], d2ParentInner[n + elementsCountAround], vector.normalise(d3ParentInner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3ParentInner[n]
                innerx = xParentInner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xParentList.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1ParentInner[n]]
                d1ParentList.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xParentInner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2ParentInner[n]]
                d2ParentList.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3ParentList.append(d3)

    xOuter = []
    #DAUGHTER1
    for n2 in range(elementsCountAlong + 1):
        print('coord inner for daughter=',n2)
        wallThickness = wallThicknessList[n2]
        print('wall thickness for daughter1 = ', wallThickness)
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Daugh1Inner[n]
            # Calculate outer coordinates
            print('print normals for d1=',norm[0],norm[1],norm[2])
            x = [xDaugh1Inner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)

            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[prevIdx], d1Daugh1Inner[prevIdx], xDaugh1Inner[n], d1Daugh1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[n], d1Daugh1Inner[n], xDaugh1Inner[nextIdx], d1Daugh1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n], d2Daugh1Inner[n], xDaugh1Inner[n + elementsCountAround], d2Daugh1Inner[n + elementsCountAround], vector.normalise(d3Daugh1Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n - elementsCountAround], d2Daugh1Inner[n - elementsCountAround], xDaugh1Inner[n], d2Daugh1Inner[n], vector.normalise(d3Daugh1Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n - elementsCountAround], d2Daugh1Inner[n - elementsCountAround], xDaugh1Inner[n], d2Daugh1Inner[n], vector.normalise(d3Daugh1Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n], d2Daugh1Inner[n], xDaugh1Inner[n + elementsCountAround], d2Daugh1Inner[n + elementsCountAround], vector.normalise(d3Daugh1Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                print('coord inner for daughter=', n3, n1)
                n = n2*elementsCountAround + n1
                norm = d3Daugh1Inner[n]
                innerx = xDaugh1Inner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xDaughter1List.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Daugh1Inner[n]]
                d1Daughter1List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xDaugh1Inner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Daugh1Inner[n]]
                d2Daughter1List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3Daughter1List.append(d3)

    xOuter = []
    #DAUGHTER2
    ###########
    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Daugh2Inner[n]
            # Calculate outer coordinates
            x = [xDaugh2Inner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)
            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[prevIdx], d1Daugh1Inner[prevIdx], xDaugh1Inner[n], d1Daugh1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[n], d1Daugh1Inner[n], xDaugh1Inner[nextIdx], d1Daugh1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n], d2Daugh2Inner[n], xDaugh2Inner[n + elementsCountAround], d2Daugh2Inner[n + elementsCountAround], vector.normalise(d3Daugh2Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n - elementsCountAround], d2Daugh2Inner[n - elementsCountAround], xDaugh2Inner[n], d2Daugh2Inner[n], vector.normalise(d3Daugh2Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n - elementsCountAround], d2Daugh2Inner[n - elementsCountAround], xDaugh2Inner[n], d2Daugh2Inner[n], vector.normalise(d3Daugh2Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n], d2Daugh1Inner[n], xDaugh2Inner[n + elementsCountAround], d2Daugh2Inner[n + elementsCountAround], vector.normalise(d3Daugh2Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3Daugh2Inner[n]
                innerx = xDaugh2Inner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xDaughter2List.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Daugh1Inner[n]]
                d1Daughter2List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xDaugh1Inner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Daugh1Inner[n]]
                d2Daughter2List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3Daughter2List.append(d3)

    return xParentList, d1ParentList, d2ParentList, d3ParentList, \
           xDaughter1List, d1Daughter1List, d2Daughter1List, d3Daughter1List,\
           xDaughter2List, d1Daughter2List, d2Daughter2List, d3Daughter2List,\
           curvatureList


def createAirwaySegmentNodesAndElements\
                (region,
                 xParent, d1Parent, d2Parent, d3Parent,
                 xDaughter1, d1Daughter1, d2Daughter1, d3Daughter1,
                 xDaughter2, d1Daughter2, d2Daughter2, d3Daughter2,
                 xjunctionOuter, xjunctionInner, d1junction, d2junction,
                 elementsCountAround, elementsCountAlong,
                 firstNodeIdentifier, firstElementIdentifier,
                 useCubicHermiteThroughWall,
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
    zero = [0.0, 0.0, 0.0]

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()

    # Coordinates field
    coordinates = findOrCreateFieldCoordinates(fm)
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    if useCubicHermiteThroughWall:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

    mesh = fm.findMeshByDimension(3)

    eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)


    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)

    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    ###################
    # Create nodes
    ##################

    # Create nodes for Parent
    # Coordinates field
    print('length of parent with thickness = ',len(xParent))
    for n in range(len(xParent)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xParent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Parent[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes for Daughter1
    # Coordinates field
    print('length of parent with thickness = ',len(xDaughter1))

    for n in range(len(xDaughter1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Daughter1[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes for Daughter2
    # Coordinates field
    for n in range(len(xDaughter2)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Daughter2[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    fm.endChange()

    return nodeIdentifier, elementIdentifier



def createAirwaySegmentSurfaceNodesAndElements(region,
                                  xParent, d1Parent, d2Parent,
                                  xDaugh1, d1Daugh1, d2Daugh1,
                                  xDaugh2, d1Daugh2, d2Daugh2,
                                  xjunctionOuter, xjunctionInner, d1junction, d2junction,
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
    zero = [0.0, 0.0, 0.0]

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
    print('coming to write nodes in create surface node')

    #####################
    # CREATE NODES
    #####################
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
        print('daughter2 nodes list', nodeIdentifier)

    # Create nodes - JUNCTION
    nodeIdentifierOuter = []
    nodeIdentifierInner = []

    for n in range(len(xjunctionOuter)):
        nodeIdentifierOuter.append(nodeIdentifier)

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    for n in range(len(xjunctionInner)):
        nodeIdentifierInner.append(nodeIdentifier)

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1
        print('adding junction inner to a list', nodeIdentifier)

    #####################
    # CREATE ELEMENTS
    #####################
    # # create elements - Parent
    # ##################
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + 1
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter1
    # ##############################
    node0 = (elementsCountAround) * (elementsCountAlong + 1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + node0
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter2
    # ##########################3##
    node0 = 2 * (elementsCountAround) * (elementsCountAlong + 1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + node0
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    # # create elements - junction
    # ##########################3##
    # for e2 in range(len(xjunctionInner)):
    #     nodeinner = nodeIdentifierInner[e2]
    #     bnp1 = (elementsCountAround)*(elementsCountAlong)+1
    #     bnd1 = (elementsCountAround)*(elementsCountAlong+1)+1
    #     bnd2 = 2*(elementsCountAround)*(elementsCountAlong+1)+1
    #
    #     bndparentouter = bnp1 if e2==0 else bnp1+2
    #     bnddaugh1outer = bnd1+1 if e2==0 else bnd1+3
    #     bnddaugh2outer = bnd2+1 if e2==0 else bnd2+3
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [bndparentouter, 10, nodeinner, 37]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [39, nodeinner, 13, bnddaugh1outer]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [12, bndparentouter, 38, nodeinner]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [nodeinner, 38, bnddaugh2outer, 27]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [39, nodeinner, 25, bnddaugh2outer]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [nodeinner, 37, bnddaugh1outer, 15]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [9, 10, 40, 37]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [12, 9, 38, 40]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [39, 40, 13, 14]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [40, 37, 14, 15]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [40, 38, 26, 27]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [39, 40, 25, 26]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        #JUNCTION BACK
        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [11, 12, 41, 38]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [10, 11, 37, 41]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [41, 39, 28, 25]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [38, 41, 27, 28]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [41, 39, 16, 13]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [37, 41, 15, 16]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1


    # node0 = 3 * (elementsCountAround) * (elementsCountAlong + 1) + 1
    # for e2 in range(elementsCountAlong):
    #     for e1 in range(elementsCountAround):
    #         element = mesh.createElement(elementIdentifier, elementtemplate)
    #         bni1 = e2 * elementsCountAround + e1 + node0
    #         bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
    #         nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
    #         result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #         elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return nodeIdentifier, elementIdentifier




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


    # # create elements - Parent
    # ##################
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + 1
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    # # create elements - Daughter1
    # ##############################
    node0 = (elementsCountAround)*(elementsCountAlong+1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + node0
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + node0
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter2
    # ##########################3##
    node0 = 2*(elementsCountAround)*(elementsCountAlong+1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + node0
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + node0
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    fm.endChange()

    return nodeIdentifier,elementIdentifier

