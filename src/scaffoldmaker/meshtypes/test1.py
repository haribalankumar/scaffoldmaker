if nodeIdentifier == 109:
    interpolatedNodes = []
    interpolatedNodes_d1 = []
    interpolatedNodes_d2 = []
    for n2 in range(elementsCountUpNeck + 1):
        xi = 1.0 - (ellipsoidal_x[m1][2] - z_bottom) / delta_z
        for n1 in range(elementsCountAround):
            phi_inner, _, phi_outer, _ = getCubicHermiteBasis(xi)
            x = [(phi_inner * tube_x[m1][c] + phi_outer * ellipsoidal_x[m1][c]) for c in range(3)]
            d1 = [(phi_inner * tube_d1[m1][c] + phi_outer * ellipsoidal_d1[m1][c]) for c in
                  range(3)]
            d2 = [(phi_inner * tube_d2[m1][c] + phi_outer * ellipsoidal_d2[m1][c]) for c in
                  range(3)]
            interpolatedNodes.append(x)
            interpolatedNodes_d1.append(d1)
            interpolatedNodes_d2.append(d2)
            m1 += 1

    # for iann in range(3):
    #    startPointsx[0][0][iann] = x[iann]
    #    endPointsx[0][0][iann] = x[iann] + axis1[iann]*daughterlength
    #    startPointsd1[0][0][iann] = dx_ds1[iann]
    #    startPointsd2[0][0][iann] = dx_ds2[iann]
    #    startPointsd3[0][0][iann] = dx_ds3[iann]
    startPointsx[0][0].append(x)
    endPointsx[0][0].append(x + axis1 * daughterlength)
    startNodeId[0][0] = nodeIdentifier
if nodeIdentifier == 37:
    # for iann in range(3):
    #    startPointsx[1][0][iann] = x[iann]
    #    endPointsx[1][0][iann] = endPointsx[0][0][iann] + delX
    #    startPointsd1[1][0][iann] = dx_ds1[iann]
    #    startPointsd2[1][0][iann] = dx_ds2[iann]
    #    startPointsd3[1][0][iann] = dx_ds3[iann]
    startPointsx[1][0].append(x)
    endPointsx[1][0].append(x + axis1 * daughterlength + delX)
    startNodeId[1][0] = nodeIdentifier
if nodeIdentifier == 110:
    # for iann in range(3):
    #    startPointsx[0][1][iann] = x[iann]
    #    endPointsx[0][1][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][1][iann] = dx_ds1[iann]
    #    startPointsd2[0][1][iann] = dx_ds2[iann]
    #    startPointsd3[0][1][iann] = dx_ds3[iann]
    startPointsx[0][1].append(x)
    endPointsx[0][1].append(x + axis1 * daughterlength)
    startNodeId[0][1] = nodeIdentifier
if nodeIdentifier == 38:
    # for iann in range(3):
    #    startPointsx[1][1][iann] = x[iann]
    #    endPointsx[1][1][iann] = endPointsx[0][1][iann] + delX
    #    startPointsd1[1][1][iann] = dx_ds1[iann]
    #    startPointsd2[1][1][iann] = dx_ds2[iann]
    #    startPointsd3[1][1][iann] = dx_ds3[iann]
    startPointsx[1][1].append(x)
    endPointsx[1][1].append(x + axis1 * daughterlength + delX)
    startNodeId[1][1] = nodeIdentifier
if nodeIdentifier == 111:
    # for iann in range(3):
    #    startPointsx[0][2][iann] = x[iann]
    #    endPointsx[0][2][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][2][iann] = dx_ds1[iann]
    #    startPointsd2[0][2][iann] = dx_ds2[iann]
    #    startPointsd3[0][2][iann] = dx_ds3[iann]
    startPointsx[0][2].append(x)
    endPointsx[0][2].append(x + axis1 * daughterlength)
    startNodeId[0][2] = nodeIdentifier
if nodeIdentifier == 39:
    # for iann in range(3):
    #    startPointsx[1][2][iann] = x[iann]
    #    endPointsx[1][2][iann] = endPointsx[0][2][iann] + delX
    #    startPointsd1[1][2][iann] = dx_ds1[iann]
    #    startPointsd2[1][2][iann] = dx_ds2[iann]
    #    startPointsd3[1][2][iann] = dx_ds3[iann]
    startPointsx[1][2].append(x)
    endPointsx[1][2].append(x + axis1 * daughterlength + delX)
    startNodeId[1][2] = nodeIdentifier
if nodeIdentifier == 112:
    # for iann in range(3):
    #    startPointsx[0][3][iann] = x[iann]
    #    endPointsx[0][3][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][3][iann] = dx_ds1[iann]
    #    startPointsd2[0][3][iann] = dx_ds2[iann]
    #    startPointsd3[0][3][iann] = dx_ds3[iann]
    startPointsx[0][3].append(x)
    endPointsx[0][3].append(x + axis1 * daughterlength)
    startNodeId[0][3] = nodeIdentifier
if nodeIdentifier == 40:
    # for iann in range(3):
    #    startPointsx[1][3][iann] = x[iann]
    #    endPointsx[1][3][iann] = endPointsx[0][3][iann] + delX
    #    startPointsd1[1][3][iann] = dx_ds1[iann]
    #    startPointsd2[1][3][iann] = dx_ds2[iann]
    #    startPointsd3[1][3][iann] = dx_ds3[iann]
    startPointsx[1][3].append(x)
    endPointsx[1][3].append(x + axis1 * daughterlength + delX)
    startNodeId[1][3] = nodeIdentifier
if nodeIdentifier == 120:
    # for iann in range(3):
    #    startPointsx[0][4][iann] = x[iann]
    #    endPointsx[0][4][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][4][iann] = dx_ds1[iann]
    #    startPointsd2[0][4][iann] = dx_ds2[iann]
    #    startPointsd3[0][4][iann] = dx_ds3[iann]
    startNodeId[0][4] = nodeIdentifier
if nodeIdentifier == 48:
    # for iann in range(3):
    #    startPointsx[1][4][iann] = x[iann]
    #    endPointsx[1][4][iann] = endPointsx[1][4][iann] + delX
    #    startPointsd1[1][4][iann] = dx_ds1[iann]
    #    startPointsd2[1][4][iann] = dx_ds2[iann]
    #    startPointsd3[1][4][iann] = dx_ds3[iann]
    startNodeId[1][4] = nodeIdentifier
if nodeIdentifier == 128:
    # for iann in range(3):
    #    startPointsx[0][5][iann] = x[iann]
    #    endPointsx[0][5][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][5][iann] = dx_ds1[iann]
    #    startPointsd2[0][5][iann] = dx_ds2[iann]
    #    startPointsd3[0][5][iann] = dx_ds3[iann]
    startNodeId[0][5] = nodeIdentifier
if nodeIdentifier == 56:
    # for iann in range(3):
    #    startPointsx[1][5][iann] = x[iann]
    #    endPointsx[1][5][iann] = endPointsx[0][5][iann] + delX
    #    startPointsd1[1][5][iann] = dx_ds1[iann]
    #    startPointsd2[1][5][iann] = dx_ds2[iann]
    #    startPointsd3[1][5][iann] = dx_ds3[iann]
    startNodeId[1][5] = nodeIdentifier
if nodeIdentifier == 127:
    # for iann in range(3):
    #    startPointsx[0][6][iann] = x[iann]
    #    endPointsx[0][6][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][6][iann] = dx_ds1[iann]
    #    startPointsd2[0][6][iann] = dx_ds2[iann]
    #    startPointsd3[0][6][iann] = dx_ds3[iann]
    startNodeId[0][6] = nodeIdentifier
if nodeIdentifier == 55:
    # for iann in range(3):
    #    startPointsx[1][6][iann] = x[iann]
    #    endPointsx[1][6][iann] = endPointsx[0][6][iann] + delX
    #    startPointsd1[1][6][iann] = dx_ds1[iann]
    #    startPointsd2[1][6][iann] = dx_ds2[iann]
    #    startPointsd3[1][6][iann] = dx_ds3[iann]
    startNodeId[1][6] = nodeIdentifier
if nodeIdentifier == 126:
    # for iann in range(3):
    #    startPointsx[0][7][iann] = x[iann]
    #    endPointsx[0][7][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][7][iann] = dx_ds1[iann]
    #    startPointsd2[0][7][iann] = dx_ds2[iann]
    #    startPointsd3[0][7][iann] = dx_ds3[iann]
    startNodeId[0][7] = nodeIdentifier
if nodeIdentifier == 54:
    # for iann in range(3):
    #    startPointsx[1][7][iann] = x[iann]
    #    endPointsx[1][7][iann] = endPointsx[0][7][iann] + delX
    #    startPointsd1[1][7][iann] = dx_ds1[iann]
    #    startPointsd2[1][7][iann] = dx_ds2[iann]
    #    startPointsd3[1][7][iann] = dx_ds3[iann]
    startNodeId[1][7] = nodeIdentifier
if nodeIdentifier == 125:
    # for iann in range(3):
    #    startPointsx[0][8][iann] = x[iann]
    #    endPointsx[0][8][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][8][iann] = dx_ds1[iann]
    #    startPointsd2[0][8][iann] = dx_ds2[iann]
    #    startPointsd3[0][8][iann] = dx_ds3[iann]
    startNodeId[0][8] = nodeIdentifier
if nodeIdentifier == 53:
    # for iann in range(3):
    #    startPointsx[1][8][iann] = x[iann]
    #    endPointsx[1][8][iann] = endPointsx[0][8][iann]  + delX
    #    startPointsd1[1][8][iann] = dx_ds1[iann]
    #    startPointsd2[1][8][iann] = dx_ds2[iann]
    #    startPointsd3[1][8][iann] = dx_ds3[iann]
    startNodeId[1][8] = nodeIdentifier
if nodeIdentifier == 117:
    # for iann in range(3):
    #    startPointsx[0][9][iann] = x[iann]
    #    endPointsx[0][9][iann] = x[iann] + axis1[iann] * daughterlength
    #    startPointsd1[0][9][iann] = dx_ds1[iann]
    #    startPointsd2[0][9][iann] = dx_ds2[iann]
    #    startPointsd3[0][9][iann] = dx_ds3[iann]
    startNodeId[0][9] = nodeIdentifier
if nodeIdentifier == 45:
    # for iann in range(3):
    #    startPointsx[1][9][iann] = x[iann]
    #    endPointsx[1][9][iann] = endPointsx[0][9][iann] + delX
    #    startPointsd1[1][9][iann] = dx_ds1[iann]
    #    startPointsd2[1][9][iann] = dx_ds3[iann]
    #    startPointsd3[1][9][iann] = dx_ds4[iann]
    startNodeId[1][9] = nodeIdentifier


============================================================
        # create annulus mesh around ostium
        #####==============================
        endPoints1_x = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d1 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endPoints1_d2 = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endNode1_Id = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]
        endDerivativesMap = [[None] * elementsCountAroundOstium, [None] * elementsCountAroundOstium]

#       nodeCountsEachWallLayer = (elementsCountUpNeck + elementsCountUpBody) * elementsCountAround - 1

        endNode1_Id[0][0] = 109
        endNode1_Id[0][1] = 110
        endNode1_Id[0][2] = 111
        endNode1_Id[0][3] = 112
        endNode1_Id[0][4] = 117
        endNode1_Id[0][5] = 120
        endNode1_Id[0][6] = 125
        endNode1_Id[0][7] = 126
        endNode1_Id[0][8] = 127
        endNode1_Id[0][9] = 128

        endNode1_Id[1][0] = 37
        endNode1_Id[1][1] = 38
        endNode1_Id[1][2] = 39
        endNode1_Id[1][3] = 40
        endNode1_Id[1][4] = 45
        endNode1_Id[1][5] = 48
        endNode1_Id[1][6] = 53
        endNode1_Id[1][7] = 54
        endNode1_Id[1][8] = 55
        endNode1_Id[1][9] = 56

        for nc1 in range(9):
            endPoints1_x[0][nc1] = innerLayer_x[nc1]
            endPoints1_d1[0][nc1] = innerLayer_d1[nc1]
            endPoints1_d2[0][nc1] = innerLayer_d2[nc1]

        for nc1 in range(9):
            endPoints1_x[1][nc1] = outerLayer_x[nc1]
            endPoints1_d1[1][nc1] = outerLayer_d1[nc1]
            endPoints1_d2[1][nc1] = outerLayer_d2[nc1]


        #for n1 in range(elementsCountAroundOstium):
        #    if n1 == 0:
        #        endDerivativesMap[0][n1] = ((-1, 0, 0), (-1, -1, 0), None, (0, 1, 0))
        #        endDerivativesMap[1][n1] = ((-1, 0, 0), (-1, -1, 0), None, (0, 1, 0))
        #    elif n1 == 1:
        #        endDerivativesMap[0][n1] = ((0, 1, 0), (-1, 0, 0), None)
        #        endDerivativesMap[1][n1] = ((0, 1, 0), (-1, 0, 0), None)
        #    elif n1 == 2:
        #        endDerivativesMap[0][n1] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
        #        endDerivativesMap[1][n1] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
        #    elif n1 == 3:
        #        endDerivativesMap[0][n1] = ((1, 0, 0), (0, 1, 0), None)
        #        endDerivativesMap[1][n1] = ((1, 0, 0), (0, 1, 0), None)
        #    elif n1 == 4:
        #        endDerivativesMap[0][n1] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
        #        endDerivativesMap[1][n1] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
        #    elif n1 == 5:
        #        endDerivativesMap[0][n1] = ((0, -1, 0), (1, 0, 0), None)
        #        endDerivativesMap[1][n1] = ((0, -1, 0), (1, 0, 0), None)
        #    elif n1 == 6:
        #        endDerivativesMap[0][n1] = ((0, -1, 0), (1, -1, 0), None, (-1, 0, 0))
        #        endDerivativesMap[1][n1] = ((0, -1, 0), (1, -1, 0), None, (-1, 0, 0))
        #    else:
        #        endDerivativesMap[0][n1] = ((-1, 0, 0), (0, -1, 0), None)
        #        endDerivativesMap[1][n1] = ((-1, 0, 0), (0, -1, 0), None)

        nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
            nodes, mesh, nodeIdentifier, elementIdentifier,
            o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
            endPoints1_x, endPoints1_d1, endPoints1_d2, None, endNode1_Id, None,
            elementsCountRadial=elementsCountAnnulusRadially, meshGroups=[])
