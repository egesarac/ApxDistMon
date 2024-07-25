'''
@author: ANIK
'''

from multiprocessing import freeze_support
import multiprocessing
import time

from z3 import *
from z3.z3 import ForAll

# set_param('parallel.enable', True)
set_option(rational_to_decimal = True)
mult = 20

def prog_tp2(eps, segCount, data_0, data_1):

    # pad = min(2 * eps, len(data_0) - 1)
    pad = 2
    peps = pad * eps

    # initialize signal duration and segment duration
    sigDur = len(data_0) - 1
    segDur = sigDur / segCount

    # check for edge cases
    t0 = data_0[0][0]
    t1 = data_0[1][0]

    if sigDur / segCount < t1:

        segCount = sigDur / t1

    if t0 != 0:
        return

    var_0 = {0 : True}
    for i in range(1, len(data_0)):
        if (data_0[i][1] != data_0[i - 1][1]):
            var_0[int(data_0[i][0])] = True
    var_0[int(data_0[-1][0])] = True

    var_1 = {0 : True}
    for i in range(1, len(data_1)):
        if (data_1[i][1] != data_1[i - 1][1]):
            var_1[int(data_1[i][0])] = True
    var_1[int(data_1[-1][0])] = True

    # encoding
    i = 0
    entryFound = True
    flag = True

    while entryFound:

        # Terminate if segment count is reached
        if i == segCount:
            return flag

        # flag to be set True if at least one entry is found in the current iteration
        entryFound = False

        # initialize solver
        s = Solver()

        # calculate upper and lower time bound for current segment
        segmentLowerBound = max(0, int((i * segDur) - eps))
        # segmentUpperBound = int((i + 1) * segDur)
        segmentUpperBound = min(sigDur, int((i + 1) * segDur + eps))

        timestamps_0 = []
        segvar_0 = []
        sig0 = Function("sig0", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_0) - 1:

                s.add(sig0(j) == data_0[j][1])
                for k in range(pad):
                    timestamps_0.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_0[j][0] in var_0 and data_0[j][0] > segmentLowerBound)):
                    segvar_0.append(int(data_0[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_0) > 0 and segvar_0[-1] < segmentUpperBound):
            segvar_0.append(segmentUpperBound)

        timestamps_1 = []
        segvar_1 = []
        sig1 = Function("sig1", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_1) - 1:

                s.add(sig1(j) == data_1[j][1])
                for k in range(pad):
                    timestamps_1.append(j * pad + k)

                # if (j > 0 and data_1[j][1] != data_1[j - 1][1]):
                #     var_1.append(j)

                if (j == segmentLowerBound or (data_1[j][0] in var_1 and data_1[j][0] > segmentLowerBound)):
                    segvar_1.append(int(data_1[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_1) > 0 and segvar_1[-1] < segmentUpperBound):
            segvar_1.append(segmentUpperBound)              

        if not entryFound:

            ##print()
            pass

        i += 1

        # force terminate after one loop
        # entryFound = False

        # global clock to local clock mappings (implicit piece-wise constant interpolation w.r.t. local clock values)
        c0 = Function("c0", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c0(i) >= 0, c0(i) < data_0[-1][0])
                    # Or([c0(i) == (min(timestamps_0[-1]//pad, (max(0, ((i - peps) + j)//pad)))) for j in range(1, 2 * int(peps))])
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c0(min(timestamps_0[-1], max(0, i * pad - peps + j))) == i for j in range(1, 2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                    # for i in range(int(data_0[-1][0]))
                ]
            )
        )

        c1 = Function("c1", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c1(i) >= 0, c1(i) < data_1[-1][0])
                    # Or([c1(i) == (min(timestamps_1[-1]//pad, (max(0, ((i - peps) + j)//pad)))) for j in range(1, 2 * int(peps))])
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c1(min(timestamps_1[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                    # for i in range(int(data_1[-1][0]))
                ]
            )
        )
        
        # retiming preserves signal variability
        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_0[i] <= c0(j), c0(j) < segvar_0[i + 1])
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_0) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_1[i] <= c1(j), c1(j) < segvar_1[i + 1])
                            for j in range(timestamps_1[0], timestamps_1[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_1) - 1)
                ]
            )
        )

        # local clocks are bound by epsilon
        s.add(
            And(
                [
                    And(c0(i) - c1(i) <= eps, c0(i) - c1(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # non-concurrent edges stay non-concurrent ???
        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)])
                        ),
                        And(c0(i) - c1(i) < eps, c0(i) - c1(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        # global clock to local clock mappings are ordered
        s.add(
            And(
                [
                    And(
                        [
                            Implies(i <= j, And(c0(i) <= c0(j), c1(i) <= c1(j)))
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # consistent cut flow
        # _flow = Function("c_flow", IntSort(), RealSort())
        c_flow = Function("c_flow", IntSort(), RealSort())

        # addition
        s.add(
            And(
                [
                    c_flow(i) == (sig0(c0(i)) + sig1(c1(i)))
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # violation check
        v = Int("v")
        s.add(And(v >= timestamps_0[0], v <= timestamps_0[-1]))

        s.add(
            Implies(
                And(v >= timestamps_0[0], v <= timestamps_0[-1]),
                z3Interpolate(c_flow, v) < 10
            )
        )

        # print(s.assertions())
        if s.check() == sat:
            m = s.model()
            # out = "%s %s" % (m[test], m[test2])
            # print(m)
            ##print("unsat in segment", i)
            flag = False
            # terminate after unsat
            return flag

        elif i <= segCount:
            ##print("sat in segment", i)
            flag = True

        s.reset()

    return flag

def prog_tp3(eps, segCount, data_0, data_1, data_2):

    # pad = min(2 * eps, len(data_0) - 1)
    pad = 2
    peps = pad * eps

    # initialize signal duration and segment duration
    sigDur = len(data_0) - 1
    segDur = sigDur / segCount

    # check for edge cases
    t0 = data_0[0][0]
    t1 = data_0[1][0]

    if sigDur / segCount < t1:

        segCount = sigDur / t1

    if t0 != 0:
        return

    var_0 = {0 : True}
    for i in range(1, len(data_0)):
        if (data_0[i][1] != data_0[i - 1][1]):
            var_0[int(data_0[i][0])] = True
    var_0[int(data_0[-1][0])] = True

    var_1 = {0 : True}
    for i in range(1, len(data_1)):
        if (data_1[i][1] != data_1[i - 1][1]):
            var_1[int(data_1[i][0])] = True
    var_1[int(data_1[-1][0])] = True

    var_2 = {0 : True}
    for i in range(1, len(data_2)):
        if (data_2[i][1] != data_2[i - 1][1]):
            var_2[int(data_2[i][0])] = True
    var_2[int(data_2[-1][0])] = True

    # encoding
    i = 0
    entryFound = True
    flag = True

    while entryFound:

        # Terminate if segment count is reached
        if i == segCount:
            return flag

        # flag to be set True if at least one entry is found in the current iteration
        entryFound = False

        # initialize solver
        s = Solver()

        # calculate upper and lower time bound for current segment
        segmentLowerBound = max(0, int((i * segDur) - eps))
        # segmentUpperBound = int((i + 1) * segDur)
        segmentUpperBound = min(sigDur, int((i + 1) * segDur + eps))

        timestamps_0 = []
        segvar_0 = []
        sig0 = Function("sig0", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_0) - 1:

                s.add(sig0(j) == data_0[j][1])
                for k in range(pad):
                    timestamps_0.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_0[j][0] in var_0 and data_0[j][0] > segmentLowerBound)):
                    segvar_0.append(int(data_0[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_0) > 0 and segvar_0[-1] < segmentUpperBound):
            segvar_0.append(segmentUpperBound)

        timestamps_1 = []
        segvar_1 = []
        sig1 = Function("sig1", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_1) - 1:

                s.add(sig1(j) == data_1[j][1])
                for k in range(pad):
                    timestamps_1.append(j * pad + k)

                # if (j > 0 and data_1[j][1] != data_1[j - 1][1]):
                #     var_1.append(j)

                if (j == segmentLowerBound or (data_1[j][0] in var_1 and data_1[j][0] > segmentLowerBound)):
                    segvar_1.append(int(data_1[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_1) > 0 and segvar_1[-1] < segmentUpperBound):
            segvar_1.append(segmentUpperBound)              


        timestamps_2 = []
        segvar_2 = []
        sig2 = Function("sig2", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_2) - 1:

                s.add(sig2(j) == data_2[j][1])
                for k in range(pad):
                    timestamps_2.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_2[j][0] in var_2 and data_2[j][0] > segmentLowerBound)):
                    segvar_2.append(int(data_2[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_2) > 0 and segvar_2[-1] < segmentUpperBound):
            segvar_2.append(segmentUpperBound)

        if not entryFound:

            ##print()
            pass

        i += 1

        # force terminate after one loop
        # entryFound = False

        # global clock to local clock mappings (implicit piece-wise constant interpolation w.r.t. local clock values)
        c0 = Function("c0", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c0(i) >= 0, c0(i) < data_0[-1][0])
                    # Or([c0(i) == (min(timestamps_0[-1]//pad, (max(0, ((i - peps) + j)//pad)))) for j in range(1, 2 * int(peps))])
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c0(min(timestamps_0[-1], max(0, i * pad - peps + j))) == i for j in range(1, 2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                    # for i in range(int(data_0[-1][0]))
                ]
            )
        )

        c1 = Function("c1", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c1(i) >= 0, c1(i) < data_1[-1][0])
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c1(min(timestamps_1[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                ]
            )
        )

        c2 = Function("c2", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c2(i) >= 0, c2(i) < data_2[-1][0])
                    for i in range(timestamps_2[0], timestamps_2[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c2(min(timestamps_2[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                ]
            )
        )
        
        # retiming preserves signal variability
        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_0[i] <= c0(j), c0(j) < segvar_0[i + 1])
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_0) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_1[i] <= c1(j), c1(j) < segvar_1[i + 1])
                            for j in range(timestamps_1[0], timestamps_1[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_1) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_2[i] <= c2(j), c2(j) < segvar_2[i + 1])
                            for j in range(timestamps_2[0], timestamps_2[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_2) - 1)
                ]
            )
        )

        # local clocks are bound by epsilon
        s.add(
            And(
                [
                    And(c0(i) - c1(i) <= eps, c0(i) - c1(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    And(c0(i) - c2(i) <= eps, c0(i) - c2(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    And(c1(i) - c2(i) <= eps, c1(i) - c2(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # non-concurrent edges stay non-concurrent ???
        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)])
                        ),
                        And(c0(i) - c1(i) < eps, c0(i) - c1(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c2(i) == segvar_2[j], c2(i) != c2(i - 1)) for j in range(1, len(segvar_2) - 1)])
                        ),
                        And(c0(i) - c2(i) < eps, c0(i) - c2(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)]),
                            Or([And(c2(i) == segvar_2[j], c2(i) != c2(i - 1)) for j in range(1, len(segvar_2) - 1)])
                        ),
                        And(c1(i) - c2(i) < eps, c1(i) - c2(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        # global clock to local clock mappings are ordered
        s.add(
            And(
                [
                    And(
                        [
                            Implies(i <= j, And(c0(i) <= c0(j), c1(i) <= c1(j), c2(i) <= c2(j)))
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # consistent cut flow
        # _flow = Function("c_flow", IntSort(), RealSort())
        c_flow = Function("c_flow", IntSort(), RealSort())

        # addition
        s.add(
            And(
                [
                    c_flow(i) == (sig0(c0(i)) + sig1(c1(i)) + sig2(c2(i)))
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # violation check
        v = Int("v")
        s.add(And(v >= timestamps_0[0], v <= timestamps_0[-1]))

        s.add(
            Implies(
                And(v >= timestamps_0[0], v <= timestamps_0[-1]),
                z3Interpolate(c_flow, v) < 10
            )
        )

        # print(s.assertions())
        if s.check() == sat:
            m = s.model()
            # out = "%s %s" % (m[test], m[test2])
            # print(m)
            ##print("unsat in segment", i)
            flag = False
            # terminate after unsat
            return flag

        elif i <= segCount:
            ##print("sat in segment", i)
            flag = True

        s.reset()

    return flag

def prog_tp4(eps, segCount, data_0, data_1, data_2, data_3):

    # pad = min(2 * eps, len(data_0) - 1)
    pad = 2
    peps = pad * eps

    # initialize signal duration and segment duration
    sigDur = len(data_0) - 1
    segDur = sigDur / segCount

    # check for edge cases
    t0 = data_0[0][0]
    t1 = data_0[1][0]

    if sigDur / segCount < t1:

        segCount = sigDur / t1

    if t0 != 0:
        return

    var_0 = {0 : True}
    for i in range(1, len(data_0)):
        if (data_0[i][1] != data_0[i - 1][1]):
            var_0[int(data_0[i][0])] = True
    var_0[int(data_0[-1][0])] = True

    var_1 = {0 : True}
    for i in range(1, len(data_1)):
        if (data_1[i][1] != data_1[i - 1][1]):
            var_1[int(data_1[i][0])] = True
    var_1[int(data_1[-1][0])] = True

    var_2 = {0 : True}
    for i in range(1, len(data_2)):
        if (data_2[i][1] != data_2[i - 1][1]):
            var_2[int(data_2[i][0])] = True
    var_2[int(data_2[-1][0])] = True

    var_3 = {0 : True}
    for i in range(1, len(data_3)):
        if (data_3[i][1] != data_3[i - 1][1]):
            var_2[int(data_3[i][0])] = True
    var_3[int(data_3[-1][0])] = True

    # encoding
    i = 0
    entryFound = True
    flag = True

    while entryFound:

        # Terminate if segment count is reached
        if i == segCount:
            return flag

        # flag to be set True if at least one entry is found in the current iteration
        entryFound = False

        # initialize solver
        s = Solver()

        # calculate upper and lower time bound for current segment
        segmentLowerBound = max(0, int((i * segDur) - eps))
        # segmentUpperBound = int((i + 1) * segDur)
        segmentUpperBound = min(sigDur, int((i + 1) * segDur + eps))

        timestamps_0 = []
        segvar_0 = []
        sig0 = Function("sig0", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_0) - 1:

                s.add(sig0(j) == data_0[j][1])
                for k in range(pad):
                    timestamps_0.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_0[j][0] in var_0 and data_0[j][0] > segmentLowerBound)):
                    segvar_0.append(int(data_0[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_0) > 0 and segvar_0[-1] < segmentUpperBound):
            segvar_0.append(segmentUpperBound)

        timestamps_1 = []
        segvar_1 = []
        sig1 = Function("sig1", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_1) - 1:

                s.add(sig1(j) == data_1[j][1])
                for k in range(pad):
                    timestamps_1.append(j * pad + k)

                # if (j > 0 and data_1[j][1] != data_1[j - 1][1]):
                #     var_1.append(j)

                if (j == segmentLowerBound or (data_1[j][0] in var_1 and data_1[j][0] > segmentLowerBound)):
                    segvar_1.append(int(data_1[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_1) > 0 and segvar_1[-1] < segmentUpperBound):
            segvar_1.append(segmentUpperBound)              


        timestamps_2 = []
        segvar_2 = []
        sig2 = Function("sig2", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_2) - 1:

                s.add(sig2(j) == data_2[j][1])
                for k in range(pad):
                    timestamps_2.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_2[j][0] in var_2 and data_2[j][0] > segmentLowerBound)):
                    segvar_2.append(int(data_2[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_2) > 0 and segvar_2[-1] < segmentUpperBound):
            segvar_2.append(segmentUpperBound)


        timestamps_3 = []
        segvar_3 = []
        sig3 = Function("sig3", IntSort(), RealSort())

        for j in range(segmentLowerBound, segmentUpperBound + 1):

            if j >= 0 and j < len(data_3) - 1:

                s.add(sig3(j) == data_3[j][1])
                for k in range(pad):
                    timestamps_3.append(j * pad + k)

                # if (j > 0 and data_0[j][1] != data_0[j - 1][1]):
                #     var_0.append(j)

                if (j == segmentLowerBound or (data_3[j][0] in var_3 and data_3[j][0] > segmentLowerBound)):
                    segvar_3.append(int(data_3[j][0]))

                if j > (i * segDur):

                    entryFound = True

                if j == sigDur:

                    entryFound = False

        if (len(segvar_3) > 0 and segvar_3[-1] < segmentUpperBound):
            segvar_3.append(segmentUpperBound)

        if not entryFound:

            ##print()
            pass

        i += 1

        # force terminate after one loop
        # entryFound = False

        # global clock to local clock mappings (implicit piece-wise constant interpolation w.r.t. local clock values)
        c0 = Function("c0", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c0(i) >= 0, c0(i) < data_0[-1][0])
                    # Or([c0(i) == (min(timestamps_0[-1]//pad, (max(0, ((i - peps) + j)//pad)))) for j in range(1, 2 * int(peps))])
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c0(min(timestamps_0[-1], max(0, i * pad - peps + j))) == i for j in range(1, 2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                    # for i in range(int(data_0[-1][0]))
                ]
            )
        )

        c1 = Function("c1", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c1(i) >= 0, c1(i) < data_1[-1][0])
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c1(min(timestamps_1[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                ]
            )
        )

        c2 = Function("c2", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c2(i) >= 0, c2(i) < data_2[-1][0])
                    for i in range(timestamps_2[0], timestamps_2[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c2(min(timestamps_2[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                ]
            )
        )

        c3 = Function("c3", IntSort(), IntSort())
        s.add(
            And(
                [
                    And(c3(i) >= 0, c3(i) < data_3[-1][0])
                    for i in range(timestamps_3[0], timestamps_3[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or([c3(min(timestamps_3[-1], max(0, i * pad - peps + j))) == i for j in range(2 * int(peps) + pad)])
                    for i in range(segmentLowerBound, segmentUpperBound)
                ]
            )
        )
        
        # retiming preserves signal variability
        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_0[i] <= c0(j), c0(j) < segvar_0[i + 1])
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_0) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_1[i] <= c1(j), c1(j) < segvar_1[i + 1])
                            for j in range(timestamps_1[0], timestamps_1[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_1) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_2[i] <= c2(j), c2(j) < segvar_2[i + 1])
                            for j in range(timestamps_2[0], timestamps_2[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_2) - 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Or(
                        [
                            And(segvar_3[i] <= c3(j), c3(j) < segvar_3[i + 1])
                            for j in range(timestamps_3[0], timestamps_3[-1] + 1)
                        ]
                    )
                    for i in range(len(segvar_3) - 1)
                ]
            )
        )

        # local clocks are bound by epsilon
        s.add(
            And(
                [
                    And(c0(i) - c1(i) <= eps, c0(i) - c1(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )
        
        s.add(
            And(
                [
                    And(c0(i) - c2(i) <= eps, c0(i) - c2(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    And(c1(i) - c2(i) <= eps, c1(i) - c2(i) >= -eps)
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        s.add(And([And(c0(i) - c3(i) <= eps, c0(i) - c3(i) >= -eps) for i in range(timestamps_3[0], timestamps_3[-1] + 1)]))
        s.add(And([And(c1(i) - c3(i) <= eps, c1(i) - c3(i) >= -eps) for i in range(timestamps_3[0], timestamps_3[-1] + 1)]))
        s.add(And([And(c2(i) - c3(i) <= eps, c2(i) - c3(i) >= -eps) for i in range(timestamps_3[0], timestamps_3[-1] + 1)]))

        # non-concurrent edges stay non-concurrent ???
        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)])
                        ),
                        And(c0(i) - c1(i) < eps, c0(i) - c1(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c2(i) == segvar_2[j], c2(i) != c2(i - 1)) for j in range(1, len(segvar_2) - 1)])
                        ),
                        And(c0(i) - c2(i) < eps, c0(i) - c2(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)]),
                            Or([And(c2(i) == segvar_2[j], c2(i) != c2(i - 1)) for j in range(1, len(segvar_2) - 1)])
                        ),
                        And(c1(i) - c2(i) < eps, c1(i) - c2(i) > -eps)
                    )
                    for i in range(timestamps_1[0], timestamps_1[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c0(i) == segvar_0[j], c0(i) != c0(i - 1)) for j in range(1, len(segvar_0) - 1)]),
                            Or([And(c3(i) == segvar_3[j], c3(i) != c3(i - 1)) for j in range(1, len(segvar_3) - 1)])
                        ),
                        And(c0(i) - c3(i) < eps, c0(i) - c3(i) > -eps)
                    )
                    for i in range(timestamps_3[0], timestamps_3[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c1(i) == segvar_1[j], c1(i) != c1(i - 1)) for j in range(1, len(segvar_1) - 1)]),
                            Or([And(c3(i) == segvar_3[j], c3(i) != c3(i - 1)) for j in range(1, len(segvar_3) - 1)])
                        ),
                        And(c1(i) - c3(i) < eps, c1(i) - c3(i) > -eps)
                    )
                    for i in range(timestamps_3[0], timestamps_3[-1] + 1)
                ]
            )
        )

        s.add(
            And(
                [
                    Implies(
                        And(
                            Or([And(c2(i) == segvar_2[j], c2(i) != c2(i - 1)) for j in range(1, len(segvar_2) - 1)]),
                            Or([And(c3(i) == segvar_3[j], c3(i) != c3(i - 1)) for j in range(1, len(segvar_3) - 1)])
                        ),
                        And(c2(i) - c3(i) < eps, c2(i) - c3(i) > -eps)
                    )
                    for i in range(timestamps_3[0], timestamps_3[-1] + 1)
                ]
            )
        )

        # global clock to local clock mappings are ordered
        s.add(
            And(
                [
                    And(
                        [
                            Implies(i <= j, And(c0(i) <= c0(j), c1(i) <= c1(j), c2(i) <= c2(j), c3(i) <= c3(j)))
                            for j in range(timestamps_0[0], timestamps_0[-1] + 1)
                        ]
                    )
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # consistent cut flow
        # _flow = Function("c_flow", IntSort(), RealSort())
        c_flow = Function("c_flow", IntSort(), RealSort())

        # addition
        s.add(
            And(
                [
                    c_flow(i) == (sig0(c0(i)) + sig1(c1(i)) + sig2(c2(i)) + sig3(c3(i)))
                    for i in range(timestamps_0[0], timestamps_0[-1] + 1)
                ]
            )
        )

        # violation check
        v = Int("v")
        s.add(And(v >= timestamps_0[0], v <= timestamps_0[-1]))

        s.add(
            Implies(
                And(v >= timestamps_0[0], v <= timestamps_0[-1]),
                z3Interpolate(c_flow, v) < 10
            )
        )

        # print(s.assertions())
        if s.check() == sat:
            m = s.model()
            # out = "%s %s" % (m[test], m[test2])
            # print(m)
            ##print("unsat in segment", i)
            flag = False
            # terminate after unsat
            return flag

        elif i <= segCount:
            ##print("sat in segment", i)
            flag = True

        s.reset()

    return flag

def z3Interpolate(f, p):

    # return (f(ToInt(p)) + ((f(ToInt(p) + 1) - f(ToInt(p))) * (p - ToInt(p))))
    return f(p)

def z3Abs(x):

    return If(x >= 0, x, -x)


def z3SqDist3d(x1, x2, y1, y2, z1, z2):

    return (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)))


def z3SqDist2d(x1, x2, y1, y2):

    return (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)))


def z3SqDist1d(x1, x2):

    return z3Abs(x2 - x1)


def getDataTank(agent_ID):
    file = open('{}/../../data/wt/s5_tank_{}'.format(os.path.dirname(os.path.abspath(__file__)), agent_ID))
    line = file.readline()

    data = []

    while line:

        param = line.split('\t')

        values = []

        for i in range(2):

            values.append(float(param[i].strip()))

        # data.append(values)
        data.append([values[0] * mult, values[1]])

        line = file.readline()

    file.close()

    # print(data)

    return data

def main():
    # set repeat count for confidence interval
    repeat = 1
    d = 1

    for n in (2, 3, 4):
        for eps in (1, 2, 4, 8):
            # eps = eps * 0.05

            if (n == 2):
                out = "0"
                flag = False
                data_0 = getDataTank(0)
                data_1 = getDataTank(1)
                t0 = time.time()
                for i in range(repeat):
                    flag = prog_tp2(eps, 1, data_0, data_1)
                t1 = time.time()
                if (flag):
                    out = "1"
                line = str(d) + " " + str(eps) + " " + str((t1 - t0) / repeat) + " " + out
                print(line)
                results = open("results_wt2.txt", "a")
                results.write(line + "\n")
                results.close()

            elif (n == 3):
                out = "0"
                flag = False
                data_0 = getDataTank(0)
                data_1 = getDataTank(1)
                data_2 = getDataTank(2)
                t0 = time.time()
                for i in range(repeat):
                    flag = prog_tp3(eps, 1, data_0, data_1, data_2)
                t1 = time.time()
                if (flag):
                    out = "1"
                line = str(d) + " " + str(eps) + " " + str((t1 - t0) / repeat) + " " + out
                print(line)
                results = open("results_wt3.txt", "a")
                results.write(line + "\n")
                results.close()

            elif (n == 4):
                out = "0"
                flag = False
                data_0 = getDataTank(0)
                data_1 = getDataTank(1)
                data_2 = getDataTank(2)
                data_3 = getDataTank(3)
                t0 = time.time()
                for i in range(repeat):
                    flag = prog_tp4(eps, 1, data_0, data_1, data_2, data_3)
                t1 = time.time()
                if (flag):
                    out = "1"
                line = str(d) + " " + str(eps) + " " + str((t1 - t0) / repeat) + " " + out
                print(line)
                results = open("results_wt4.txt", "a")
                results.write(line + "\n")
                results.close()          

if __name__ == "__main__":
    main()
    pass
