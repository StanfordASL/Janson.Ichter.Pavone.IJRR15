def fill_box(box, d):
    return [box[0]+[0]*d, box[1]+[1]*d]

def box_vol(box):
    return reduce(lambda x, y: x*y, [abs(a-b) for a,b in zip(*box)])

def barriers(d):    # d >= 2
    if d == 1:
        return []
    dboxes = []
    for i in range(d-2):
        a = [0]*(d-1) + [1./3]
        a[i] = 1./3
        b = [1./3]*i + [1]*(d-1-i) + [2./3]
        dboxes.append([a, b])
    dboxes.append([[0]*(d-1) + [1./3], [1./3]*(d-2) + [2./3]*2])
    # print d, sum([box_vol(box) for box in dboxes])*3
    # return sum([[fill_box(box, d-c) for box in barriers(c)] for c in range(2, d)], []) + dboxes
    return [fill_box(box, 1) for box in barriers(d-1)] + dboxes

n = 10
f = open('maze%d'%n, 'w')
f.write(str(n))
for box in barriers(n):
    f.write('\n' + ' '.join([str(x) for x in box[0]]))
    f.write('\n' + ' '.join([str(x) for x in box[1]]))
# for x in barriers(4):
#     print x
# print sum([range(i) for i in range(5)], [])