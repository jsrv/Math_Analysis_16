def brute(l):
    leastdis = (l[0][0]-l[1][0])**2 + (l[0][1]-l[1][1])**2
    for i in xrange(len(l)):
        for j in xrange(len(l)):
            if j > i:
                dis = (l[i][0]-l[j][0])**2 + (l[i][1]-l[j][1])**2
                if dis <= leastdis:
                    leastdis = dis
                    index = [i,j]
    return l[index[0]],l[index[1]]

if __name__ == "__main__":
    l = [[3,2],[3,3],[4,5],[-2,4],[5,7],[-3,4]]
    print brute(l)
