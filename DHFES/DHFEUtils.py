import copy
import DHFES.DHFEs
from DHFES.DHFEs import *

def DHFEs_distance(d1,d2,lam):
    ## 任意两个对偶犹豫模糊元素的距离
    # d1 表示第一个对偶犹豫模糊元素，d2 表示第二个犹豫模糊元素
    # q 表示 q 级正交模糊，q=1 为对偶犹豫直觉模糊，q=2 为对偶犹豫毕达哥拉斯模糊，q=3 位对偶犹豫费马模糊
    # lam 表示广义对偶犹豫模糊元素距离参数计算，lam=1 为 Hamming distance，lam=2 为 Euclidean distance.
    # 该距离为对偶犹豫费马模糊元素距离公式
    try:
        assert d1.getQRung() == d2.getQRung() and d1.isEmpty()==False and d2.isEmpty()==False         ## 判断是否是同一种对偶犹豫模糊集
        q = d1.getQRung()
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!!',es)
        return None
    
    def normalized(D1,D2):
        ## 标准化，若 len(d1.MD)≠len(d2.MD) 则将隶属度的值数量匹配相等，采用乐观匹配 
        d1 = copy.deepcopy(D1)
        d2 = copy.deepcopy(D2)
        
        lmd = len(d1.MD)-len(d2.MD)
        lnmd = len(d1.NMD)-len(d2.NMD)
        i=0
        if lmd>0:
            while i<lmd:
                d2.MD.append(max(d2.MD))
                i += 1
        else:
            while i<(-lmd):
                d1.MD.append(max(d1.MD))
                i += 1
        j=0
        if lnmd>0:
            while j<lnmd:
                d2.NMD.append(max(d2.NMD))
                j += 1
        else:
            while j<(-lnmd):
                d1.NMD.append(max(d1.NMD))
                j += 1
        # print(d1.MD,d1.NMD)
        # print(d2.MD,d2.NMD)
        # print(len(d1.MD),len(d1.NMD))
        return d1,d2
    
    ## 先对对偶犹豫模糊元素进行排序，打印排序后的元素
    
    d1 = d1.DHFEs_Qsort()
    d2 = d2.DHFEs_Qsort()
    # print("Sorted MD & NMD of DHFE-d1 & DHFE-d2:")
    # print(d1.MD,d1.NMD)
    # print(d2.MD,d2.NMD)
    
    ## 标准化 d1 和 d2，采用乐观匹配
    d1,d2 = normalized(d1,d2)
    # print("Normalized DHFE-d1 & DHFE-d2:")
    # print(d1.MD,d1.NMD)
    # print(d2.MD,d2.NMD)
   
    mds,nmds = 0,0
    m = 1/(len(d1.MD)+len(d1.NMD))
    
    for x in range(len(d1.MD)):
        mds += pow(fabs(pow(d1.MD[x],q)-pow(d2.MD[x],q)),lam)
    for y in range(len(d1.NMD)):
        nmds += pow(fabs(pow(d1.NMD[y],q)-pow(d2.NMD[y],q)),lam)
    distance = format(pow(m*(mds+nmds),1/lam),'.4f')
    
    if lam == 1:
        print('The Hamming distance is '+distance)
    elif lam == 2:
        print('The Euclidean distance is '+distance)
    return float(distance)
    
    

## 对偶犹豫模糊元素支持度
def DHFEs_support(d1,d2,lam):
    # d1 表示第一个对偶犹豫模糊元素,d2 表示第二个犹豫模糊元素,该函数表示 d1 对 d2 的支持度
    # q 表示对偶犹豫 q 级正交模糊
    # lam 表示广义对偶犹豫模糊元素距离参数计算，lam=1 为 Hamming distance，lam=2 为 Euclidean distance.
    return 1-DHFEs_distance(d1,d2,q,lam)



