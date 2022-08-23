import copy
import DHFEs
from DHFEs import *

def DHFEs_distance(d1,d2,q,lam):
    ## 任意两个对偶犹豫模糊元素的距离
    # d1 表示第一个对偶犹豫模糊元素，d2 表示第二个犹豫模糊元素
    # q 表示 q 级正交模糊，q=1 为对偶犹豫直觉模糊，q=2 为对偶犹豫毕达哥拉斯模糊，q=3 位对偶犹豫费马模糊
    # lam 表示广义对偶犹豫模糊元素距离参数计算，lam=1 为 Hamming distance，lam=2 为 Euclidean distance.
    # 该距离为对偶犹豫费马模糊元素距离公式
    
    def DHFEs_Qsort(dhfes,q):
        ## 对偶犹豫模糊元素隶属度非隶属度快速排序算法
        
        def quick_sort(lists,i,j):
            if i >= j:
                return list
            pivot = lists[i]
            low = i
            high = j
            while i < j:
                while i < j and lists[j] >= pivot:
                    j -= 1
                lists[i]=lists[j]
                while i < j and lists[i] <=pivot:
                    i += 1
                lists[j]=lists[i]
            lists[j] = pivot
            quick_sort(lists,low,i-1)
            quick_sort(lists,i+1,high)
            return lists
        
        if q == 1:
            x = DHIFE([],[])
        elif q==2:
            x = DHPFE([],[])
        elif q== 3:
            x = DHFFE([],[])
        
        if len(dhfes.MD)>1:
            md = quick_sort(dhfes.MD,0,len(dhfes.MD)-1)
        else:
            md = dhfes.MD
        if len(dhfes.NMD)>1:
            nmd = quick_sort(dhfes.NMD,0,len(dhfes.NMD)-1)
        else:
            nmd = dhfes.NMD
        # print(md,nmd)
        for i in md:
            x.MD.append(i)
        for j in nmd:
            x.NMD.append(j)
        return x
    
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
    
    d1 = DHFEs_Qsort(d1,q)
    d2 = DHFEs_Qsort(d2,q)
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
    return float(format(pow(m*(mds+nmds),1/lam),'.4f'))

## 对偶犹豫模糊元素支持度
def DHFEs_support(d1,d2,q,lam):
    # d1 表示第一个对偶犹豫模糊元素,d2 表示第二个犹豫模糊元素,该函数表示 d1 对 d2 的支持度
    # q 表示对偶犹豫 q 级正交模糊
    # lam 表示广义对偶犹豫模糊元素距离参数计算，lam=1 为 Hamming distance，lam=2 为 Euclidean distance.
    return 1-DHFEs_distance(d1,d2,q,lam)

def DHFFEs_Qsort(dhfes,q):
	## 对偶犹豫模糊元素隶属度与非隶属度排序
    ## q 表示 q-rung
    ## 快速排序
    def quick_sort(lists,i,j):
        if i >= j:
            return list
        pivot = lists[i]
        low = i
        high = j
        while i < j:
            while i < j and lists[j] >= pivot:
                j -= 1
            lists[i]=lists[j]
            while i < j and lists[i] <=pivot:
                i += 1
            lists[j]=lists[i]
        lists[j] = pivot
        quick_sort(lists,low,i-1)
        quick_sort(lists,i+1,high)
        return lists

    if q==1:
        x = DHIFE([],[])
    elif q==2:
        x = DHPFE([], [])
    elif q==3:
        x = DHFFE([],[])

    for i in quick_sort(dhfes.MD,0,len(dhfes.MD)-1):
        x.MD.append(i)
    for j in quick_sort(dhfes.NMD,0,len(dhfes.NMD)-1):
        x.NMD.append(j)
    return x

