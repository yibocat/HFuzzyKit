import copy

import numpy as np
np.set_printoptions(suppress=True) ## 非科学计数法
# np.set_printoptions(precision=4)

################# Normalization ##################
def normalization(D1, D2, tao):
    """
        输入：
            D1: 待标准化的犹豫模糊元素
            D2: 待标准化的犹豫模糊元素
            tao:参数，决策者的风险系数，当 tao=1 时为乐观标准化；tao=0 时为悲观标准化
    """

    def adj(d, tao):
        """
            内置函数，用于计算风险因素得到的标准化元素
            输入：
                d: array 类型
        """
        assert 0 <= tao <= 1, 'The parameter must be in the interval 0-1.'
        return tao * d.max() + (1 - tao) * d.min()

    d1 = copy.deepcopy(D1)
    d2 = copy.deepcopy(D2)

    md_len = len(d1.MD) - len(d2.MD)
    nmd_len = len(d1.NMD) - len(d2.NMD)
    if md_len > 0:  # 说明 d1 隶属度元素个数大于 d2, 需要增加 d2 隶属度的元素
        i = 0
        m = d2.MD
        while i < md_len:
            d2.MD = np.append(d2.MD, adj(m, tao))
            i += 1
    else:  # 说明 d1 隶属度元素个数小于 d2， 需要增加 d1 隶属度的元素
        i = 0
        m = d1.MD
        while i < (-md_len):
            d1.MD = np.append(d1.MD, adj(m, tao))
            i += 1

    if nmd_len > 0:
        i = 0
        m = d2.NMD
        while i < nmd_len:
            d2.NMD = np.append(d2.NMD, adj(m, tao))
            i += 1
    else:
        i = 0
        m = d1.NMD
        while i < (-nmd_len):
            d1.NMD = np.append(d1.NMD, adj(m, tao))
            i += 1

    d1 = d1.DHFEs_Qsort()
    d2 = d2.DHFEs_Qsort()

    return d1, d2


############ Destance function ############
def stand_distance(d1, d2, lam=1, tao=1):
    """
        犹豫模糊元素距离方法
            d1 表示第一个犹豫模糊元素，d2 表示第二个犹豫模糊元素
            lam 表示广义犹豫模糊元素距离计算参数，lam=1 表示 Hamming distance, lam=2 表示 Euclidean distance
            tao 表示标准化风险系数，默认为 1。tao=0 为悲观标准化，tao=1 表示乐观标准化，且 tao 的值在区间[0,1]
        return 距离
    """

    def _normalization(D1, D2, tao):
        def adj(d, tao):
            """
                内置函数，用于计算风险因素得到的标准化元素
                输入：
                    d: array 类型
            """
            assert 0 <= tao <= 1, 'The parameter must be in the interval 0-1.'
            return tao * d.max() + (1 - tao) * d.min()

        d1 = copy.deepcopy(D1)
        d2 = copy.deepcopy(D2)

        md_len = len(d1.MD) - len(d2.MD)
        nmd_len = len(d1.NMD) - len(d2.NMD)
        if md_len > 0:  # 说明 d1 隶属度元素个数大于 d2, 需要增加 d2 隶属度的元素
            i = 0
            m = d2.MD
            while i < md_len:
                d2.MD = np.append(d2.MD, adj(m, tao))
                i += 1
        else:  # 说明 d1 隶属度元素个数小于 d2， 需要增加 d1 隶属度的元素
            i = 0
            m = d1.MD
            while i < (-md_len):
                d1.MD = np.append(d1.MD, adj(m, tao))
                i += 1

        if nmd_len > 0:
            i = 0
            m = d2.NMD
            while i < nmd_len:
                d2.NMD = np.append(d2.NMD, adj(m, tao))
                i += 1
        else:
            i = 0
            m = d1.NMD
            while i < (-nmd_len):
                d1.NMD = np.append(d1.NMD, adj(m, tao))
                i += 1

        d1 = d1.DHFEs_Qsort()
        d2 = d2.DHFEs_Qsort()

        return d1, d2

    assert 0 <= tao <= 1, 'ERROR: The parameter \'tao\' must be in the interval 0-1'
    assert d1.qrung == d2.qrung and d1.isEmpty() == False and d2.isEmpty() == False, 'ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!'
    q = d1.qrung

    # 标准化
    d1, d2 = _normalization(d1, d2, tao)

    # 不确定度
    pi1 = (1 - (d1.MD ** q).sum() / len(d1.MD) - (d1.NMD ** q).sum() / len(d1.NMD)) ** (1 / q)
    pi2 = (1 - (d2.MD ** q).sum() / len(d2.MD) - (d2.NMD ** q).sum() / len(d2.NMD)) ** (1 / q)
    pi = np.fabs(pi1 ** q - pi2 ** q) ** lam

    mds, nmds = 0, 0
    # m = 1/(len(d1.MD)+len(d1.NMD))

    for x in range(len(d1.MD)):
        mds += np.fabs(d1.MD[x] ** q - d2.MD[x] ** q) ** lam
    for y in range(len(d1.NMD)):
        nmds += np.fabs(d1.NMD[y] ** q - d2.NMD[y] ** q) ** lam

    mds = mds / len(d1.MD)
    nmds = nmds / len(d1.NMD)

    distance = 0.5 * (mds + nmds + pi) ** (1 / lam)

    return distance

def hausdorff_distance(d1,d2,lam=1,tao=1):
    ## 广义 Hausdorff 距离，lam为参数， method 表示乐观标准化和悲观标准化，默认为乐观标准化
    assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False,'ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!'         ## 判断是否是同一种对偶犹豫模糊集
    q = d1.qrung
    
    def _normalization(D1, D2, tao):
        def adj(d, tao):
            """
                内置函数，用于计算风险因素得到的标准化元素
                输入：
                    d: array 类型
            """
            assert 0 <= tao <= 1, 'The parameter must be in the interval 0-1.'
            return tao * d.max() + (1 - tao) * d.min()

        d1 = copy.deepcopy(D1)
        d2 = copy.deepcopy(D2)

        md_len = len(d1.MD) - len(d2.MD)
        nmd_len = len(d1.NMD) - len(d2.NMD)
        if md_len > 0:  # 说明 d1 隶属度元素个数大于 d2, 需要增加 d2 隶属度的元素
            i = 0
            m = d2.MD
            while i < md_len:
                d2.MD = np.append(d2.MD, adj(m, tao))
                i += 1
        else:  # 说明 d1 隶属度元素个数小于 d2， 需要增加 d1 隶属度的元素
            i = 0
            m = d1.MD
            while i < (-md_len):
                d1.MD = np.append(d1.MD, adj(m, tao))
                i += 1

        if nmd_len > 0:
            i = 0
            m = d2.NMD
            while i < nmd_len:
                d2.NMD = np.append(d2.NMD, adj(m, tao))
                i += 1
        else:
            i = 0
            m = d1.NMD
            while i < (-nmd_len):
                d1.NMD = np.append(d1.NMD, adj(m, tao))
                i += 1

        d1 = d1.DHFEs_Qsort()
        d2 = d2.DHFEs_Qsort()

        return d1, d2

    ## 先对对偶犹豫模糊元素进行排序，打印排序后的元素
    d1 = d1.DHFEs_Qsort()
    d2 = d2.DHFEs_Qsort()
    
    # 标准化
    d1, d2 = _normalization(d1, d2, tao)
    
    mds,nmds=0,0
    for x in range(len(d1.MD)):
        mds = mds if mds>np.fabs(d1.MD[x]-d2.MD[x])**lam else np.fabs(d1.MD[x]-d2.MD[x])**lam
    for y in range(len(d1.NMD)):
        nmds = nmds if nmds>np.fabs(d1.NMD[y]-d2.NMD[y])**lam else np.fabs(d1.NMD[y]-d2.NMD[y])**lam
    distance = max(mds,nmds)
    
    # if lam == 1:
    #     print('The '+method+' Hamming-Hausdorff distance is '+format(distance,'.4f'))
    # elif lam == 2:
    #     print('The '+method+' Euclidean-Hausdorff distance is '+format(distance,'.4f'))
    return distance

def DHFEs_Distance(d1,d2,distance='Standard',lam=1,tao=1):
    ## 对偶犹豫模糊距离公式
    ## d1,d2表示两个对偶犹豫模糊元素,distance 为距离方法,可选'Hausdorff',默认'Standard',lam 为距离参数,1:Hamming，2:Euclidean
    ## method 表示使用哪种标准化方法，'opt'为乐观标准化，'pess'悲观标准化

    # 标准化
    d1, d2 = _normalization(d1, d2, tao)
    
    if distance=='Standard':
        return stand_distance(d1,d2,lam,tao)
    elif distance=='Hausdorff':
        return hausdorff_distance(d1,d2,lam,tao)
    else:
        print('ERROR distance, please select the right distance!')
        return -1

## 对偶犹豫模糊元素支持度
def DHFEs_support(d1,d2,distance='Standard',lam=1,tao=1):
    # d1 表示第一个对偶犹豫模糊元素,d2 表示第二个犹豫模糊元素,该函数表示 d1 对 d2 的支持度
    # q 表示对偶犹豫 q 级正交模糊
    # lam 表示广义对偶犹豫模糊元素距离参数计算，lam=1 为 Hamming distance，lam=2 为 Euclidean distance.
    return 1-DHFEs_distance(d1,d2,distance,lam,tao)


################# Correlation coefficient  ##################
## Correlation coefficient
def Corr_coefficient_1(d1,d2,tao=1):
    '''
        对偶犹豫模糊元素相关系数 C1
        这里是对偶犹豫模糊元素的相关系数，与对偶犹豫模糊集是不一样的
        
        reference:  
            1.Ye, J. Correlation coefficient of dual hesitant 
        fuzzy sets and its application to multiple attribute decision 
        making. Appl Math Model 38, 659–666 (2014).
            2.Xu,Z., Zhao,H. ,《犹豫模糊集理论及应用》,160-161(2021)
            3.Wang, L., Ni, M. & Zhu, L. Correlation Measures of Dual 
        Hesitant Fuzzy Sets. J Appl Math 2013, 1–12 (2013).
    '''
    assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False,'ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!'         ## 判断是否是同一种对偶犹豫模糊集
    q = d1.qrung
    
    def corr(d1,d2):
        ## d1和d2的相关性
        ## 这里隶属度和非隶属度采用直接加的形式，有待商榷
        assert len(d1.MD)==len(d2.MD) and len(d1.NMD)==len(d2.NMD)
        sum_md,sum_nmd = 0,0
        for i in range(len(d1.MD)):
            sum_md += d1.MD[i]*d2.MD[i]
        for j in range(len(d1.NMD)):
            sum_nmd += d1.NMD[j]*d2.NMD[j]
        return sum_md/len(d1.MD)+sum_nmd/len(d1.NMD)                 ## 这里隶属度和非隶属度采用直接加的形式，有待商榷
    
    d1 = d1.DHFEs_Qsort()
    d2 = d2.DHFEs_Qsort()
    
    # 标准化
    d1, d2 = _normalization(d1, d2, tao)
    
    corr_coe = corr(d1,d2)/(corr(d1,d1)*corr(d2,d2)**(1/2))         ## 相关系数公式
    return corr_coe

def Corr_coefficient_2(d1,d2,tao=1):
    '''
        对偶犹豫模糊元素相关系数 C2
        分子采用 C1 同样的相关性，分母采用最大值形式
    '''
    assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False,'ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!'         ## 判断是否是同一种对偶犹豫模糊集
    q = d1.qrung
    
    def corr(d1,d2):
        ## d1和d2的相关性
        assert len(d1.MD)==len(d2.MD) and len(d1.NMD)==len(d2.NMD)
        sum_md,sum_nmd = 0,0
        for i in range(len(d1.MD)):
            sum_md += d1.MD[i]*d2.MD[i]
        for j in range(len(d1.NMD)):
            sum_nmd += d1.NMD[j]*d2.NMD[j]
        return sum_md/len(d1.MD)+sum_nmd/len(d1.NMD)
    
    d1 = d1.DHFEs_Qsort()
    d2 = d2.DHFEs_Qsort()
    
    # 标准化
    d1, d2 = _normalization(d1, d2, tao)
    
    corr_coe = corr(d1,d2)/max(corr(d1,d1),corr(d2,d2))          ## 相关系数公式
    return corr_coe

## Some bugs no fix ##
# def Corr_coefficient_3(d1,d2,tao=1):
#     '''
#         对偶犹豫模糊元素相关系数 C3
#         该相关系数相关性采用方差形式
#         该相关系数形式有较大问题，当隶属度和非隶属度的数量较少时，根据乐观与悲观匹配规则，
#             有可能造成相关性为 0 的情况，导致分母为 0，舍弃
#     '''
#     try:
#         assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False         ## 判断是否是同一种对偶犹豫模糊集
#         q = d1.qrung
#     except AssertionError as es:
#         print('ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!!',es)
#         return None
    
#     def corr(d1,d2):
#         assert len(d1.MD)==len(d2.MD) and len(d1.NMD)==len(d2.NMD)
#         ## 求平均
#         ave_md1,ave_md2,ave_nmd1,ave_nmd2 = 0,0,0,0
#         for x1 in d1.MD:
#             ave_md1 += x1
#         ave_md1 = ave_md1/len(d1.MD)
#         for x2 in d2.MD:
#             ave_md2 += x2
#         ave_md2 = ave_md2/len(d2.MD)
#         for x3 in d1.NMD:
#             ave_nmd1 += x3
#         ave_nmd1 = ave_nmd1/len(d1.NMD)
#         for x4 in d2.NMD:
#             ave_nmd2 += x4
#         ave_nmd2 = ave_nmd2/len(d2.NMD)
        
#         sum_md,sum_nmd=0,0
#         for i in range(len(d1.MD)):
#             sum_md += (d1.MD[i]-ave_md1)*(d2.MD[i]-ave_md2)
#         for j in range(len(d1.NMD)):
#             sum_nmd += (d1.NMD[j]-ave_nmd1)*(d2.NMD[j]-ave_nmd2)
#         return sum_md/len(d1.MD)+sum_nmd/len(d1.NMD)
    
#     d1 = d1.DHFEs_Qsort()
#     d2 = d2.DHFEs_Qsort()
    
#     ## 标准化 d1 和 d2，默认采用乐观匹配
#     if method=='opt':
#         d1,d2 = opt_normalized(d1,d2)
#     elif method=='pess':
#         d1,d2 = pess_normalized(d1,d2)
#     else:
#         print('ERROR method, please select \'opt\' or \'pess\'!')
#         return -1
    
#     corr_coe = corr(d1,d2)/(corr(d1,d1)*corr(d2,d2)**(1/2)         ## 相关系数公式
#     return corr_coe

# def Corr_coefficient_4(d1,d2,tao=1):
#     '''
#         对偶犹豫模糊元素相关系数 C4
#         该相关系数相关性采用方差形式
#         与相关系数 C3 同样的问题
#         该相关系数形式有较大问题，当隶属度和非隶属度的数量较少时，根据乐观与悲观匹配规则，
#             有可能造成相关性为 0 的情况，导致分母为 0，舍弃
#     '''
#     try:
#         assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False         ## 判断是否是同一种对偶犹豫模糊集
#         q = d1.qrung
#     except AssertionError as es:
#         print('ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!!',es)
#         return None
    
#     def corr(d1,d2):
#         assert len(d1.MD)==len(d2.MD) and len(d1.NMD)==len(d2.NMD)
#         ## 求平均
#         ave_md1,ave_md2,ave_nmd1,ave_nmd2 = 0,0,0,0
#         for x1 in d1.MD:
#             ave_md1 += x1
#         ave_md1 = ave_md1/len(d1.MD)
#         for x2 in d2.MD:
#             ave_md2 += x2
#         ave_md2 = ave_md2/len(d2.MD)
#         for x3 in d1.NMD:
#             ave_nmd1 += x3
#         ave_nmd1 = ave_nmd1/len(d1.NMD)
#         for x4 in d2.NMD:
#             ave_nmd2 += x4
#         ave_nmd2 = ave_nmd2/len(d2.NMD)
        
#         sum_md,sum_nmd=0,0
#         for i in range(len(d1.MD)):
#             sum_md += (d1.MD[i]-ave_md1)*(d2.MD[i]-ave_md2)
#         for j in range(len(d1.NMD)):
#             sum_nmd += (d1.NMD[j]-ave_nmd1)*(d2.NMD[j]-ave_nmd2)
#         return sum_md/len(d1.MD)+sum_nmd/len(d1.NMD)
    
#     d1 = d1.DHFEs_Qsort()
#     d2 = d2.DHFEs_Qsort()
    
#     ## 标准化 d1 和 d2，默认采用乐观匹配
#     if method=='opt':
#         d1,d2 = opt_normalized(d1,d2)
#     elif method=='pess':
#         d1,d2 = pess_normalized(d1,d2)
#     else:
#         print('ERROR method, please select \'opt\' or \'pess\'!')
#         return -1
    
#     corr_coe = corr(d1,d2)/max(corr(d1,d1),corr(d2,d2))         ## 相关系数公式
#     return corr_coe

# def Corr_coefficient_5(d1,d2,tao=1):
#     '''
#         对偶犹豫模糊元素相关系数 C4
#         Reference: 
#             1.Wang, L., Ni, M. & Zhu, L. Correlation Measures of Dual Hesitant Fuzzy Sets. J Appl Math 2013, 1–12 (2013).
#         同样有问题，会出现分母为 0 的情况，比如 d1=[[0.9],[0.3]],d2=[[0.9],[0.2]]
#     '''
#     try:
#         assert d1.qrung == d2.qrung and d1.isEmpty()==False and d2.isEmpty()==False         ## 判断是否是同一种对偶犹豫模糊集
#         q = d1.qrung
#     except AssertionError as es:
#         print('ERROR! The two DHFEs are not the same DHFE or one of them is a empty DHFE!!',es)
#         return None
    
#     d1 = d1.DHFEs_Qsort()
#     d2 = d2.DHFEs_Qsort()
    
#     ## 标准化 d1 和 d2，默认采用乐观匹配
#     if method=='opt':
#         d1,d2 = opt_normalized(d1,d2)
#     elif method=='pess':
#         d1,d2 = pess_normalized(d1,d2)
#     else:
#         print('ERROR method, please select \'opt\' or \'pess\'!')
#         return -1
#     # print(d1.show())
#     # print(d2.show())
    
#     ## membership degree
#     max_md,min_md = 0,0
#     for i in range(len(d1.MD)):
#         max_md = max_md if max_md>=fabs(d1.MD[i]-d2.MD[i]) else fabs(d1.MD[i]-d2.MD[i])
#         min_md = min_md if min_md<=fabs(d1.MD[i]-d2.MD[i]) else fabs(d1.MD[i]-d2.MD[i])
    
#     max_nmd,min_nmd = 0,0
#     for j in range(len(d1.NMD)):
#         max_nmd = max_nmd if max_nmd>=fabs(d1.NMD[j]-d2.NMD[j]) else fabs(d1.NMD[j]-d2.NMD[j])
#         min_nmd = min_nmd if min_nmd<=fabs(d1.NMD[j]-d2.NMD[j]) else fabs(d1.NMD[j]-d2.NMD[j])
        
#     # print(max_md,min_md)
#     # print(max_nmd,min_nmd)
    
#     sum_md,sum_nmd = 0,0
#     for i in range(len(d1.MD)):
#         sum_md += (min_md+max_md)/(fabs(d1.MD[i]-d2.MD[i])+max_md)
#     for j in range(len(d1.NMD)):
#         sum_nmd += (min_nmd+max_nmd)/(fabs(d1.NMD[j]-d2.NMD[j])+max_nmd)
#     return ((1/len(d1.MD))*sum_md+(1/len(d1.NMD))*sum_nmd)/2

###########  DHFEs Entropy and similarity measure ###############
def DHFE_Entropy(d,distance='Standard',lam=1,tao=1):
    ## 初等对偶犹豫模糊熵
    ## lam 为距离参数,当 lam=1,距离使用 Hamming 距离;lam=2,为 Euclidean 距离
    ## 和距离公式一样，distance 和 method 表示不同的距离公式，从而生成不同的对偶犹豫模糊熵
    return 1-DHFEs_Distance(d,d.comp(),distance,lam,tao)


########### generate a random DHFE ##############
from .DHIFE import DHIFE
from .DHPFE import DHPFE
from .DHFFE import DHFFE
def randomDHFE(q,n=5):
    ## q 级正交对偶犹豫模糊随机数生成
    ## n 表示生成的数的最大随机个数，默认为 5
    md = np.random.rand(np.random.randint(1,n))
    nmd = np.random.rand(np.random.randint(1,n))

    if q==1:
        newDHFE = DHIFE([],[])
    elif q==2:
        newDHFE = DHPFE([],[])
    elif q==3:
        newDHFE = DHFFE([],[])

    if max(md)**q+max(nmd)**q <= 1:
        newDHFE.MD = md
        newDHFE.NMD = nmd
        return newDHFE
    else:
        return randomDHFE(q,n)