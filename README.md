# Fuzzy-Number-4Operations





### The definition of basic fuzzy numbers and basic operations are included.

1. Definition of Intuitionistic fuzzy numbers and the four basic operations.
2. Definition of interval-valued Intuitionistic fuzzy numbers and the four basic operations.
3. Definition of Pythagorean fuzzy numbers and the four basic operations.
4. Definition of interval-valued Pythagorean fuzzy numbers and the four basic operations.
5. Definition of Fermatean fuzzy numbers and the four basic operations.
6. Definition of interval-valued Fermatean fuzzy numbers and the four basic operations.
7. Interaction Intuitionistic fuzzy basic operations.
8. Interaction interval-valued Intuitionistic fuzzy basic operations.
9. Interaction Pythagorean fuzzy basic operations.
10. Interaction interval-valued Pythagorean fuzzy basic operations.
11. Interaction Fermatean fuzzy basic operations.
12. Interaction interval-valued Fermatean fuzzy basic operations.

operators based on Archimedes' t-parametrization, including algebraic t-parametrization and t-remaining parametrization and Einstein t-parametrization and t-remaining parametrization



### Example of the Definition and Operations of Fuzzy Number

PFN,IVPFN,FFN,IVFFN are same as IFN and IVIFN.

```python
ifn = IFN(0.5,0.3)							#the definition of Intuitionistic fuzzy number
ivifn = IVIFN(0.5,0.7,0.1,0.3)				#the definition of Interval-Valued Intuitionistic fuzzy number

ifn.show()									#show the Intuitionistic fuzzy number
>>(0.5, 0.3)
ivifn.show()								#show the Interval-Valued Intuitionistic fuzzy number
>>([0.5, 0.7], [0.1, 0.3])

# Intuitionistic fuzzy number Power and Times operations, show() means show the IFN.
ifn.Algebraic_Power(3).show()				#3 powers of the IFN in Algebraic t-norm & t-conorm
ifn.Algebraic_Times(3).show()				#3 times of the IFN in Algebraic t-norm & t-conorm
ifn.Einstein_Power(3).show()				#3 powers of the IFN in Einstein r-norm & t-conorm
ifn.Einstein_Times(3).show()				#3 times of the IFN in Einstein r-norm & t-conorm
>>(0.125, 0.657)
>>(0.875, 0.026999999999999993)
>>(0.07142857142857145, 0.7299212598425198)
>>(0.9285714285714286, 0.010931174089068817)

# Intered-valued Intuitionistic fuzzy number Power and Times operations, show() means show the IVIFN.
ivifn.Algebraic_Power(3).show()				#3 powers of the IVIFN in Algebraic t-norm & t-conorm
ivifn.Algebraic_Times(3).show()				#3 times of the IVIFN in Algebraic t-norm & t-conorm
ivifn.Einstein_Power(3).show()				#3 powers of the IVIFN in Einstein r-norm & t-conorm
ivifn.Einstein_Times(3).show()				#3 times of the IVIFN in Einstein r-norm & t-conorm
>>([0.125, 0.34299999999999997], [0.2709999999999999, 0.657])
>>([0.875, 0.973], [0.001, 0.026999999999999993])
>>([0.07142857142857145, 0.2700787401574802], [0.29223300970873795, 0.7299212598425198])
>>([0.9285714285714286, 0.9890688259109311], [0.00029154518950437334, 0.010931174089068817])
```

Interaction operations are similar to the above.



### Next

List all method descriptions.



### Updating...

------



### 包括了基本模糊数的定义和基本运算：

1. 直觉模糊数的定义和四项基本运算；
2. 区间值直觉模糊数的定义和四项基本运算；
3. 毕达哥拉斯模糊数的定义和四项基本运算；
4. 区间值毕达哥拉斯模糊数的定义和四项基本运算；
5. 费马模糊数的定义和四项基本运算；
6. 区间值费马模糊数的定义和四项基本运算。
7. 交互直觉模糊基本运算
8. 交互区间值直觉模糊基本运算
9. 交互毕达哥拉斯模糊基本运算
10. 交互区间值毕达哥拉斯模糊基本运算
11. 交互费马模糊基本运算
12. 交互区间值费马模糊基本运算

基于阿基米德 t-范数下的运算法则，包括代数 t-范数和 t-余范数以及爱因斯坦 t-范数和 t-余范数。

### 模糊数定义和运算举例

毕达哥拉斯模糊数，区间值毕达哥拉斯模糊数，费马模糊数，区间值费马模糊数运算与直觉模糊，区间值直觉模糊数一样。

```python
ifn = IFN(0.5,0.3)							#直觉模糊数定义
ivifn = IVIFN(0.5,0.7,0.1,0.3)				#区间值直觉模糊数定义

ifn.show()									#显示直觉模糊数
>>(0.5, 0.3)
ivifn.show()								#显示区间值直觉模糊数
>>([0.5, 0.7], [0.1, 0.3])

# 直觉模糊数的幂和倍运算，show()函数表示显示结果
ifn.Algebraic_Power(3).show()				# IFN 在代数范数下的 3 次幂运算
ifn.Algebraic_Times(3).show()				# IFN 在代数范数下的 3 倍运算
ifn.Einstein_Power(3).show()				# IFN 在爱因斯坦范数下的 3 次幂运算
ifn.Einstein_Times(3).show()				# IFN 在爱因斯坦范数下的 3 倍运算
>>(0.125, 0.657)
>>(0.875, 0.026999999999999993)
>>(0.07142857142857145, 0.7299212598425198)
>>(0.9285714285714286, 0.010931174089068817)

# 区间值直觉模糊数的幂和倍运算，show()函数表示显示结果
ivifn.Algebraic_Power(3).show()				# IVIFN 在代数范数下的 3 次幂运算
ivifn.Algebraic_Times(3).show()				# IVIFN 在代数范数下的 3 倍运算
ivifn.Einstein_Power(3).show()				# IVIFN 在爱因斯坦范数下的 3 次幂运算
ivifn.Einstein_Times(3).show()				# IVIFN 在爱因斯坦范数下的 3 倍运算
>>([0.125, 0.34299999999999997], [0.2709999999999999, 0.657])
>>([0.875, 0.973], [0.001, 0.026999999999999993])
>>([0.07142857142857145, 0.2700787401574802], [0.29223300970873795, 0.7299212598425198])
>>([0.9285714285714286, 0.9890688259109311], [0.00029154518950437334, 0.010931174089068817])
```

交互运算和以上运算一样。

### 随后

列出所有的方法说明。

### 更新中...
